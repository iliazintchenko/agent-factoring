/*
 * Parallel factoring: spawns multiple ECM worker processes.
 * Each worker runs ECM with different random seeds.
 * First worker to find a factor reports it and kills siblings.
 * Usage: ./pfactor <N> [timeout_secs] [num_workers]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <signal.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <gmp.h>
#include <ecm.h>

static struct timespec g_start;
static double elapsed_sec(void) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec - g_start.tv_sec) + (now.tv_nsec - g_start.tv_nsec) / 1e9;
}

/* Trial division up to limit */
static int try_trial(mpz_t factor, const mpz_t n, unsigned long limit) {
    if (mpz_divisible_ui_p(n, 2)) { mpz_set_ui(factor, 2); return 1; }
    for (unsigned long p = 3; p <= limit; p += 2) {
        if (mpz_divisible_ui_p(n, p)) { mpz_set_ui(factor, p); return 1; }
    }
    return 0;
}

/* Pollard's rho with Brent's improvement */
static int try_rho(mpz_t factor, const mpz_t n, unsigned long max_iters) {
    mpz_t y, x, ys, q, tmp;
    mpz_inits(y, x, ys, q, tmp, NULL);
    int found = 0;

    for (unsigned long c = 1; c <= 20 && !found; c++) {
        mpz_set_ui(y, c + 1);
        mpz_set(x, y);
        mpz_set_ui(q, 1);
        unsigned long r = 1, m = 256;

        for (unsigned long total = 0; total < max_iters; ) {
            mpz_set(x, y);
            for (unsigned long i = 0; i < r; i++) {
                mpz_mul(y, y, y); mpz_add_ui(y, y, c); mpz_mod(y, y, n);
            }
            unsigned long k = 0;
            while (k < r && !found) {
                mpz_set(ys, y);
                unsigned long batch = (m < r - k) ? m : (r - k);
                for (unsigned long i = 0; i < batch; i++) {
                    mpz_mul(y, y, y); mpz_add_ui(y, y, c); mpz_mod(y, y, n);
                    mpz_sub(tmp, x, y); mpz_abs(tmp, tmp);
                    mpz_mul(q, q, tmp); mpz_mod(q, q, n);
                }
                mpz_gcd(factor, q, n);
                if (mpz_cmp_ui(factor, 1) > 0) found = 1;
                k += batch; total += batch;
            }
            if (found && mpz_cmp(factor, n) == 0) {
                found = 0;
                mpz_set(y, ys);
                for (unsigned long i = 0; i < m + 10; i++) {
                    mpz_mul(y, y, y); mpz_add_ui(y, y, c); mpz_mod(y, y, n);
                    mpz_sub(tmp, x, y); mpz_abs(tmp, tmp);
                    mpz_gcd(factor, tmp, n);
                    if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, n) != 0) {
                        found = 1; break;
                    }
                }
                if (!found) break;
            }
            if (found) break;
            r *= 2;
        }
    }
    mpz_clears(y, x, ys, q, tmp, NULL);
    return found;
}

/* ECM worker: runs ECM curves at given B1, writes factor to pipe if found */
static void ecm_worker(const char *n_str, double b1, int worker_id, int write_fd, int timeout) {
    signal(SIGALRM, SIG_DFL);
    alarm(timeout);

    mpz_t n, factor;
    mpz_inits(n, factor, NULL);
    mpz_set_str(n, n_str, 10);

    ecm_params params;
    ecm_init(params);

    /* Use worker_id to seed differently */
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, worker_id * 999983 + time(NULL));

    for (int i = 0; ; i++) {
        mpz_set(factor, n);
        params->B1done = 1.0;
        /* Set a random sigma for variety */
        mpz_t sigma;
        mpz_init(sigma);
        mpz_urandomb(sigma, rstate, 64);
        mpz_add_ui(sigma, sigma, 6); /* sigma must be >= 6 */
        mpz_set(params->sigma, sigma);
        params->param = ECM_PARAM_SUYAMA;
        mpz_clear(sigma);

        int ret = ecm_factor(factor, n, b1, params);
        if (ret > 0 && mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, n) < 0) {
            /* Found! Write factor string to pipe */
            char *fstr = mpz_get_str(NULL, 10, factor);
            write(write_fd, fstr, strlen(fstr));
            write(write_fd, "\n", 1);
            free(fstr);
            ecm_clear(params);
            gmp_randclear(rstate);
            mpz_clears(n, factor, NULL);
            _exit(0);
        }
    }
    /* Should not reach here */
    ecm_clear(params);
    gmp_randclear(rstate);
    mpz_clears(n, factor, NULL);
    _exit(1);
}

static void print_result(mpz_t factor, const mpz_t n, const char *method) {
    mpz_t other;
    mpz_init(other);
    mpz_divexact(other, n, factor);
    if (mpz_cmp(factor, other) > 0) mpz_swap(factor, other);
    gmp_printf("%Zd\n", factor);
    fprintf(stderr, "Found by %s in %.3f seconds\n", method, elapsed_sec());
    mpz_clear(other);
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N> [timeout_secs] [num_workers]\n", argv[0]);
        return 1;
    }

    int timeout = 290;
    int num_workers = 8;
    if (argc >= 3) timeout = atoi(argv[2]);
    if (argc >= 4) num_workers = atoi(argv[3]);

    clock_gettime(CLOCK_MONOTONIC, &g_start);

    mpz_t n, factor;
    mpz_inits(n, factor, NULL);
    if (mpz_set_str(n, argv[1], 10) != 0) {
        fprintf(stderr, "Invalid number\n");
        return 1;
    }

    size_t digits = mpz_sizeinbase(n, 10);
    int target_fd = ((int)digits + 1) / 2;
    fprintf(stderr, "Factoring %zu digits (target factor ~%d digits), %d workers\n",
            digits, target_fd, num_workers);

    /* Stage 1: Trial division */
    if (try_trial(factor, n, 1000000)) {
        print_result(factor, n, "trial division");
        mpz_clears(n, factor, NULL);
        return 0;
    }

    /* Stage 2: Brief Pollard's rho */
    {
        unsigned long iters = (digits <= 40) ? 2000000 : 500000;
        if (try_rho(factor, n, iters)) {
            print_result(factor, n, "Pollard rho");
            mpz_clears(n, factor, NULL);
            return 0;
        }
        fprintf(stderr, "Rho done (%.3fs)\n", elapsed_sec());
    }

    /* Stage 3: Parallel ECM */
    /* Choose B1 based on expected factor size */
    struct { int fdigits; double b1; } b1_table[] = {
        {15, 2000},
        {20, 11000},
        {25, 50000},
        {30, 250000},
        {35, 1000000},
        {40, 3000000},
        {45, 11000000},
        {50, 43000000},
        {55, 110000000},
    };
    int ntable = sizeof(b1_table) / sizeof(b1_table[0]);

    /* Find appropriate B1 for target factor size */
    double b1 = 50000; /* default */
    for (int i = 0; i < ntable; i++) {
        if (b1_table[i].fdigits >= target_fd) {
            b1 = b1_table[i].b1;
            break;
        }
        b1 = b1_table[i].b1;
    }

    /* Also try a smaller B1 with some workers (might get lucky faster) */
    double b1_small = b1;
    for (int i = 0; i < ntable; i++) {
        if (b1_table[i].fdigits >= target_fd - 5) {
            b1_small = b1_table[i].b1;
            break;
        }
    }

    fprintf(stderr, "Launching %d ECM workers (B1=%.0f/%.0f)\n", num_workers, b1_small, b1);

    int remaining_timeout = timeout - (int)elapsed_sec() - 2;
    if (remaining_timeout < 5) {
        fprintf(stderr, "Not enough time left for ECM\n");
        mpz_clears(n, factor, NULL);
        return 1;
    }

    /* Create pipes and fork workers */
    pid_t *pids = calloc(num_workers, sizeof(pid_t));
    int *pipes = calloc(num_workers * 2, sizeof(int));

    for (int w = 0; w < num_workers; w++) {
        int pfd[2];
        if (pipe(pfd) < 0) { perror("pipe"); return 1; }
        pipes[w * 2] = pfd[0];     /* read end */
        pipes[w * 2 + 1] = pfd[1]; /* write end */

        pid_t pid = fork();
        if (pid == 0) {
            /* Child: close read ends */
            for (int i = 0; i <= w; i++) close(pipes[i * 2]);
            /* Half workers use smaller B1, half use target B1 */
            double worker_b1 = (w < num_workers / 3) ? b1_small : b1;
            ecm_worker(argv[1], worker_b1, w, pfd[1], remaining_timeout);
            _exit(1);
        }
        close(pfd[1]); /* Parent closes write end */
        pids[w] = pid;
    }

    /* Parent: wait for any child to finish with factor */
    int found = 0;
    while (!found) {
        int status;
        pid_t p = wait(&status);
        if (p < 0) break;

        /* Find which worker finished */
        for (int w = 0; w < num_workers; w++) {
            if (pids[w] == p) {
                if (WIFEXITED(status) && WEXITSTATUS(status) == 0) {
                    /* Read factor from pipe */
                    char buf[256] = {0};
                    int n_read = read(pipes[w * 2], buf, sizeof(buf) - 1);
                    if (n_read > 0) {
                        /* Remove newline */
                        buf[strcspn(buf, "\n")] = 0;
                        mpz_set_str(factor, buf, 10);
                        found = 1;
                    }
                }
                pids[w] = 0;
                close(pipes[w * 2]);
                break;
            }
        }
    }

    /* Kill remaining workers */
    for (int w = 0; w < num_workers; w++) {
        if (pids[w] > 0) {
            kill(pids[w], SIGKILL);
            waitpid(pids[w], NULL, 0);
            close(pipes[w * 2]);
        }
    }

    free(pids);
    free(pipes);

    if (found) {
        char method[64];
        snprintf(method, 64, "parallel ECM (B1=%.0f/%.0f)", b1_small, b1);
        print_result(factor, n, method);
        mpz_clears(n, factor, NULL);
        return 0;
    }

    fprintf(stderr, "FAILED after %.3f seconds\n", elapsed_sec());
    mpz_clears(n, factor, NULL);
    return 1;
}
