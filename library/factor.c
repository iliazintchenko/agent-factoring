/*
 * Parallel factoring: trial division + parallel ECM workers
 * Usage: ./factor <N>
 * Forks NCORES workers, each running ECM with random curves.
 * First worker to find a factor writes it to a pipe and all exit.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <signal.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <gmp.h>
#include <ecm.h>

/* Number of parallel workers */
static int get_ncores(void) {
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    return (n > 0) ? (int)n : 4;
}

int trial_div(mpz_t factor, const mpz_t n) {
    if (mpz_even_p(n)) { mpz_set_ui(factor, 2); return 1; }
    for (unsigned long p = 3; p <= 1000000; p += 2) {
        if (mpz_divisible_ui_p(n, p)) { mpz_set_ui(factor, p); return 1; }
    }
    return 0;
}

/*
 * ECM worker function. Runs in a forked child.
 * Writes factor string to fd on success.
 */
void ecm_worker(const char *n_str, double b1, int curves_per_batch, int fd) {
    mpz_t n, f;
    mpz_inits(n, f, NULL);
    mpz_set_str(n, n_str, 10);

    ecm_params params;
    ecm_init(params);
    params->verbose = 0;

    for (int i = 0; i < curves_per_batch; i++) {
        mpz_set(f, n);
        int ret = ecm_factor(f, n, b1, params);
        if (ret > 0 && mpz_cmp_ui(f, 1) > 0 && mpz_cmp(f, n) < 0) {
            /* Found factor - write to pipe */
            char buf[512];
            gmp_sprintf(buf, "%Zd", f);
            write(fd, buf, strlen(buf) + 1);
            ecm_clear(params);
            mpz_clears(n, f, NULL);
            _exit(0);
        }
        mpz_set_ui(params->sigma, 0);
    }
    ecm_clear(params);
    mpz_clears(n, f, NULL);
    _exit(1);
}

int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    mpz_t n, factor;
    mpz_inits(n, factor, NULL);
    if (mpz_set_str(n, argv[1], 10) != 0) {
        fprintf(stderr, "Invalid number: %s\n", argv[1]); return 1;
    }
    int digits = (int)mpz_sizeinbase(n, 10);
    int found = 0;

    /* Trial division */
    found = trial_div(factor, n);
    if (found) goto done;

    /* Parallel ECM */
    {
        int ncores = get_ncores();
        int half = digits / 2; /* expected factor digit count */

        /* B1 levels: each level targets factors of increasing size */
        struct { double b1; int curves_per_worker; } levels[] = {
            {500,       10},    /* 12d */
            {2000,      10},    /* 15d */
            {11000,     15},    /* 20d */
            {50000,     25},    /* 25d */
            {250000,    40},    /* 30d */
            {1000000,   60},    /* 35d */
            {3000000,   80},    /* 40d */
            {11000000,  120},   /* 45d */
            {43000000,  200},   /* 50d */
            {110000000, 300},   /* 55d */
        };
        int n_levels = sizeof(levels) / sizeof(levels[0]);

        /* Start from appropriate level */
        int start = 0;
        if (half > 27) start = 3;
        else if (half > 22) start = 2;
        else if (half > 17) start = 1;

        for (int lvl = start; lvl < n_levels && !found; lvl++) {
            fprintf(stderr, "ECM B1=%.0f x%d workers x%d curves\n",
                    levels[lvl].b1, ncores, levels[lvl].curves_per_worker);

            int pipefd[2];
            if (pipe(pipefd) < 0) { perror("pipe"); return 1; }

            pid_t children[ncores];
            for (int i = 0; i < ncores; i++) {
                pid_t pid = fork();
                if (pid == 0) {
                    close(pipefd[0]);
                    ecm_worker(argv[1], levels[lvl].b1, levels[lvl].curves_per_worker, pipefd[1]);
                    _exit(1); /* should not reach here */
                }
                children[i] = pid;
            }
            close(pipefd[1]);

            /* Read from pipe - blocks until a child writes or all close */
            char buf[512] = {0};
            ssize_t nread = read(pipefd[0], buf, sizeof(buf) - 1);
            close(pipefd[0]);

            if (nread > 0) {
                mpz_set_str(factor, buf, 10);
                found = 1;
            }

            /* Kill all children */
            for (int i = 0; i < ncores; i++) {
                kill(children[i], SIGKILL);
            }
            /* Reap all children */
            for (int i = 0; i < ncores; i++) {
                waitpid(children[i], NULL, 0);
            }
        }
    }

done:
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;

    if (found) {
        mpz_t other; mpz_init(other);
        mpz_divexact(other, n, factor);
        if (mpz_cmp(factor, other) > 0) mpz_swap(factor, other);
        gmp_printf("%Zd\n", factor);
        fprintf(stderr, "Time: %.3f seconds (%d digits)\n", elapsed, digits);
        mpz_clear(other);
    } else {
        fprintf(stderr, "FAILED after %.3f seconds (%d digits)\n", elapsed, digits);
        return 1;
    }
    mpz_clears(n, factor, NULL);
    return 0;
}
