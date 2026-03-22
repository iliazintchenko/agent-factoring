/*
 * Combined factoring: trial division + Pollard's rho (Brent) + GMP-ECM
 * Tuned for balanced semiprimes (p*q where p,q have ~equal digit count).
 * Usage: ./factor <N> [timeout_secs]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <signal.h>
#include <unistd.h>
#include <gmp.h>
#include <ecm.h>

static volatile int timed_out = 0;
static void alarm_handler(int sig) { (void)sig; timed_out = 1; }

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

/* Pollard's rho with Brent's cycle detection, batch GCD */
static int try_rho(mpz_t factor, const mpz_t n, unsigned long max_iters) {
    mpz_t y, x, ys, q, tmp;
    mpz_inits(y, x, ys, q, tmp, NULL);
    int found = 0;

    for (unsigned long c = 1; c <= 20 && !found && !timed_out; c++) {
        mpz_set_ui(y, c + 1);
        mpz_set(x, y);
        mpz_set_ui(q, 1);
        unsigned long r = 1, m = 256;

        for (unsigned long total = 0; total < max_iters && !timed_out; ) {
            mpz_set(x, y);
            for (unsigned long i = 0; i < r; i++) {
                mpz_mul(y, y, y); mpz_add_ui(y, y, c); mpz_mod(y, y, n);
            }
            unsigned long k = 0;
            while (k < r && !found && !timed_out) {
                mpz_set(ys, y);
                unsigned long batch = (m < r - k) ? m : (r - k);
                for (unsigned long i = 0; i < batch; i++) {
                    mpz_mul(y, y, y); mpz_add_ui(y, y, c); mpz_mod(y, y, n);
                    mpz_sub(tmp, x, y); mpz_abs(tmp, tmp);
                    mpz_mul(q, q, tmp); mpz_mod(q, q, n);
                }
                mpz_gcd(factor, q, n);
                if (mpz_cmp_ui(factor, 1) > 0) found = 1;
                k += batch;
                total += batch;
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

/* ECM with specified B1, running up to max_curves */
static int try_ecm(mpz_t factor, const mpz_t n, double b1, unsigned long max_curves) {
    ecm_params params;
    ecm_init(params);
    int found = 0;

    for (unsigned long i = 0; i < max_curves && !found && !timed_out; i++) {
        mpz_set(factor, n);
        params->B1done = 1.0;
        mpz_set_ui(params->sigma, 0);
        int ret = ecm_factor(factor, n, b1, params);
        if (ret > 0 && mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, n) < 0) {
            found = 1;
        }
    }
    ecm_clear(params);
    return found;
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
        fprintf(stderr, "Usage: %s <N> [timeout_secs]\n", argv[0]);
        return 1;
    }

    int timeout = 290;
    if (argc >= 3) timeout = atoi(argv[2]);
    signal(SIGALRM, alarm_handler);
    alarm(timeout);
    clock_gettime(CLOCK_MONOTONIC, &g_start);

    mpz_t n, factor;
    mpz_inits(n, factor, NULL);
    if (mpz_set_str(n, argv[1], 10) != 0) {
        fprintf(stderr, "Invalid number\n");
        return 1;
    }

    size_t digits = mpz_sizeinbase(n, 10);
    int target_fd = ((int)digits + 1) / 2; /* expected factor digits for balanced semiprime */
    fprintf(stderr, "Factoring %zu digits (target factor ~%d digits)\n", digits, target_fd);

    /* Stage 1: Trial division to 1M */
    if (try_trial(factor, n, 1000000)) {
        print_result(factor, n, "trial division");
        mpz_clears(n, factor, NULL);
        return 0;
    }

    /* Stage 2: Brief Pollard's rho - only useful for small factors or lucky hits */
    {
        unsigned long iters = (digits <= 40) ? 2000000 : 500000;
        if (try_rho(factor, n, iters)) {
            print_result(factor, n, "Pollard rho");
            mpz_clears(n, factor, NULL);
            return 0;
        }
        fprintf(stderr, "Rho done (%.3fs)\n", elapsed_sec());
    }

    /* Stage 3: ECM - start at appropriate B1 for expected factor size */
    /* B1 table: factor_digits -> (B1, expected_curves) */
    struct { int fdigits; double b1; unsigned long curves; } ecm_plan[] = {
        {15, 2000,       25},
        {20, 11000,      90},
        {25, 50000,      300},
        {30, 250000,     700},
        {35, 1000000,    1700},
        {40, 3000000,    5100},
        {45, 11000000,   10600},
        {50, 43000000,   19300},
        {55, 110000000,  49000},
    };
    int nstages = sizeof(ecm_plan) / sizeof(ecm_plan[0]);

    /* Find starting stage based on target factor size */
    int start_stage = 0;
    for (int i = 0; i < nstages; i++) {
        if (ecm_plan[i].fdigits >= target_fd - 5) {
            start_stage = (i > 0) ? i - 1 : 0;
            break;
        }
        start_stage = i;
    }

    /* Run a few curves at each smaller B1 first (cheap and might get lucky) */
    for (int s = 0; s < start_stage && !timed_out; s++) {
        unsigned long quick_curves = ecm_plan[s].curves / 4;
        if (quick_curves < 5) quick_curves = 5;
        if (try_ecm(factor, n, ecm_plan[s].b1, quick_curves)) {
            char m[64]; snprintf(m, 64, "ECM B1=%.0f", ecm_plan[s].b1);
            print_result(factor, n, m);
            mpz_clears(n, factor, NULL);
            return 0;
        }
    }

    /* Full curves at target and above */
    for (int s = start_stage; s < nstages && !timed_out; s++) {
        fprintf(stderr, "ECM B1=%.0f (%lu curves, %.1fs elapsed)\n",
                ecm_plan[s].b1, ecm_plan[s].curves, elapsed_sec());
        if (try_ecm(factor, n, ecm_plan[s].b1, ecm_plan[s].curves)) {
            char m[64]; snprintf(m, 64, "ECM B1=%.0f", ecm_plan[s].b1);
            print_result(factor, n, m);
            mpz_clears(n, factor, NULL);
            return 0;
        }
    }

    /* If all ECM stages exhausted, loop with highest B1 until timeout */
    while (!timed_out) {
        fprintf(stderr, "ECM extended (B1=110000000, %.1fs)\n", elapsed_sec());
        if (try_ecm(factor, n, 110000000, 1000)) {
            print_result(factor, n, "ECM extended");
            mpz_clears(n, factor, NULL);
            return 0;
        }
    }

    fprintf(stderr, "FAILED after %.3f seconds\n", elapsed_sec());
    mpz_clears(n, factor, NULL);
    return 1;
}
