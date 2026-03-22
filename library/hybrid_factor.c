/*
 * hybrid_factor.c - Hybrid factoring for 90-digit balanced semiprimes
 *
 * Strategy:
 * 1. Quick ECM shots (5 seconds) - might get lucky if factors aren't perfectly balanced
 * 2. If no factor found, invoke YAFU SIQS with optimal parameters
 *
 * Compile: gcc -O3 -o hybrid_factor library/hybrid_factor.c -lgmp -lecm
 * Usage: ./hybrid_factor <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/wait.h>
#include <unistd.h>
#include <gmp.h>
#include <ecm.h>

#define ECM_TIMEOUT 5   /* seconds for ECM phase */
#define SIQS_TIMEOUT 285 /* seconds for SIQS phase (total 290) */

/* Run ECM with increasing B1 for ECM_TIMEOUT seconds */
int try_ecm(mpz_t n, mpz_t factor) {
    ecm_params params;
    int ret;
    time_t start = time(NULL);

    ecm_init(params);

    /* Seed = 42 as required */
    mpz_set_ui(params->sigma, 42);
    params->method = ECM_ECM;

    /* Try increasing B1 values */
    double b1_values[] = {1000, 5000, 11000, 50000, 250000, 1000000};
    int num_b1 = sizeof(b1_values) / sizeof(b1_values[0]);

    for (int i = 0; i < num_b1; i++) {
        if (difftime(time(NULL), start) >= ECM_TIMEOUT) break;

        params->B1done = 1.0;
        mpz_set_ui(params->sigma, 42 + i);

        ret = ecm_factor(factor, n, b1_values[i], params);

        if (ret > 0 && mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, n) < 0) {
            printf("ECM found factor at B1=%g: ", b1_values[i]);
            mpz_out_str(stdout, 10, factor);
            printf("\n");
            ecm_clear(params);
            return 1;
        }
    }

    ecm_clear(params);
    return 0;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    mpz_t n, factor;
    mpz_init(n);
    mpz_init(factor);

    mpz_set_str(n, argv[1], 10);

    int digits = mpz_sizeinbase(n, 10);
    printf("Factoring %d-digit number\n", digits);

    struct timespec ts_start, ts_now;
    clock_gettime(CLOCK_MONOTONIC, &ts_start);

    /* Phase 1: Quick ECM */
    printf("Phase 1: ECM (up to %ds)...\n", ECM_TIMEOUT);
    if (try_ecm(n, factor)) {
        mpz_t cofactor;
        mpz_init(cofactor);
        mpz_divexact(cofactor, n, factor);

        gmp_printf("P%d = %Zd\n", (int)mpz_sizeinbase(factor, 10), factor);
        gmp_printf("P%d = %Zd\n", (int)mpz_sizeinbase(cofactor, 10), cofactor);

        mpz_clear(cofactor);
        mpz_clear(n);
        mpz_clear(factor);
        return 0;
    }

    clock_gettime(CLOCK_MONOTONIC, &ts_now);
    double ecm_time = (ts_now.tv_sec - ts_start.tv_sec) +
                      (ts_now.tv_nsec - ts_start.tv_nsec) / 1e9;
    printf("ECM: no factor found in %.1fs\n", ecm_time);

    /* Phase 2: YAFU SIQS */
    int remaining = 295 - (int)(ecm_time + 0.5) - 2; /* 2s safety margin */
    if (remaining < 200) remaining = 200;

    printf("Phase 2: YAFU SIQS (timeout %ds)...\n", remaining);

    /* Determine optimal NB/B based on digit size */
    int nb = 20, b = 120000;
    if (digits <= 85) { nb = 18; b = 70000; }
    else if (digits <= 86) { nb = 18; b = 80000; }
    else if (digits <= 87) { nb = 18; b = 90000; }
    else if (digits <= 88) { nb = 18; b = 80000; }
    else if (digits <= 89) { nb = 18; b = 100000; }

    char cmd[4096];
    snprintf(cmd, sizeof(cmd),
        "WORKDIR=$(mktemp -d /tmp/yafu_XXXXXX) && cd $WORKDIR && "
        "echo \"siqs(%s)\" | LD_LIBRARY_PATH=/usr/local/lib "
        "timeout %d /tmp/agent-factoring-9/yafu/yafu "
        "-threads 1 -seed 42 -siqsNB %d -siqsB %d 2>&1 && "
        "rm -rf $WORKDIR",
        argv[1], remaining, nb, b);

    int ret = system(cmd);

    clock_gettime(CLOCK_MONOTONIC, &ts_now);
    double total_time = (ts_now.tv_sec - ts_start.tv_sec) +
                        (ts_now.tv_nsec - ts_start.tv_nsec) / 1e9;
    printf("Total time: %.1fs\n", total_time);

    mpz_clear(n);
    mpz_clear(factor);
    return WEXITSTATUS(ret);
}
