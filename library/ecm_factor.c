/*
 * ecm_factor.c — Minimal ECM-based factoring utility
 * Usage: ./ecm_factor <N>
 * Output: FACTOR: <p> (smallest prime factor)
 * Uses GMP-ECM library with escalating bounds.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <ecm.h>
#include <time.h>

/* Try trial division up to limit */
static int trial_divide(mpz_t n, mpz_t factor, unsigned long limit) {
    if (mpz_even_p(n)) {
        mpz_set_ui(factor, 2);
        return 1;
    }
    for (unsigned long p = 3; p <= limit; p += 2) {
        if (mpz_divisible_ui_p(n, p)) {
            mpz_set_ui(factor, p);
            return 1;
        }
    }
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

    if (mpz_cmp_ui(n, 1) <= 0 || mpz_probab_prime_p(n, 25)) {
        fprintf(stderr, "N is trivial or prime\n");
        mpz_clear(n); mpz_clear(factor);
        return 1;
    }

    /* Phase 1: Trial division up to 10^6 */
    if (trial_divide(n, factor, 1000000)) {
        gmp_printf("FACTOR: %Zd\n", factor);
        mpz_clear(n); mpz_clear(factor);
        return 0;
    }

    /* Phase 2: ECM with escalating bounds */
    ecm_params params;
    ecm_init(params);
    /* Use parameterization that allows arbitrary sigma, seed deterministically */
    params->param = ECM_PARAM_SUYAMA;

    /* ECM bound schedule: (B1, curves) pairs */
    struct { double b1; int curves; } schedule[] = {
        {     2000,    20 },
        {    11000,    50 },
        {    50000,   100 },
        {   250000,   200 },
        {  1000000,   400 },
        {  3000000,   700 },
        { 11000000,  1200 },
        { 44000000,  2400 },
        {110000000,  4800 },
        {260000000,  8000 },
    };
    int nlevels = sizeof(schedule) / sizeof(schedule[0]);

    for (int level = 0; level < nlevels; level++) {
        double b1 = schedule[level].b1;
        int ncurves = schedule[level].curves;

        for (int curve = 0; curve < ncurves; curve++) {
            mpz_set_ui(factor, 0);
            params->B1done = 1.0;
            /* Deterministic sigma from seed 42 + curve index */
            mpz_set_ui(params->sigma, 42 + (unsigned long)level * 10000 + (unsigned long)curve);
            int ret = ecm_factor(factor, n, b1, params);

            if (ret > 0 && mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, n) < 0) {
                /* Found a non-trivial factor. Return the smaller one. */
                mpz_t cofactor;
                mpz_init(cofactor);
                mpz_divexact(cofactor, n, factor);
                if (mpz_cmp(factor, cofactor) > 0)
                    mpz_set(factor, cofactor);
                gmp_printf("FACTOR: %Zd\n", factor);
                mpz_clear(cofactor);
                ecm_clear(params);
                mpz_clear(n); mpz_clear(factor);
                return 0;
            }
        }
    }

    fprintf(stderr, "ECM failed to factor after all levels\n");
    ecm_clear(params);
    mpz_clear(n); mpz_clear(factor);
    return 1;
}
