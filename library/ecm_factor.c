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

    int digits = mpz_sizeinbase(n, 10);

    /* Phase 2: ECM with escalating bounds */
    ecm_params params;
    ecm_init(params);
    /* Use parameterization that allows arbitrary sigma, seed deterministically */
    params->param = ECM_PARAM_SUYAMA;

    /* ECM bound schedule: (B1, curves) pairs
     * Skip early levels for large N (balanced semiprime has factors ≈ √N) */
    struct { double b1; int curves; int min_digits; } schedule[] = {
        {     2000,    30,  0 },
        {    11000,    60,  0 },
        {    50000,   120, 20 },
        {   250000,   250, 28 },
        {  1000000,   500, 34 },
        {  3000000,   900, 38 },
        { 11000000,  1600, 44 },
        { 44000000,  3000, 50 },
        {110000000,  5000, 56 },
        {260000000, 10000, 60 },
        {850000000, 20000, 66 },
    };
    int nlevels = sizeof(schedule) / sizeof(schedule[0]);

    for (int level = 0; level < nlevels; level++) {
        if (digits < schedule[level].min_digits) continue; /* skip low B1 for large N */
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
