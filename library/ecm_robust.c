/*
 * ecm_robust.c — Robust ECM factoring with escalating bounds
 *
 * Automatically selects B1 based on N size and escalates if needed.
 * Uses GMP-ECM library with Suyama parametrization.
 *
 * For balanced semiprimes N=pq with p~q~√N:
 * - Factor size: ~n/2 digits
 * - Optimal B1: exp(sqrt(2 * ln(p) * ln(ln(p)))) ≈ L_p[1/2, √2]
 * - Expected curves: L_p[1/2, 1/√2]
 *
 * This implementation tries multiple B1 levels, spending more curves
 * at each level before escalating.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <ecm.h>

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    mpz_t N, f;
    mpz_inits(N, f, NULL);
    mpz_set_str(N, argv[1], 10);

    int nbits = mpz_sizeinbase(N, 2);
    int ndig = mpz_sizeinbase(N, 10);

    /* Check if N is trivially factorable */
    if (mpz_probab_prime_p(N, 25)) {
        gmp_printf("%Zd 1\n", N);
        mpz_clears(N, f, NULL);
        return 0;
    }

    /* Small trial division */
    mpz_t small_f;
    mpz_init(small_f);
    for (unsigned long p = 2; p < 100000; p++) {
        if (mpz_divisible_ui_p(N, p)) {
            mpz_set_ui(f, p);
            mpz_t cof;
            mpz_init(cof);
            mpz_divexact(cof, N, f);
            gmp_printf("%Zd %Zd\n", f, cof);
            mpz_clears(N, f, small_f, cof, NULL);
            return 0;
        }
    }
    mpz_clear(small_f);

    /* Estimate optimal B1 for p ~ N^{1/2} */
    double ln_p = nbits * log(2.0) / 2.0;
    double ln_ln_p = log(ln_p);

    /* ECM B1 levels: start lower, escalate */
    double base_B1 = exp(0.6 * sqrt(2 * ln_p * ln_ln_p));
    if (base_B1 < 2000) base_B1 = 2000;

    /* B1 escalation schedule */
    double B1_levels[] = {base_B1 * 0.3, base_B1, base_B1 * 3, base_B1 * 10,
                          base_B1 * 30, base_B1 * 100};
    int curves_per_level[] = {100, 500, 1000, 2000, 5000, 10000};
    int n_levels = 6;

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    int total_curves = 0;
    int found = 0;

    for (int level = 0; level < n_levels && !found; level++) {
        double B1 = B1_levels[level];
        int n_curves = curves_per_level[level];

        for (int c = 0; c < n_curves && !found; c++) {
            ecm_params params;
            ecm_init(params);
            params->B1done = 1.0;
            params->param = ECM_PARAM_SUYAMA;
            mpz_set_ui(params->sigma, 6 + total_curves);

            int ret = ecm_factor(f, N, B1, params);
            ecm_clear(params);
            total_curves++;

            if (ret > 0 && mpz_cmp_ui(f, 1) > 0 && mpz_cmp(f, N) < 0) {
                found = 1;
                break;
            }
        }

        if (!found) {
            clock_gettime(CLOCK_MONOTONIC, &t1);
            double elapsed = (t1.tv_sec - t0.tv_sec) +
                           (t1.tv_nsec - t0.tv_nsec) / 1e9;
            fprintf(stderr, "ECM: Level %d (B1=%.0f, %d curves) done. "
                    "Total: %d curves, %.1fs\n",
                    level, B1, curves_per_level[level], total_curves, elapsed);
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;

    if (found) {
        mpz_t cof;
        mpz_init(cof);
        mpz_divexact(cof, N, f);
        if (mpz_cmp(f, cof) > 0) mpz_swap(f, cof);
        gmp_printf("%Zd %Zd\n", f, cof);
        fprintf(stderr, "ECM: Factor found after %d curves in %.3fs\n",
                total_curves, elapsed);
        mpz_clear(cof);
    } else {
        fprintf(stderr, "ECM: FAIL after %d curves (%.1fs)\n",
                total_curves, elapsed);
    }

    mpz_clears(N, f, NULL);
    return found ? 0 : 1;
}
