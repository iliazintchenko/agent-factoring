/*
 * ecm_ultra.c — Ultra-aggressive ECM for 60-70 digit semiprimes
 *
 * Uses GMP-ECM with:
 * 1. Large B1 values (up to 10^8) with Stage 2
 * 2. Stage 2 bounds automatically set by GMP-ECM
 * 3. Many curves per B1 level
 * 4. Deterministic sigma sequence for reproducibility (seed 42)
 *
 * For balanced 62-digit semiprimes (factor ~31 digits = ~103 bits):
 * Optimal B1 ≈ 10^6, expected curves ≈ 2000-5000
 * With Stage 2: effective bound ~100*B1, significant speedup
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

    /* Small trial division first */
    for (unsigned long p = 2; p < 1000000; p++) {
        if (mpz_divisible_ui_p(N, p)) {
            mpz_set_ui(f, p);
            mpz_t cof; mpz_init(cof);
            mpz_divexact(cof, N, f);
            gmp_printf("%Zd %Zd\n", f, cof);
            mpz_clears(N, f, cof, NULL);
            return 0;
        }
    }

    /* Optimal B1 for balanced semiprimes */
    double ln_p = nbits * log(2.0) / 2.0;
    double ln_ln_p = log(ln_p);
    double opt_B1 = exp(sqrt(2.0 * ln_p * ln_ln_p));

    /* B1 schedule: start at opt_B1/10, escalate */
    double B1_schedule[] = {
        opt_B1 * 0.1,
        opt_B1 * 0.3,
        opt_B1,
        opt_B1 * 3,
        opt_B1 * 10,
        opt_B1 * 30,
        opt_B1 * 100,
    };
    int curves_schedule[] = {50, 100, 300, 500, 1000, 2000, 5000};
    int n_levels = 7;

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    int total_curves = 0;
    int found = 0;

    for (int level = 0; level < n_levels && !found; level++) {
        double B1 = B1_schedule[level];
        if (B1 < 1000) B1 = 1000;
        if (B1 > 1e10) B1 = 1e10;

        int n_curves = curves_schedule[level];

        for (int c = 0; c < n_curves && !found; c++) {
            ecm_params params;
            ecm_init(params);
            params->B1done = 1.0;
            params->param = ECM_PARAM_SUYAMA;
            mpz_set_ui(params->sigma, 42 + total_curves);
            /* Let ECM choose stage 2 bound automatically */
            /* params->B2 is set automatically when < 0 */

            int ret = ecm_factor(f, N, B1, params);
            ecm_clear(params);
            total_curves++;

            if (ret > 0 && mpz_cmp_ui(f, 1) > 0 && mpz_cmp(f, N) < 0) {
                found = 1;
                break;
            }

            /* Time check: bail if approaching 290s */
            if (total_curves % 100 == 0) {
                clock_gettime(CLOCK_MONOTONIC, &t1);
                double elapsed = (t1.tv_sec - t0.tv_sec) +
                               (t1.tv_nsec - t0.tv_nsec) / 1e9;
                if (elapsed > 280) {
                    fprintf(stderr, "ECM: Time limit approaching, stopping\n");
                    goto done;
                }
            }
        }

        clock_gettime(CLOCK_MONOTONIC, &t1);
        double elapsed = (t1.tv_sec - t0.tv_sec) +
                       (t1.tv_nsec - t0.tv_nsec) / 1e9;
        if (!found)
            fprintf(stderr, "ECM: Level %d (B1=%.0f, %d curves) done. "
                    "Total: %d curves, %.1fs\n",
                    level, B1, curves_schedule[level], total_curves, elapsed);
    }

done:
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;

    if (found) {
        mpz_t cof; mpz_init(cof);
        mpz_divexact(cof, N, f);
        if (mpz_cmp(f, cof) > 0) mpz_swap(f, cof);
        gmp_printf("%Zd %Zd\n", f, cof);
        fprintf(stderr, "ECM: Factor found after %d curves in %.3fs\n",
                total_curves, elapsed);
        mpz_clear(cof);
    } else {
        fprintf(stderr, "ECM: FAIL after %d curves (%.1fs)\n", total_curves, elapsed);
    }

    mpz_clears(N, f, NULL);
    return found ? 0 : 1;
}
