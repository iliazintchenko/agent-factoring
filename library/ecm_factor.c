/*
 * ecm_factor.c - ECM-based factoring wrapper using GMP-ECM library.
 * Used as baseline for scaling comparison.
 * Usage: ./ecm_factor <N>
 * Output: "p q" where N = p * q
 * Seed: deterministic starting from sigma=42.
 */
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <ecm.h>

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    mpz_t n, f;
    mpz_init(n);
    mpz_init(f);

    if (mpz_set_str(n, argv[1], 10) != 0) {
        fprintf(stderr, "Invalid number: %s\n", argv[1]);
        return 1;
    }

    /* Trial division up to 10^6 */
    for (unsigned long p = 2; p < 1000000; p++) {
        if (mpz_divisible_ui_p(n, p)) {
            mpz_t cofactor;
            mpz_init(cofactor);
            mpz_set_ui(f, p);
            mpz_divexact(cofactor, n, f);
            if (mpz_cmp_ui(cofactor, 1) > 0) {
                gmp_printf("%Zd %Zd\n", f, cofactor);
                mpz_clear(cofactor);
                mpz_clear(n); mpz_clear(f);
                return 0;
            }
            mpz_clear(cofactor);
        }
    }

    /* ECM with increasing B1 bounds.
     * Must re-init ecm_params for each curve to reset internal state. */
    double b1_levels[] = {2e3, 1.1e4, 5e4, 2.5e5, 1e6, 3e6, 1.1e7, 4.3e7, 1.1e8, 0};
    int curves_per_level[] = {25, 90, 300, 700, 1800, 5100, 10600, 19300, 49000};

    unsigned long sigma = 42;

    for (int level = 0; b1_levels[level] > 0; level++) {
        double b1 = b1_levels[level];
        int ncurves = curves_per_level[level];

        for (int curve = 0; curve < ncurves; curve++) {
            sigma++;
            if (sigma < 6) sigma = 6;

            ecm_params params;
            ecm_init(params);
            mpz_set_ui(params->sigma, sigma);
            params->sigma_is_A = 0;
            params->method = ECM_ECM;

            int ret = ecm_factor(f, n, b1, params);
            ecm_clear(params);

            if (ret > 0 && mpz_cmp_ui(f, 1) > 0 && mpz_cmp(f, n) != 0) {
                mpz_t cofactor;
                mpz_init(cofactor);
                mpz_divexact(cofactor, n, f);
                if (mpz_cmp(f, cofactor) > 0)
                    mpz_swap(f, cofactor);
                gmp_printf("%Zd %Zd\n", f, cofactor);
                mpz_clear(cofactor);
                mpz_clear(n); mpz_clear(f);
                return 0;
            }
        }
    }

    fprintf(stderr, "ECM failed to factor\n");
    mpz_clear(n); mpz_clear(f);
    return 1;
}
