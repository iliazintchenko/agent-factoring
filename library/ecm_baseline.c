/*
 * ecm_baseline.c — ECM baseline factoring with timing
 *
 * Uses GMP-ECM to factor semiprimes. This establishes baseline timing
 * for comparison with novel approaches.
 *
 * ECM has complexity L_p[1/2, sqrt(2)] where p is the smallest factor.
 * For balanced semiprimes (p ~ q ~ N^{1/2}), this is L[1/4, sqrt(2)] in N.
 * Wait — actually for balanced N=pq, ECM is L_p[1/2] = L[1/4] in terms of N!
 * But the constant is bad. NFS is L[1/3] which beats ECM for large N.
 *
 * OUTPUT FORMAT: prints "factor1 factor2" on success, or "FAIL" on failure.
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
        fprintf(stderr, "Usage: %s <N_decimal> [B1] [max_curves]\n", argv[0]);
        return 1;
    }

    mpz_t N, f;
    mpz_inits(N, f, NULL);
    mpz_set_str(N, argv[1], 10);

    int n_bits = mpz_sizeinbase(N, 2);
    int n_digits = mpz_sizeinbase(N, 10);

    /* Default B1: heuristic for balanced semiprimes */
    /* For ECM finding p ~ N^{1/2}, optimal B1 ~ L_p[1/2, 1/sqrt(2)] */
    double ln_p = n_bits * log(2.0) / 2.0; /* log of expected factor size */
    double ln_ln_p = log(ln_p);
    double default_B1 = exp(sqrt(2.0 * ln_p * ln_ln_p) / 2.0);
    if (default_B1 < 1000) default_B1 = 1000;
    if (default_B1 > 1e12) default_B1 = 1e12;

    double B1 = argc > 2 ? atof(argv[2]) : default_B1;
    int max_curves = argc > 3 ? atoi(argv[3]) : 10000;

    fprintf(stderr, "ECM: N=%d digits (%d bits), B1=%.0f, max_curves=%d\n",
            n_digits, n_bits, B1, max_curves);

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    int found = 0;
    for (int curve = 0; curve < max_curves; curve++) {
        ecm_params params;
        ecm_init(params);
        params->B1done = 1.0;
        params->param = ECM_PARAM_SUYAMA;
        /* sigma must be >= 6 for Suyama parametrization */
        mpz_set_ui(params->sigma, 6 + curve);

        int ret = ecm_factor(f, N, B1, params);

        if (ret > 0) {
            /* Verify it's a non-trivial factor */
            if (mpz_cmp_ui(f, 1) > 0 && mpz_cmp(f, N) < 0) {
                found = 1;

                clock_gettime(CLOCK_MONOTONIC, &t1);
                double elapsed = (t1.tv_sec - t0.tv_sec) +
                               (t1.tv_nsec - t0.tv_nsec) / 1e9;

                mpz_t cofactor;
                mpz_init(cofactor);
                mpz_divexact(cofactor, N, f);

                /* Ensure f < cofactor for consistent output */
                if (mpz_cmp(f, cofactor) > 0)
                    mpz_swap(f, cofactor);

                gmp_printf("%Zd %Zd\n", f, cofactor);
                fprintf(stderr, "ECM: Found factor on curve %d in %.3fs (B1=%.0f)\n",
                        curve + 1, elapsed, B1);

                mpz_clear(cofactor);
                ecm_clear(params);
                break;
            }
        }

        ecm_clear(params);

        if ((curve + 1) % 500 == 0) {
            clock_gettime(CLOCK_MONOTONIC, &t1);
            double elapsed = (t1.tv_sec - t0.tv_sec) +
                           (t1.tv_nsec - t0.tv_nsec) / 1e9;
            fprintf(stderr, "ECM: %d curves done (%.1fs)\n", curve + 1, elapsed);
        }
    }

    if (!found) {
        clock_gettime(CLOCK_MONOTONIC, &t1);
        double elapsed = (t1.tv_sec - t0.tv_sec) +
                       (t1.tv_nsec - t0.tv_nsec) / 1e9;
        fprintf(stderr, "ECM: FAIL after %d curves (%.1fs, B1=%.0f)\n",
                max_curves, elapsed, B1);
        printf("FAIL\n");
    }

    mpz_clears(N, f, NULL);
    return found ? 0 : 1;
}
