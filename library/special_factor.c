/*
 * special_factor.c - Pollard p-1, Williams p+1, and ECM factoring
 *
 * Tries special-purpose methods that are O(p^epsilon) rather than
 * O(N^epsilon). For balanced semiprimes where p ≈ sqrt(N), these
 * methods are only useful if p-1 or p+1 has specific smoothness
 * properties. But for random (non-crypto) semiprimes, this is worth
 * checking.
 *
 * Compile: gcc -O3 -march=native -o special_factor library/special_factor.c -lgmp -lecm -lm
 * Usage: ./special_factor <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#include <ecm.h>

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    mpz_t N, factor;
    mpz_init(N);
    mpz_init(factor);
    mpz_set_str(N, argv[1], 10);

    int digits = mpz_sizeinbase(N, 10);

    /* ECM parameters based on expected factor size */
    /* For balanced semiprimes: factors are ~digits/2 digits */
    int factor_digits = digits / 2;

    /* P-1 method: works if p-1 is B-smooth */
    ecm_params params;
    ecm_init(params);

    /* Set seed */
    mpz_set_ui(params->sigma, 42);

    /* Try p-1 with escalating bounds */
    double B1_pm1[] = {1e4, 1e5, 1e6, 5e6, 1e7, 5e7, 0};
    for (int i = 0; B1_pm1[i] > 0; i++) {
        clock_gettime(CLOCK_MONOTONIC, &t1);
        double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
        if (elapsed > 280.0) break;

        params->method = ECM_PM1;
        int ret = ecm_factor(factor, N, B1_pm1[i], params);

        if (ret > 0 && mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0) {
            clock_gettime(CLOCK_MONOTONIC, &t1);
            elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;

            mpz_t cofactor;
            mpz_init(cofactor);
            mpz_divexact(cofactor, N, factor);
            gmp_printf("%Zd = %Zd * %Zd\n", N, factor, cofactor);
            fprintf(stderr, "P-1 with B1=%.0f found factor in %.3fs\n", B1_pm1[i], elapsed);
            mpz_clear(cofactor);
            ecm_clear(params);
            mpz_clear(N); mpz_clear(factor);
            return 0;
        }
    }

    /* Try p+1 with escalating bounds */
    double B1_pp1[] = {1e4, 1e5, 1e6, 5e6, 1e7, 5e7, 0};
    for (int i = 0; B1_pp1[i] > 0; i++) {
        clock_gettime(CLOCK_MONOTONIC, &t1);
        double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
        if (elapsed > 280.0) break;

        params->method = ECM_PP1;
        int ret = ecm_factor(factor, N, B1_pp1[i], params);

        if (ret > 0 && mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0) {
            clock_gettime(CLOCK_MONOTONIC, &t1);
            elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;

            mpz_t cofactor;
            mpz_init(cofactor);
            mpz_divexact(cofactor, N, factor);
            gmp_printf("%Zd = %Zd * %Zd\n", N, factor, cofactor);
            fprintf(stderr, "P+1 with B1=%.0f found factor in %.3fs\n", B1_pp1[i], elapsed);
            mpz_clear(cofactor);
            ecm_clear(params);
            mpz_clear(N); mpz_clear(factor);
            return 0;
        }
    }

    /* ECM with increasing curves */
    /* For balanced semiprimes, ECM is generally not competitive with QS,
     * but might get lucky on some numbers */
    int ecm_curves[] = {10, 20, 50, 100, 0};
    double B1_ecm[] = {1e4, 5e4, 1e5, 5e5, 0};

    for (int b = 0; B1_ecm[b] > 0; b++) {
        for (int c = 0; ecm_curves[c] > 0; c++) {
            clock_gettime(CLOCK_MONOTONIC, &t1);
            double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
            if (elapsed > 280.0) goto done;

            for (int curve = 0; curve < ecm_curves[c]; curve++) {
                params->method = ECM_ECM;
                mpz_set_ui(params->sigma, 42 + b * 1000 + c * 100 + curve);

                int ret = ecm_factor(factor, N, B1_ecm[b], params);
                if (ret > 0 && mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0) {
                    clock_gettime(CLOCK_MONOTONIC, &t1);
                    elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;

                    mpz_t cofactor;
                    mpz_init(cofactor);
                    mpz_divexact(cofactor, N, factor);
                    gmp_printf("%Zd = %Zd * %Zd\n", N, factor, cofactor);
                    fprintf(stderr, "ECM (B1=%.0f, curve %d) found factor in %.3fs\n",
                            B1_ecm[b], curve, elapsed);
                    mpz_clear(cofactor);
                    ecm_clear(params);
                    mpz_clear(N); mpz_clear(factor);
                    return 0;
                }
            }
        }
    }

done:
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
    fprintf(stderr, "No factor found after %.3fs\n", elapsed);
    gmp_printf("FAIL %Zd\n", N);

    ecm_clear(params);
    mpz_clear(N); mpz_clear(factor);
    return 0;
}
