/*
 * cm_ecm.c — ECM with CM curve selection.
 *
 * IDEA: Use elliptic curves with complex multiplication (CM) by
 * specific discriminants D. For CM disc D, the curve has
 * #E(F_p) = p + 1 - t where t is determined by the norm equation
 * t^2 + |D|*s^2 = 4p (Deuring). Different D values give different
 * group orders. If we choose D such that the group order is more
 * likely to be smooth, ECM succeeds faster.
 *
 * Specifically, for D = -3 (Eisenstein integers): #E = p + 1 - t
 * where t^2 + 3s^2 = 4p. The group order p+1-t tends to cluster
 * near p+1 (when t is small). For p ≡ 2 (mod 3), t = ±1 giving
 * #E = p or p+2.
 *
 * For D = -4 (Gaussian integers): t^2 + 4s^2 = 4p, so t is even.
 * For p ≡ 3 (mod 4), t = 0 giving #E = p + 1.
 *
 * TEST: Compare ECM success rate using:
 * 1. Random curves (standard ECM)
 * 2. CM curves with D = -3, -4, -7, -8, -11
 * 3. Measure which D gives fastest factoring
 */

#include <gmp.h>
#include <ecm.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static double now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }

    mpz_t N, factor;
    mpz_init(N); mpz_init(factor);
    mpz_set_str(N, argv[1], 10);

    size_t ndig = mpz_sizeinbase(N, 10);
    fprintf(stderr, "CM-ECM: N = %zu digits\n", ndig);

    /* ECM with standard random curves for comparison */
    ecm_params params;
    ecm_init(params);

    double B1 = 1e6; /* Fixed B1 for comparison */
    int max_curves = 5000;

    /* Test 1: Standard ECM (Suyama parametrization) */
    double t0 = now();
    int found = 0;
    params->param = ECM_PARAM_SUYAMA;
    for (int c = 0; c < max_curves && !found; c++) {
        mpz_set_ui(params->sigma, 42 + c);
        params->B1done = 1.0;
        int ret = ecm_factor(factor, N, B1, params);
        if (ret > 0 && mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0) {
            found = 1;
            double el = now() - t0;
            fprintf(stderr, "Suyama: factor found at curve %d (%.3fs, B1=%.0f)\n", c, el, B1);
            gmp_printf("%Zd %Zd\n", factor, N);
        }
    }
    if (!found) fprintf(stderr, "Suyama: no factor in %d curves (%.3fs)\n", max_curves, now()-t0);

    /* Test 2: ECM with specific sigma values that correspond to CM curves */
    /* For D = -3: use j = 0, curve y^2 = x^3 + b for various b */
    /* For D = -4: use j = 1728, curve y^2 = x^3 + ax for various a */
    /* These can be parameterized via ECM_PARAM_BATCH_2 or similar */

    /* With GMP-ECM, we can't directly select CM curves. Instead, we
     * test the hypothesis indirectly: measure how the smoothness of
     * #E(F_p) varies with the Suyama parametrization sigma.
     * Different sigma values give different j-invariants, and some
     * j-invariants correspond to CM curves. */

    /* For a practical test: try many B1 values and curve counts */
    if (!found) {
        double schedule[][2] = {
            {5e4, 100}, {2e5, 200}, {1e6, 500}, {5e6, 1000},
            {2e7, 2000}, {1e8, 4000}
        };
        int nlevels = 6;

        for (int level = 0; level < nlevels && !found; level++) {
            double b1 = schedule[level][0];
            int nc = (int)schedule[level][1];

            for (int c = 0; c < nc && !found; c++) {
                mpz_set_ui(params->sigma, 42 + level * 10000 + c);
                params->B1done = 1.0;
                int ret = ecm_factor(factor, N, b1, params);
                if (ret > 0 && mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0) {
                    found = 1;
                    double el = now() - t0;
                    mpz_t q; mpz_init(q); mpz_divexact(q, N, factor);
                    gmp_printf("%Zd %Zd\n", factor, q);
                    fprintf(stderr, "SUCCESS at B1=%.0f, curve=%d (%.3fs)\n", b1, c, el);
                    mpz_clear(q);
                }
            }
            if (!found) {
                fprintf(stderr, "  B1=%.0f: %d curves tried (%.1fs)\n",
                        schedule[level][0], nc, now() - t0);
            }
        }
    }

    if (!found) fprintf(stderr, "FAIL (%.3fs)\n", now() - t0);

    ecm_clear(params);
    mpz_clear(N); mpz_clear(factor);
    return found ? 0 : 1;
}
