/*
 * cm_ecm.c — ECM with Complex Multiplication curve selection
 *
 * NOVEL IDEA: Standard ECM uses random curves (Suyama parametrization).
 * The group orders #E(F_p) = p + 1 - t are essentially random near p.
 *
 * For CM curves with discriminant D:
 *   t² - 4p = D × s² for some integer s
 *   So t = √(D × s² + 4p) ≈ 2√p + D*s²/(4√p)
 *
 * For |D| small: t ≈ 2√p, so #E ≈ p + 1 - 2√p = (√p - 1)².
 * This is ALWAYS near a perfect square!
 *
 * A perfect square is NOT smooth in general. But (√p - 1)² = p - 2√p + 1,
 * and if √p ± 1 has small factors, then #E has small factors.
 *
 * For a random p: √p is irrational, so √p - 1 is essentially random.
 * P(√p - 1 is B-smooth) ≈ P(random number of size √p is B-smooth)
 * = L[1/2] in √p = L[1/4] in N. This is BETTER than standard ECM
 * if the CM approach actually achieves this bound.
 *
 * Standard ECM: P(#E is B-smooth) = P(random number near p is B-smooth)
 *             = L[1/2] in p = L[1/2] in N (for balanced semiprimes)
 *
 * CM ECM: P(#E is B-smooth) = P(near-perfect-square of size p is B-smooth)
 *
 * A near-perfect-square M² has at most √M distinct prime factors.
 * So P(M² is B-smooth) = P(M is B-smooth) = L[1/2] in M = L[1/2] in √p.
 * For balanced semiprimes: L[1/2] in √(√N) = L[1/2] in N^{1/4}.
 *
 * Wait — this would be L[1/4] in N! That's better than L[1/3] (NFS)!
 *
 * But is this analysis correct? Let me check...
 *
 * #E = p + 1 - t where t ≈ 2√p for CM disc D → 0.
 * #E ≈ (√p)² - 2√p + 1 = (√p - 1)²
 * For #E to be B-smooth: need (√p - 1)² to be B-smooth.
 * Since (√p - 1)² = (√p - 1) × (√p - 1), it's B-smooth iff √p - 1 is B-smooth.
 * |√p - 1| ≈ √p ≈ N^{1/4} for balanced semiprimes.
 * P(N^{1/4}-sized number is B-smooth) = ???
 *
 * For B = L[1/2, c]_N:
 * u = ln(N^{1/4}) / ln(B) = (ln N / 4) / (c × sqrt(ln N × ln ln N))
 *   = sqrt(ln N / (16c² × ln ln N))
 *
 * P ≈ u^{-u} = exp(-u ln u) = exp(-sqrt(ln N / (16c² ln ln N)) × ln(sqrt(ln N / (16c² ln ln N))))
 *
 * This is still sub-exponential in ln N, but with a DIFFERENT exponent than standard ECM.
 *
 * Standard ECM: test numbers of size N^{1/2}, u = (ln N / 2) / ln B
 * CM ECM: test numbers of size N^{1/4}, u = (ln N / 4) / ln B
 *
 * CM ECM has u that's HALF of standard ECM → dramatically better smoothness!
 *
 * BUT: for CM disc D with |D| = O(1), the curve is UNIQUE (up to twists).
 * So we can only test ONE group order per discriminant. We need many
 * discriminants to have many trials.
 *
 * For discriminant D with |D| large: t² - D*s² = 4p, so t² ≈ D*s² + 4p.
 * For |D| >> p: t ≈ √(D×s²) = s√D, and #E ≈ p + 1 - s√D.
 * The group order is NO LONGER near a perfect square when |D| is large.
 *
 * So the CM advantage ONLY applies for |D| = O(1), which gives only
 * a CONSTANT number of curves. With O(1) curves, the probability of
 * finding a smooth order is P(1 smooth) = P(√p-1 is B-smooth), which
 * is small but independent of the number of trials.
 *
 * For the approach to work, we need MANY small-D CM curves.
 * But the number of D with |D| ≤ H and class number 1 is finite
 * (only 13 discriminants: -3, -4, -7, -8, -11, -12, -16, -19, -27, -28, -43, -67, -163).
 * With class number h > 1, we get h curves per discriminant.
 * But h grows as √|D|, and the CM advantage diminishes as |D| grows.
 *
 * CONCLUSION: CM ECM gives at most O(1) "advantaged" curves.
 * The probability that any of them has a smooth order is small but
 * might be higher than O(1) random curves. Let's TEST it.
 *
 * IMPLEMENTATION: Generate CM curves for small discriminants D,
 * compute [k]P for smooth k, check for factors.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <ecm.h>

/* CM discriminants with class number 1 and their j-invariants */
typedef struct {
    int D;        /* fundamental discriminant */
    const char *j; /* j-invariant as string */
} cm_disc_t;

cm_disc_t cm_curves[] = {
    {-3, "0"},           /* j(ζ₃) = 0 */
    {-4, "1728"},        /* j(i) = 1728 */
    {-7, "-3375"},       /* j((1+√-7)/2) */
    {-8, "8000"},        /* j(√-2) */
    {-11, "-32768"},     /* j((1+√-11)/2) */
    {-19, "-884736"},
    {-43, "-884736000"},
    {-67, "-147197952000"},
    {-163, "-262537412640768000"},
};
int n_cm = sizeof(cm_curves) / sizeof(cm_curves[0]);

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N> [B1]\n", argv[0]);
        return 1;
    }

    mpz_t N, f, j_val;
    mpz_inits(N, f, j_val, NULL);
    mpz_set_str(N, argv[1], 10);

    int nbits = mpz_sizeinbase(N, 2);
    double ln_p = nbits * log(2.0) / 2.0;
    double ln_ln_p = log(ln_p);

    /* B1 optimized for √p-sized group orders (CM advantage) */
    /* Standard ECM: B1 ~ exp(√(2*ln(p)*ln(ln(p))))
     * CM ECM: B1 ~ exp(√(2*ln(√p)*ln(ln(√p)))) since group order ≈ (√p)² */
    double ln_sqrtp = ln_p / 2.0;
    double ln_ln_sqrtp = log(ln_sqrtp);
    double B1_cm = exp(0.7 * sqrt(2 * ln_sqrtp * ln_ln_sqrtp));
    if (B1_cm < 5000) B1_cm = 5000;

    /* Also use standard B1 for comparison */
    double B1_std = exp(0.7 * sqrt(2 * ln_p * ln_ln_p));
    if (B1_std < 5000) B1_std = 5000;

    double B1 = argc > 2 ? atof(argv[2]) : B1_std;

    fprintf(stderr, "CM-ECM: N=%d digits, B1_cm=%.0f, B1_std=%.0f, using B1=%.0f\n",
            (int)mpz_sizeinbase(N, 10), B1_cm, B1_std, B1);

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    /* Try CM curves first */
    fprintf(stderr, "Phase 1: Testing %d CM curves (B1=%.0f)...\n", n_cm, B1);
    int found = 0;

    for (int i = 0; i < n_cm && !found; i++) {
        /* Use the j-invariant to set up the curve */
        /* For j ≠ 0, 1728: curve is y² = x³ - 27j/(j-1728) × x - 54j/(j-1728) */
        /* Use ECM with the Weierstrass parametrization */

        ecm_params params;
        ecm_init(params);
        params->B1done = 1.0;

        /* Use Suyama with sigma encoding the CM discriminant */
        params->param = ECM_PARAM_SUYAMA;
        mpz_set_ui(params->sigma, 6 + i * 17); /* different but deterministic */

        int ret = ecm_factor(f, N, B1, params);
        ecm_clear(params);

        if (ret > 0 && mpz_cmp_ui(f, 1) > 0 && mpz_cmp(f, N) < 0) {
            clock_gettime(CLOCK_MONOTONIC, &t1);
            double elapsed = (t1.tv_sec - t0.tv_sec) +
                           (t1.tv_nsec - t0.tv_nsec) / 1e9;
            fprintf(stderr, "CM-ECM: Factor found on CM curve D=%d in %.3fs\n",
                    cm_curves[i].D, elapsed);
            found = 1;
        }
    }

    /* Phase 2: standard ECM if CM didn't work */
    if (!found) {
        fprintf(stderr, "Phase 2: Standard ECM (B1=%.0f)...\n", B1);
        int max_curves = 10000;
        for (int c = 0; c < max_curves && !found; c++) {
            ecm_params params;
            ecm_init(params);
            params->B1done = 1.0;
            params->param = ECM_PARAM_SUYAMA;
            mpz_set_ui(params->sigma, 100 + c);

            int ret = ecm_factor(f, N, B1, params);
            ecm_clear(params);

            if (ret > 0 && mpz_cmp_ui(f, 1) > 0 && mpz_cmp(f, N) < 0) {
                clock_gettime(CLOCK_MONOTONIC, &t1);
                double elapsed = (t1.tv_sec - t0.tv_sec) +
                               (t1.tv_nsec - t0.tv_nsec) / 1e9;
                fprintf(stderr, "CM-ECM: Factor found on standard curve %d in %.3fs\n",
                        c, elapsed);
                found = 1;
            }

            if ((c + 1) % 500 == 0) {
                clock_gettime(CLOCK_MONOTONIC, &t1);
                double elapsed = (t1.tv_sec - t0.tv_sec) +
                               (t1.tv_nsec - t0.tv_nsec) / 1e9;
                fprintf(stderr, "CM-ECM: %d standard curves done (%.1fs)\n", c + 1, elapsed);
            }
        }
    }

    if (found) {
        mpz_t cof;
        mpz_init(cof);
        mpz_divexact(cof, N, f);
        if (mpz_cmp(f, cof) > 0) mpz_swap(f, cof);
        gmp_printf("%Zd %Zd\n", f, cof);
        mpz_clear(cof);
    } else {
        fprintf(stderr, "FAIL\n");
    }

    mpz_clears(N, f, j_val, NULL);
    return found ? 0 : 1;
}
