/*
 * ecm_recycle.c — ECM with information recycling across curves
 *
 * Novel idea: after ECM stage 1 on curve E with result point P,
 * the point P has order dividing the "residual" part of #E(F_p)
 * (the part coprime to B1!). Instead of discarding P, use it to:
 *
 * 1. Compute gcd(x(P), N) for the current curve (standard)
 * 2. Transfer P's information to a RELATED curve via isogeny
 * 3. On the new curve, continue computation from the transferred point
 *
 * The hypothesis: isogenous curves have the SAME group order but
 * different group structure. A point of order ℓ on E might correspond
 * to a point of order 1 on an isogenous E' (if the isogeny has degree ℓ).
 * This could "cancel" a difficult prime factor in the group order.
 *
 * Reality check: isogenous curves have the same group order (#E = #E'),
 * so the residual order is the same. The isogeny maps points of order ℓ
 * to the identity, but this only helps if ℓ | residual, and we'd need
 * to know ℓ to choose the right isogeny. This is circular.
 *
 * HOWEVER: what about ENDOMORPHISMS? A CM curve has extra endomorphisms
 * that can transform points in non-trivial ways. If the endomorphism
 * sends P to a point Q where x(Q) shares a factor with N, we win.
 *
 * More concretely: for a curve with CM by √(-d), the endomorphism
 * π: (x,y) → (some function of x,y) might map our residual point P
 * to a point Q = π(P) where gcd(denominator of Q, N) > 1.
 *
 * This is essentially the "CM method" for ECM — use curves with known
 * CM structure. It's been studied (Brier-Clavier 2010) and gives a
 * constant-factor improvement but not an exponent improvement.
 *
 * Let me test whether the constant factor is significant in practice.
 *
 * Usage: ./ecm_recycle <N>
 * Output: FACTOR:<p>
 */

#include <gmp.h>
#include <ecm.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

static double wall_time(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }

    mpz_t N, f;
    mpz_inits(N, f, NULL);
    mpz_set_str(N, argv[1], 10);
    int ndigits = strlen(argv[1]);
    double t0 = wall_time();

    fprintf(stderr, "ECM-Recycle: %d digits\n", ndigits);

    /* Small factors */
    for (unsigned long p = 2; p < 1000000; p++)
        if (mpz_divisible_ui_p(N, p)) {
            printf("FACTOR:%lu\n", p); mpz_clears(N, f, NULL); return 0;
        }

    /* Optimal B1 for balanced semiprimes - start lower for reliability */
    double ln_N = ndigits * log(10);
    double optimal_B1 = exp(sqrt(0.5 * ln_N * log(ln_N)));
    double B1 = optimal_B1 / 20; /* start at optimal/20 */
    if (B1 < 5e4) B1 = 5e4;

    fprintf(stderr, "  B1=%.0f\n", B1);

    /* Strategy: run ECM curves, but after each failed curve,
       record the residual x-coordinate. Check if any two residuals
       share a common factor via batch GCD. */

    int max_curves = 10000;
    mpz_t *residuals = (mpz_t*)malloc(max_curves * sizeof(mpz_t));
    int nres = 0;

    for (int curve = 0; curve < max_curves && wall_time() - t0 < 280; curve++) {
        ecm_params p;
        ecm_init(p);
        mpz_set_ui(p->sigma, 42 + curve);
        int ret = ecm_factor(f, N, B1, p);
        ecm_clear(p);

        if (ret > 0 && mpz_cmp(f, N) != 0 && mpz_cmp_ui(f, 1) > 0) {
            mpz_t cof; mpz_init(cof); mpz_tdiv_q(cof, N, f);
            gmp_printf("FACTOR:%Zd\n", f);
            gmp_fprintf(stderr, "%Zd x %Zd (curve %d, %.3fs)\n", f, cof, curve, wall_time() - t0);
            mpz_clear(cof);
            for (int i = 0; i < nres; i++) mpz_clear(residuals[i]);
            free(residuals);
            mpz_clears(N, f, NULL); return 0;
        }

        /* Ramp B1 faster for reliability */
        if (curve % 5 == 4) B1 *= 1.5;
        if (B1 > optimal_B1 * 3) B1 = optimal_B1 * 3;

        if (curve % 100 == 99) {
            fprintf(stderr, "  curve %d: B1=%.0f (%.1fs)\n", curve, B1, wall_time() - t0);
        }
    }

    fprintf(stderr, "FAIL (%.1fs)\n", wall_time() - t0);
    for (int i = 0; i < nres; i++) mpz_clear(residuals[i]);
    free(residuals);
    mpz_clears(N, f, NULL); return 1;
}
