/*
 * divpoly_factor.c — Factor via division polynomial GCD.
 *
 * KEY INSIGHT from the Schoof analysis: the monic modulus prevented
 * factor-finding during polynomial reduction. Division polynomials
 * have NON-UNIT leading coefficients, creating GCD opportunities.
 *
 * For E: y^2 = x^3 + ax + b, the ℓ-th division polynomial ψ_ℓ(x)
 * has degree (ℓ^2-1)/2 and leading coefficient ℓ (for odd ℓ).
 * When computing x^N mod ψ_ℓ(x) mod N, the reduction step requires
 * dividing by ℓ. If gcd(ℓ, N) > 1, we find a factor.
 *
 * More usefully: the RESULTANT or GCD of ψ_ℓ values at different
 * points, computed mod N, creates many opportunities for non-invertible
 * elements.
 *
 * APPROACH: For the ℓ-th division polynomial ψ_ℓ(x):
 * 1. Compute ψ_ℓ(r) for random x = r mod N
 * 2. gcd(ψ_ℓ(r), N) might be non-trivial
 * 3. If #E(F_p) = ℓ*k, then ψ_ℓ vanishes on the ℓ-torsion points of E(F_p)
 * 4. Testing many r values is like searching for ℓ-torsion points
 *
 * This IS essentially ECM — but the division polynomial approach
 * explicitly constructs the ℓ-torsion structure, potentially allowing
 * more targeted search.
 *
 * NOVEL VARIANT: Instead of computing ψ_ℓ for a fixed ℓ, compute
 * the PRODUCT ψ_2 * ψ_3 * ψ_5 * ... * ψ_L, which vanishes at
 * ALL torsion points of order ≤ L. This is analogous to ECM stage 1
 * but working with division polynomials directly.
 */

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* Compute ψ_n(x) for the elliptic curve y^2 = x^3 + ax + b
 * Division polynomials satisfy the recurrence:
 *   ψ_1 = 1
 *   ψ_2 = 2y = 2*(x^3+ax+b)^{1/2} [we use ψ_2^2 = 4*(x^3+ax+b)]
 *   ψ_3 = 3x^4 + 6ax^2 + 12bx - a^2
 *   ψ_4 = 4y*(2x^6 + 10ax^4 + 40bx^3 - 10a^2x^2 - 8abx - 2a^3 - 16b^2)
 *   ...
 * For our purposes, we evaluate ψ_n at a specific point (x0, y0) mod N.
 * We work with the "numerator" version that avoids y (works with ψ^2 for even n).
 */

/* Evaluate ψ_n(P) where P = (x0, y0) on E: y^2 = x^3+ax+b mod N.
 * Returns ψ_n(P) mod N in result. Returns 1 if factor found. */
int eval_divpoly(mpz_t result, int n, const mpz_t x0, const mpz_t y0,
                 const mpz_t a, const mpz_t b, const mpz_t N, mpz_t factor) {
    /* Use the recurrence for ψ_n evaluated at (x0, y0):
     * ψ_1 = 1
     * ψ_2 = 2*y0
     * ψ_3 = 3*x0^4 + 6*a*x0^2 + 12*b*x0 - a^2
     * ψ_{2m+1} = ψ_{m+2}*ψ_m^3 - ψ_{m-1}*ψ_{m+1}^3
     * ψ_{2m} = (ψ_{m+2}*ψ_{m-1}^2 - ψ_{m-2}*ψ_{m+1}^2) * ψ_m / (2*y0)
     * (Division by 2y0 can reveal factors!)
     */

    if (n <= 0) { mpz_set_ui(result, 0); return 0; }

    /* Compute ψ values bottom-up */
    int sz = n + 3;
    mpz_t *psi = malloc(sz * sizeof(mpz_t));
    for (int i = 0; i < sz; i++) mpz_init(psi[i]);

    mpz_t x2, x3, x4, tmp1, tmp2, two_y;
    mpz_init(x2); mpz_init(x3); mpz_init(x4);
    mpz_init(tmp1); mpz_init(tmp2); mpz_init(two_y);

    mpz_mul(x2, x0, x0); mpz_mod(x2, x2, N);
    mpz_mul(x3, x2, x0); mpz_mod(x3, x3, N);
    mpz_mul(x4, x3, x0); mpz_mod(x4, x4, N);
    mpz_mul_ui(two_y, y0, 2); mpz_mod(two_y, two_y, N);

    /* ψ_0 = 0 */
    mpz_set_ui(psi[0], 0);
    /* ψ_1 = 1 */
    mpz_set_ui(psi[1], 1);
    /* ψ_2 = 2*y0 */
    mpz_set(psi[2], two_y);
    /* ψ_3 = 3x0^4 + 6ax0^2 + 12bx0 - a^2 */
    mpz_mul_ui(psi[3], x4, 3);
    mpz_mul(tmp1, a, x2); mpz_mul_ui(tmp1, tmp1, 6); mpz_add(psi[3], psi[3], tmp1);
    mpz_mul(tmp1, b, x0); mpz_mul_ui(tmp1, tmp1, 12); mpz_add(psi[3], psi[3], tmp1);
    mpz_mul(tmp1, a, a); mpz_sub(psi[3], psi[3], tmp1);
    mpz_mod(psi[3], psi[3], N);
    /* ψ_4 = 2y0*(2x0^6 + 10ax0^4 + 40bx0^3 - 10a^2*x0^2 - 8abx0 - 2a^3 - 16b^2) */
    {
        mpz_t x6; mpz_init(x6);
        mpz_mul(x6, x3, x3); mpz_mod(x6, x6, N);

        mpz_mul_ui(psi[4], x6, 2);
        mpz_mul(tmp1, a, x4); mpz_mul_ui(tmp1, tmp1, 10); mpz_add(psi[4], psi[4], tmp1);
        mpz_mul(tmp1, b, x3); mpz_mul_ui(tmp1, tmp1, 40); mpz_add(psi[4], psi[4], tmp1);
        mpz_mul(tmp1, a, a); mpz_mul(tmp1, tmp1, x2); mpz_mul_ui(tmp1, tmp1, 10);
        mpz_sub(psi[4], psi[4], tmp1);
        mpz_mul(tmp1, a, b); mpz_mul(tmp1, tmp1, x0); mpz_mul_ui(tmp1, tmp1, 8);
        mpz_sub(psi[4], psi[4], tmp1);
        mpz_mul(tmp1, a, a); mpz_mul(tmp1, tmp1, a); mpz_mul_ui(tmp1, tmp1, 2);
        mpz_sub(psi[4], psi[4], tmp1);
        mpz_mul(tmp1, b, b); mpz_mul_ui(tmp1, tmp1, 16);
        mpz_sub(psi[4], psi[4], tmp1);
        mpz_mul(psi[4], psi[4], two_y);
        mpz_mod(psi[4], psi[4], N);
        mpz_clear(x6);
    }

    /* Compute ψ_5, ..., ψ_n using the recurrence */
    for (int m = 3; m <= n; m++) {
        if (m + 1 >= sz || m - 1 < 0) break;
        /* Already have ψ_0 through ψ_{m-1} and ψ_m if m ≤ 4 */
        if (m >= 5) {
            int half = m / 2;
            if (m % 2 == 1) {
                /* ψ_m = ψ_{half+2} * ψ_{half}^3 - ψ_{half-1} * ψ_{half+1}^3
                 * where half = (m-1)/2 */
                half = (m - 1) / 2;
                if (half + 2 >= sz || half - 1 < 0) break;

                mpz_mul(tmp1, psi[half], psi[half]); mpz_mod(tmp1, tmp1, N);
                mpz_mul(tmp1, tmp1, psi[half]); mpz_mod(tmp1, tmp1, N);
                mpz_mul(tmp1, tmp1, psi[half + 2]); mpz_mod(tmp1, tmp1, N);

                mpz_mul(tmp2, psi[half + 1], psi[half + 1]); mpz_mod(tmp2, tmp2, N);
                mpz_mul(tmp2, tmp2, psi[half + 1]); mpz_mod(tmp2, tmp2, N);
                mpz_mul(tmp2, tmp2, psi[half - 1]); mpz_mod(tmp2, tmp2, N);

                mpz_sub(psi[m], tmp1, tmp2);
                mpz_mod(psi[m], psi[m], N);
            } else {
                /* ψ_m = ψ_{half} * (ψ_{half+2}*ψ_{half-1}^2 - ψ_{half-2}*ψ_{half+1}^2) / (2y0)
                 * The division by 2y0 can reveal a factor! */
                half = m / 2;
                if (half + 2 >= sz || half - 2 < 0) break;

                mpz_mul(tmp1, psi[half - 1], psi[half - 1]); mpz_mod(tmp1, tmp1, N);
                mpz_mul(tmp1, tmp1, psi[half + 2]); mpz_mod(tmp1, tmp1, N);

                mpz_mul(tmp2, psi[half + 1], psi[half + 1]); mpz_mod(tmp2, tmp2, N);
                mpz_mul(tmp2, tmp2, psi[half - 2]); mpz_mod(tmp2, tmp2, N);

                mpz_sub(tmp1, tmp1, tmp2);
                mpz_mul(tmp1, tmp1, psi[half]);
                mpz_mod(tmp1, tmp1, N);

                /* Divide by 2*y0 */
                if (!mpz_invert(tmp2, two_y, N)) {
                    /* 2*y0 not invertible — factor found! */
                    mpz_gcd(factor, two_y, N);
                    if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0) {
                        for (int i = 0; i < sz; i++) mpz_clear(psi[i]);
                        free(psi);
                        mpz_clear(x2);mpz_clear(x3);mpz_clear(x4);
                        mpz_clear(tmp1);mpz_clear(tmp2);mpz_clear(two_y);
                        return 1;
                    }
                }
                mpz_mul(psi[m], tmp1, tmp2);
                mpz_mod(psi[m], psi[m], N);
            }
        }

        /* Check if ψ_m shares a factor with N */
        mpz_gcd(tmp1, psi[m], N);
        if (mpz_cmp_ui(tmp1, 1) > 0 && mpz_cmp(tmp1, N) < 0) {
            mpz_set(factor, tmp1);
            for (int i = 0; i < sz; i++) mpz_clear(psi[i]);
            free(psi);
            mpz_clear(x2);mpz_clear(x3);mpz_clear(x4);
            mpz_clear(tmp1);mpz_clear(tmp2);mpz_clear(two_y);
            return 1;
        }
    }

    mpz_set(result, psi[n]);

    for (int i = 0; i < sz; i++) mpz_clear(psi[i]);
    free(psi);
    mpz_clear(x2);mpz_clear(x3);mpz_clear(x4);
    mpz_clear(tmp1);mpz_clear(tmp2);mpz_clear(two_y);
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N> [max_ℓ]\n", argv[0]); return 1; }

    mpz_t N; mpz_init(N); mpz_set_str(N, argv[1], 10);
    int max_l = (argc >= 3) ? atoi(argv[2]) : 100;

    size_t ndig = mpz_sizeinbase(N, 10);
    fprintf(stderr, "DivPoly-factor: N = %zu digits, max ℓ = %d\n", ndig, max_l);

    struct timespec t0; clock_gettime(CLOCK_MONOTONIC, &t0);

    gmp_randstate_t rng;
    gmp_randinit_mt(rng);
    gmp_randseed_ui(rng, 42);

    int found = 0;
    mpz_t factor, a, b, x0, y0, y0_sq, result;
    mpz_init(factor); mpz_init(a); mpz_init(b);
    mpz_init(x0); mpz_init(y0); mpz_init(y0_sq); mpz_init(result);

    for (int curve = 0; curve < 10000 && !found; curve++) {
        /* Random curve and point */
        mpz_urandomm(a, rng, N);
        mpz_urandomm(x0, rng, N);

        /* y0^2 = x0^3 + a*x0 + b for random b, or compute b = y0^2 - x0^3 - a*x0 */
        mpz_urandomm(y0, rng, N);
        mpz_mul(y0_sq, y0, y0); mpz_mod(y0_sq, y0_sq, N);
        /* b = y0^2 - x0^3 - a*x0 */
        mpz_t x0_3; mpz_init(x0_3);
        mpz_mul(x0_3, x0, x0); mpz_mul(x0_3, x0_3, x0); mpz_mod(x0_3, x0_3, N);
        mpz_t ax0; mpz_init(ax0);
        mpz_mul(ax0, a, x0); mpz_mod(ax0, ax0, N);
        mpz_sub(b, y0_sq, x0_3); mpz_sub(b, b, ax0); mpz_mod(b, b, N);
        mpz_clear(x0_3); mpz_clear(ax0);

        /* For each prime ℓ up to max_l, evaluate ψ_ℓ(x0, y0) */
        for (int l = 2; l <= max_l && !found; l++) {
            /* Quick primality check */
            int is_prime = 1;
            for (int d = 2; d * d <= l; d++) if (l % d == 0) { is_prime = 0; break; }
            if (!is_prime && l != 4) continue; /* Also try ℓ=4 for 2-power torsion */

            int ret = eval_divpoly(result, l, x0, y0, a, b, N, factor);
            if (ret) {
                found = 1;
                fprintf(stderr, "Factor from ψ_%d inversion! curve=%d\n", l, curve);
                break;
            }

            /* Also check gcd(result, N) */
            mpz_t g; mpz_init(g);
            mpz_gcd(g, result, N);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                mpz_set(factor, g);
                found = 1;
                fprintf(stderr, "Factor from ψ_%d value! curve=%d\n", l, curve);
            }
            mpz_clear(g);
        }

        struct timespec tn; clock_gettime(CLOCK_MONOTONIC, &tn);
        double el = (tn.tv_sec-t0.tv_sec)+(tn.tv_nsec-t0.tv_nsec)/1e9;
        if (el > 280) { fprintf(stderr, "TIMEOUT at curve %d\n", curve); break; }
        if ((curve + 1) % 500 == 0) fprintf(stderr, "  %d curves (%.1fs)\n", curve + 1, el);
    }

    struct timespec t1; clock_gettime(CLOCK_MONOTONIC, &t1);
    double total = (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)/1e9;

    if (found) {
        mpz_t q; mpz_init(q); mpz_divexact(q, N, factor);
        gmp_printf("%Zd %Zd\n", factor, q);
        fprintf(stderr, "SUCCESS in %.3fs\n", total);
        mpz_clear(q);
    } else {
        fprintf(stderr, "FAIL after %.3fs\n", total);
    }

    mpz_clear(factor); mpz_clear(a); mpz_clear(b);
    mpz_clear(x0); mpz_clear(y0); mpz_clear(y0_sq); mpz_clear(result);
    gmp_randclear(rng); mpz_clear(N);
    return found ? 0 : 1;
}
