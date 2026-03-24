/*
 * surface_factor.c — Factoring via smooth values on algebraic surfaces
 *
 * NOVEL APPROACH: Instead of sieving a 1D polynomial Q(x) = (x+m)² - N,
 * we search for PAIRS (x, y) where x² + y² - N has a large smooth part.
 *
 * The value V(x,y) = x² + y² - N is a degree-2 form in TWO variables.
 * For fixed N, the set of (x,y) with V(x,y) ≡ 0 (mod p) is a curve
 * in the (x,y) plane — specifically, the circle x² + y² ≡ N (mod p).
 *
 * KEY OBSERVATION: For each prime p where N is a sum of two squares mod p,
 * the number of points on x² + y² ≡ N (mod p) is exactly p (for p > 2).
 * This means the "sieve density" is p/p² = 1/p per prime — SAME as 1D QS.
 *
 * BUT: the VALUES V(x,y) can be much SMALLER than QS values:
 * - QS: Q(x) = (x+m)² - N ≈ 2√N·x for small x. Values ~ √N × M.
 * - Surface: V(x,y) = x² + y² - N. For x,y ≈ √(N/2), V ≈ 0.
 *   Deviations from √(N/2) give V ≈ 2√(N/2) × |x - √(N/2)| ≈ √N × δ.
 *
 * So V(x,y) values are ~√N × δ where δ is the deviation. For δ ~ M,
 * values are the same size as QS. NO improvement.
 *
 * UNLESS we use the 2D structure to find (x,y) pairs where V is
 * ACCIDENTALLY small. On the variety V = 0, x² + y² = N exactly.
 * Near this variety, V is small. The "width" of the near-zero band
 * is δ ~ 1 (integer precision), giving V ~ √N. Same as QS minimum.
 *
 * ALTERNATIVE: use V(x,y) = x·y - N with x ≈ √N, y ≈ √N.
 * Then V = x·y - N ≈ √N(x - √N) + √N(y - √N) = √N(δx + δy).
 * Values are √N × (δx + δy). For |δx|, |δy| ≤ M: values ≈ √N × M.
 * Same again.
 *
 * The fundamental limit: for ANY degree-2 form F(x₁,...,x_k) with
 * F(x*) = N, the values near x* are F ≈ √N × ||x - x*||.
 * The smooth probability depends on |F|, which is always ≈ √N × M
 * for a sieve region of radius M. This gives L[1/2] regardless
 * of the number of variables or the specific form.
 *
 * CONCLUSION: Multi-variable quadratic forms don't help.
 * The only way to get smaller values is via HIGHER-DEGREE forms
 * in number fields (NFS).
 *
 * HOWEVER: this implementation tests whether the 2D geometry gives
 * any PRACTICAL advantage (better constants, more candidates per
 * unit of computation, etc.) even if the asymptotic scaling is the same.
 *
 * We use x² - N (standard QS polynomial) and compare with
 * (x+y)² - N = x² + 2xy + y² - N for various y values simultaneously.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

/* Test whether multi-polynomial evaluation + batch smooth detection
 * gives better smoothness yield than single-polynomial sieving */
void measure_multi_yield(mpz_t N, unsigned long B, int M, int n_polys) {
    int n_digits = mpz_sizeinbase(N, 10);
    int n_bits = mpz_sizeinbase(N, 2);

    mpz_t sqrtN, m, qx, abs_qx, primorial, g, rem;
    mpz_inits(sqrtN, m, qx, abs_qx, primorial, g, rem, NULL);
    mpz_sqrt(sqrtN, N);
    mpz_add_ui(m, sqrtN, 1);

    /* Build primorial */
    mpz_t bound;
    mpz_init(bound);
    mpz_mul_ui(bound, sqrtN, 2 * M);
    mpz_set_ui(primorial, 1);
    mpz_t pp, pk;
    mpz_inits(pp, pk, NULL);
    mpz_set_ui(pp, 2);
    while (mpz_cmp_ui(pp, B) <= 0) {
        mpz_set(pk, pp);
        while (1) {
            mpz_t t; mpz_init(t);
            mpz_mul(t, pk, pp);
            if (mpz_cmp(t, bound) > 0) { mpz_clear(t); break; }
            mpz_set(pk, t); mpz_clear(t);
        }
        mpz_mul(primorial, primorial, pk);
        mpz_nextprime(pp, pp);
    }
    mpz_clears(pp, pk, bound, NULL);

    /* Strategy A: Single polynomial, M candidates */
    int smooth_A = 0;
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    for (int x = 1; x <= M; x++) {
        mpz_add_ui(qx, m, x);
        mpz_mul(qx, qx, qx);
        mpz_sub(qx, qx, N);
        mpz_abs(abs_qx, qx);

        mpz_set(rem, abs_qx);
        while (1) {
            mpz_gcd(g, rem, primorial);
            if (mpz_cmp_ui(g, 1) == 0) break;
            mpz_divexact(rem, rem, g);
        }
        if (mpz_cmp_ui(rem, 1) == 0) smooth_A++;
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double time_A = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;

    /* Strategy B: n_polys polynomials, M/n_polys candidates each */
    /* Use Q_k(x) = (x + m + k*step)² - N for k = 0..n_polys-1 */
    /* Pick x with smallest |Q_k(x)| across all k polynomials */
    int smooth_B = 0;
    int M_per_poly = M / n_polys;
    int step = M_per_poly;

    clock_gettime(CLOCK_MONOTONIC, &t0);

    for (int x = 1; x <= M_per_poly; x++) {
        /* Find polynomial k that gives smallest value */
        mpz_t best_val, test_val, best_xm;
        mpz_inits(best_val, test_val, best_xm, NULL);
        int best_k = 0;

        for (int k = 0; k < n_polys; k++) {
            mpz_add_ui(qx, m, x + k * step);
            mpz_mul(qx, qx, qx);
            mpz_sub(qx, qx, N);
            mpz_abs(test_val, qx);

            if (k == 0 || mpz_cmp(test_val, best_val) < 0) {
                mpz_set(best_val, test_val);
                mpz_add_ui(best_xm, m, x + k * step);
                best_k = k;
            }
        }

        /* Test smoothness of best value */
        mpz_set(rem, best_val);
        while (1) {
            mpz_gcd(g, rem, primorial);
            if (mpz_cmp_ui(g, 1) == 0) break;
            mpz_divexact(rem, rem, g);
        }
        if (mpz_cmp_ui(rem, 1) == 0) smooth_B++;

        mpz_clears(best_val, test_val, best_xm, NULL);
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double time_B = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;

    printf("N: %d digits, B: %lu, M: %d, n_polys: %d\n", n_digits, B, M, n_polys);
    printf("Strategy A (single poly, %d candidates): %d smooth (%.4f%%) in %.2fs\n",
           M, smooth_A, 100.0 * smooth_A / M, time_A);
    printf("Strategy B (best-of-%d, %d candidates each): %d smooth (%.4f%%) in %.2fs\n",
           n_polys, M_per_poly, smooth_B, 100.0 * smooth_B / M_per_poly, time_B);
    printf("Smoothness rate ratio B/A: %.2f\n",
           (double)smooth_B / M_per_poly / ((double)smooth_A / M));
    printf("Work-adjusted ratio (accounts for n_polys evaluations per step): %.2f\n",
           (double)smooth_B / smooth_A * M / M_per_poly / n_polys);

    mpz_clears(sqrtN, m, qx, abs_qx, primorial, g, rem, NULL);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N> [B] [M] [n_polys]\n", argv[0]);
        return 1;
    }

    mpz_t N;
    mpz_init_set_str(N, argv[1], 10);

    int nbits = mpz_sizeinbase(N, 2);
    double lnN = nbits * log(2.0);
    double lnlnN = log(lnN);

    unsigned long B = argc > 2 ? atol(argv[2]) :
        (unsigned long)(exp(0.5 * sqrt(lnN * lnlnN)));
    int M = argc > 3 ? atoi(argv[3]) : 500000;
    int n_polys = argc > 4 ? atoi(argv[4]) : 10;

    if (B < 2000) B = 2000;

    measure_multi_yield(N, B, M, n_polys);

    mpz_clear(N);
    return 0;
}
