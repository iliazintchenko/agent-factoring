/*
 * pric.c — Polynomial Ring Index Calculus / Berlekamp factoring
 *
 * Correct Berlekamp approach:
 * 1. Choose random degree-d polynomial f(x) over Z/NZ
 * 2. Compute g(x) = x^N mod f(x) mod N
 * 3. Compute gcd(g(x) - x, f(x)) over Z/NZ
 * 4. During the GCD computation, non-invertible leading coefficients reveal factors of N
 *
 * Over F_p: gcd(x^N - x, f) = product of linear factors of f
 * Over F_q: different set of linear factors
 * If f splits differently mod p vs q, the GCD over Z/NZ reveals the factorization
 *
 * Complexity: O(d^2 * log(N)) polynomial multiplications, each O(M(n)) for n-bit N
 * For d=2: O(n * M(n)) = O(n^2 log n) with FFT multiplication
 *
 * This is fundamentally O(N^epsilon) — exponential! — because it only works
 * when f has different splitting patterns mod p vs q, which requires
 * trying O(1) random polynomials on average.
 * The real question: can we use this structure to build a sub-exponential algorithm?
 *
 * Usage: ./pric <N>
 * Output: FACTOR:<p>
 */

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

static double wall_time(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* Polynomial over Z/NZ represented as (c0, c1) for degree ≤ 1 when working mod degree-2 f */
/* f(x) = x^2 + a*x + b (monic degree 2) */
/* Element in ring: h0 + h1*x */

/* Multiply two elements: (h0 + h1*x)(g0 + g1*x) = h0*g0 + (h0*g1 + h1*g0)*x + h1*g1*x^2 */
/* x^2 = -a*x - b, so h1*g1*x^2 = h1*g1*(-a*x - b) */
/* Result: (h0*g0 - h1*g1*b) + (h0*g1 + h1*g0 - h1*g1*a)*x */

static void ring_mul(mpz_t r0, mpz_t r1,
                     const mpz_t h0, const mpz_t h1,
                     const mpz_t g0, const mpz_t g1,
                     const mpz_t a, const mpz_t b, const mpz_t N) {
    mpz_t t0, t1, t2;
    mpz_inits(t0, t1, t2, NULL);

    /* t2 = h1*g1 */
    mpz_mul(t2, h1, g1); mpz_mod(t2, t2, N);
    /* r0 = h0*g0 - t2*b */
    mpz_mul(t0, h0, g0); mpz_mul(t1, t2, b); mpz_sub(t0, t0, t1); mpz_mod(r0, t0, N);
    /* r1 = h0*g1 + h1*g0 - t2*a */
    mpz_mul(t0, h0, g1); mpz_mul(t1, h1, g0); mpz_add(t0, t0, t1);
    mpz_mul(t1, t2, a); mpz_sub(t0, t0, t1); mpz_mod(r1, t0, N);

    mpz_clears(t0, t1, t2, NULL);
}

/* Compute x^exp mod (x^2 + a*x + b) mod N */
static void ring_powmod(mpz_t r0, mpz_t r1,
                        const mpz_t exp,
                        const mpz_t a, const mpz_t b, const mpz_t N) {
    /* result = 1 */
    mpz_set_ui(r0, 1); mpz_set_ui(r1, 0);
    /* base = x = (0, 1) */
    mpz_t ba0, ba1, t0, t1;
    mpz_inits(ba0, ba1, t0, t1, NULL);
    mpz_set_ui(ba0, 0); mpz_set_ui(ba1, 1);

    int nbits = mpz_sizeinbase(exp, 2);
    for (int i = nbits - 1; i >= 0; i--) {
        /* Square result */
        ring_mul(t0, t1, r0, r1, r0, r1, a, b, N);
        mpz_set(r0, t0); mpz_set(r1, t1);

        if (mpz_tstbit(exp, i)) {
            /* Multiply by x = (0, 1) */
            ring_mul(t0, t1, r0, r1, ba0, ba1, a, b, N);
            mpz_set(r0, t0); mpz_set(r1, t1);
        }
    }

    mpz_clears(ba0, ba1, t0, t1, NULL);
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }
    mpz_t N, factor;
    mpz_inits(N, factor, NULL);
    mpz_set_str(N, argv[1], 10);
    int ndigits = strlen(argv[1]);
    double t0 = wall_time();

    fprintf(stderr, "PRIC: factoring %d-digit number\n", ndigits);

    /* Small factors */
    for (unsigned long p = 2; p < 1000000; p++)
        if (mpz_divisible_ui_p(N, p)) {
            printf("FACTOR:%lu\n", p); mpz_clears(N, factor, NULL); return 0; }

    gmp_randstate_t rs;
    gmp_randinit_default(rs);
    gmp_randseed_ui(rs, 42);

    /* Try random degree-2 polynomials */
    for (int trial = 0; trial < 100000 && wall_time() - t0 < 280; trial++) {
        mpz_t a_coef, b_coef;
        mpz_inits(a_coef, b_coef, NULL);
        mpz_urandomm(a_coef, rs, N);
        mpz_urandomm(b_coef, rs, N);

        /* Compute h(x) = x^N mod (x^2 + a*x + b) mod N */
        mpz_t h0, h1;
        mpz_inits(h0, h1, NULL);
        ring_powmod(h0, h1, N, a_coef, b_coef, N);

        /* g(x) = h(x) - x = h0 + (h1 - 1)*x */
        mpz_sub_ui(h1, h1, 1); mpz_mod(h1, h1, N);

        /* Now gcd(g(x), f(x)) over Z/NZ */
        /* g has degree ≤ 1, f has degree 2 (monic) */
        /* If g = 0, f splits completely — try gcd(coefficients, N) */
        /* If g is constant (h1 = 0), gcd = gcd(h0, N) */
        /* If g is degree 1: compute f mod g. f = x^2+a*x+b, g = h0 + h1_adj*x */
        /* We need to divide x^2+a*x+b by h0+h1_adj*x */

        mpz_t g;
        mpz_init(g);

        if (mpz_cmp_ui(h1, 0) == 0) {
            /* g is constant: gcd = gcd(h0, f) which is gcd(h0, N) if h0 != 0 */
            mpz_gcd(g, h0, N);
        } else {
            /* g(x) = h0 + h1*x. Need gcd(g, f) over Z/NZ. */
            /* f = x^2 + a*x + b. Root of g: x = -h0/h1. f at that root: */
            /* f(-h0/h1) = h0^2/h1^2 - a*h0/h1 + b = (h0^2 - a*h0*h1 + b*h1^2) / h1^2 */
            /* gcd(numerator, N) might give a factor */

            mpz_t num, h1_sq;
            mpz_inits(num, h1_sq, NULL);
            mpz_mul(num, h0, h0); /* h0^2 */
            mpz_t tmp; mpz_init(tmp);
            mpz_mul(tmp, a_coef, h0); mpz_mul(tmp, tmp, h1); /* a*h0*h1 */
            mpz_sub(num, num, tmp); /* h0^2 - a*h0*h1 */
            mpz_mul(h1_sq, h1, h1); mpz_mul(tmp, b_coef, h1_sq); /* b*h1^2 */
            mpz_add(num, num, tmp); /* h0^2 - a*h0*h1 + b*h1^2 */
            mpz_mod(num, num, N);

            mpz_gcd(g, num, N);

            /* Also check h1 itself */
            if (mpz_cmp_ui(g, 1) == 0) {
                mpz_gcd(g, h1, N);
            }

            mpz_clears(num, h1_sq, tmp, NULL);
        }

        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_set(factor, g);
            mpz_t cof; mpz_init(cof); mpz_tdiv_q(cof, N, factor);
            gmp_printf("FACTOR:%Zd\n", factor);
            gmp_fprintf(stderr, "Found (trial %d): %Zd x %Zd (%.3fs)\n",
                        trial, factor, cof, wall_time() - t0);
            mpz_clear(cof); mpz_clear(g);
            mpz_clears(a_coef, b_coef, h0, h1, NULL);
            gmp_randclear(rs); mpz_clears(N, factor, NULL); return 0;
        }

        mpz_clear(g);
        mpz_clears(a_coef, b_coef, h0, h1, NULL);
    }

    fprintf(stderr, "FAIL (%.1fs)\n", wall_time() - t0);
    gmp_randclear(rs); mpz_clears(N, factor, NULL); return 1;
}
