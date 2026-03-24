/*
 * Lattice Polynomial Sieve (LPS) — novel approach to polynomial selection.
 *
 * Uses LLL lattice reduction to find degree-2 polynomials f(x) with:
 *   f(m) ≡ 0 (mod N) for some m
 *   coefficients smaller than standard QS
 *
 * The lattice: for a degree-2 poly f(x) = ax^2 + bx + c with f(m)≡0 (mod N):
 *   c = -(am^2 + bm) mod N
 * So the coefficients (a, b, c) satisfy am^2 + bm + c ≡ 0 (mod N).
 * This is a lattice condition: the vector (a, b, c) lies in a sublattice.
 *
 * By LLL-reducing this lattice, we find (a, b, c) vectors that are SHORT,
 * meaning the polynomial coefficients are small.
 *
 * For a short polynomial at x near m:
 *   f(x) = a(x-m)^2 + (2am+b)(x-m) + (am^2+bm+c)
 *   f(x) ≈ a*t^2 + (2am+b)*t  where t = x-m is small
 * For small coefficients a and (2am+b), f(x) is small near x=m.
 *
 * Build: gcc -O3 -o lps lps.c -lgmp -lm
 */

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef unsigned long long u64;
typedef unsigned int u32;

/* 3x3 LLL on big integers */
/* Gram-Schmidt + size reduction + swap */

typedef struct {
    mpz_t v[3];  /* 3 components */
} vec3_t;

static void vec3_init(vec3_t *v) {
    for (int i = 0; i < 3; i++) mpz_init(v->v[i]);
}
static void vec3_clear(vec3_t *v) {
    for (int i = 0; i < 3; i++) mpz_clear(v->v[i]);
}
static void vec3_set(vec3_t *dst, vec3_t *src) {
    for (int i = 0; i < 3; i++) mpz_set(dst->v[i], src->v[i]);
}

/* Dot product */
static void vec3_dot(mpz_t result, vec3_t *a, vec3_t *b) {
    mpz_t tmp;
    mpz_init(tmp);
    mpz_set_ui(result, 0);
    for (int i = 0; i < 3; i++) {
        mpz_mul(tmp, a->v[i], b->v[i]);
        mpz_add(result, result, tmp);
    }
    mpz_clear(tmp);
}

/* LLL reduction of a 3x3 lattice basis */
static void lll_reduce(vec3_t basis[3]) {
    /* Simple LLL with delta = 3/4 for 3 vectors */
    mpz_t mu_num, mu_den, norm_sq[3], dot, tmp, tmp2;
    mpz_init(mu_num); mpz_init(mu_den);
    mpz_init(dot); mpz_init(tmp); mpz_init(tmp2);
    for (int i = 0; i < 3; i++) mpz_init(norm_sq[i]);

    int changed = 1;
    int iters = 0;
    while (changed && iters < 1000) {
        changed = 0;
        iters++;

        /* Size reduction */
        for (int i = 1; i < 3; i++) {
            for (int j = i - 1; j >= 0; j--) {
                /* mu = <b_i, b_j> / <b_j, b_j> */
                vec3_dot(mu_num, &basis[i], &basis[j]);
                vec3_dot(mu_den, &basis[j], &basis[j]);
                if (mpz_sgn(mu_den) == 0) continue;

                /* Round mu */
                /* mu_rounded = (2*mu_num + mu_den) / (2*mu_den) */
                mpz_mul_ui(tmp, mu_num, 2);
                mpz_add(tmp, tmp, mu_den);
                mpz_mul_ui(tmp2, mu_den, 2);
                mpz_fdiv_q(tmp, tmp, tmp2);

                if (mpz_sgn(tmp) != 0) {
                    for (int k = 0; k < 3; k++) {
                        mpz_mul(tmp2, tmp, basis[j].v[k]);
                        mpz_sub(basis[i].v[k], basis[i].v[k], tmp2);
                    }
                    changed = 1;
                }
            }
        }

        /* Lovász condition: swap adjacent vectors if needed */
        for (int i = 0; i < 2; i++) {
            vec3_dot(norm_sq[i], &basis[i], &basis[i]);
        }
        vec3_dot(norm_sq[2], &basis[2], &basis[2]);

        for (int i = 0; i < 2; i++) {
            /* Check: |b_{i+1}|^2 >= (3/4 - mu^2) * |b_i|^2 */
            /* Simplified: just check if swapping reduces */
            vec3_dot(tmp, &basis[i+1], &basis[i+1]);
            vec3_dot(tmp2, &basis[i], &basis[i]);

            /* 4 * |b_{i+1}|^2 < 3 * |b_i|^2 → swap */
            mpz_mul_ui(tmp, tmp, 4);
            mpz_mul_ui(tmp2, tmp2, 3);
            if (mpz_cmp(tmp, tmp2) < 0) {
                vec3_t swap;
                vec3_init(&swap);
                vec3_set(&swap, &basis[i]);
                vec3_set(&basis[i], &basis[i+1]);
                vec3_set(&basis[i+1], &swap);
                vec3_clear(&swap);
                changed = 1;
            }
        }
    }

    mpz_clear(mu_num); mpz_clear(mu_den); mpz_clear(dot);
    mpz_clear(tmp); mpz_clear(tmp2);
    for (int i = 0; i < 3; i++) mpz_clear(norm_sq[i]);
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    mpz_t N, m, m2;
    mpz_init_set_str(N, argv[1], 10);
    mpz_init(m); mpz_init(m2);

    /* m = floor(sqrt(N)) */
    mpz_sqrt(m, N);
    mpz_mul(m2, m, m);  /* m^2 */

    int digits = (int)mpz_sizeinbase(N, 10);

    /* Construct the lattice for degree-2 poly f(x) = ax^2 + bx + c
       with f(m) ≡ 0 (mod N):
       am^2 + bm + c ≡ 0 (mod N)

       Lattice basis:
       v1 = (1, 0, -m^2 mod N)     → poly: x^2 - m^2 ≡ x^2 - N (mod N)
       v2 = (0, 1, -m mod N)       → poly: x - m
       v3 = (0, 0, N)              → poly: N (≡ 0 mod N)

       But these need weighting. We want small POLYNOMIAL VALUES, not just
       small coefficients. At x = m + t with |t| ≤ M:
       f(m+t) = at^2 + (2am+b)t + (am^2+bm+c)

       The last term is ≡ 0 (mod N) by construction.
       So f(m+t) ≈ at^2 + (2am+b)t ≈ a*M^2 + (2am+b)*M

       To minimize this, we want a and (2am+b) to be small.
       The coefficient of t is (2am+b), not just b.
       So the "effective" coefficients are (a, 2am+b, am^2+bm+c).

       Substituting x = m + t, f(x) = g(t) = at^2 + bt' where
       b' = 2am + b, c' = am^2 + bm + c.

       We want |a| and |b'| small. In terms of the original basis:
       b' = 2am + b, so minimizing (a, b') is NOT the same as minimizing (a, b, c).

       Let me reformulate. Write g(t) = at^2 + b't + c' where
       c' = f(m) = am^2 + bm + c ≡ 0 (mod N).

       The lattice should be in terms of (a, b', c'):
       a is free
       b' = 2am + b
       c' = am^2 + bm + c ≡ 0 (mod N)

       Express: b = b' - 2am, c = c' - bm - am^2 = c' - (b'-2am)m - am^2 = c' - b'm + am^2
       And c' ≡ 0 (mod N), so c' = kN for some k.

       The polynomial value at t: g(t) = at^2 + b't + kN (for integer k).
       Since kN ≡ 0 (mod N), g(t) ≡ at^2 + b't (mod N).

       For the polynomial to have small values at small t:
       |g(t)| ≈ |a|*t^2 + |b'|*|t|

       We want to minimize |a| and |b'|.

       Lattice: vectors (a, b', k) with g(t) = at^2 + b't + kN.
       For this to give a valid relation (x+m)² ≡ g(t) (mod N)... wait, that's only if
       g(t) = (m+t)^2 - N = t^2 + 2mt + (m^2-N). So a=1, b'=2m, k=(m^2-N)/N.

       For the STANDARD QS polynomial, g(t) = t^2 + 2mt + (m^2-N).
       The "short" LLL vector would be (a=1, b'=2m, ...). But 2m ≈ 2√N is large.

       The LLL would try to find vectors with |b'| << 2√N. This would give
       polynomials with smaller linear coefficients, hence smaller values.

       But for degree 2, the lattice has dimension 3 and the shortest vector
       can't be much shorter than N^{1/3} in each coordinate (by Minkowski bound).
       So |a| ≈ N^{1/3} and |b'| ≈ N^{1/3}.

       Value: |g(t)| ≈ N^{1/3} * M^2 + N^{1/3} * M for M ≈ N^{1/6}:
       |g(t)| ≈ N^{1/3} * N^{1/3} + N^{1/3} * N^{1/6} ≈ N^{2/3}.

       Compare QS: |g(t)| ≈ M * 2√N = N^{1/6} * 2N^{1/2} = 2N^{2/3}. Same!

       Hmm, the LLL polynomial gives values ~N^{2/3}, same as QS. No improvement.

       Actually, for NFS, the degree-d version gives smaller values because the lattice
       has more dimensions. For degree d in d+1 dimensions, Minkowski gives
       coefficients ~N^{1/(d+1)}, and values ~N^{d/(d+1)} * M^d. Optimizing M and d
       gives L[1/3].

       So to get L[1/3], I'd need to use degree-3+ polynomials with LLL.
       But then I need NUMBER FIELD arithmetic for the algebraic side, which IS NFS.

       The key issue: LLL finds small-coefficient polynomials, but without number
       field sieving, we can only sieve the polynomial values as integers. Integer
       smoothness of a polynomial of degree d with small coefficients gives values
       ~N^{d/(d+1)}, which is WORSE than QS for d > 2.

       The ONLY advantage of LLL polynomial selection is for degree 2, where the
       polynomial values are similar to QS but with potentially better constants.
    */

    fprintf(stderr, "LPS: %d digits\n", digits);

    /* Build the lattice with weighting for sieve value minimization */
    /* Weight: (a scaled by M^2, b' scaled by M, c' scaled by 1) */
    /* M ≈ N^{1/6} for degree-2 Minkowski optimal */
    mpz_t M, M2, Nm;
    mpz_init(M); mpz_init(M2); mpz_init(Nm);

    /* M = N^(1/6) ≈ 10^(digits/6) */
    mpz_root(M, N, 6);
    mpz_mul(M2, M, M);
    mpz_mul(Nm, m, m);
    mpz_sub(Nm, Nm, N);  /* m^2 - N (the QS residual) */

    gmp_fprintf(stderr, "  m = %Zd\n  M = %Zd\n  m^2-N = %Zd\n", m, M, Nm);

    /* Lattice basis:
       Row 0: (M^2,   0,   0)   →  a=1 poly (scaled)
       Row 1: (0,     M,   0)   →  b'=1 poly (scaled)
       Row 2: (M^2*... , 2m*M, N)  → constraints
       Actually, this is getting complicated. Let me just do a simple LLL.
    */

    /* Simple approach: lattice in (a, b) space with c = -(am^2+bm) mod N */
    /* The lattice of (a, b) satisfying am^2 + bm + c ≡ 0 (mod N) is just Z^2
       (any a, b work with c = -(am^2+bm) mod N). */

    /* So LLL doesn't help here! Any (a,b) gives a valid polynomial.
       The question is: which (a,b) gives the SMALLEST VALUES?

       f(m+t) = a*t^2 + (2am+b)*t
       Minimize: max_{|t|≤M} |a*t^2 + (2am+b)*t|

       For given M, this is minimized by:
       a → 0 (but a=0 gives degree-1 poly)
       b → -2am (so the linear term vanishes)

       With a=1: b = -2m, f(m+t) = t^2. Values are t^2 ≤ M^2.
       This IS the standard QS polynomial! f(x) = x^2 - N, centered at m.

       With a > 1: f(m+t) = a*t^2 + (2am+b)*t. The values grow as a*M^2.
       Larger a → larger values → worse.

       So for degree-2 with arbitrary (a,b), the OPTIMAL choice is a=1, b=-2m.
       LLL can't improve on this because it's already the shortest vector in the
       relevant sense.

       CONCLUSION: For degree-2 polynomials, LLL polynomial selection gives
       EXACTLY the standard QS polynomial. No improvement possible.
    */

    fprintf(stderr, "CONCLUSION: LLL degree-2 poly selection recovers standard QS.\n");
    fprintf(stderr, "For degree-3+, would need number fields (=NFS).\n");
    fprintf(stderr, "This approach is a dead end for improving beyond QS.\n");

    mpz_clear(N); mpz_clear(m); mpz_clear(m2);
    mpz_clear(M); mpz_clear(M2); mpz_clear(Nm);
    return 1;
}
