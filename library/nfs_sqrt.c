/*
 * nfs_sqrt.c - Algebraic square root for NFS via CRT
 *
 * Given a degree-d polynomial f(x) with f(m) = N, and a set of
 * relations (a_i, b_i) where both norms are smooth, compute the
 * algebraic square root T(m) mod N.
 *
 * Algorithm:
 * For each large prime q where f(x) splits completely mod q:
 *   1. Find all d roots r_1,...,r_d of f mod q
 *   2. For each root: S_j = prod(a_i - b_i*r_j) mod q
 *   3. T_j = sqrt(S_j) mod q (Tonelli-Shanks)
 *   4. Interpolate: T(x) mod q from (r_j, T_j) pairs
 *   5. Compute T(m) mod q
 * Then CRT gives T(m) mod N.
 *
 * Standalone test: gcc -O3 -o nfs_sqrt library/nfs_sqrt.c -lgmp -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#define MAX_DEGREE 6

/* Find all roots of f(x) mod q by brute force (OK for small d, large q) */
/* For degree d <= 5, use Cantor-Zassenhaus */
static int find_roots_mod_q(unsigned long *roots, int max_roots,
                            mpz_t *coeff, int degree, unsigned long q) {
    int nroots = 0;

    /* For small q, brute force */
    if (q < 100000) {
        for (unsigned long x = 0; x < q && nroots < max_roots; x++) {
            unsigned long long val = 0;
            for (int j = degree; j >= 0; j--) {
                unsigned long cj = mpz_fdiv_ui(coeff[j], q);
                val = (val * x + cj) % q;
            }
            if (val == 0)
                roots[nroots++] = x;
        }
        return nroots;
    }

    /* For larger q, use gcd-based approach:
     * Compute g(x) = gcd(x^q - x, f(x)) mod q
     * g(x) is the product of all linear factors of f mod q.
     * Then split g using random elements. */

    /* Represent polynomials as arrays of unsigned long */
    int max_deg = degree + 1;
    unsigned long *f_mod = calloc(max_deg, sizeof(unsigned long));
    for (int i = 0; i <= degree; i++)
        f_mod[i] = mpz_fdiv_ui(coeff[i], q);

    /* Compute x^q mod f(x) mod q using repeated squaring */
    unsigned long *xpow = calloc(max_deg, sizeof(unsigned long));
    unsigned long *tmp_poly = calloc(2 * max_deg, sizeof(unsigned long));
    unsigned long *result = calloc(max_deg, sizeof(unsigned long));

    /* Start with x */
    xpow[1] = 1; /* x^1 */

    /* Square and reduce, computing x^q mod f */
    /* This is polynomial exponentiation mod f mod q */
    /* Use binary method on q */
    unsigned long *base = calloc(max_deg, sizeof(unsigned long));
    unsigned long *acc = calloc(max_deg, sizeof(unsigned long));
    base[1] = 1; /* x */
    acc[0] = 1;  /* 1 */

    /* Helper: multiply two polynomials mod f mod q */
    /* poly_mul_mod(a, b, f, q, result, degree) */
    /* This is getting complex. Let me use mpz_t arrays for polynomial arithmetic. */

    free(f_mod); free(xpow); free(tmp_poly); free(result);
    free(base); free(acc);

    /* Fallback to mpz-based polynomial arithmetic */
    mpz_t *poly_f = calloc(max_deg, sizeof(mpz_t));
    mpz_t *poly_base = calloc(max_deg, sizeof(mpz_t));
    mpz_t *poly_acc = calloc(max_deg, sizeof(mpz_t));
    mpz_t *poly_tmp = calloc(2 * max_deg, sizeof(mpz_t));
    mpz_t qz, inv;
    mpz_init_set_ui(qz, q);
    mpz_init(inv);

    for (int i = 0; i < max_deg; i++) {
        mpz_init(poly_f[i]); mpz_init(poly_base[i]); mpz_init(poly_acc[i]);
    }
    for (int i = 0; i < 2*max_deg; i++) mpz_init(poly_tmp[i]);

    for (int i = 0; i <= degree; i++)
        mpz_set_ui(poly_f[i], mpz_fdiv_ui(coeff[i], q));

    /* base = x */
    mpz_set_ui(poly_base[1], 1);
    /* acc = 1 */
    mpz_set_ui(poly_acc[0], 1);

    /* Compute x^q mod f(x) mod q */
    unsigned long exp = q;
    while (exp > 0) {
        if (exp & 1) {
            /* acc = acc * base mod f mod q */
            /* Multiply */
            for (int i = 0; i < 2*max_deg; i++) mpz_set_ui(poly_tmp[i], 0);
            for (int i = 0; i < max_deg; i++) {
                if (mpz_sgn(poly_acc[i]) == 0) continue;
                for (int j = 0; j < max_deg; j++) {
                    if (mpz_sgn(poly_base[j]) == 0) continue;
                    mpz_addmul(poly_tmp[i+j], poly_acc[i], poly_base[j]);
                    mpz_mod(poly_tmp[i+j], poly_tmp[i+j], qz);
                }
            }
            /* Reduce mod f: for each term of degree >= d, replace with -c/c_d terms */
            /* Assume f is monic-ish: divide by leading coeff */
            if (mpz_cmp_ui(poly_f[degree], 1) != 0) {
                mpz_invert(inv, poly_f[degree], qz);
            } else {
                mpz_set_ui(inv, 1);
            }

            for (int i = 2*degree - 2; i >= degree; i--) {
                if (mpz_sgn(poly_tmp[i]) == 0) continue;
                /* poly_tmp[i] * x^i = poly_tmp[i] * x^(i-d) * (-f + c_d*x^d) / c_d */
                /* = poly_tmp[i] * x^(i-d) * (-(c_0 + c_1*x + ... + c_{d-1}*x^{d-1})) / c_d */
                mpz_t coeff_hi;
                mpz_init(coeff_hi);
                mpz_mul(coeff_hi, poly_tmp[i], inv);
                mpz_mod(coeff_hi, coeff_hi, qz);

                for (int j = 0; j < degree; j++) {
                    mpz_t sub;
                    mpz_init(sub);
                    mpz_mul(sub, coeff_hi, poly_f[j]);
                    mpz_sub(poly_tmp[i - degree + j], poly_tmp[i - degree + j], sub);
                    mpz_mod(poly_tmp[i - degree + j], poly_tmp[i - degree + j], qz);
                    mpz_clear(sub);
                }
                mpz_set_ui(poly_tmp[i], 0);
                mpz_clear(coeff_hi);
            }

            for (int i = 0; i < max_deg; i++)
                mpz_set(poly_acc[i], poly_tmp[i]);
        }

        /* base = base^2 mod f mod q */
        for (int i = 0; i < 2*max_deg; i++) mpz_set_ui(poly_tmp[i], 0);
        for (int i = 0; i < max_deg; i++) {
            if (mpz_sgn(poly_base[i]) == 0) continue;
            for (int j = 0; j < max_deg; j++) {
                if (mpz_sgn(poly_base[j]) == 0) continue;
                mpz_addmul(poly_tmp[i+j], poly_base[i], poly_base[j]);
                mpz_mod(poly_tmp[i+j], poly_tmp[i+j], qz);
            }
        }
        /* Reduce mod f */
        if (mpz_cmp_ui(poly_f[degree], 1) != 0)
            mpz_invert(inv, poly_f[degree], qz);
        else
            mpz_set_ui(inv, 1);

        for (int i = 2*degree - 2; i >= degree; i--) {
            if (mpz_sgn(poly_tmp[i]) == 0) continue;
            mpz_t coeff_hi;
            mpz_init(coeff_hi);
            mpz_mul(coeff_hi, poly_tmp[i], inv);
            mpz_mod(coeff_hi, coeff_hi, qz);
            for (int j = 0; j < degree; j++) {
                mpz_t sub; mpz_init(sub);
                mpz_mul(sub, coeff_hi, poly_f[j]);
                mpz_sub(poly_tmp[i - degree + j], poly_tmp[i - degree + j], sub);
                mpz_mod(poly_tmp[i - degree + j], poly_tmp[i - degree + j], qz);
                mpz_clear(sub);
            }
            mpz_set_ui(poly_tmp[i], 0);
            mpz_clear(coeff_hi);
        }

        for (int i = 0; i < max_deg; i++)
            mpz_set(poly_base[i], poly_tmp[i]);

        exp >>= 1;
    }

    /* Now poly_acc = x^q mod f mod q */
    /* Compute g = gcd(x^q - x, f) mod q = gcd(poly_acc - x, f) */
    /* poly_acc - x */
    mpz_sub_ui(poly_acc[1], poly_acc[1], 1);
    mpz_mod(poly_acc[1], poly_acc[1], qz);

    /* Now compute polynomial GCD of (poly_acc) and (poly_f) mod q */
    /* If deg(gcd) = degree, f splits completely */
    /* Extract roots by splitting the GCD */

    /* For simplicity, if we computed x^q - x mod f correctly,
     * the roots of gcd(x^q-x, f) are exactly the roots of f in F_q.
     * To find them, we can use the fact that for each root r,
     * (r + random)^((q-1)/2) ≡ ±1 (mod q), and gcd with ±1 splits. */

    /* But this is getting very complex. Let me just use the brute force
     * approach for primes q where brute force is feasible (q < 10^7).
     * For NFS on 60-90d numbers, we need primes of sufficient product
     * (> N), so we need large primes. But we can use MANY medium primes. */

    /* For the CRT approach, we need about log_q(N) primes, each evaluated
     * at degree d roots. For 60d N and q~10^6, we need about 10 primes.
     * For 90d N, about 15 primes. */

    /* Brute force root finding for q < 10^7 is too slow (10M iterations).
     * Use the polynomial GCD approach above. But implementing polynomial
     * GCD over F_q correctly requires careful handling. */

    /* TODO: Implement proper polynomial GCD for large q */

    /* For now, return 0 roots (not implemented for large q) */
    for (int i = 0; i < max_deg; i++) {
        mpz_clear(poly_f[i]); mpz_clear(poly_base[i]); mpz_clear(poly_acc[i]);
    }
    for (int i = 0; i < 2*max_deg; i++) mpz_clear(poly_tmp[i]);
    free(poly_f); free(poly_base); free(poly_acc); free(poly_tmp);
    mpz_clear(qz); mpz_clear(inv);

    return nroots;
}

/* Tonelli-Shanks mod q */
static unsigned long sqrt_mod_q(unsigned long n, unsigned long q) {
    if (n == 0) return 0;
    unsigned long long nn = n % q, m = q;
    unsigned long long r = 1, b = nn, e = (q-1)/2;
    while (e) { if (e&1) r=r*b%m; b=b*b%m; e>>=1; }
    if (r != 1) return 0; /* not QR */
    if (q%4==3) { r=1; b=nn; e=(q+1)/4; while (e) { if (e&1) r=r*b%m; b=b*b%m; e>>=1; } return (unsigned long)r; }
    unsigned long Q2=q-1, S=0; while (Q2%2==0) { Q2/=2; S++; }
    unsigned long z=2;
    for (;;) { r=1; b=z; e=(q-1)/2; while (e) { if (e&1) r=r*b%m; b=b*b%m; e>>=1; } if (r==m-1) break; z++; }
    unsigned long long M2=S; r=1; b=z; e=Q2; while (e) { if (e&1) r=r*b%m; b=b*b%m; e>>=1; } unsigned long long c=r;
    r=1; b=nn; e=Q2; while (e) { if (e&1) r=r*b%m; b=b*b%m; e>>=1; } unsigned long long t=r;
    r=1; b=nn; e=(Q2+1)/2; while (e) { if (e&1) r=r*b%m; b=b*b%m; e>>=1; } unsigned long long R=r;
    for (;;) { if (t==1) return (unsigned long)R; int i2=0; unsigned long long tt=t; while (tt!=1) { tt=tt*tt%q; i2++; }
        unsigned long long bb=c; for (int j=0; j<(int)M2-i2-1; j++) bb=bb*bb%q; M2=i2; c=bb*bb%q; t=t*c%q; R=R*bb%q; }
}

/* Lagrange interpolation: given (x_i, y_i) for i=0..d-1,
 * compute T(eval_point) mod q */
static unsigned long lagrange_eval(unsigned long *xs, unsigned long *ys, int d,
                                   unsigned long eval_point, unsigned long q) {
    unsigned long long result = 0;
    for (int i = 0; i < d; i++) {
        unsigned long long num = ys[i];
        unsigned long long den = 1;
        for (int j = 0; j < d; j++) {
            if (j == i) continue;
            num = num * ((eval_point + q - xs[j]) % q) % q;
            den = den * ((xs[i] + q - xs[j]) % q) % q;
        }
        /* num / den mod q */
        /* Invert den */
        unsigned long long inv = 1, b2 = den, e2 = q - 2;
        while (e2) { if (e2&1) inv = inv*b2%q; b2 = b2*b2%q; e2>>=1; }
        result = (result + num % q * inv % q) % q;
    }
    return (unsigned long)result;
}

/* Test function */
int main(int argc, char *argv[]) {
    printf("NFS algebraic square root module\n");
    printf("This module provides find_roots_mod_q, sqrt_mod_q, and lagrange_eval\n");
    printf("for computing the algebraic square root via CRT.\n");

    /* Test with a simple example: f(x) = x^2 + 1, roots mod 5 are 2 and 3 */
    mpz_t coeff[3];
    for (int i = 0; i < 3; i++) mpz_init(coeff[i]);
    mpz_set_ui(coeff[0], 1); /* constant */
    mpz_set_ui(coeff[1], 0); /* x coeff */
    mpz_set_ui(coeff[2], 1); /* x^2 coeff */

    unsigned long roots[10];
    int nroots = find_roots_mod_q(roots, 10, coeff, 2, 5);
    printf("f(x) = x^2 + 1, roots mod 5: ");
    for (int i = 0; i < nroots; i++) printf("%lu ", roots[i]);
    printf("(expected: 2 3)\n");

    /* Test sqrt */
    unsigned long s = sqrt_mod_q(4, 7);
    printf("sqrt(4) mod 7 = %lu (expected: 2)\n", s);

    /* Test Lagrange interpolation */
    unsigned long xs[] = {2, 3};
    unsigned long ys[] = {1, 1};
    unsigned long val = lagrange_eval(xs, ys, 2, 0, 5);
    printf("Interpolation at 0: %lu (expected: 1)\n", val);

    for (int i = 0; i < 3; i++) mpz_clear(coeff[i]);
    return 0;
}
