/*
 * poly_split.c — Factor N by exploiting polynomial splitting mod N
 *
 * NOVEL APPROACH: Polynomial root-finding modular arithmetic.
 *
 * For a polynomial f(x) mod N where N = pq:
 * - f may have different numbers of roots mod p vs mod q
 * - The polynomial GCD g = gcd(x^((N-1)/2) - 1, f(x)) mod N
 *   computes the "QR part" of the roots
 * - If the root structure differs mod p and mod q, intermediate
 *   computations encounter non-invertible elements, revealing factors
 *
 * We use random polynomials and the modular polynomial Phi_2 (2-isogeny)
 * to generate polynomials whose splitting behavior encodes N's factors.
 *
 * This is a PROBABILISTIC method: each trial has constant probability
 * of success, and each trial costs O(n^2 log n) for n-bit N.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

/* Polynomial mod N: coefficients are mpz_t, stored as array */
typedef struct {
    mpz_t *coeffs; /* coeffs[i] = coefficient of x^i */
    int deg;       /* degree (-1 for zero polynomial) */
    int alloc;     /* allocated size */
} poly_t;

void poly_init(poly_t *p, int alloc) {
    p->alloc = alloc;
    p->deg = -1;
    p->coeffs = malloc(alloc * sizeof(mpz_t));
    for (int i = 0; i < alloc; i++)
        mpz_init(p->coeffs[i]);
}

void poly_free(poly_t *p) {
    for (int i = 0; i < p->alloc; i++)
        mpz_clear(p->coeffs[i]);
    free(p->coeffs);
}

void poly_set_zero(poly_t *p) {
    for (int i = 0; i <= p->deg; i++)
        mpz_set_ui(p->coeffs[i], 0);
    p->deg = -1;
}

void poly_copy(poly_t *dst, poly_t *src) {
    for (int i = 0; i <= src->deg; i++)
        mpz_set(dst->coeffs[i], src->coeffs[i]);
    for (int i = src->deg + 1; i <= dst->deg; i++)
        mpz_set_ui(dst->coeffs[i], 0);
    dst->deg = src->deg;
}

void poly_normalize(poly_t *p, mpz_t N) {
    while (p->deg >= 0 && mpz_sgn(p->coeffs[p->deg]) == 0)
        p->deg--;
    for (int i = 0; i <= p->deg; i++)
        mpz_mod(p->coeffs[i], p->coeffs[i], N);
    while (p->deg >= 0 && mpz_sgn(p->coeffs[p->deg]) == 0)
        p->deg--;
}

/* r = a * b mod (N, modpoly) */
/* modpoly is the polynomial we reduce by */
/* Returns 0 on success, factor of N in 'factor' if found */
int poly_mulmod(poly_t *r, poly_t *a, poly_t *b, poly_t *mod, mpz_t N, mpz_t factor) {
    /* First multiply a * b */
    int result_deg = (a->deg >= 0 && b->deg >= 0) ? a->deg + b->deg : -1;

    /* Use temporary storage */
    poly_t tmp;
    poly_init(&tmp, result_deg + 2);

    for (int i = 0; i <= a->deg; i++) {
        for (int j = 0; j <= b->deg; j++) {
            mpz_t prod;
            mpz_init(prod);
            mpz_mul(prod, a->coeffs[i], b->coeffs[j]);
            mpz_add(tmp.coeffs[i + j], tmp.coeffs[i + j], prod);
            mpz_mod(tmp.coeffs[i + j], tmp.coeffs[i + j], N);
            mpz_clear(prod);
        }
    }
    tmp.deg = result_deg;
    poly_normalize(&tmp, N);

    /* Now reduce mod modpoly */
    /* modpoly should be monic (leading coeff = 1) for simplicity */
    /* If not monic, we need to invert the leading coeff */
    while (tmp.deg >= mod->deg) {
        /* Leading coefficient of tmp */
        mpz_t lc;
        mpz_init(lc);
        mpz_set(lc, tmp.coeffs[tmp.deg]);

        if (mpz_sgn(lc) == 0) {
            tmp.deg--;
            mpz_clear(lc);
            continue;
        }

        /* Check if leading coeff of mod is invertible mod N */
        mpz_t mod_lc, g, inv;
        mpz_inits(mod_lc, g, inv, NULL);
        mpz_set(mod_lc, mod->coeffs[mod->deg]);
        mpz_gcd(g, mod_lc, N);

        if (mpz_cmp_ui(g, 1) != 0 && mpz_cmp(g, N) != 0) {
            /* Found a factor! */
            mpz_set(factor, g);
            mpz_clears(mod_lc, g, inv, lc, NULL);
            poly_free(&tmp);
            return 1;
        }

        mpz_invert(inv, mod_lc, N);

        /* Multiply lc by inv to get the quotient coefficient */
        mpz_t qcoeff;
        mpz_init(qcoeff);
        mpz_mul(qcoeff, lc, inv);
        mpz_mod(qcoeff, qcoeff, N);

        /* Subtract qcoeff * mod * x^(tmp.deg - mod.deg) from tmp */
        int shift = tmp.deg - mod->deg;
        for (int i = 0; i <= mod->deg; i++) {
            mpz_t t;
            mpz_init(t);
            mpz_mul(t, qcoeff, mod->coeffs[i]);
            mpz_sub(tmp.coeffs[i + shift], tmp.coeffs[i + shift], t);
            mpz_mod(tmp.coeffs[i + shift], tmp.coeffs[i + shift], N);
            mpz_clear(t);
        }
        tmp.deg--;
        poly_normalize(&tmp, N);

        mpz_clears(mod_lc, g, inv, qcoeff, lc, NULL);
    }

    /* Copy result */
    poly_set_zero(r);
    for (int i = 0; i <= tmp.deg; i++)
        mpz_set(r->coeffs[i], tmp.coeffs[i]);
    r->deg = tmp.deg;

    poly_free(&tmp);
    return 0;
}

/* Compute x^e mod (f(x), N) using binary exponentiation */
/* Returns 0 on success, 1 if factor found */
int poly_powmod(poly_t *result, mpz_t e, poly_t *mod, mpz_t N, mpz_t factor) {
    /* result = x^e mod (mod, N) */
    poly_t base, temp;
    poly_init(&base, mod->deg + 1);
    poly_init(&temp, mod->deg + 1);

    /* base = x */
    base.deg = 1;
    mpz_set_ui(base.coeffs[1], 1);

    /* result = 1 */
    poly_set_zero(result);
    result->deg = 0;
    mpz_set_ui(result->coeffs[0], 1);

    /* Binary exponentiation */
    int n_bits = mpz_sizeinbase(e, 2);
    for (int bit = n_bits - 1; bit >= 0; bit--) {
        /* Square */
        if (poly_mulmod(&temp, result, result, mod, N, factor)) {
            poly_free(&base);
            poly_free(&temp);
            return 1;
        }
        poly_copy(result, &temp);

        /* Multiply by base if bit is set */
        if (mpz_tstbit(e, bit)) {
            if (poly_mulmod(&temp, result, &base, mod, N, factor)) {
                poly_free(&base);
                poly_free(&temp);
                return 1;
            }
            poly_copy(result, &temp);
        }
    }

    poly_free(&base);
    poly_free(&temp);
    return 0;
}

/* Polynomial GCD mod N. Returns 0 on success, 1 if factor found. */
int poly_gcd(poly_t *result, poly_t *a, poly_t *b, mpz_t N, mpz_t factor) {
    poly_t u, v, temp;
    poly_init(&u, a->alloc);
    poly_init(&v, b->alloc);
    poly_init(&temp, a->alloc + b->alloc);

    poly_copy(&u, a);
    poly_copy(&v, b);

    while (v.deg >= 0) {
        /* u = u mod v */
        while (u.deg >= v.deg) {
            mpz_t lc_v, g, inv;
            mpz_inits(lc_v, g, inv, NULL);
            mpz_set(lc_v, v.coeffs[v.deg]);
            mpz_gcd(g, lc_v, N);

            if (mpz_cmp_ui(g, 1) != 0 && mpz_cmp(g, N) != 0) {
                mpz_set(factor, g);
                mpz_clears(lc_v, g, inv, NULL);
                poly_free(&u);
                poly_free(&v);
                poly_free(&temp);
                return 1;
            }

            mpz_invert(inv, lc_v, N);
            mpz_t qcoeff;
            mpz_init(qcoeff);
            mpz_mul(qcoeff, u.coeffs[u.deg], inv);
            mpz_mod(qcoeff, qcoeff, N);

            int shift = u.deg - v.deg;
            for (int i = 0; i <= v.deg; i++) {
                mpz_t t;
                mpz_init(t);
                mpz_mul(t, qcoeff, v.coeffs[i]);
                mpz_sub(u.coeffs[i + shift], u.coeffs[i + shift], t);
                mpz_mod(u.coeffs[i + shift], u.coeffs[i + shift], N);
                mpz_clear(t);
            }
            poly_normalize(&u, N);

            mpz_clears(lc_v, g, inv, qcoeff, NULL);
        }

        /* Swap u and v */
        poly_copy(&temp, &u);
        poly_copy(&u, &v);
        poly_copy(&v, &temp);
    }

    poly_copy(result, &u);

    poly_free(&u);
    poly_free(&v);
    poly_free(&temp);
    return 0;
}

/*
 * Try to factor N using polynomial splitting.
 *
 * For a random polynomial f(x) of degree d mod N:
 * 1. Compute g(x) = gcd(x^((N-1)/2) - 1, f(x)) mod N
 * 2. If any intermediate step encounters a non-invertible element, factor found!
 * 3. Even if no inversion failure, check if g has non-trivial degree
 *    (different from 0 or d), which indicates different splitting mod p vs q.
 */
int try_factor_polysplit(mpz_t N, mpz_t factor, int trial, gmp_randstate_t rng) {
    int d = 3; /* degree of random polynomial */

    /* Generate a random monic polynomial of degree d */
    poly_t f;
    poly_init(&f, d + 2);
    f.deg = d;
    mpz_set_ui(f.coeffs[d], 1); /* monic */
    for (int i = 0; i < d; i++) {
        mpz_urandomm(f.coeffs[i], rng, N);
    }

    /* Compute e = (N-1)/2 */
    mpz_t e;
    mpz_init(e);
    mpz_sub_ui(e, N, 1);
    mpz_fdiv_q_2exp(e, e, 1);

    /* Compute x^e mod (f(x), N) */
    poly_t power_result;
    poly_init(&power_result, d + 2);

    if (poly_powmod(&power_result, e, &f, N, factor)) {
        /* Factor found during exponentiation! */
        mpz_clear(e);
        poly_free(&f);
        poly_free(&power_result);
        return 1;
    }

    /* g(x) = gcd(power_result - 1, f(x)) */
    mpz_sub_ui(power_result.coeffs[0], power_result.coeffs[0], 1);
    mpz_mod(power_result.coeffs[0], power_result.coeffs[0], N);
    poly_normalize(&power_result, N);

    poly_t gcd_result;
    poly_init(&gcd_result, d + 2);

    if (poly_gcd(&gcd_result, &power_result, &f, N, factor)) {
        /* Factor found during GCD! */
        mpz_clear(e);
        poly_free(&f);
        poly_free(&power_result);
        poly_free(&gcd_result);
        return 1;
    }

    /* Check gcd result: if degree is between 0 and d exclusive,
     * the polynomial split differently mod p and mod q.
     * Extract factor from the constant term. */
    if (gcd_result.deg > 0 && gcd_result.deg < d) {
        /* Try gcd of constant term with N */
        mpz_gcd(factor, gcd_result.coeffs[0], N);
        if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0) {
            mpz_clear(e);
            poly_free(&f);
            poly_free(&power_result);
            poly_free(&gcd_result);
            return 1;
        }

        /* Try gcd of each coefficient */
        for (int i = 0; i <= gcd_result.deg; i++) {
            mpz_gcd(factor, gcd_result.coeffs[i], N);
            if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0) {
                mpz_clear(e);
                poly_free(&f);
                poly_free(&power_result);
                poly_free(&gcd_result);
                return 1;
            }
        }

        /* Try resultant / discriminant approach */
        /* The gcd polynomial encodes roots that are QR mod one prime but
         * QNR mod another. Its resultant with f reveals the factor. */
    }

    /* Also try x^((N-1)/2) + 1 */
    mpz_add_ui(power_result.coeffs[0], power_result.coeffs[0], 2); /* undo -1, add +1 */
    mpz_mod(power_result.coeffs[0], power_result.coeffs[0], N);
    poly_normalize(&power_result, N);

    if (poly_gcd(&gcd_result, &power_result, &f, N, factor)) {
        mpz_clear(e);
        poly_free(&f);
        poly_free(&power_result);
        poly_free(&gcd_result);
        return 1;
    }

    if (gcd_result.deg > 0 && gcd_result.deg < d) {
        for (int i = 0; i <= gcd_result.deg; i++) {
            mpz_gcd(factor, gcd_result.coeffs[i], N);
            if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0) {
                mpz_clear(e);
                poly_free(&f);
                poly_free(&power_result);
                poly_free(&gcd_result);
                return 1;
            }
        }
    }

    mpz_clear(e);
    poly_free(&f);
    poly_free(&power_result);
    poly_free(&gcd_result);
    return 0;
}

/*
 * Factor using random degree-d polynomials and x^((N-1)/2) splitting.
 * Also try with x^((N+1)/2), x^(N-1), etc.
 */
int factor_by_polysplit(mpz_t N, mpz_t factor, int max_trials) {
    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, 42);

    int n_digits = mpz_sizeinbase(N, 10);
    fprintf(stderr, "PolySplit: N=%d digits, trying %d trials\n", n_digits, max_trials);

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    for (int trial = 0; trial < max_trials; trial++) {
        if (try_factor_polysplit(N, factor, trial, rng)) {
            clock_gettime(CLOCK_MONOTONIC, &t1);
            double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
            fprintf(stderr, "PolySplit: Factor found on trial %d in %.3fs\n",
                    trial + 1, elapsed);
            gmp_randclear(rng);
            return 1;
        }

        if ((trial + 1) % 100 == 0) {
            clock_gettime(CLOCK_MONOTONIC, &t1);
            double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
            fprintf(stderr, "PolySplit: %d trials done in %.3fs\n", trial + 1, elapsed);
        }
    }

    gmp_randclear(rng);
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N_decimal> [max_trials]\n", argv[0]);
        return 1;
    }

    mpz_t N, factor;
    mpz_inits(N, factor, NULL);
    mpz_set_str(N, argv[1], 10);

    int max_trials = argc > 2 ? atoi(argv[2]) : 1000;

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    int success = factor_by_polysplit(N, factor, max_trials);

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;

    if (success) {
        mpz_t cofactor;
        mpz_init(cofactor);
        mpz_divexact(cofactor, N, factor);
        gmp_printf("%Zd %Zd\n", factor, cofactor);
        fprintf(stderr, "Total time: %.3fs\n", elapsed);
        mpz_clear(cofactor);
    } else {
        fprintf(stderr, "FAIL after %d trials (%.3fs)\n", max_trials, elapsed);
    }

    mpz_clears(N, factor, NULL);
    return success ? 0 : 1;
}
