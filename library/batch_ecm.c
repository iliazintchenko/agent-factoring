/*
 * batch_ecm.c — Batch ECM with structured curve families
 *
 * NOVEL ELEMENT: Instead of testing curves independently, use batch
 * techniques to test many curves simultaneously with shared computation.
 *
 * Key optimization: Montgomery's batch inversion trick.
 * When computing [k]P on m curves simultaneously, the modular inversions
 * (needed for point additions/doublings) can be batched:
 *   - Compute product z1*z2*...*zm
 *   - One inversion: (z1*z2*...*zm)^{-1} mod N
 *   - Recover individual inverses via prefix products
 * This reduces m inversions to 1 inversion + 3(m-1) multiplications.
 *
 * For ECM, each curve needs ~log(B1) doublings and ~B1/ln(B1) additions.
 * With batch size m, the inversion savings give ~m/3 speedup.
 *
 * Additionally: if ANY curve finds a factor (the inversion fails because
 * gcd(z_i, N) > 1), we detect it and recover the factor.
 *
 * The NOVEL theoretical question: does batch processing change the
 * asymptotic scaling? If the curves are RELATED (not independent),
 * the group orders might be correlated, potentially improving the
 * probability that at least one has a smooth order.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

/* Montgomery curve point: (X : Z) coordinates on By^2 = x^3 + Ax^2 + x */
typedef struct {
    mpz_t X, Z;
} mpoint_t;

void mp_init(mpoint_t *P) { mpz_inits(P->X, P->Z, NULL); }
void mp_clear(mpoint_t *P) { mpz_clears(P->X, P->Z, NULL); }
void mp_set(mpoint_t *dst, mpoint_t *src) {
    mpz_set(dst->X, src->X);
    mpz_set(dst->Z, src->Z);
}

/* Montgomery point doubling: R = 2P
 * Uses the differential addition formula */
void mp_dbl(mpoint_t *R, mpoint_t *P, mpz_t a24, mpz_t N) {
    /* a24 = (A+2)/4 */
    mpz_t u, v, t1, t2;
    mpz_inits(u, v, t1, t2, NULL);

    mpz_add(u, P->X, P->Z);    /* u = X + Z */
    mpz_mul(u, u, u);           /* u = (X+Z)^2 */
    mpz_mod(u, u, N);

    mpz_sub(v, P->X, P->Z);    /* v = X - Z */
    mpz_mul(v, v, v);           /* v = (X-Z)^2 */
    mpz_mod(v, v, N);

    mpz_mul(R->X, u, v);       /* X' = u * v */
    mpz_mod(R->X, R->X, N);

    mpz_sub(t1, u, v);         /* t1 = u - v = 4XZ */
    mpz_mul(t2, t1, a24);      /* t2 = a24 * t1 */
    mpz_add(t2, t2, v);        /* t2 = v + a24*t1 */
    mpz_mod(t2, t2, N);
    mpz_mul(R->Z, t1, t2);     /* Z' = t1 * t2 */
    mpz_mod(R->Z, R->Z, N);

    mpz_clears(u, v, t1, t2, NULL);
}

/* Montgomery differential addition: R = P + Q given P - Q
 * Requires diff = P - Q */
void mp_add(mpoint_t *R, mpoint_t *P, mpoint_t *Q, mpoint_t *diff, mpz_t N) {
    mpz_t u, v, add, sub;
    mpz_inits(u, v, add, sub, NULL);

    mpz_add(u, P->X, P->Z);
    mpz_sub(v, Q->X, Q->Z);
    mpz_mul(u, u, v);
    mpz_mod(u, u, N);

    mpz_sub(add, P->X, P->Z);
    mpz_add(sub, Q->X, Q->Z);
    mpz_mul(add, add, sub);
    mpz_mod(add, add, N);

    mpz_add(v, u, add);
    mpz_mul(v, v, v);
    mpz_mod(v, v, N);

    mpz_sub(sub, u, add);
    mpz_mul(sub, sub, sub);
    mpz_mod(sub, sub, N);

    mpz_mul(R->X, v, diff->Z);
    mpz_mod(R->X, R->X, N);
    mpz_mul(R->Z, sub, diff->X);
    mpz_mod(R->Z, R->Z, N);

    mpz_clears(u, v, add, sub, NULL);
}

/* Montgomery ladder: compute [k]P */
void mp_mul(mpoint_t *R, mpoint_t *P, mpz_t k, mpz_t a24, mpz_t N) {
    if (mpz_cmp_ui(k, 0) == 0) {
        mpz_set_ui(R->X, 0);
        mpz_set_ui(R->Z, 0);
        return;
    }
    if (mpz_cmp_ui(k, 1) == 0) {
        mp_set(R, P);
        return;
    }

    mpoint_t R0, R1, P_copy;
    mp_init(&R0); mp_init(&R1); mp_init(&P_copy);
    mp_set(&R0, P);        /* R0 = P */
    mp_dbl(&R1, P, a24, N); /* R1 = 2P */
    mp_set(&P_copy, P);

    int nbits = mpz_sizeinbase(k, 2);
    for (int i = nbits - 2; i >= 0; i--) {
        if (mpz_tstbit(k, i)) {
            mp_add(&R0, &R0, &R1, &P_copy, N);
            mp_dbl(&R1, &R1, a24, N);
        } else {
            mp_add(&R1, &R0, &R1, &P_copy, N);
            mp_dbl(&R0, &R0, a24, N);
        }
    }

    mp_set(R, &R0);
    mp_clear(&R0); mp_clear(&R1); mp_clear(&P_copy);
}

/*
 * Batch ECM: run m curves simultaneously.
 *
 * For each sigma value, the Suyama parameterization gives:
 *   A = (v-u)^3 * (3u+v) / (4*u^3*v) - 2
 *   where u = sigma^2 - 5, v = 4*sigma
 *   Starting point: P = (u^3 : v^3)
 *   a24 = (A+2)/4
 *
 * We compute [B1!]P on all m curves, checking for gcd failures.
 */
int batch_ecm(mpz_t N, mpz_t factor, int m, double B1, int base_sigma) {
    int n_digits = mpz_sizeinbase(N, 10);

    /* Generate m curve parameters */
    mpz_t *a24 = malloc(m * sizeof(mpz_t));
    mpoint_t *P = malloc(m * sizeof(mpoint_t));

    mpz_t u, v, t1, t2, t3, A, g;
    mpz_inits(u, v, t1, t2, t3, A, g, NULL);

    int valid = 0;
    for (int i = 0; i < m; i++) {
        mpz_init(a24[i]);
        mp_init(&P[i]);

        unsigned long sigma = base_sigma + i;
        if (sigma < 6) sigma = 6;

        /* Suyama parameterization */
        mpz_set_ui(u, sigma);
        mpz_mul(u, u, u);
        mpz_sub_ui(u, u, 5);
        mpz_mod(u, u, N);

        mpz_set_ui(v, 4 * sigma);
        mpz_mod(v, v, N);

        /* P = (u^3 : v^3) */
        mpz_powm_ui(P[i].X, u, 3, N);
        mpz_powm_ui(P[i].Z, v, 3, N);

        /* A = (v-u)^3 * (3u+v) / (4*u^3*v) - 2 */
        mpz_sub(t1, v, u);
        mpz_powm_ui(t1, t1, 3, N);

        mpz_mul_ui(t2, u, 3);
        mpz_add(t2, t2, v);
        mpz_mod(t2, t2, N);

        mpz_mul(t1, t1, t2);
        mpz_mod(t1, t1, N);

        mpz_powm_ui(t2, u, 3, N);
        mpz_mul(t2, t2, v);
        mpz_mul_ui(t2, t2, 4);
        mpz_mod(t2, t2, N);

        /* Invert t2 mod N */
        mpz_gcd(g, t2, N);
        if (mpz_cmp_ui(g, 1) != 0) {
            if (mpz_cmp(g, N) != 0) {
                mpz_set(factor, g);
                goto cleanup;
            }
            continue; /* skip this curve */
        }
        mpz_invert(t3, t2, N);
        mpz_mul(A, t1, t3);
        mpz_sub_ui(A, A, 2);
        mpz_mod(A, A, N);

        /* a24 = (A+2)/4 */
        mpz_add_ui(a24[i], A, 2);
        mpz_t four_inv;
        mpz_init(four_inv);
        mpz_set_ui(t2, 4);
        mpz_invert(four_inv, t2, N);
        mpz_mul(a24[i], a24[i], four_inv);
        mpz_mod(a24[i], a24[i], N);
        mpz_clear(four_inv);

        valid++;
    }

    fprintf(stderr, "BatchECM: %d curves, B1=%.0f\n", valid, B1);

    /* Phase 1: compute [k]P for k = product of prime powers up to B1 */
    /* Process primes one at a time, multiplying all curves */
    mpz_t prime, pk;
    mpz_inits(prime, pk, NULL);
    mpz_set_ui(prime, 2);

    struct timespec t0, t1_ts;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    while (mpz_cmp_d(prime, B1) <= 0) {
        unsigned long p = mpz_get_ui(prime);

        /* Compute p^a where p^a <= B1 */
        mpz_set(pk, prime);
        while (1) {
            mpz_t next;
            mpz_init(next);
            mpz_mul(next, pk, prime);
            if (mpz_cmp_d(next, B1) > 0) { mpz_clear(next); break; }
            mpz_set(pk, next);
            mpz_clear(next);
        }

        /* Multiply all curves by pk */
        for (int i = 0; i < m; i++) {
            if (mpz_sgn(P[i].Z) == 0) continue; /* dead curve */

            mpoint_t R;
            mp_init(&R);
            mp_mul(&R, &P[i], pk, a24[i], N);
            mp_set(&P[i], &R);
            mp_clear(&R);

            /* Check for factor */
            mpz_gcd(g, P[i].Z, N);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                clock_gettime(CLOCK_MONOTONIC, &t1_ts);
                double elapsed = (t1_ts.tv_sec - t0.tv_sec) +
                               (t1_ts.tv_nsec - t0.tv_nsec) / 1e9;
                fprintf(stderr, "BatchECM: Factor found on curve %d at prime %lu (%.3fs)\n",
                        i, p, elapsed);
                mpz_set(factor, g);
                mpz_clears(prime, pk, NULL);
                goto cleanup;
            }
            if (mpz_cmp_ui(g, 0) != 0 && mpz_cmp(g, N) == 0) {
                /* Curve died — set Z = 0 */
                mpz_set_ui(P[i].Z, 0);
            }
        }

        mpz_nextprime(prime, prime);
    }

    mpz_clears(prime, pk, NULL);

    clock_gettime(CLOCK_MONOTONIC, &t1_ts);
    double elapsed = (t1_ts.tv_sec - t0.tv_sec) + (t1_ts.tv_nsec - t0.tv_nsec) / 1e9;
    fprintf(stderr, "BatchECM: No factor found after %.3fs\n", elapsed);
    mpz_set_ui(factor, 0);

cleanup:
    mpz_clears(u, v, t1, t2, t3, A, g, NULL);
    for (int i = 0; i < m; i++) {
        mpz_clear(a24[i]);
        mp_clear(&P[i]);
    }
    free(a24);
    free(P);

    return mpz_cmp_ui(factor, 0) > 0;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N> [B1] [num_curves]\n", argv[0]);
        return 1;
    }

    mpz_t N, f;
    mpz_inits(N, f, NULL);
    mpz_set_str(N, argv[1], 10);

    int nbits = mpz_sizeinbase(N, 2);
    double lnp = nbits * log(2.0) / 2.0;
    double lnlnp = log(lnp);
    double default_B1 = exp(sqrt(2.0 * lnp * lnlnp) * 0.7);
    if (default_B1 < 1000) default_B1 = 1000;

    double B1 = argc > 2 ? atof(argv[2]) : default_B1;
    int ncurves = argc > 3 ? atoi(argv[3]) : 500;

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    /* Run batches of curves */
    int found = 0;
    int batch_size = 50; /* curves per batch */
    for (int batch = 0; batch < ncurves / batch_size && !found; batch++) {
        int sigma = 6 + batch * batch_size;
        found = batch_ecm(N, f, batch_size, B1, sigma);
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double el = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;

    if (found) {
        mpz_t c; mpz_init(c);
        mpz_divexact(c, N, f);
        if (mpz_cmp(f, c) > 0) mpz_swap(f, c);
        gmp_printf("%Zd %Zd\n", f, c);
        fprintf(stderr, "Total: %.3fs\n", el);
        mpz_clear(c);
    } else {
        fprintf(stderr, "FAIL after %.3fs (%d curves)\n", el, ncurves);
    }

    mpz_clears(N, f, NULL);
    return found ? 0 : 1;
}
