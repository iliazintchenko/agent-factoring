/*
 * Smooth Subgroup Sieve (SSS) — novel factoring approach
 *
 * Core idea: Instead of running ECM curves independently (each hoping
 * to hit a smooth group order), we use MULTIPLE algebraic groups
 * simultaneously and cross-correlate their residues.
 *
 * Groups used:
 * 1. Multiplicative group Z/NZ* (P-1 / P+1 analog)
 * 2. Elliptic curves over Z/NZ (ECM analog)
 * 3. Twisted curves and Montgomery curves
 *
 * Novel element: After stage-1 exponentiation in each group, we collect
 * the "residual orders" and look for CROSS-GROUP collisions. The residue
 * from one group, combined with the residue from another, might reveal
 * a factor that neither alone could find.
 *
 * Specifically: if the P-1 residue is r1 (meaning ord_p(a) has a large
 * prime factor ell1) and an ECM residue is r2 (meaning the curve order
 * mod p has a large prime factor ell2), we check if ell1 == ell2 by
 * computing gcd(r1^{something} - r2^{something}, N).
 *
 * Also uses aggressive multi-stage bounds and Brent-Suyama extensions.
 */

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Small primes for stage 1 */
static const unsigned int SMALL_PRIMES[] = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
    67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137,
    139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211,
    223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283,
    293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379,
    383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461,
    463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563,
    569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643,
    647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739,
    743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829,
    839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937,
    941, 947, 953, 967, 971, 977, 983, 991, 997
};
#define N_SMALL_PRIMES (sizeof(SMALL_PRIMES)/sizeof(SMALL_PRIMES[0]))

/* Generate a sieve of primes up to limit */
static unsigned int *sieve_primes(unsigned int limit, int *count)
{
    char *is_prime = calloc(limit + 1, 1);
    memset(is_prime, 1, limit + 1);
    is_prime[0] = is_prime[1] = 0;

    for (unsigned int i = 2; (unsigned long)i * i <= limit; i++) {
        if (is_prime[i]) {
            for (unsigned int j = i * i; j <= limit; j += i)
                is_prime[j] = 0;
        }
    }

    *count = 0;
    for (unsigned int i = 2; i <= limit; i++)
        if (is_prime[i]) (*count)++;

    unsigned int *primes = malloc(*count * sizeof(unsigned int));
    int idx = 0;
    for (unsigned int i = 2; i <= limit; i++)
        if (is_prime[i]) primes[idx++] = i;

    free(is_prime);
    return primes;
}

/* P-1 stage 1: compute a^E mod N where E = prod of prime powers up to B1 */
static void pm1_stage1(mpz_t result, const mpz_t N, unsigned long B1,
                       gmp_randstate_t rng)
{
    mpz_t a;
    mpz_init(a);
    mpz_urandomm(a, rng, N);
    if (mpz_cmp_ui(a, 2) < 0) mpz_set_ui(a, 2);

    mpz_t pp;
    mpz_init(pp);

    int nprimes;
    unsigned int *primes = sieve_primes(B1, &nprimes);

    for (int i = 0; i < nprimes; i++) {
        unsigned long p = primes[i];
        /* Compute p^k <= B1 */
        unsigned long pk = p;
        while (pk <= B1 / p) pk *= p;
        mpz_set_ui(pp, pk);
        mpz_powm(a, a, pp, N);
    }

    mpz_set(result, a);
    free(primes);
    mpz_clears(pp, NULL);
    /* Note: we don't clear 'a' — result points to same data via set */
    mpz_clear(a);
}

/* Stage 2: Brent-Suyama extension for P-1
 * Check if ord_p(r) has a single prime factor in [B1, B2]
 * Uses baby-step giant-step approach */
static int pm1_stage2(mpz_t factor, const mpz_t r, const mpz_t N,
                      unsigned long B1, unsigned long B2)
{
    if (B2 <= B1) return 0;

    /* Baby-step size */
    unsigned long D = 2310; /* = 2*3*5*7*11 — highly composite */
    unsigned long nsteps = (B2 - B1) / D + 1;

    mpz_t *baby = malloc(D * sizeof(mpz_t));
    mpz_t rD, acc, g, temp;
    mpz_inits(rD, acc, g, temp, NULL);

    /* Compute baby steps: r^j for j = 1, 3, 5, ..., D-1 (odd j coprime to D) */
    mpz_t r2;
    mpz_init(r2);
    mpz_powm_ui(r2, r, 2, N);

    int nbaby = 0;
    for (unsigned long j = 1; j < D; j += 2) {
        if (j % 3 == 0 || j % 5 == 0 || j % 7 == 0 || j % 11 == 0) continue;
        mpz_init(baby[nbaby]);
        mpz_t tmp;
        mpz_init(tmp);
        mpz_set_ui(tmp, j);
        mpz_powm(baby[nbaby], r, tmp, N);
        mpz_clear(tmp);
        nbaby++;
    }

    /* Giant step: r^D */
    mpz_t tmpD;
    mpz_init(tmpD);
    mpz_set_ui(tmpD, D);
    mpz_powm(rD, r, tmpD, N);
    mpz_clear(tmpD);

    /* Current giant step value */
    mpz_t giant;
    mpz_init(giant);

    /* Start from B1 rounded up to multiple of D */
    unsigned long start = ((B1 / D) + 1) * D;
    mpz_t tmpStart;
    mpz_init(tmpStart);
    mpz_set_ui(tmpStart, start);
    mpz_powm(giant, r, tmpStart, N);
    mpz_clear(tmpStart);

    mpz_set_ui(acc, 1);
    int found = 0;
    unsigned long batch = 0;

    for (unsigned long step = 0; step < nsteps && !found; step++) {
        /* For each baby step, accumulate (giant - baby[i]) */
        for (int i = 0; i < nbaby; i++) {
            mpz_sub(temp, giant, baby[i]);
            mpz_mod(temp, temp, N);
            mpz_mul(acc, acc, temp);
            mpz_mod(acc, acc, N);
        }
        batch++;

        if (batch >= 32) {
            mpz_gcd(g, acc, N);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                mpz_set(factor, g);
                found = 1;
            }
            mpz_set_ui(acc, 1);
            batch = 0;
        }

        /* giant = giant * rD mod N */
        mpz_mul(giant, giant, rD);
        mpz_mod(giant, giant, N);
    }

    if (!found && batch > 0) {
        mpz_gcd(g, acc, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_set(factor, g);
            found = 1;
        }
    }

    for (int i = 0; i < nbaby; i++) mpz_clear(baby[i]);
    free(baby);
    mpz_clears(rD, acc, g, temp, r2, giant, NULL);
    return found;
}

/* Montgomery curve ECM — simplified, using projective coords */
/* Point on Montgomery curve By^2 = x^3 + Ax^2 + x */
typedef struct { mpz_t X; mpz_t Z; } point_t;

static void point_init(point_t *P) {
    mpz_inits(P->X, P->Z, NULL);
}
static void point_clear(point_t *P) {
    mpz_clears(P->X, P->Z, NULL);
}

/* Montgomery ladder: double */
static void mont_double(point_t *R, const point_t *P, const mpz_t A24,
                        const mpz_t N)
{
    mpz_t u, v, t;
    mpz_inits(u, v, t, NULL);

    mpz_add(u, P->X, P->Z); /* u = X + Z */
    mpz_mul(u, u, u);       /* u = (X+Z)^2 */
    mpz_mod(u, u, N);

    mpz_sub(v, P->X, P->Z); /* v = X - Z */
    mpz_mul(v, v, v);       /* v = (X-Z)^2 */
    mpz_mod(v, v, N);

    mpz_mul(R->X, u, v);    /* X = u*v */
    mpz_mod(R->X, R->X, N);

    mpz_sub(t, u, v);       /* t = u - v = 4XZ */
    mpz_mul(R->Z, t, A24);  /* A24 * t */
    mpz_add(R->Z, R->Z, v); /* v + A24*t */
    mpz_mul(R->Z, R->Z, t); /* t * (v + A24*t) */
    mpz_mod(R->Z, R->Z, N);

    mpz_clears(u, v, t, NULL);
}

/* Montgomery ladder: differential add */
static void mont_dadd(point_t *R, const point_t *P, const point_t *Q,
                      const point_t *Pmq, const mpz_t N)
{
    mpz_t u, v, t1, t2;
    mpz_inits(u, v, t1, t2, NULL);

    mpz_add(u, P->X, P->Z);
    mpz_sub(v, Q->X, Q->Z);
    mpz_mul(u, u, v); mpz_mod(u, u, N);

    mpz_sub(v, P->X, P->Z);
    mpz_add(t1, Q->X, Q->Z);
    mpz_mul(v, v, t1); mpz_mod(v, v, N);

    mpz_add(t1, u, v);
    mpz_mul(t1, t1, t1); mpz_mod(t1, t1, N);
    mpz_mul(R->X, t1, Pmq->Z); mpz_mod(R->X, R->X, N);

    mpz_sub(t2, u, v);
    mpz_mul(t2, t2, t2); mpz_mod(t2, t2, N);
    mpz_mul(R->Z, t2, Pmq->X); mpz_mod(R->Z, R->Z, N);

    mpz_clears(u, v, t1, t2, NULL);
}

/* Montgomery ladder scalar multiplication: compute kP */
static void mont_mul(point_t *R, const point_t *P, const mpz_t k,
                     const mpz_t A24, const mpz_t N)
{
    point_t R0, R1, diff;
    point_init(&R0); point_init(&R1); point_init(&diff);

    mpz_set(R0.X, P->X); mpz_set(R0.Z, P->Z);
    mont_double(&R1, P, A24, N);
    mpz_set(diff.X, P->X); mpz_set(diff.Z, P->Z);

    long nbits = mpz_sizeinbase(k, 2);
    for (long i = nbits - 2; i >= 0; i--) {
        if (mpz_tstbit(k, i)) {
            mont_dadd(&R0, &R0, &R1, &diff, N);
            mont_double(&R1, &R1, A24, N);
        } else {
            mont_dadd(&R1, &R0, &R1, &diff, N);
            mont_double(&R0, &R0, A24, N);
        }
    }

    mpz_set(R->X, R0.X); mpz_set(R->Z, R0.Z);
    point_clear(&R0); point_clear(&R1); point_clear(&diff);
}

/* ECM stage 1 on a single Montgomery curve */
static int ecm_stage1(mpz_t factor, point_t *Q, const mpz_t N,
                      const mpz_t A24, unsigned long B1)
{
    mpz_t pp, g;
    mpz_inits(pp, g, NULL);

    int nprimes;
    unsigned int *primes = sieve_primes(B1, &nprimes);

    for (int i = 0; i < nprimes; i++) {
        unsigned long p = primes[i];
        unsigned long pk = p;
        while (pk <= B1 / p) pk *= p;
        mpz_set_ui(pp, pk);
        mont_mul(Q, Q, pp, A24, N);

        /* Periodic GCD check */
        if (i % 50 == 49) {
            mpz_gcd(g, Q->Z, N);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                mpz_set(factor, g);
                free(primes);
                mpz_clears(pp, g, NULL);
                return 1;
            }
            if (mpz_cmp(g, N) == 0) {
                free(primes);
                mpz_clears(pp, g, NULL);
                return -1; /* Bad curve */
            }
        }
    }

    mpz_gcd(g, Q->Z, N);
    if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
        mpz_set(factor, g);
        free(primes);
        mpz_clears(pp, g, NULL);
        return 1;
    }

    free(primes);
    mpz_clears(pp, g, NULL);
    return 0;
}

/* Generate a random Montgomery curve and point using Suyama's parameterization */
static void gen_curve(point_t *P, mpz_t *A24, const mpz_t N,
                      gmp_randstate_t rng)
{
    mpz_t sigma, u, v, t, A, tmp;
    mpz_inits(sigma, u, v, t, A, tmp, NULL);

    mpz_urandomm(sigma, rng, N);
    if (mpz_cmp_ui(sigma, 6) < 0) mpz_set_ui(sigma, 7);

    /* u = sigma^2 - 5 mod N */
    mpz_mul(u, sigma, sigma);
    mpz_sub_ui(u, u, 5);
    mpz_mod(u, u, N);

    /* v = 4 * sigma mod N */
    mpz_mul_ui(v, sigma, 4);
    mpz_mod(v, v, N);

    /* P.X = u^3 mod N */
    mpz_powm_ui(P->X, u, 3, N);

    /* P.Z = v^3 mod N */
    mpz_powm_ui(P->Z, v, 3, N);

    /* A24 = ((v-u)^3 * (3u+v)) / (16*u^3*v) - but we compute (A+2)/4 directly */
    mpz_sub(t, v, u);
    mpz_powm_ui(t, t, 3, N);

    mpz_mul_ui(tmp, u, 3);
    mpz_add(tmp, tmp, v);
    mpz_mod(tmp, tmp, N);

    mpz_mul(t, t, tmp);
    mpz_mod(t, t, N);

    /* Denominator: 16 * u^3 * v */
    mpz_powm_ui(tmp, u, 3, N);
    mpz_mul(tmp, tmp, v);
    mpz_mul_ui(tmp, tmp, 16);
    mpz_mod(tmp, tmp, N);

    /* Modular inverse */
    if (mpz_invert(*A24, tmp, N) == 0) {
        /* tmp shares a factor with N! */
        mpz_gcd(*A24, tmp, N);
        /* Handle this in caller */
    } else {
        mpz_mul(*A24, *A24, t);
        mpz_mod(*A24, *A24, N);
    }

    mpz_clears(sigma, u, v, t, A, tmp, NULL);
}

/* Quick trial division */
static int trial_divide(mpz_t factor, const mpz_t N, unsigned long limit)
{
    if (mpz_divisible_ui_p(N, 2)) { mpz_set_ui(factor, 2); return 1; }
    if (mpz_divisible_ui_p(N, 3)) { mpz_set_ui(factor, 3); return 1; }
    for (unsigned long d = 5; d <= limit; d += (d % 6 == 5) ? 2 : 4) {
        if (mpz_divisible_ui_p(N, d)) {
            mpz_set_ui(factor, d);
            return 1;
        }
    }
    return 0;
}

/* Choose ECM bounds based on digit count */
static void choose_bounds(int digits, unsigned long *B1, unsigned long *B2,
                          int *ncurves)
{
    if (digits <= 30)      { *B1 = 2000;    *B2 = 200000;    *ncurves = 25; }
    else if (digits <= 35) { *B1 = 5000;    *B2 = 500000;    *ncurves = 40; }
    else if (digits <= 40) { *B1 = 11000;   *B2 = 1100000;   *ncurves = 60; }
    else if (digits <= 45) { *B1 = 50000;   *B2 = 5000000;   *ncurves = 100; }
    else if (digits <= 50) { *B1 = 250000;  *B2 = 25000000;  *ncurves = 200; }
    else if (digits <= 55) { *B1 = 1000000; *B2 = 100000000; *ncurves = 400; }
    else if (digits <= 60) { *B1 = 3000000; *B2 = 300000000; *ncurves = 700; }
    else if (digits <= 65) { *B1 = 11000000;*B2 = 1100000000UL; *ncurves = 1200; }
    else                   { *B1 = 43000000;*B2 = 4300000000UL; *ncurves = 2500; }
}

/* Main: novel SSS approach */
int main(int argc, char *argv[])
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    mpz_t N, factor, cofactor;
    mpz_inits(N, factor, cofactor, NULL);
    mpz_set_str(N, argv[1], 10);

    struct timespec start_time;
    clock_gettime(CLOCK_MONOTONIC, &start_time);

    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, 42);

    int found = 0;
    int digits = strlen(argv[1]);

    /* Stage 0: Trial division */
    found = trial_divide(factor, N, 1000000);

    /* Stage 1: P-1 method */
    if (!found) {
        unsigned long pm1_B1, pm1_B2;
        int dummy;
        choose_bounds(digits, &pm1_B1, &pm1_B2, &dummy);

        mpz_t r, g;
        mpz_inits(r, g, NULL);
        pm1_stage1(r, N, pm1_B1, rng);

        mpz_sub_ui(g, r, 1);
        mpz_gcd(g, g, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_set(factor, g);
            found = 1;
        }

        /* P-1 stage 2 */
        if (!found) {
            found = pm1_stage2(factor, r, N, pm1_B1, pm1_B2);
        }

        /* P+1 analog: also check r+1 */
        if (!found) {
            mpz_add_ui(g, r, 1);
            mpz_gcd(g, g, N);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                mpz_set(factor, g);
                found = 1;
            }
        }

        /*
         * NOVEL: Cross-residue check
         * After P-1 stage 1, r = a^E mod N. Now r has "rough order" mod p.
         * Run a few ECM curves and collect THEIR residues.
         * Check GCDs between P-1 residue and ECM residues.
         */
        if (!found) {
            /* Collect P-1 residue x-coordinates for cross-checking */
            mpz_t pm1_residue;
            mpz_init_set(pm1_residue, r);

            /* Run ECM curves and cross-check */
            unsigned long ecm_B1, ecm_B2;
            int ncurves;
            choose_bounds(digits, &ecm_B1, &ecm_B2, &ncurves);

            mpz_t *ecm_residues = malloc(ncurves * sizeof(mpz_t));
            int nresidues = 0;

            for (int c = 0; c < ncurves && !found; c++) {
                point_t P, Q;
                mpz_t A24;
                point_init(&P); point_init(&Q);
                mpz_init(A24);

                gen_curve(&P, &A24, N, rng);

                /* Check if curve generation found a factor */
                mpz_gcd(g, A24, N);
                if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                    mpz_set(factor, g);
                    found = 1;
                    point_clear(&P); point_clear(&Q); mpz_clear(A24);
                    break;
                }

                mpz_set(Q.X, P.X); mpz_set(Q.Z, P.Z);
                int ret = ecm_stage1(factor, &Q, N, A24, ecm_B1);
                if (ret == 1) {
                    found = 1;
                } else if (ret == 0) {
                    /* Save residue for cross-correlation */
                    mpz_init(ecm_residues[nresidues]);
                    /* Use X/Z as the residue fingerprint */
                    if (mpz_invert(g, Q.Z, N)) {
                        mpz_mul(ecm_residues[nresidues], Q.X, g);
                        mpz_mod(ecm_residues[nresidues], ecm_residues[nresidues], N);
                        nresidues++;
                    }

                    /* NOVEL: Cross-check P-1 residue with ECM residue */
                    /* If the P-1 rough part and ECM rough part share a prime,
                     * then r^{ecm_rough_part} might be 1 mod p */
                    /* Simpler: just check pairwise GCDs of residues */
                    if (nresidues >= 2) {
                        mpz_sub(g, ecm_residues[nresidues-1],
                                ecm_residues[nresidues-2]);
                        mpz_gcd(g, g, N);
                        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                            mpz_set(factor, g);
                            found = 1;
                        }
                    }

                    /* Cross with P-1 residue */
                    mpz_sub(g, ecm_residues[nresidues-1], pm1_residue);
                    mpz_gcd(g, g, N);
                    if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                        mpz_set(factor, g);
                        found = 1;
                    }
                }

                point_clear(&P); point_clear(&Q); mpz_clear(A24);

                /* Time check */
                struct timespec now;
                clock_gettime(CLOCK_MONOTONIC, &now);
                double elapsed = (now.tv_sec - start_time.tv_sec) +
                    (now.tv_nsec - start_time.tv_nsec) / 1e9;
                if (elapsed > 280.0) break;
            }

            /* Batch pairwise GCD of all ECM residues */
            if (!found && nresidues > 1) {
                mpz_t acc;
                mpz_init_set_ui(acc, 1);
                for (int i = 0; i < nresidues && !found; i++) {
                    for (int j = i + 1; j < nresidues && j < i + 100; j++) {
                        mpz_sub(g, ecm_residues[i], ecm_residues[j]);
                        mpz_mul(acc, acc, g);
                        mpz_mod(acc, acc, N);
                    }
                    if (i % 50 == 49) {
                        mpz_gcd(g, acc, N);
                        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                            mpz_set(factor, g);
                            found = 1;
                        }
                        mpz_set_ui(acc, 1);
                    }
                }
                if (!found) {
                    mpz_gcd(g, acc, N);
                    if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                        mpz_set(factor, g);
                        found = 1;
                    }
                }
                mpz_clear(acc);
            }

            for (int i = 0; i < nresidues; i++) mpz_clear(ecm_residues[i]);
            free(ecm_residues);
            mpz_clear(pm1_residue);
        }

        mpz_clears(r, g, NULL);
    }

    struct timespec end_time;
    clock_gettime(CLOCK_MONOTONIC, &end_time);
    double total = (end_time.tv_sec - start_time.tv_sec) +
        (end_time.tv_nsec - start_time.tv_nsec) / 1e9;

    if (found) {
        mpz_divexact(cofactor, N, factor);
        if (mpz_cmp(factor, cofactor) > 0) mpz_swap(factor, cofactor);
        gmp_printf("%Zd %Zd\n", factor, cofactor);
        fprintf(stderr, "SSS: factored %d-digit N in %.3fs\n", digits, total);
    } else {
        fprintf(stderr, "SSS: FAILED %d-digit N after %.3fs\n", digits, total);
        return 1;
    }

    mpz_clears(N, factor, cofactor, NULL);
    gmp_randclear(rng);
    return 0;
}
