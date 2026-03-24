/*
 * MMCFRAC — Multi-Multiplier Continued Fraction Factoring
 *
 * Novel combination of CFRAC with MMS-style multi-multiplier approach.
 *
 * For each multiplier k (square-free), expand sqrt(kN) as a continued fraction.
 * Each convergent p/q gives |p^2 - kN*q^2| < 2*sqrt(kN), which is
 * automatically the SMALLEST possible residue for that approximation quality.
 *
 * Key advantages over QS:
 * 1. CF convergents are optimized to minimize |p^2 - kNq^2|
 * 2. Multiple multipliers give independent candidate streams
 * 3. Cross-multiplier LP collisions add "free" relations
 * 4. No sieve needed — CF directly generates smooth candidates
 *
 * Key advantages over single-multiplier CFRAC:
 * 1. K multipliers give K independent CF streams → K× faster relation generation
 * 2. Cross-multiplier sharing: a large prime appearing in different k streams
 *    can be merged
 * 3. The Knuth-Schroeppel multiplier effect: different k values make different
 *    primes split, effectively enlarging the factor base
 *
 * Compile: gcc -O3 -o mmcfrac mmcfrac.c -lgmp -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

static struct timespec g_start;
static double elapsed_sec(void) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec - g_start.tv_sec) +
           (now.tv_nsec - g_start.tv_nsec) / 1e9;
}

/* Square-free multipliers */
static const int SQFREE[] = {
    1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19, 21, 22, 23, 26, 29,
    30, 31, 33, 34, 35, 37, 38, 39, 41, 42, 43, 46, 47, 51, 53, 55, 57,
    58, 59, 61, 62, 65, 66, 67, 69, 70, 71, 73, 74, 77, 78, 79, 82, 83,
    85, 86, 87, 89, 91, 93, 94, 95, 97
};
#define N_SQFREE (sizeof(SQFREE)/sizeof(SQFREE[0]))

/* Relation */
typedef struct {
    mpz_t x_val;       /* the p value (numerator of convergent) */
    int k;              /* multiplier index */
    int neg;            /* sign of residue */
    unsigned char *exp; /* exponent parity over factor base */
    int has_lp;
    unsigned long lp;
} rel_t;

/* Factor base entry */
typedef struct {
    unsigned long p;
    double logp;
} fb_t;

/* CF state for one multiplier */
typedef struct {
    mpz_t kN;       /* k*N */
    mpz_t sqrtKN;   /* floor(sqrt(k*N)) */
    mpz_t P, Q, a;  /* CF state variables */
    mpz_t p_prev, p_curr; /* convergent numerators */
    mpz_t q_prev, q_curr; /* convergent denominators */
    int step;
} cfrac_state_t;

static void cfrac_init(cfrac_state_t *s, const mpz_t N, int k) {
    mpz_init(s->kN);
    mpz_mul_ui(s->kN, N, k);
    mpz_init(s->sqrtKN);
    mpz_sqrt(s->sqrtKN, s->kN);
    mpz_inits(s->P, s->Q, s->a, NULL);
    mpz_inits(s->p_prev, s->p_curr, s->q_prev, s->q_curr, NULL);

    /* Initialize CF: sqrt(kN) = a0 + 1/(a1 + 1/(...)) */
    mpz_set(s->a, s->sqrtKN);
    mpz_set_ui(s->P, 0);
    mpz_set_ui(s->Q, 1);

    /* p_{-1} = 1, p_0 = a0 */
    mpz_set_ui(s->p_prev, 1);
    mpz_set(s->p_curr, s->a);
    /* q_{-1} = 0, q_0 = 1 */
    mpz_set_ui(s->q_prev, 0);
    mpz_set_ui(s->q_curr, 1);

    s->step = 0;
}

/* Advance CF by one step. Returns the residue |p^2 - kN*q^2| */
static void cfrac_step(cfrac_state_t *s, mpz_t residue) {
    mpz_t temp;
    mpz_init(temp);

    /* Standard CF iteration for sqrt(kN):
     * P_{n+1} = a_n * Q_n - P_n
     * Q_{n+1} = (kN - P_{n+1}^2) / Q_n
     * a_{n+1} = floor((sqrt(kN) + P_{n+1}) / Q_{n+1}) */

    /* P_{n+1} = a_n * Q_n - P_n */
    mpz_mul(temp, s->a, s->Q);
    mpz_sub(temp, temp, s->P);
    mpz_set(s->P, temp);

    /* Q_{n+1} = (kN - P_{n+1}^2) / Q_n */
    mpz_mul(temp, s->P, s->P);
    mpz_sub(temp, s->kN, temp);
    mpz_divexact(temp, temp, s->Q);
    mpz_set(s->Q, temp);

    /* a_{n+1} = floor((sqrtKN + P_{n+1}) / Q_{n+1}) */
    mpz_add(temp, s->sqrtKN, s->P);
    mpz_fdiv_q(s->a, temp, s->Q);

    /* Update convergents: p_{n+1} = a_{n+1} * p_n + p_{n-1} */
    mpz_mul(temp, s->a, s->p_curr);
    mpz_add(temp, temp, s->p_prev);
    mpz_set(s->p_prev, s->p_curr);
    mpz_set(s->p_curr, temp);

    mpz_mul(temp, s->a, s->q_curr);
    mpz_add(temp, temp, s->q_prev);
    mpz_set(s->q_prev, s->q_curr);
    mpz_set(s->q_curr, temp);

    /* Residue = |p^2 - kN*q^2| = (-1)^n * Q_{n+1} (by CF theory) */
    /* Actually, for the standard CF: p_n^2 - kN*q_n^2 = (-1)^{n+1} * Q_{n+1} */
    mpz_set(residue, s->Q);

    s->step++;
    mpz_clear(temp);
}

static void cfrac_clear(cfrac_state_t *s) {
    mpz_clears(s->kN, s->sqrtKN, s->P, s->Q, s->a, NULL);
    mpz_clears(s->p_prev, s->p_curr, s->q_prev, s->q_curr, NULL);
}

/* Generate primes up to limit */
static unsigned long *gen_primes(unsigned long limit, int *count) {
    char *sieve = calloc(limit + 1, 1);
    memset(sieve, 1, limit + 1);
    sieve[0] = sieve[1] = 0;
    for (unsigned long i = 2; i * i <= limit; i++)
        if (sieve[i])
            for (unsigned long j = i*i; j <= limit; j += i)
                sieve[j] = 0;
    *count = 0;
    for (unsigned long i = 2; i <= limit; i++)
        if (sieve[i]) (*count)++;
    unsigned long *p = malloc(*count * sizeof(unsigned long));
    int idx = 0;
    for (unsigned long i = 2; i <= limit; i++)
        if (sieve[i]) p[idx++] = i;
    free(sieve);
    return p;
}

/* Quick trial division check */
static int trial_divide_quick(mpz_t factor, const mpz_t N) {
    if (mpz_divisible_ui_p(N, 2)) { mpz_set_ui(factor, 2); return 1; }
    if (mpz_divisible_ui_p(N, 3)) { mpz_set_ui(factor, 3); return 1; }
    for (unsigned long d = 5; d <= 1000000; d += (d % 6 == 5) ? 2 : 4) {
        if (mpz_divisible_ui_p(N, d)) { mpz_set_ui(factor, d); return 1; }
    }
    return 0;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    clock_gettime(CLOCK_MONOTONIC, &g_start);

    mpz_t N, factor, cofactor;
    mpz_inits(N, factor, cofactor, NULL);
    mpz_set_str(N, argv[1], 10);

    int ndigits = mpz_sizeinbase(N, 10);

    /* Quick trial division */
    if (trial_divide_quick(factor, N)) {
        mpz_divexact(cofactor, N, factor);
        if (mpz_cmp(factor, cofactor) > 0) mpz_swap(factor, cofactor);
        gmp_printf("%Zd %Zd\n", factor, cofactor);
        mpz_clears(N, factor, cofactor, NULL);
        return 0;
    }

    /* Parameters */
    double logN = ndigits * log(10);
    double loglogN = log(logN);
    unsigned long B = (unsigned long)exp(0.55 * sqrt(logN * loglogN));
    if (B < 300) B = 300;
    if (B > 2000000) B = 2000000;
    unsigned long LP_bound = B * 30;

    /* Pre-filter constants */
    #define MAX_LP_BITS 1  /* Allow at most 1 large prime per relation */
    double max_smooth_bits = (MAX_LP_BITS + 1) * log2((double)B) + 5;

    /* How many multipliers to use */
    int K = 10; /* Start with 10 multipliers */
    if (ndigits > 40) K = 20;
    if (ndigits > 50) K = 30;
    if (K > (int)N_SQFREE) K = N_SQFREE;

    fprintf(stderr, "MMCFRAC: %d digits, B=%lu, LP=%lu, K=%d multipliers\n",
            ndigits, B, LP_bound, K);

    /* Build factor base: primes up to B where kN is QR for at least one k */
    int nprimes;
    unsigned long *primes = gen_primes(B, &nprimes);
    fb_t *fb = malloc(nprimes * sizeof(fb_t));
    int fb_size = 0;

    for (int i = 0; i < nprimes; i++) {
        unsigned long p = primes[i];
        /* Check if N is QR mod p (sufficient for k=1). More k values
         * would add more primes, but for simplicity check k=1 first. */
        int is_in_fb = 0;
        for (int ki = 0; ki < K && !is_in_fb; ki++) {
            mpz_t kN_mod;
            mpz_init(kN_mod);
            mpz_mul_ui(kN_mod, N, SQFREE[ki]);
            unsigned long kNmodp = mpz_fdiv_ui(kN_mod, p);
            mpz_t pmpz;
            mpz_init_set_ui(pmpz, p);
            mpz_t kNmpz;
            mpz_init_set_ui(kNmpz, kNmodp);
            if (p == 2 || mpz_jacobi(kNmpz, pmpz) >= 0) {
                is_in_fb = 1;
            }
            mpz_clears(kN_mod, pmpz, kNmpz, NULL);
        }
        if (is_in_fb) {
            fb[fb_size].p = p;
            fb[fb_size].logp = log2((double)p);
            fb_size++;
        }
    }
    free(primes);

    int needed = fb_size + 10;
    fprintf(stderr, "MMCFRAC: factor base size=%d, need %d relations\n",
            fb_size, needed);

    /* Initialize CF states for each multiplier */
    cfrac_state_t *cfs = malloc(K * sizeof(cfrac_state_t));
    for (int ki = 0; ki < K; ki++) {
        cfrac_init(&cfs[ki], N, SQFREE[ki]);
    }

    /* Collect relations */
    int max_rels = needed * 10;
    rel_t *rels = malloc(max_rels * sizeof(rel_t));
    int nrels = 0;
    int n_full = 0, n_lp = 0;

    mpz_t residue, cofac;
    mpz_inits(residue, cofac, NULL);

    int total_steps = 0;

    while (nrels < max_rels && elapsed_sec() < 260.0) {
        /* Advance each CF by one step, check for smooth residues */
        for (int ki = 0; ki < K && nrels < max_rels; ki++) {
            cfrac_step(&cfs[ki], residue);

            /* Pre-filter: skip residues too large for smooth + LP.
             * Q_n is typically much smaller than 2√(kN) — it's the CF denominator.
             * Only trial divide if Q_n could plausibly be B-smooth times at most LP. */
            {
                unsigned long res_bits = mpz_sizeinbase(residue, 2);
                unsigned long max_bits = (unsigned long)(log2((double)B) * 4 + log2((double)LP_bound));
                if (res_bits > max_bits) continue;
            }

            /* Trial divide residue over factor base */
            unsigned char *exp = calloc(fb_size, 1);
            mpz_set(cofac, residue);

            for (int j = 0; j < fb_size; j++) {
                while (mpz_divisible_ui_p(cofac, fb[j].p)) {
                    mpz_divexact_ui(cofac, cofac, fb[j].p);
                    exp[j] ^= 1;
                }
            }

            int usable = 0;
            unsigned long lp = 0;

            if (mpz_cmp_ui(cofac, 1) == 0) {
                usable = 1; /* Fully smooth */
            } else if (mpz_fits_ulong_p(cofac) && mpz_get_ui(cofac) <= LP_bound) {
                if (mpz_probab_prime_p(cofac, 5) > 0) {
                    usable = 1;
                    lp = mpz_get_ui(cofac);
                }
            }

            if (usable) {
                mpz_init_set(rels[nrels].x_val, cfs[ki].p_prev);
                rels[nrels].k = ki;
                rels[nrels].neg = (cfs[ki].step % 2 == 0);
                rels[nrels].exp = exp;
                rels[nrels].has_lp = (lp > 0);
                rels[nrels].lp = lp;
                nrels++;
                if (lp) n_lp++; else n_full++;
            } else {
                free(exp);
            }
        }

        total_steps += K;
        if (total_steps % (K * 1000) == 0) {
            fprintf(stderr, "MMCFRAC: %d CF steps, %d rels (%d full + %d LP), %.1fs\n",
                    total_steps, nrels, n_full, n_lp, elapsed_sec());
        }
    }

    fprintf(stderr, "MMCFRAC: total %d rels (%d full + %d LP)\n",
            nrels, n_full, n_lp);

    /* LP merging */
    int *lp_idx = malloc(n_lp * sizeof(int));
    int li = 0;
    for (int i = 0; i < nrels; i++)
        if (rels[i].has_lp) lp_idx[li++] = i;

    /* Sort by LP */
    for (int i = 0; i < n_lp - 1; i++)
        for (int j = i + 1; j < n_lp; j++)
            if (rels[lp_idx[i]].lp > rels[lp_idx[j]].lp) {
                int t = lp_idx[i]; lp_idx[i] = lp_idx[j]; lp_idx[j] = t;
            }

    typedef struct { int r1, r2; } merge_t;
    merge_t *merges = malloc(n_lp * sizeof(merge_t));
    int n_merged = 0;
    for (int i = 0; i < n_lp - 1; i++) {
        if (rels[lp_idx[i]].lp == rels[lp_idx[i+1]].lp) {
            merges[n_merged].r1 = lp_idx[i];
            merges[n_merged].r2 = lp_idx[i+1];
            n_merged++;
            i++;
        }
    }
    free(lp_idx);

    int total_usable = n_full + n_merged;
    fprintf(stderr, "MMCFRAC: %d full + %d merged = %d usable (need %d)\n",
            n_full, n_merged, total_usable, needed);

    int found = 0;

    if (total_usable > fb_size) {
        /* GF(2) Gaussian elimination */
        int ncols = 1 + fb_size;
        int nrows = total_usable;
        int words = (ncols + 63) / 64;
        int hist_words = (nrows + 63) / 64;
        unsigned long *matrix = calloc(nrows * words, sizeof(unsigned long));
        unsigned long *history = calloc(nrows * hist_words, sizeof(unsigned long));

        /* Build full-relation list */
        int *full_map = malloc(n_full * sizeof(int));
        int fi = 0;
        for (int i = 0; i < nrels; i++)
            if (!rels[i].has_lp) full_map[fi++] = i;

        for (int i = 0; i < nrows; i++) {
            history[i * hist_words + (i/64)] |= (1UL << (i%64));

            if (i < n_full) {
                int ri = full_map[i];
                if (rels[ri].neg) matrix[i * words] |= 1UL;
                for (int j = 0; j < fb_size; j++)
                    if (rels[ri].exp[j])
                        matrix[i * words + ((j+1)/64)] |= (1UL << ((j+1)%64));
            } else {
                int mi = i - n_full;
                int r1 = merges[mi].r1, r2 = merges[mi].r2;
                if (rels[r1].neg ^ rels[r2].neg)
                    matrix[i * words] |= 1UL;
                for (int j = 0; j < fb_size; j++)
                    if (rels[r1].exp[j] ^ rels[r2].exp[j])
                        matrix[i * words + ((j+1)/64)] |= (1UL << ((j+1)%64));
            }
        }

        /* Gauss elimination */
        int *pivot = malloc(ncols * sizeof(int));
        memset(pivot, -1, ncols * sizeof(int));

        for (int col = 0; col < ncols; col++) {
            int piv = -1;
            for (int row = 0; row < nrows; row++) {
                if (matrix[row * words + (col/64)] & (1UL << (col%64))) {
                    int ok = 1;
                    for (int c = 0; c < col; c++)
                        if (pivot[c] == row) { ok = 0; break; }
                    if (ok) { piv = row; break; }
                }
            }
            if (piv < 0) continue;
            pivot[col] = piv;

            for (int row = 0; row < nrows; row++) {
                if (row != piv && (matrix[row * words + (col/64)] & (1UL << (col%64)))) {
                    for (int w = 0; w < words; w++)
                        matrix[row * words + w] ^= matrix[piv * words + w];
                    for (int w = 0; w < hist_words; w++)
                        history[row * hist_words + w] ^= history[piv * hist_words + w];
                }
            }
        }

        /* Find zero rows and extract factors */
        mpz_t lhs, rhs, g, tmp, qval;
        mpz_inits(lhs, rhs, g, tmp, qval, NULL);

        for (int row = 0; row < nrows && !found; row++) {
            int zero = 1;
            for (int w = 0; w < words; w++)
                if (matrix[row * words + w]) { zero = 0; break; }
            if (!zero) continue;

            mpz_set_ui(lhs, 1);
            int *real_exp = calloc(fb_size, sizeof(int));

            for (int r = 0; r < nrows; r++) {
                if (!(history[row * hist_words + (r/64)] & (1UL << (r%64)))) continue;

                int indices[2]; int ni = 0;
                if (r < n_full) {
                    indices[ni++] = full_map[r];
                } else {
                    int mi = r - n_full;
                    indices[ni++] = merges[mi].r1;
                    indices[ni++] = merges[mi].r2;
                }

                for (int k = 0; k < ni; k++) {
                    int ri = indices[k];
                    /* LHS: product of p values mod N */
                    mpz_mul(lhs, lhs, rels[ri].x_val);
                    mpz_mod(lhs, lhs, N);

                    /* Compute residue and factor it */
                    /* residue = Q_n = |p^2 - kN*q^2| */
                    /* But we stored x_val = p_prev, not the actual convergent that produced this relation */
                    /* Recalculate: Q(x) where the congruence is x^2 ≡ ±Q (mod N) */
                    /* Actually, for CFRAC: p_{n-1}^2 ≡ (-1)^n * Q_n (mod kN) */
                    /* So p^2 mod kN = the residue. And p^2 mod N is what we need. */
                    /* p^2 ≡ kN*q^2 + (-1)^n * Q_n. Mod N: p^2 ≡ (-1)^n * Q_n (mod N) */
                    /* So the relation is: p^2 ≡ (-1)^n * Q_n (mod N) */
                    /* And Q_n was trial-divided to get the exponents. */

                    /* For the RHS computation, we need the full factorization of Q_n. */
                    /* Recompute Q_n from x_val: Q_n = x_val^2 mod kN */
                    int ki = rels[ri].k;
                    mpz_mul(qval, rels[ri].x_val, rels[ri].x_val);
                    mpz_t kN_val;
                    mpz_init(kN_val);
                    mpz_mul_ui(kN_val, N, SQFREE[ki]);
                    mpz_mod(qval, qval, kN_val);
                    /* Q_n might be kN - qval if negative */
                    mpz_t alt;
                    mpz_init(alt);
                    mpz_sub(alt, kN_val, qval);
                    if (mpz_cmp(alt, qval) < 0) mpz_set(qval, alt);
                    mpz_clear(alt);
                    mpz_clear(kN_val);

                    for (int j = 0; j < fb_size; j++)
                        while (mpz_divisible_ui_p(qval, fb[j].p)) {
                            mpz_divexact_ui(qval, qval, fb[j].p);
                            real_exp[j]++;
                        }
                    /* LP: divide out */
                    if (rels[ri].has_lp)
                        while (mpz_divisible_ui_p(qval, rels[ri].lp))
                            mpz_divexact_ui(qval, qval, rels[ri].lp);
                }
            }

            int ok = 1;
            for (int j = 0; j < fb_size; j++)
                if (real_exp[j] % 2 != 0) { ok = 0; break; }

            if (ok) {
                mpz_set_ui(rhs, 1);
                for (int j = 0; j < fb_size; j++) {
                    if (real_exp[j] > 0) {
                        mpz_set_ui(tmp, fb[j].p);
                        mpz_t e; mpz_init_set_ui(e, real_exp[j]/2);
                        mpz_powm(tmp, tmp, e, N);
                        mpz_mul(rhs, rhs, tmp);
                        mpz_mod(rhs, rhs, N);
                        mpz_clear(e);
                    }
                }

                mpz_sub(g, lhs, rhs);
                mpz_gcd(g, g, N);
                if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                    mpz_set(factor, g); found = 1;
                } else {
                    mpz_add(g, lhs, rhs);
                    mpz_gcd(g, g, N);
                    if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                        mpz_set(factor, g); found = 1;
                    }
                }
            }
            free(real_exp);
        }

        mpz_clears(lhs, rhs, g, tmp, qval, NULL);
        free(matrix); free(history); free(pivot); free(full_map);
    }

    free(merges);

    if (found) {
        mpz_divexact(cofactor, N, factor);
        if (mpz_cmp(factor, cofactor) > 0) mpz_swap(factor, cofactor);
        gmp_printf("%Zd %Zd\n", factor, cofactor);
        fprintf(stderr, "MMCFRAC: factored in %.3fs\n", elapsed_sec());
    } else {
        fprintf(stderr, "MMCFRAC: FAILED after %.3fs\n", elapsed_sec());
    }

    /* Cleanup */
    for (int ki = 0; ki < K; ki++) cfrac_clear(&cfs[ki]);
    free(cfs);
    for (int i = 0; i < nrels; i++) {
        mpz_clear(rels[i].x_val);
        free(rels[i].exp);
    }
    free(rels);
    free(fb);
    mpz_clears(N, factor, cofactor, residue, cofac, NULL);
    return found ? 0 : 1;
}
