/*
 * BatchQS — Novel Quadratic Sieve variant using batch smoothness detection
 *
 * Key insight from Bernstein (2004): replace the sieve step with
 * batch GCD using product/remainder trees. This changes the
 * per-candidate cost from O(B/log B) (sieve) to O(log²B) (batch),
 * allowing LARGER smoothness bounds B for the same total cost.
 *
 * Larger B → higher smoothness probability → fewer candidates needed
 * → smaller sieve range → smaller sieve values → even higher smoothness
 * → potentially better L[1/2] constant
 *
 * Algorithm:
 * 1. Generate M candidates: Q(x) = (x+m)^2 - N for x in sieve range
 * 2. Compute primorial P = ∏_{p≤B} p
 * 3. Use product tree to compute all Q(x) values
 * 4. Compute P^k mod Q(x) for large k (extracting smooth part)
 * 5. The smooth part is gcd(P^k mod Q(x), Q(x))
 * 6. Candidates where Q(x) / smooth_part is 1 or a single LP are relations
 * 7. GF(2) linear algebra for congruence of squares
 *
 * Compile: gcc -O3 -o batchqs batchqs.c -lgmp -lm
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

/* Compute primorial: product of primes up to B */
static void compute_primorial(mpz_t result, unsigned long B) {
    mpz_set_ui(result, 1);
    int np;
    unsigned long *primes = gen_primes(B, &np);
    for (int i = 0; i < np; i++) {
        unsigned long pk = primes[i];
        /* Include prime powers up to B */
        while (pk <= B / primes[i]) pk *= primes[i];
        mpz_mul_ui(result, result, pk);
    }
    free(primes);
}

/* Product tree: compute product of values[0..n-1] */
/* Returns a tree where tree[level][i] = product of values in range */
typedef struct {
    mpz_t **levels;
    int *level_sizes;
    int n_levels;
} prod_tree_t;

static void build_product_tree(prod_tree_t *tree, mpz_t *values, int n) {
    /* Count levels */
    int levels = 1;
    int sz = n;
    while (sz > 1) { sz = (sz + 1) / 2; levels++; }

    tree->n_levels = levels;
    tree->levels = malloc(levels * sizeof(mpz_t *));
    tree->level_sizes = malloc(levels * sizeof(int));

    /* Level 0: the input values */
    tree->levels[0] = malloc(n * sizeof(mpz_t));
    tree->level_sizes[0] = n;
    for (int i = 0; i < n; i++) {
        mpz_init_set(tree->levels[0][i], values[i]);
    }

    /* Build upward */
    for (int lev = 1; lev < levels; lev++) {
        int prev_sz = tree->level_sizes[lev - 1];
        int cur_sz = (prev_sz + 1) / 2;
        tree->levels[lev] = malloc(cur_sz * sizeof(mpz_t));
        tree->level_sizes[lev] = cur_sz;

        for (int i = 0; i < cur_sz; i++) {
            mpz_init(tree->levels[lev][i]);
            if (2 * i + 1 < prev_sz) {
                mpz_mul(tree->levels[lev][i],
                        tree->levels[lev-1][2*i],
                        tree->levels[lev-1][2*i+1]);
            } else {
                mpz_set(tree->levels[lev][i], tree->levels[lev-1][2*i]);
            }
        }
    }
}

/* Remainder tree: given product tree of values and a number P,
 * compute P mod values[i] for each i.
 * Actually, we compute P^k mod values[i] for large k to extract smooth parts.
 * Uses the remainder tree approach: compute P mod (product tree node) top-down. */
static void compute_remainders(mpz_t *remainders, const prod_tree_t *tree,
                               const mpz_t P) {
    /* Top-down: compute P mod each node */
    int top = tree->n_levels - 1;

    /* Allocate remainder arrays for each level */
    mpz_t **rem = malloc(tree->n_levels * sizeof(mpz_t *));
    for (int lev = 0; lev < tree->n_levels; lev++) {
        rem[lev] = malloc(tree->level_sizes[lev] * sizeof(mpz_t));
        for (int i = 0; i < tree->level_sizes[lev]; i++)
            mpz_init(rem[lev][i]);
    }

    /* Top level: P mod top_product */
    mpz_mod(rem[top][0], P, tree->levels[top][0]);

    /* Descend */
    for (int lev = top - 1; lev >= 0; lev--) {
        for (int i = 0; i < tree->level_sizes[lev]; i++) {
            int parent = i / 2;
            mpz_mod(rem[lev][i], rem[lev+1][parent], tree->levels[lev][i]);
        }
    }

    /* Copy level 0 remainders to output */
    for (int i = 0; i < tree->level_sizes[0]; i++) {
        mpz_set(remainders[i], rem[0][i]);
    }

    /* Cleanup */
    for (int lev = 0; lev < tree->n_levels; lev++) {
        for (int i = 0; i < tree->level_sizes[lev]; i++)
            mpz_clear(rem[lev][i]);
        free(rem[lev]);
    }
    free(rem);
}

static void free_product_tree(prod_tree_t *tree) {
    for (int lev = 0; lev < tree->n_levels; lev++) {
        for (int i = 0; i < tree->level_sizes[lev]; i++)
            mpz_clear(tree->levels[lev][i]);
        free(tree->levels[lev]);
    }
    free(tree->levels);
    free(tree->level_sizes);
}

/* Extract smooth part of value using batch remainder:
 * smooth_part = gcd(P^k mod value, value) where P is the primorial
 * and k is large enough that P^k contains all prime power factors up to value's size */
static void extract_smooth_parts(mpz_t *smooth_parts, mpz_t *values, int n,
                                 const mpz_t primorial) {
    /* Compute primorial^2 (or higher power) to handle prime powers */
    mpz_t P2;
    mpz_init(P2);
    mpz_mul(P2, primorial, primorial);

    /* Build product tree of values */
    prod_tree_t tree;
    build_product_tree(&tree, values, n);

    /* Compute P^2 mod each value using remainder tree */
    mpz_t *remainders = malloc(n * sizeof(mpz_t));
    for (int i = 0; i < n; i++) mpz_init(remainders[i]);

    compute_remainders(remainders, &tree, P2);

    /* For each value, iteratively compute gcd to extract smooth part */
    for (int i = 0; i < n; i++) {
        /* Repeated squaring of remainder to extract all prime power factors */
        mpz_t r, g;
        mpz_inits(r, g, NULL);
        mpz_set(r, remainders[i]);

        /* Iterate: r = r^2 mod value, g = gcd(r, value) */
        /* This converges to the smooth part */
        mpz_gcd(smooth_parts[i], r, values[i]);

        /* Iterate a few times for prime powers */
        for (int iter = 0; iter < 10; iter++) {
            mpz_mul(r, r, r);
            mpz_mod(r, r, values[i]);
            mpz_gcd(g, r, values[i]);
            if (mpz_cmp(g, smooth_parts[i]) > 0) {
                mpz_set(smooth_parts[i], g);
            }
            if (mpz_cmp(smooth_parts[i], values[i]) == 0) break;
        }

        mpz_clears(r, g, NULL);
    }

    /* Cleanup */
    for (int i = 0; i < n; i++) mpz_clear(remainders[i]);
    free(remainders);
    free_product_tree(&tree);
    mpz_clear(P2);
}

/* Factor base entry */
typedef struct {
    unsigned long p;
    double logp;
} fb_entry_t;

/* Relation */
typedef struct {
    int x_off;
    int neg;
    unsigned char *exp; /* exponent parity */
    int has_lp;
    unsigned long lp;
} rel_t;

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    clock_gettime(CLOCK_MONOTONIC, &g_start);

    mpz_t N, sqrtN, factor, cofactor;
    mpz_inits(N, sqrtN, factor, cofactor, NULL);
    mpz_set_str(N, argv[1], 10);

    int ndigits = mpz_sizeinbase(N, 10);
    mpz_sqrt(sqrtN, N);

    /* Quick trial division */
    for (unsigned long d = 2; d <= 100000; d++) {
        if (mpz_divisible_ui_p(N, d)) {
            mpz_set_ui(factor, d);
            mpz_divexact(cofactor, N, factor);
            if (mpz_cmp(factor, cofactor) > 0) mpz_swap(factor, cofactor);
            gmp_printf("%Zd %Zd\n", factor, cofactor);
            fprintf(stderr, "BatchQS: trial division %.3fs\n", elapsed_sec());
            mpz_clears(N, sqrtN, factor, cofactor, NULL);
            return 0;
        }
    }

    /* Choose parameters */
    double logN = ndigits * log(10);
    double loglogN = log(logN);
    /* Start with moderate B, tuned for batch GCD */
    unsigned long B = (unsigned long)exp(0.5 * sqrt(logN * loglogN));
    if (B < 300) B = 300;
    if (B > 2000000) B = 2000000;

    unsigned long LP_bound = B * 20;

    fprintf(stderr, "BatchQS: %d digits, B=%lu, LP=%lu\n", ndigits, B, LP_bound);

    /* Build factor base */
    int nprimes;
    unsigned long *primes = gen_primes(B, &nprimes);
    fb_entry_t *fb = malloc(nprimes * sizeof(fb_entry_t));
    int fb_size = 0;
    mpz_t Nmod, pmpz;
    mpz_inits(Nmod, pmpz, NULL);

    for (int i = 0; i < nprimes; i++) {
        mpz_set_ui(pmpz, primes[i]);
        mpz_mod(Nmod, N, pmpz);
        if (primes[i] == 2 || mpz_jacobi(Nmod, pmpz) >= 0) {
            fb[fb_size].p = primes[i];
            fb[fb_size].logp = log2((double)primes[i]);
            fb_size++;
        }
    }
    free(primes);
    mpz_clears(Nmod, pmpz, NULL);

    fprintf(stderr, "BatchQS: factor base size=%d\n", fb_size);

    /* Compute primorial */
    mpz_t primorial;
    mpz_init(primorial);
    compute_primorial(primorial, B);
    fprintf(stderr, "BatchQS: primorial computed (%.1f bits), %.1fs\n",
            (double)mpz_sizeinbase(primorial, 2), elapsed_sec());

    /* Generate candidates and batch-test for smoothness */
    int needed = fb_size + 20;
    rel_t *rels = malloc(needed * 5 * sizeof(rel_t));
    int nrels = 0;

    int BATCH_SIZE = 4096;
    mpz_t *candidates = malloc(BATCH_SIZE * sizeof(mpz_t));
    mpz_t *smooth_parts = malloc(BATCH_SIZE * sizeof(mpz_t));
    int *x_offsets = malloc(BATCH_SIZE * sizeof(int));
    for (int i = 0; i < BATCH_SIZE; i++) {
        mpz_init(candidates[i]);
        mpz_init(smooth_parts[i]);
    }

    int batch_num = 0;
    int sieve_pos = 1;

    while (nrels < needed * 2 && elapsed_sec() < 260.0) {
        /* Generate a batch of candidates */
        int actual = 0;
        for (int i = 0; i < BATCH_SIZE && sieve_pos < 100000000; i++) {
            int x = (sieve_pos % 2 == 0) ? sieve_pos / 2 : -(sieve_pos / 2 + 1);
            sieve_pos++;

            /* Q(x) = (sqrtN + x)^2 - N */
            mpz_set(candidates[actual], sqrtN);
            if (x >= 0) mpz_add_ui(candidates[actual], candidates[actual], x);
            else mpz_sub_ui(candidates[actual], candidates[actual], -x);
            mpz_mul(candidates[actual], candidates[actual], candidates[actual]);
            mpz_sub(candidates[actual], candidates[actual], N);

            if (mpz_sgn(candidates[actual]) < 0)
                mpz_neg(candidates[actual], candidates[actual]);

            x_offsets[actual] = x;
            actual++;
        }

        if (actual == 0) break;

        /* Batch extract smooth parts */
        extract_smooth_parts(smooth_parts, candidates, actual, primorial);

        /* Check which candidates are (nearly) smooth */
        for (int i = 0; i < actual && nrels < needed * 5; i++) {
            mpz_t cofac;
            mpz_init(cofac);
            mpz_divexact(cofac, candidates[i], smooth_parts[i]);

            int usable = 0;
            unsigned long lp = 0;

            if (mpz_cmp_ui(cofac, 1) == 0) {
                usable = 1; /* Fully smooth */
            } else if (mpz_fits_ulong_p(cofac) && mpz_get_ui(cofac) <= LP_bound) {
                /* Single large prime */
                if (mpz_probab_prime_p(cofac, 5) > 0) {
                    usable = 1;
                    lp = mpz_get_ui(cofac);
                }
            }

            if (usable) {
                /* Factor the smooth part over the factor base to get exponents */
                unsigned char *exp = calloc(fb_size, 1);
                mpz_t tmp;
                mpz_init_set(tmp, smooth_parts[i]);

                for (int j = 0; j < fb_size; j++) {
                    while (mpz_divisible_ui_p(tmp, fb[j].p)) {
                        mpz_divexact_ui(tmp, tmp, fb[j].p);
                        exp[j] ^= 1;
                    }
                }
                mpz_clear(tmp);

                rels[nrels].x_off = x_offsets[i];
                rels[nrels].neg = (mpz_sgn(candidates[i]) < 0) ? 1 : 0;
                /* Actually, we negated already. Check original sign */
                mpz_t orig;
                mpz_init(orig);
                mpz_set(orig, sqrtN);
                int x = x_offsets[i];
                if (x >= 0) mpz_add_ui(orig, orig, x);
                else mpz_sub_ui(orig, orig, -x);
                mpz_mul(orig, orig, orig);
                mpz_sub(orig, orig, N);
                rels[nrels].neg = (mpz_sgn(orig) < 0) ? 1 : 0;
                mpz_clear(orig);

                rels[nrels].exp = exp;
                rels[nrels].has_lp = (lp > 0);
                rels[nrels].lp = lp;
                nrels++;
            }
            mpz_clear(cofac);
        }

        batch_num++;
        if (batch_num % 10 == 0) {
            fprintf(stderr, "BatchQS: batch %d, sieve_pos=%d, %d rels, %.1fs\n",
                    batch_num, sieve_pos, nrels, elapsed_sec());
        }
    }

    fprintf(stderr, "BatchQS: collected %d relations\n", nrels);

    /* Count full vs LP */
    int n_full = 0, n_lp = 0;
    for (int i = 0; i < nrels; i++) {
        if (rels[i].has_lp) n_lp++; else n_full++;
    }
    fprintf(stderr, "BatchQS: %d full + %d LP\n", n_full, n_lp);

    /* Merge LP pairs */
    /* Sort LP relations by their large prime */
    int *lp_idx = malloc(n_lp * sizeof(int));
    int li = 0;
    for (int i = 0; i < nrels; i++)
        if (rels[i].has_lp) lp_idx[li++] = i;

    /* Simple sort */
    for (int i = 0; i < n_lp - 1; i++)
        for (int j = i + 1; j < n_lp; j++)
            if (rels[lp_idx[i]].lp > rels[lp_idx[j]].lp) {
                int t = lp_idx[i]; lp_idx[i] = lp_idx[j]; lp_idx[j] = t;
            }

    /* Merge pairs */
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
    fprintf(stderr, "BatchQS: %d full + %d merged = %d usable\n",
            n_full, n_merged, total_usable);

    int found = 0;

    if (total_usable > fb_size + 1) {
        /* GF(2) Gaussian elimination */
        int ncols = 1 + fb_size; /* sign + fb */
        int nrows = total_usable;
        int words = (ncols + 63) / 64;
        int hist_words = (nrows + 63) / 64;
        unsigned long *matrix = calloc(nrows * words, sizeof(unsigned long));
        unsigned long *history = calloc(nrows * hist_words, sizeof(unsigned long));

        /* Fill: first n_full rows from full relations, then n_merged from merges */
        int *full_map = malloc(n_full * sizeof(int));
        int fi = 0;
        for (int i = 0; i < nrels; i++)
            if (!rels[i].has_lp) full_map[fi++] = i;

        for (int i = 0; i < nrows; i++) {
            history[i * hist_words + (i / 64)] |= (1UL << (i % 64));

            if (i < n_full) {
                int ri = full_map[i];
                if (rels[ri].neg)
                    matrix[i * words] |= 1UL;
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

        /* Gaussian elimination */
        int *pivot = malloc(ncols * sizeof(int));
        memset(pivot, -1, ncols * sizeof(int));

        for (int col = 0; col < ncols; col++) {
            int piv = -1;
            for (int row = 0; row < nrows; row++) {
                if (matrix[row * words + (col/64)] & (1UL << (col%64))) {
                    int ok = 1;
                    for (int c2 = 0; c2 < col; c2++)
                        if (pivot[c2] == row) { ok = 0; break; }
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

        /* Find zero rows */
        mpz_t lhs, rhs, g, qval, tmp;
        mpz_inits(lhs, rhs, g, qval, tmp, NULL);

        for (int row = 0; row < nrows && !found; row++) {
            int zero = 1;
            for (int w = 0; w < words; w++)
                if (matrix[row * words + w]) { zero = 0; break; }
            if (!zero) continue;

            /* Collect original relations */
            mpz_set_ui(lhs, 1);
            int *real_exp = calloc(fb_size, sizeof(int));

            for (int r = 0; r < nrows; r++) {
                if (!(history[row * hist_words + (r/64)] & (1UL << (r%64)))) continue;

                /* Get original relation indices */
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
                    mpz_set(tmp, sqrtN);
                    if (rels[ri].x_off >= 0) mpz_add_ui(tmp, tmp, rels[ri].x_off);
                    else mpz_sub_ui(tmp, tmp, -rels[ri].x_off);
                    mpz_mul(lhs, lhs, tmp);
                    mpz_mod(lhs, lhs, N);

                    /* Factor Q(x) for exponents */
                    mpz_mul(qval, tmp, tmp);
                    mpz_sub(qval, qval, N);
                    if (mpz_sgn(qval) < 0) mpz_neg(qval, qval);

                    for (int j = 0; j < fb_size; j++)
                        while (mpz_divisible_ui_p(qval, fb[j].p)) {
                            mpz_divexact_ui(qval, qval, fb[j].p);
                            real_exp[j]++;
                        }
                    /* LP (if any) appears an even number of times in merged pairs */
                }
            }

            /* Check even exponents */
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

        mpz_clears(lhs, rhs, g, qval, tmp, NULL);
        free(matrix); free(history); free(pivot); free(full_map);
    }

    free(merges);

    if (found) {
        mpz_divexact(cofactor, N, factor);
        if (mpz_cmp(factor, cofactor) > 0) mpz_swap(factor, cofactor);
        gmp_printf("%Zd %Zd\n", factor, cofactor);
        fprintf(stderr, "BatchQS: factored in %.3fs\n", elapsed_sec());
    } else {
        fprintf(stderr, "BatchQS: FAILED after %.3fs\n", elapsed_sec());
    }

    /* Cleanup */
    for (int i = 0; i < BATCH_SIZE; i++) {
        mpz_clear(candidates[i]);
        mpz_clear(smooth_parts[i]);
    }
    free(candidates); free(smooth_parts); free(x_offsets);
    for (int i = 0; i < nrels; i++) free(rels[i].exp);
    free(rels); free(fb);
    mpz_clears(N, sqrtN, factor, cofactor, primorial, NULL);
    return found ? 0 : 1;
}
