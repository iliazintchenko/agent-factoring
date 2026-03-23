/*
 * batch_smooth.c - Batch smoothness detection using product trees
 *
 * Implements Bernstein's "Smooth parts of integers" algorithm for
 * efficiently testing many numbers for B-smoothness simultaneously.
 *
 * The algorithm:
 * 1. Compute P = product of all primes up to B
 * 2. For each batch of M numbers v_1, ..., v_M:
 *    - Build a product tree of the v_i values
 *    - Compute remainder tree: P mod v_i for all i
 *    - Extract smooth parts via GCD
 *
 * Time complexity: O(M * log^2(B) * log(M)) vs O(M * B / ln(B)) for trial division
 *
 * This is a standalone tool and also a component for future sieve-free factoring.
 *
 * Compile: gcc -O3 -march=native -o batch_smooth library/batch_smooth.c -lgmp -lm
 * Usage: ./batch_smooth <N> [fb_size] [batch_size]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

/*
 * Product tree: given v[0..n-1], compute a binary tree where
 * level 0 = v[i], and each parent = product of two children.
 * Root = product of all v[i].
 *
 * Returns: array of levels. levels[0] = leaves, levels[k] = root.
 * Each level[j] has ceil(n / 2^j) entries.
 */
typedef struct {
    mpz_t **levels;
    int *level_sizes;
    int num_levels;
} product_tree_t;

static void build_product_tree(product_tree_t *tree, mpz_t *values, int n) {
    /* Count levels */
    int levels = 1;
    int sz = n;
    while (sz > 1) { sz = (sz + 1) / 2; levels++; }

    tree->num_levels = levels;
    tree->levels = malloc(levels * sizeof(mpz_t*));
    tree->level_sizes = malloc(levels * sizeof(int));

    /* Level 0 = input values */
    tree->level_sizes[0] = n;
    tree->levels[0] = malloc(n * sizeof(mpz_t));
    for (int i = 0; i < n; i++) {
        mpz_init_set(tree->levels[0][i], values[i]);
    }

    /* Build up */
    for (int lev = 1; lev < levels; lev++) {
        int prev_size = tree->level_sizes[lev - 1];
        int cur_size = (prev_size + 1) / 2;
        tree->level_sizes[lev] = cur_size;
        tree->levels[lev] = malloc(cur_size * sizeof(mpz_t));

        for (int i = 0; i < cur_size; i++) {
            mpz_init(tree->levels[lev][i]);
            int left = 2 * i;
            int right = 2 * i + 1;
            if (right < prev_size) {
                mpz_mul(tree->levels[lev][i],
                         tree->levels[lev-1][left],
                         tree->levels[lev-1][right]);
            } else {
                mpz_set(tree->levels[lev][i], tree->levels[lev-1][left]);
            }
        }
    }
}

/*
 * Remainder tree: given the product tree and a value P,
 * compute P mod v[i] for all i, using the product tree structure.
 *
 * Start at root: r = P mod (product of all v[i])
 * At each level: left_child_remainder = parent_remainder mod left_product
 *                right_child_remainder = parent_remainder mod right_product
 */
static void remainder_tree(product_tree_t *tree, mpz_t P, mpz_t *remainders) {
    int levels = tree->num_levels;

    /* Compute remainders top-down */
    mpz_t *cur_rems = malloc(tree->level_sizes[levels - 1] * sizeof(mpz_t));
    for (int i = 0; i < tree->level_sizes[levels - 1]; i++) {
        mpz_init(cur_rems[i]);
    }
    /* Root remainder = P mod root_product */
    mpz_mod(cur_rems[0], P, tree->levels[levels - 1][0]);

    for (int lev = levels - 2; lev >= 0; lev--) {
        int cur_size = tree->level_sizes[lev];
        int parent_size = tree->level_sizes[lev + 1];
        mpz_t *new_rems = malloc(cur_size * sizeof(mpz_t));

        for (int i = 0; i < cur_size; i++) {
            mpz_init(new_rems[i]);
            int parent = i / 2;
            /* remainder = parent_rem mod v[i]_product */
            mpz_mod(new_rems[i], cur_rems[parent], tree->levels[lev][i]);
        }

        /* Free old remainders */
        for (int i = 0; i < parent_size; i++) mpz_clear(cur_rems[i]);
        free(cur_rems);
        cur_rems = new_rems;
    }

    /* Level 0 remainders are what we want: P mod v[i] */
    for (int i = 0; i < tree->level_sizes[0]; i++) {
        mpz_set(remainders[i], cur_rems[i]);
    }

    for (int i = 0; i < tree->level_sizes[0]; i++) mpz_clear(cur_rems[i]);
    free(cur_rems);
}

static void free_product_tree(product_tree_t *tree) {
    for (int lev = 0; lev < tree->num_levels; lev++) {
        for (int i = 0; i < tree->level_sizes[lev]; i++) {
            mpz_clear(tree->levels[lev][i]);
        }
        free(tree->levels[lev]);
    }
    free(tree->levels);
    free(tree->level_sizes);
}

/*
 * Compute the smooth part of a number v with respect to primes up to B.
 * Uses P = product of primes up to B.
 * smooth_part(v) = gcd(P^k, v) for large enough k.
 * Efficient: compute gcd(P mod v, v), then iterate.
 */
static void smooth_part(mpz_t result, mpz_t v, mpz_t P_mod_v) {
    /* g = gcd(P mod v, v) */
    mpz_gcd(result, P_mod_v, v);

    /* Iterate: keep dividing v by g until g = 1 */
    mpz_t cofactor, g;
    mpz_inits(cofactor, g, NULL);
    mpz_divexact(cofactor, v, result);

    while (1) {
        mpz_gcd(g, result, cofactor);
        if (mpz_cmp_ui(g, 1) == 0) break;
        mpz_mul(result, result, g);
        mpz_divexact(cofactor, cofactor, g);
    }

    /* Wait, this isn't right. We need to repeatedly extract smooth part. */
    /* Let me redo: smooth_part = v / (cofactor that's coprime to P) */
    mpz_set(result, v);
    while (1) {
        mpz_gcd(g, P_mod_v, result);
        if (mpz_cmp_ui(g, 1) == 0) break;
        /* Divide out g repeatedly */
        while (mpz_divisible_p(result, g)) {
            mpz_divexact(result, result, g);
        }
        /* Recompute P mod result */
        mpz_mod(P_mod_v, P_mod_v, result);
    }
    /* result now contains the B-rough part. Smooth part = v / result */
    mpz_divexact(result, v, result);

    mpz_clears(cofactor, g, NULL);
}

/*
 * Batch test M numbers for B-smoothness.
 * Returns count of smooth numbers, fills smooth_flags[i] = 1 if values[i] is B-smooth.
 */
static int batch_smooth_test(mpz_t *values, int M, mpz_t P, int *smooth_flags) {
    /* Build product tree of values */
    product_tree_t tree;
    build_product_tree(&tree, values, M);

    /* Compute P^2 for better smooth extraction */
    mpz_t P_sq;
    mpz_init(P_sq);
    mpz_mul(P_sq, P, P);

    /* Compute remainder tree: P^2 mod v[i] for all i */
    mpz_t *remainders = malloc(M * sizeof(mpz_t));
    for (int i = 0; i < M; i++) mpz_init(remainders[i]);

    remainder_tree(&tree, P_sq, remainders);

    /* For each value, check if smooth */
    int count = 0;
    for (int i = 0; i < M; i++) {
        /* g = gcd(P^2 mod v[i], v[i]) */
        mpz_t g, cofactor;
        mpz_inits(g, cofactor, NULL);

        mpz_gcd(g, remainders[i], values[i]);
        mpz_divexact(cofactor, values[i], g);

        /* Remove all B-smooth factors from cofactor */
        while (mpz_cmp_ui(cofactor, 1) > 0) {
            mpz_gcd(g, remainders[i], cofactor);
            if (mpz_cmp_ui(g, 1) == 0) break;
            mpz_divexact(cofactor, cofactor, g);
        }

        smooth_flags[i] = (mpz_cmp_ui(cofactor, 1) == 0) ? 1 : 0;
        if (smooth_flags[i]) count++;

        mpz_clears(g, cofactor, NULL);
    }

    /* Cleanup */
    for (int i = 0; i < M; i++) mpz_clear(remainders[i]);
    free(remainders);
    mpz_clear(P_sq);
    free_product_tree(&tree);

    return count;
}

/*
 * Dixon's random squares factoring with batch smoothness detection.
 *
 * Strategy:
 * 1. Choose factor base of primes up to B where N is QR
 * 2. Generate random x values, compute x^2 mod N
 * 3. Batch-test x^2 mod N for B-smoothness
 * 4. Collect smooth relations
 * 5. Linear algebra over GF(2) to find square product
 * 6. Compute GCD to find factor
 */
static int dixon_factor(mpz_t N, mpz_t factor, int B_bound) {
    fprintf(stderr, "Dixon's method with batch smoothness detection\n");
    fprintf(stderr, "N = "); mpz_out_str(stderr, 10, N); fprintf(stderr, "\n");
    fprintf(stderr, "B = %d\n", B_bound);

    /* Generate primes up to B */
    int limit = B_bound;
    char *sieve = calloc(limit + 1, 1);
    for (int i = 2; i <= limit; i++) sieve[i] = 1;
    for (int i = 2; (long)i * i <= limit; i++)
        if (sieve[i]) for (int j = i * i; j <= limit; j += i) sieve[j] = 0;

    int *primes = NULL;
    int nprimes = 0;
    for (int i = 2; i <= limit; i++) {
        if (sieve[i] && mpz_kronecker_ui(N, i) >= 0) {
            primes = realloc(primes, (nprimes + 1) * sizeof(int));
            primes[nprimes++] = i;
        }
    }
    free(sieve);
    fprintf(stderr, "Factor base: %d primes\n", nprimes);

    /* Compute P = product of primes in factor base */
    mpz_t P;
    mpz_init_set_ui(P, 1);
    for (int i = 0; i < nprimes; i++) {
        mpz_mul_ui(P, P, primes[i]);
    }
    fprintf(stderr, "P has %zu bits\n", mpz_sizeinbase(P, 2));

    /* Generate batches of x^2 mod N and test for smoothness */
    mpz_t sqrt_N;
    mpz_init(sqrt_N);
    mpz_sqrt(sqrt_N, N);

    gmp_randstate_t rng;
    gmp_randinit_mt(rng);
    gmp_randseed_ui(rng, 42);

    int BATCH = 1024;
    mpz_t *x_vals = malloc(BATCH * sizeof(mpz_t));
    mpz_t *qx_vals = malloc(BATCH * sizeof(mpz_t));
    int *smooth_flags = calloc(BATCH, sizeof(int));

    for (int i = 0; i < BATCH; i++) {
        mpz_init(x_vals[i]);
        mpz_init(qx_vals[i]);
    }

    int total_tested = 0;
    int total_smooth = 0;
    int target = nprimes + 20;

    struct timespec start, now;
    clock_gettime(CLOCK_MONOTONIC, &start);

    while (total_smooth < target) {
        /* Generate batch of x values near sqrt(N) */
        for (int i = 0; i < BATCH; i++) {
            /* x = sqrt(N) + random offset */
            mpz_urandomm(x_vals[i], rng, N);
            mpz_add(x_vals[i], x_vals[i], sqrt_N);
            mpz_mod(x_vals[i], x_vals[i], N);

            /* Compute x^2 mod N */
            mpz_mul(qx_vals[i], x_vals[i], x_vals[i]);
            mpz_mod(qx_vals[i], qx_vals[i], N);
        }

        /* Batch smoothness test */
        int found = batch_smooth_test(qx_vals, BATCH, P, smooth_flags);
        total_smooth += found;
        total_tested += BATCH;

        clock_gettime(CLOCK_MONOTONIC, &now);
        double elapsed = (now.tv_sec - start.tv_sec) + (now.tv_nsec - start.tv_nsec) / 1e9;

        if (total_tested % (BATCH * 10) == 0 || found > 0) {
            fprintf(stderr, "Tested %d, found %d smooth (%.6f%%), %.1f tests/sec\n",
                    total_tested, total_smooth,
                    100.0 * total_smooth / total_tested,
                    total_tested / elapsed);
        }

        if (elapsed > 280) {
            fprintf(stderr, "Timeout\n");
            break;
        }
    }

    /* TODO: linear algebra step to combine smooth relations */
    fprintf(stderr, "Total: %d smooth out of %d tested (%.6f%%)\n",
            total_smooth, total_tested, 100.0 * total_smooth / total_tested);

    /* Cleanup */
    for (int i = 0; i < BATCH; i++) {
        mpz_clear(x_vals[i]);
        mpz_clear(qx_vals[i]);
    }
    free(x_vals);
    free(qx_vals);
    free(smooth_flags);
    free(primes);
    mpz_clear(P);
    mpz_clear(sqrt_N);
    gmp_randclear(rng);

    return 0;
}

/*
 * QS-style batch smooth: Instead of random x, use QS polynomial values
 * Q(x) = (x + ceil(sqrt(N)))^2 - N for consecutive x.
 * These values are much smaller than N, so smoothness probability is much higher.
 */
static int qs_batch_factor(mpz_t N, mpz_t factor, int B_bound, int M) {
    fprintf(stderr, "QS + batch smoothness detection\n");

    /* Generate primes up to B where N is QR */
    int limit = B_bound;
    char *sieve = calloc(limit + 1, 1);
    for (int i = 2; i <= limit; i++) sieve[i] = 1;
    for (int i = 2; (long)i * i <= limit; i++)
        if (sieve[i]) for (int j = i * i; j <= limit; j += i) sieve[j] = 0;

    int *primes = NULL;
    int nprimes = 0;
    for (int i = 2; i <= limit; i++) {
        if (sieve[i] && mpz_kronecker_ui(N, i) >= 0) {
            primes = realloc(primes, (nprimes + 1) * sizeof(int));
            primes[nprimes++] = i;
        }
    }
    free(sieve);
    fprintf(stderr, "Factor base: %d primes, B=%d\n", nprimes, B_bound);

    /* Compute P_star = product of p^floor(log_p(V)) for all FB primes
     * where V is the max value we'll test. This handles prime powers. */
    mpz_t P;
    mpz_init_set_ui(P, 1);
    /* For QS, Q(x) ≈ 2*sqrt(N)*M, so log2(V) ≈ bits/2 + log2(M) */
    int log2_V = (int)(mpz_sizeinbase(N, 2) / 2 + log2((double)M) + 1);
    for (int i = 0; i < nprimes; i++) {
        int k = (int)(log2_V / log2((double)primes[i]));
        if (k < 1) k = 1;
        mpz_t pk;
        mpz_init(pk);
        mpz_ui_pow_ui(pk, primes[i], k);
        mpz_mul(P, P, pk);
        mpz_clear(pk);
    }
    fprintf(stderr, "P_star has %zu bits (handles prime powers up to 2^%d)\n",
            mpz_sizeinbase(P, 2), log2_V);

    /* sqrt(N) */
    mpz_t sqrtN;
    mpz_init(sqrtN);
    mpz_sqrt(sqrtN, N);
    mpz_add_ui(sqrtN, sqrtN, 1); /* ceil */

    int BATCH = 4096;
    if (BATCH > M) BATCH = M;

    mpz_t *qx_vals = malloc(BATCH * sizeof(mpz_t));
    int *smooth_flags = calloc(BATCH, sizeof(int));
    for (int i = 0; i < BATCH; i++) mpz_init(qx_vals[i]);

    int total_tested = 0;
    int total_smooth = 0;
    int target = nprimes + 20;
    int x_offset = 1;

    struct timespec start, now;
    clock_gettime(CLOCK_MONOTONIC, &start);

    while (total_smooth < target) {
        /* Generate batch of Q(x) values */
        for (int i = 0; i < BATCH; i++) {
            /* Q(x) = (x + sqrtN)^2 - N = x^2 + 2*sqrtN*x + sqrtN^2 - N */
            mpz_set_si(qx_vals[i], x_offset + i);
            mpz_add(qx_vals[i], qx_vals[i], sqrtN);
            mpz_mul(qx_vals[i], qx_vals[i], qx_vals[i]);
            mpz_sub(qx_vals[i], qx_vals[i], N);
            /* Make positive */
            if (mpz_sgn(qx_vals[i]) < 0) mpz_neg(qx_vals[i], qx_vals[i]);
            if (mpz_sgn(qx_vals[i]) == 0) mpz_set_ui(qx_vals[i], 1);
        }

        /* Batch smoothness test */
        int found = batch_smooth_test(qx_vals, BATCH, P, smooth_flags);
        total_smooth += found;
        total_tested += BATCH;
        x_offset += BATCH;

        clock_gettime(CLOCK_MONOTONIC, &now);
        double elapsed = (now.tv_sec - start.tv_sec) + (now.tv_nsec - start.tv_nsec) / 1e9;

        fprintf(stderr, "x=%d: tested %d, smooth %d (%.4f%%), %.0f tests/sec, %.1f smooth/sec\n",
                x_offset, total_tested, total_smooth,
                100.0 * total_smooth / total_tested,
                total_tested / elapsed,
                total_smooth / elapsed);

        if (elapsed > 280) break;
    }

    fprintf(stderr, "\nResult: %d smooth of %d tested (%.4f%%)\n",
            total_smooth, total_tested, 100.0 * total_smooth / total_tested);
    fprintf(stderr, "Effective smooth rate: %.1f/sec\n",
            total_smooth / ((clock_gettime(CLOCK_MONOTONIC, &now),
            (now.tv_sec - start.tv_sec) + (now.tv_nsec - start.tv_nsec) / 1e9)));

    /* Cleanup */
    for (int i = 0; i < BATCH; i++) mpz_clear(qx_vals[i]);
    free(qx_vals);
    free(smooth_flags);
    free(primes);
    mpz_clear(P);
    mpz_clear(sqrtN);

    return 0;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N> [fb_bound] [sieve_size]\n", argv[0]);
        fprintf(stderr, "  fb_bound: factor base bound (default: auto)\n");
        fprintf(stderr, "  sieve_size: sieve interval size (default: auto)\n");
        return 1;
    }

    mpz_t N, factor;
    mpz_inits(N, factor, NULL);
    mpz_set_str(N, argv[1], 10);

    int digits = mpz_sizeinbase(N, 10);
    int bits = mpz_sizeinbase(N, 2);

    /* Auto-select parameters based on L(N) */
    double ln_N = bits * log(2.0);
    double ln_ln_N = log(ln_N);
    double L_N = exp(sqrt(ln_N * ln_ln_N));

    int B_bound = argc > 2 ? atoi(argv[2]) : (int)fmin(pow(L_N, 1.0/sqrt(2.0)), 500000);
    int M = argc > 3 ? atoi(argv[3]) : (int)fmin(pow(L_N, sqrt(2.0)), 10000000);

    fprintf(stderr, "N: %d digits, %d bits\n", digits, bits);
    fprintf(stderr, "L(N) = %.1f, B = %d, M = %d\n", L_N, B_bound, M);

    /* Use QS-style batch smoothness (much better than Dixon's) */
    qs_batch_factor(N, factor, B_bound, M);

    mpz_clears(N, factor, NULL);
    return 0;
}
