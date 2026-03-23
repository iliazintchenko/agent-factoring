/*
 * batch_qs.c - Quadratic Sieve with Batch Smoothness Detection
 *
 * Novel approach: Instead of traditional sieving (O(M*B) per polynomial),
 * uses Bernstein's product tree / remainder tree technique to detect
 * B-smooth numbers in batch. This is O(M * log²M * log B) which can
 * be faster for large factor bases.
 *
 * Algorithm:
 * 1. Generate many Q(x) = (x + floor(sqrt(N)))² - N values
 * 2. Compute P = product of all primes up to B
 * 3. Use remainder tree: compute P^k mod each Q(x)
 * 4. GCD(P^k mod Q(x), Q(x)) gives the B-smooth part
 * 5. If smooth part == |Q(x)|, the number is fully B-smooth
 * 6. Allow single/double large primes for partial relations
 *
 * Compile: gcc -O3 -march=native -o batch_qs library/batch_qs.c -lgmp -lm
 * Usage: ./batch_qs <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

/* Parameters - will be auto-tuned based on N size */
typedef struct {
    int fb_size;        /* Factor base size (number of primes) */
    int batch_size;     /* Number of candidates per batch */
    int lp_bits;        /* Large prime bound in bits */
    int num_extra;      /* Extra relations beyond fb_size */
    int sieve_radius;   /* How far to search from sqrt(N) */
} params_t;

/* Factor base entry */
typedef struct {
    unsigned int p;
    unsigned int r1, r2;  /* sqrt(N) mod p */
    double logp;
} fb_entry_t;

/* Relation: Q(x) = product of factor base primes (with exponents) */
typedef struct {
    int x_val;
    int x_val2;          /* second x value for merged SLP pairs (0 if single) */
    mpz_t qx;           /* Q(x) value (or product for merged) */
    unsigned char *exponents; /* exponent vector mod 2 */
    unsigned int lp1, lp2;   /* large primes (0 if none) */
    int sign;            /* sign of Q(x) */
    int is_merged;       /* 1 if this is a merged SLP pair */
} relation_t;

/* Global state */
static mpz_t N, sqrtN;
static fb_entry_t *fb;
static int fb_count;
static relation_t *relations;
static int rel_count, rel_alloc;
static int *primes;
static int prime_count;

/* Sieve of Eratosthenes up to limit */
static void generate_primes(int limit) {
    char *sieve = calloc(limit + 1, 1);
    prime_count = 0;
    primes = malloc(sizeof(int) * (limit / 2 + 100));

    for (int i = 2; i <= limit; i++) {
        if (!sieve[i]) {
            primes[prime_count++] = i;
            for (long j = (long)i * i; j <= limit; j += i)
                sieve[j] = 1;
        }
    }
    free(sieve);
}

/* Tonelli-Shanks: compute sqrt(N) mod p */
static int mod_sqrt(mpz_t n, int p) {
    if (p == 2) return mpz_odd_p(n) ? 1 : 0;

    mpz_t tmp, nn;
    mpz_init(tmp);
    mpz_init(nn);
    mpz_mod_ui(nn, n, p);

    /* Check quadratic residue */
    mpz_set_ui(tmp, p);
    mpz_powm_ui(tmp, nn, (p - 1) / 2, tmp);
    if (mpz_cmp_ui(tmp, 1) != 0 && mpz_sgn(tmp) != 0) {
        mpz_clear(tmp);
        mpz_clear(nn);
        return -1; /* Not a QR */
    }

    unsigned long n_mod = mpz_get_ui(nn);
    mpz_clear(tmp);
    mpz_clear(nn);

    if (n_mod == 0) return 0;

    /* Simple case: p ≡ 3 (mod 4) */
    if (p % 4 == 3) {
        mpz_t base, mod, exp;
        mpz_init_set_ui(base, n_mod);
        mpz_init_set_ui(mod, p);
        mpz_init_set_ui(exp, (p + 1) / 4);
        mpz_powm(base, base, exp, mod);
        int result = (int)mpz_get_ui(base);
        mpz_clear(base);
        mpz_clear(mod);
        mpz_clear(exp);
        return result;
    }

    /* Tonelli-Shanks */
    int Q = p - 1, S = 0;
    while (Q % 2 == 0) { Q /= 2; S++; }

    /* Find quadratic non-residue */
    int z = 2;
    {
        mpz_t ztmp, ptmp;
        mpz_init_set_ui(ptmp, p);
        mpz_init(ztmp);
        while (z < p) {
            mpz_set_ui(ztmp, z);
            mpz_powm_ui(ztmp, ztmp, (p - 1) / 2, ptmp);
            if (mpz_cmp_ui(ztmp, p - 1) == 0) break;
            z++;
        }
        mpz_clear(ztmp);
        mpz_clear(ptmp);
    }

    mpz_t M, c, t, R, ptmp, qtmp;
    mpz_init_set_ui(M, S);
    mpz_init(c);
    mpz_init(t);
    mpz_init(R);
    mpz_init_set_ui(ptmp, p);
    mpz_init_set_ui(qtmp, Q);

    mpz_set_ui(c, z);
    mpz_powm(c, c, qtmp, ptmp);

    mpz_set_ui(t, n_mod);
    mpz_powm(t, t, qtmp, ptmp);

    mpz_set_ui(R, n_mod);
    mpz_set_ui(qtmp, (Q + 1) / 2);
    mpz_powm(R, R, qtmp, ptmp);

    while (1) {
        if (mpz_cmp_ui(t, 0) == 0) {
            mpz_clear(M); mpz_clear(c); mpz_clear(t); mpz_clear(R); mpz_clear(ptmp); mpz_clear(qtmp);
            return 0;
        }
        if (mpz_cmp_ui(t, 1) == 0) {
            int result = (int)mpz_get_ui(R);
            mpz_clear(M); mpz_clear(c); mpz_clear(t); mpz_clear(R); mpz_clear(ptmp); mpz_clear(qtmp);
            return result;
        }

        /* Find least i such that t^(2^i) ≡ 1 (mod p) */
        int i = 0;
        mpz_set(qtmp, t);
        while (mpz_cmp_ui(qtmp, 1) != 0) {
            mpz_mul(qtmp, qtmp, qtmp);
            mpz_mod(qtmp, qtmp, ptmp);
            i++;
        }

        int m_val = (int)mpz_get_ui(M);
        if (i == m_val) {
            mpz_clear(M); mpz_clear(c); mpz_clear(t); mpz_clear(R); mpz_clear(ptmp); mpz_clear(qtmp);
            return -1;
        }

        /* b = c^(2^(M-i-1)) */
        mpz_t b;
        mpz_init_set(b, c);
        for (int j = 0; j < m_val - i - 1; j++) {
            mpz_mul(b, b, b);
            mpz_mod(b, b, ptmp);
        }

        mpz_set_ui(M, i);
        mpz_mul(c, b, b);
        mpz_mod(c, c, ptmp);
        mpz_mul(t, t, c);
        mpz_mod(t, t, ptmp);
        mpz_mul(R, R, b);
        mpz_mod(R, R, ptmp);

        mpz_clear(b);
    }
}

/* Build factor base: all primes up to smooth_bound where N is QR */
static int build_factor_base(int smooth_bound) {
    generate_primes(smooth_bound);

    fb = malloc(sizeof(fb_entry_t) * (prime_count + 2));
    fb_count = 0;

    /* p = -1 (sign) */
    fb[fb_count].p = 0; /* sentinel for sign */
    fb[fb_count].logp = 0;
    fb_count++;

    for (int i = 0; i < prime_count; i++) {
        int p = primes[i];
        int r = mod_sqrt(N, p);
        if (r < 0) continue;

        fb[fb_count].p = p;
        fb[fb_count].r1 = r;
        fb[fb_count].r2 = (p - r) % p;
        fb[fb_count].logp = log2(p);
        fb_count++;
    }

    return fb_count;
}

/* Compute Q(x) = (x + sqrtN)^2 - N */
static void compute_qx(mpz_t result, int x) {
    mpz_set(result, sqrtN);
    if (x >= 0)
        mpz_add_ui(result, result, x);
    else
        mpz_sub_ui(result, result, -x);
    mpz_mul(result, result, result);
    mpz_sub(result, result, N);
}

/*
 * BATCH SMOOTHNESS DETECTION
 *
 * Given a batch of values v[0..n-1], and a primorial P (product of
 * factor base primes raised to appropriate powers), we want to find
 * which v[i] are B-smooth.
 *
 * Method:
 * 1. Build product tree of |v[i]| values
 * 2. Compute P mod (product tree root) using remainder tree
 * 3. The GCD of the remainder with each v[i] gives smooth part
 *
 * Optimization: Instead of computing full primorial, we use iterated
 * squaring: compute P^(2^k) mod v[i] where k is chosen so smooth
 * parts are fully extracted.
 */

/* Product tree node */
typedef struct ptree_node {
    mpz_t val;
} ptree_node_t;

/* Build product tree bottom-up. Returns array of levels.
 * Level 0 = leaves (the input values)
 * Level k = pairwise products of level k-1
 * Top level has 1 element = product of all inputs
 */
static mpz_t **build_product_tree(mpz_t *vals, int n, int *num_levels) {
    /* Calculate number of levels */
    int levels = 1;
    int sz = n;
    while (sz > 1) { sz = (sz + 1) / 2; levels++; }
    *num_levels = levels;

    /* Allocate level arrays */
    mpz_t **tree = malloc(sizeof(mpz_t*) * levels);
    int *level_size = malloc(sizeof(int) * levels);
    level_size[0] = n;

    /* Level 0: copy input values */
    tree[0] = malloc(sizeof(mpz_t) * n);
    for (int i = 0; i < n; i++) {
        mpz_init_set(tree[0][i], vals[i]);
    }

    /* Build up */
    for (int lev = 1; lev < levels; lev++) {
        int prev_size = level_size[lev - 1];
        int cur_size = (prev_size + 1) / 2;
        level_size[lev] = cur_size;
        tree[lev] = malloc(sizeof(mpz_t) * cur_size);

        for (int i = 0; i < cur_size; i++) {
            mpz_init(tree[lev][i]);
            if (2 * i + 1 < prev_size) {
                mpz_mul(tree[lev][i], tree[lev-1][2*i], tree[lev-1][2*i+1]);
            } else {
                mpz_set(tree[lev][i], tree[lev-1][2*i]);
            }
        }
    }

    free(level_size);
    return tree;
}

/* Remainder tree: given product tree and a value P, compute P mod each leaf.
 * Works top-down: start with P mod root, then for each node,
 * remainder[child] = remainder[parent] mod child_product_tree_value
 */
static void remainder_tree(mpz_t *remainders, mpz_t **tree, int n,
                          int num_levels, mpz_t P) {
    int *level_size = malloc(sizeof(int) * num_levels);
    level_size[0] = n;
    for (int lev = 1; lev < num_levels; lev++)
        level_size[lev] = (level_size[lev-1] + 1) / 2;

    /* Allocate remainder at each level */
    mpz_t **rem = malloc(sizeof(mpz_t*) * num_levels);
    for (int lev = 0; lev < num_levels; lev++) {
        rem[lev] = malloc(sizeof(mpz_t) * level_size[lev]);
        for (int i = 0; i < level_size[lev]; i++)
            mpz_init(rem[lev][i]);
    }

    /* Top level: P mod root */
    int top = num_levels - 1;
    mpz_mod(rem[top][0], P, tree[top][0]);

    /* Work down */
    for (int lev = top - 1; lev >= 0; lev--) {
        for (int i = 0; i < level_size[lev]; i++) {
            int parent = i / 2;
            mpz_mod(rem[lev][i], rem[lev+1][parent], tree[lev][i]);
        }
    }

    /* Copy leaf remainders to output */
    for (int i = 0; i < n; i++)
        mpz_set(remainders[i], rem[0][i]);

    /* Cleanup */
    for (int lev = 0; lev < num_levels; lev++) {
        for (int i = 0; i < level_size[lev]; i++)
            mpz_clear(rem[lev][i]);
        free(rem[lev]);
    }
    free(rem);
    free(level_size);
}

static void free_product_tree(mpz_t **tree, int n, int num_levels) {
    int *level_size = malloc(sizeof(int) * num_levels);
    level_size[0] = n;
    for (int lev = 1; lev < num_levels; lev++)
        level_size[lev] = (level_size[lev-1] + 1) / 2;

    for (int lev = 0; lev < num_levels; lev++) {
        for (int i = 0; i < level_size[lev]; i++)
            mpz_clear(tree[lev][i]);
        free(tree[lev]);
    }
    free(tree);
    free(level_size);
}

/* Trial divide a value by factor base primes, extract exponents */
static int trial_divide(mpz_t val, unsigned char *exponents,
                       unsigned int *lp1, unsigned int *lp2,
                       unsigned long lp_bound) {
    mpz_t tmp;
    mpz_init_set(tmp, val);

    memset(exponents, 0, fb_count);
    *lp1 = 0;
    *lp2 = 0;

    /* Trial divide by factor base primes */
    for (int i = 1; i < fb_count; i++) {
        unsigned int p = fb[i].p;
        while (mpz_divisible_ui_p(tmp, p)) {
            exponents[i] ^= 1;
            mpz_divexact_ui(tmp, tmp, p);
        }
    }

    /* Check remainder */
    if (mpz_cmp_ui(tmp, 1) == 0) {
        mpz_clear(tmp);
        return 1; /* Fully smooth */
    }

    /* Check for single large prime */
    if (mpz_fits_ulong_p(tmp) && mpz_get_ui(tmp) <= lp_bound) {
        *lp1 = mpz_get_ui(tmp);
        mpz_clear(tmp);
        return 2; /* SLP */
    }

    /* DLP: cofactor might be product of two large primes - skip for now */

    mpz_clear(tmp);
    return 0; /* Not smooth enough */
}

/* Add a relation */
static void add_relation(int x, mpz_t qx, unsigned char *exponents,
                        int sign, unsigned int lp1, unsigned int lp2) {
    if (rel_count >= rel_alloc) {
        rel_alloc = rel_alloc ? rel_alloc * 2 : 1024;
        relations = realloc(relations, sizeof(relation_t) * rel_alloc);
    }
    relation_t *r = &relations[rel_count];
    r->x_val = x;
    r->x_val2 = 0;
    r->is_merged = 0;
    mpz_init_set(r->qx, qx);
    r->exponents = malloc(fb_count);
    memcpy(r->exponents, exponents, fb_count);
    r->sign = sign;
    r->lp1 = lp1;
    r->lp2 = lp2;
    rel_count++;
}

/* Gaussian elimination over GF(2) to find null space vectors */
static int gaussian_elimination(unsigned char **matrix, int rows, int cols,
                               int *pivot_row, unsigned char **history) {
    /* history[i] tracks which rows were XORed to form row i */
    for (int i = 0; i < rows; i++) {
        history[i] = calloc((rows + 7) / 8, 1);
        history[i][i / 8] |= (1 << (i % 8));
    }

    int rank = 0;
    for (int col = 0; col < cols && rank < rows; col++) {
        /* Find pivot */
        int piv = -1;
        for (int row = rank; row < rows; row++) {
            if (matrix[row][col]) { piv = row; break; }
        }
        if (piv < 0) continue;

        pivot_row[rank] = col;

        /* Swap */
        if (piv != rank) {
            unsigned char *tmp = matrix[piv];
            matrix[piv] = matrix[rank];
            matrix[rank] = tmp;
            tmp = history[piv];
            history[piv] = history[rank];
            history[rank] = tmp;
        }

        /* Eliminate */
        for (int row = 0; row < rows; row++) {
            if (row != rank && matrix[row][col]) {
                for (int c = 0; c < cols; c++)
                    matrix[row][c] ^= matrix[rank][c];
                int hbytes = (rows + 7) / 8;
                for (int b = 0; b < hbytes; b++)
                    history[row][b] ^= history[rank][b];
            }
        }
        rank++;
    }

    return rank;
}

/* Try to find factor from null space vector */
static int try_factor(int *dep_indices, int dep_count, mpz_t factor) {
    mpz_t x, y, tmp;
    mpz_init(x);
    mpz_init(y);
    mpz_init(tmp);

    mpz_set_ui(x, 1);

    /* For each relation in dependency, multiply (x_val + sqrtN) into x.
     * For merged SLP pairs (x_val == 0), we need to use the qx product.
     * Actually for merged relations, we stored the product of Q(x) values.
     * The x side needs the product of the corresponding (x_i + sqrtN).
     * But we lost that info. Let me use a different approach:
     * Compute y from the product of |Q(x)| values directly using mpz_sqrt.
     */

    /* Accumulate product of Q(x) and product of (x+sqrtN) */
    mpz_t prod_qx;
    mpz_init_set_ui(prod_qx, 1);

    for (int i = 0; i < dep_count; i++) {
        relation_t *r = &relations[dep_indices[i]];

        /* x *= (x_val + sqrtN) mod N -- for both x values if merged */
        mpz_set(tmp, sqrtN);
        if (r->x_val >= 0)
            mpz_add_ui(tmp, tmp, r->x_val);
        else
            mpz_sub_ui(tmp, tmp, -(r->x_val));
        mpz_mul(x, x, tmp);
        mpz_mod(x, x, N);

        if (r->is_merged && r->x_val2 != 0) {
            mpz_set(tmp, sqrtN);
            if (r->x_val2 >= 0)
                mpz_add_ui(tmp, tmp, r->x_val2);
            else
                mpz_sub_ui(tmp, tmp, -(r->x_val2));
            mpz_mul(x, x, tmp);
            mpz_mod(x, x, N);
        }

        /* Track exponents via actual Q(x) */
        mpz_t absqx;
        mpz_init(absqx);
        mpz_abs(absqx, r->qx);
        mpz_mul(prod_qx, prod_qx, absqx);
        mpz_clear(absqx);
    }

    /* y = sqrt(product of |Q(x)|) -- must be a perfect square */
    mpz_t sq_root, rem;
    mpz_init(sq_root);
    mpz_init(rem);
    mpz_sqrtrem(sq_root, rem, prod_qx);

    if (mpz_sgn(rem) != 0) {
        /* Not a perfect square -- dependency is wrong or SLP cofactors */
        /* Fall back to computing y via individual prime exponents */
        int *total_exp = calloc(fb_count, sizeof(int));

        for (int i = 0; i < dep_count; i++) {
            relation_t *r = &relations[dep_indices[i]];
            mpz_t val;
            mpz_init(val);
            mpz_abs(val, r->qx);

            for (int j = 1; j < fb_count; j++) {
                unsigned int p = fb[j].p;
                while (mpz_divisible_ui_p(val, p)) {
                    total_exp[j]++;
                    mpz_divexact_ui(val, val, p);
                }
            }
            if (mpz_sgn(r->qx) < 0) total_exp[0]++;
            mpz_clear(val);
        }

        /* Check all exponents are even */
        int valid = 1;
        for (int j = 0; j < fb_count; j++) {
            if (total_exp[j] % 2 != 0) { valid = 0; break; }
        }

        if (!valid) {
            free(total_exp);
            mpz_clear(x); mpz_clear(y); mpz_clear(tmp);
            mpz_clear(prod_qx); mpz_clear(sq_root); mpz_clear(rem);
            return 0;
        }

        /* Compute y = product of p^(e/2) mod N */
        mpz_set_ui(y, 1);
        for (int j = 1; j < fb_count; j++) {
            if (total_exp[j] > 0) {
                mpz_set_ui(tmp, fb[j].p);
                mpz_powm_ui(tmp, tmp, total_exp[j] / 2, N);
                mpz_mul(y, y, tmp);
                mpz_mod(y, y, N);
            }
        }
        free(total_exp);
    } else {
        /* Perfect square - use mpz_sqrt result */
        mpz_mod(y, sq_root, N);
    }

    mpz_clear(prod_qx);
    mpz_clear(sq_root);
    mpz_clear(rem);

    /* factor = gcd(x - y, N) */
    mpz_sub(tmp, x, y);
    mpz_gcd(factor, tmp, N);

    int found = (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0);

    if (!found) {
        /* Try x + y */
        mpz_add(tmp, x, y);
        mpz_gcd(factor, tmp, N);
        found = (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0);
    }

    mpz_clear(x); mpz_clear(y); mpz_clear(tmp);
    return found;
}

/* Set parameters based on digit count */
static void set_params(params_t *params, int digits) {
    params->fb_size = fb_count;

    /* Batch size - larger batches amortize tree overhead better */
    /* Product/remainder tree is O(M log^2 M), so moderate batches */
    params->batch_size = fb_count * 6;
    if (params->batch_size < 500) params->batch_size = 500;
    if (params->batch_size > 100000) params->batch_size = 100000;

    /* Sieve radius */
    params->sieve_radius = params->batch_size * 20;

    /* Large prime bound: ~30x largest FB prime */
    params->lp_bits = (int)(log2(fb[fb_count-1].p) + 5);
    if (params->lp_bits > 28) params->lp_bits = 28;

    params->num_extra = 20;

    fprintf(stderr, "Params: fb=%d, batch=%d, lp_bits=%d\n",
            params->fb_size, params->batch_size, params->lp_bits);
}

/*
 * BATCH SMOOTH DETECTION:
 * Given candidates[], compute which are B-smooth using product/remainder trees.
 *
 * For each candidate v:
 * 1. Compute primorial P = prod(p^k for p in FB, where p^k <= v)
 * 2. Compute r = P mod v using remainder tree
 * 3. smooth_part = gcd(r, v) -- this extracts all FB-prime factors
 * 4. Repeatedly: v = v / smooth_part, recompute until stable
 * 5. If remaining cofactor == 1: fully smooth
 * 6. If cofactor < LP_bound: SLP relation
 */
static int batch_smooth_detect(mpz_t *candidates, int *x_vals, int ncand,
                              unsigned long lp_bound, params_t *params) {
    if (ncand == 0) return 0;

    int found = 0;

    /* Compute primorial P = product of p^floor(log_p(2^bits)) for each FB prime */
    mpz_t primorial;
    mpz_init_set_ui(primorial, 1);

    /* Use prime powers up to a bound */
    int bits_bound = mpz_sizeinbase(candidates[0], 2) + 5;

    mpz_t pp;
    mpz_init(pp);
    for (int i = 1; i < fb_count; i++) {
        unsigned int p = fb[i].p;
        /* p^k where p^k <= 2^bits_bound */
        unsigned long pk = p;
        while (pk <= (1UL << (bits_bound < 40 ? bits_bound : 40)) / p) {
            pk *= p;
        }
        mpz_mul_ui(primorial, primorial, pk);
    }

    /* Make absolute values for tree */
    mpz_t *abs_vals = malloc(sizeof(mpz_t) * ncand);
    for (int i = 0; i < ncand; i++) {
        mpz_init(abs_vals[i]);
        mpz_abs(abs_vals[i], candidates[i]);
        if (mpz_cmp_ui(abs_vals[i], 0) == 0)
            mpz_set_ui(abs_vals[i], 1); /* avoid zero */
    }

    /* Build product tree */
    int num_levels;
    mpz_t **ptree = build_product_tree(abs_vals, ncand, &num_levels);

    /* Compute remainders: primorial mod each candidate */
    mpz_t *remainders = malloc(sizeof(mpz_t) * ncand);
    for (int i = 0; i < ncand; i++)
        mpz_init(remainders[i]);

    remainder_tree(remainders, ptree, ncand, num_levels, primorial);

    /* For each candidate, extract smooth part */
    mpz_t smooth_part, cofactor;
    mpz_init(smooth_part);
    mpz_init(cofactor);
    unsigned char *exponents = calloc(fb_count, 1);

    for (int i = 0; i < ncand; i++) {
        /* GCD(remainder, candidate) = smooth part of candidate */
        mpz_gcd(smooth_part, remainders[i], abs_vals[i]);
        mpz_divexact(cofactor, abs_vals[i], smooth_part);

        /* Iteratively extract more: sometimes GCD misses powers */
        int changed = 1;
        while (changed && mpz_cmp_ui(cofactor, 1) > 0) {
            changed = 0;
            mpz_gcd(smooth_part, cofactor, primorial);
            if (mpz_cmp_ui(smooth_part, 1) > 0) {
                mpz_divexact(cofactor, cofactor, smooth_part);
                changed = 1;
            }
        }

        /* Check if smooth or has small cofactor */
        int is_smooth = 0;
        unsigned int lp1_val = 0;

        if (mpz_cmp_ui(cofactor, 1) == 0) {
            is_smooth = 1; /* Fully B-smooth */
        } else if (mpz_fits_ulong_p(cofactor) && mpz_get_ui(cofactor) <= lp_bound) {
            is_smooth = 2; /* SLP */
            lp1_val = mpz_get_ui(cofactor);
        }

        if (is_smooth) {
            /* Full trial division to get exact exponents */
            int sign = (mpz_sgn(candidates[i]) < 0) ? 1 : 0;
            memset(exponents, 0, fb_count);
            exponents[0] = sign;

            mpz_t tdiv_val;
            mpz_init(tdiv_val);
            mpz_abs(tdiv_val, candidates[i]);

            for (int j = 1; j < fb_count; j++) {
                unsigned int p = fb[j].p;
                while (mpz_divisible_ui_p(tdiv_val, p)) {
                    exponents[j] ^= 1;
                    mpz_divexact_ui(tdiv_val, tdiv_val, p);
                }
            }

            /* Verify */
            if (is_smooth == 1 && mpz_cmp_ui(tdiv_val, 1) != 0) {
                /* Batch detection false positive - cofactor not caught */
                mpz_clear(tdiv_val);
                continue;
            }
            if (is_smooth == 2) {
                unsigned long rem = mpz_get_ui(tdiv_val);
                if (rem != lp1_val && rem > lp_bound) {
                    mpz_clear(tdiv_val);
                    continue;
                }
                lp1_val = rem;
            }

            mpz_clear(tdiv_val);

            add_relation(x_vals[i], candidates[i], exponents, sign, lp1_val, 0);
            found++;
        }
    }

    /* Cleanup */
    free(exponents);
    mpz_clear(smooth_part);
    mpz_clear(cofactor);
    for (int i = 0; i < ncand; i++) {
        mpz_clear(remainders[i]);
        mpz_clear(abs_vals[i]);
    }
    free(remainders);
    free(abs_vals);
    free_product_tree(ptree, ncand, num_levels);
    mpz_clear(primorial);
    mpz_clear(pp);

    return found;
}

/* Main factoring routine */
static int factor(mpz_t n, mpz_t result) {
    struct timespec start, now;
    clock_gettime(CLOCK_MONOTONIC, &start);

    mpz_set(N, n);

    /* Check small factors first */
    for (int p = 2; p < 1000; p++) {
        if (mpz_divisible_ui_p(N, p)) {
            mpz_set_ui(result, p);
            return 1;
        }
    }

    int digits = mpz_sizeinbase(N, 10);
    fprintf(stderr, "Factoring %d-digit number\n", digits);

    /* Compute floor(sqrt(N)) */
    mpz_init(sqrtN);
    mpz_sqrt(sqrtN, N);

    /* Build factor base */
    /* Smoothness bound B ≈ exp(0.5 * sqrt(ln N * ln ln N)) */
    int fb_target;
    int smooth_bound;
    {
        double ln_n = digits * log(10);
        double ln_ln_n = log(ln_n);
        smooth_bound = (int)exp(0.5 * sqrt(ln_n * ln_ln_n));
        if (smooth_bound < 500) smooth_bound = 500;
        if (smooth_bound > 5000000) smooth_bound = 5000000;
        /* Count primes up to smooth_bound (about half are QR) */
        fb_target = smooth_bound; /* build_factor_base will use this as prime limit */
        fprintf(stderr, "Smoothness bound B = %d\n", smooth_bound);
    }

    build_factor_base(smooth_bound);
    fprintf(stderr, "Factor base: %d primes, largest = %u\n",
            fb_count, fb[fb_count-1].p);

    /* Set parameters */
    params_t params;
    set_params(&params, digits);

    int target_relations = fb_count + params.num_extra;
    fprintf(stderr, "Need %d relations\n", target_relations);

    unsigned long lp_bound = (unsigned long)fb[fb_count-1].p * 30;
    if (lp_bound > (1UL << 30)) lp_bound = (1UL << 30);

    /* Relation collection via batch smoothness */
    rel_count = 0;
    rel_alloc = 0;
    relations = NULL;

    int batch_size = params.batch_size;
    mpz_t *batch_vals = malloc(sizeof(mpz_t) * batch_size);
    int *batch_x = malloc(sizeof(int) * batch_size);
    for (int i = 0; i < batch_size; i++)
        mpz_init(batch_vals[i]);

    int x_offset = 1;
    int total_tested = 0;
    int batch_num = 0;

    while (rel_count < target_relations) {
        clock_gettime(CLOCK_MONOTONIC, &now);
        double elapsed = (now.tv_sec - start.tv_sec) + (now.tv_nsec - start.tv_nsec) / 1e9;
        if (elapsed > 290.0) {
            fprintf(stderr, "Timeout after %.1fs, %d relations found\n", elapsed, rel_count);
            break;
        }

        /* Fill batch with Q(x) values */
        int bcnt = 0;
        for (int i = 0; i < batch_size && bcnt < batch_size; i++) {
            int x = x_offset + i;
            compute_qx(batch_vals[bcnt], x);
            batch_x[bcnt] = x;
            bcnt++;
        }
        x_offset += batch_size;
        total_tested += bcnt;

        /* Batch smoothness detection */
        int found = batch_smooth_detect(batch_vals, batch_x, bcnt, lp_bound, &params);

        batch_num++;
        if (batch_num % 10 == 0 || found > 0) {
            clock_gettime(CLOCK_MONOTONIC, &now);
            elapsed = (now.tv_sec - start.tv_sec) + (now.tv_nsec - start.tv_nsec) / 1e9;
            fprintf(stderr, "[%.1fs] Batch %d: tested %d, found %d smooth (%d/%d total, %.1f/s)\n",
                    elapsed, batch_num, total_tested, found,
                    rel_count, target_relations,
                    rel_count / elapsed);
        }
    }

    if (rel_count < target_relations) {
        fprintf(stderr, "Not enough relations (%d/%d)\n", rel_count, target_relations);
        /* Cleanup */
        for (int i = 0; i < batch_size; i++) mpz_clear(batch_vals[i]);
        free(batch_vals); free(batch_x);
        return 0;
    }

    fprintf(stderr, "Collected %d relations, starting linear algebra\n", rel_count);

    /* SLP pair merging: find pairs with same large prime */
    /* Sort SLP relations by large prime, then merge pairs */
    int slp_count = 0;
    int *slp_idx = malloc(sizeof(int) * rel_count);
    for (int i = 0; i < rel_count; i++) {
        if (relations[i].lp1 > 0)
            slp_idx[slp_count++] = i;
    }

    /* Simple O(n log n) approach: sort by LP, merge adjacent pairs */
    /* Insertion sort (fine for moderate sizes) */
    for (int i = 1; i < slp_count; i++) {
        int key = slp_idx[i];
        unsigned int kval = relations[key].lp1;
        int j = i - 1;
        while (j >= 0 && relations[slp_idx[j]].lp1 > kval) {
            slp_idx[j+1] = slp_idx[j];
            j--;
        }
        slp_idx[j+1] = key;
    }

    /* Create merged relations from SLP pairs */
    int merged_count = 0;
    typedef struct { int idx1, idx2; } slp_pair_t;
    slp_pair_t *merged = malloc(sizeof(slp_pair_t) * (slp_count / 2 + 1));

    for (int i = 0; i < slp_count - 1; i++) {
        if (relations[slp_idx[i]].lp1 == relations[slp_idx[i+1]].lp1) {
            merged[merged_count].idx1 = slp_idx[i];
            merged[merged_count].idx2 = slp_idx[i+1];
            merged_count++;
            i++; /* skip second in pair */
        }
    }

    /* Count: fully smooth + merged SLP pairs */
    int smooth_count = 0;
    int *smooth_idx = malloc(sizeof(int) * (rel_count + merged_count));
    for (int i = 0; i < rel_count; i++) {
        if (relations[i].lp1 == 0) {
            smooth_idx[smooth_count++] = i;
        }
    }
    int base_smooth = smooth_count;

    /* Add merged pairs as "virtual" relations with combined exponents */
    /* We'll store them as new relation entries */
    for (int i = 0; i < merged_count; i++) {
        int a = merged[i].idx1, b = merged[i].idx2;
        if (rel_count >= rel_alloc) {
            rel_alloc = rel_alloc ? rel_alloc * 2 : 1024;
            relations = realloc(relations, sizeof(relation_t) * rel_alloc);
        }
        relation_t *r = &relations[rel_count];
        r->x_val = relations[a].x_val;
        r->x_val2 = relations[b].x_val;
        r->is_merged = 1;
        mpz_init(r->qx);
        mpz_mul(r->qx, relations[a].qx, relations[b].qx);
        r->exponents = malloc(fb_count);
        for (int j = 0; j < fb_count; j++)
            r->exponents[j] = relations[a].exponents[j] ^ relations[b].exponents[j];
        r->sign = 0;
        r->lp1 = 0;
        r->lp2 = 0;
        smooth_idx[smooth_count++] = rel_count;
        rel_count++;
    }

    fprintf(stderr, "Fully smooth: %d, SLP pairs merged: %d, total usable: %d\n",
            base_smooth, merged_count, smooth_count);
    free(slp_idx);
    free(merged);

    if (smooth_count <= fb_count) {
        fprintf(stderr, "Not enough fully smooth relations (%d need %d)\n", smooth_count, fb_count + 1);
        /* Cleanup */
        for (int i = 0; i < batch_size; i++) mpz_clear(batch_vals[i]);
        free(batch_vals); free(batch_x); free(smooth_idx);
        return 0;
    }

    /* Build matrix from fully smooth relations only */
    int rows = smooth_count;
    int cols = fb_count;
    unsigned char **matrix = malloc(sizeof(unsigned char*) * rows);
    for (int i = 0; i < rows; i++) {
        matrix[i] = malloc(cols);
        memcpy(matrix[i], relations[smooth_idx[i]].exponents, cols);
    }

    int *pivot_row = malloc(sizeof(int) * cols);
    unsigned char **history = malloc(sizeof(unsigned char*) * rows);

    int rank = gaussian_elimination(matrix, rows, cols, pivot_row, history);
    fprintf(stderr, "Matrix rank = %d, null space dim = %d\n", rank, rows - rank);

    /* Try each dependency */
    int *dep_indices = malloc(sizeof(int) * rows);
    int factored = 0;

    for (int try_idx = rank; try_idx < rows && !factored; try_idx++) {
        int dep_count = 0;
        for (int i = 0; i < rows; i++) {
            if (history[try_idx][i / 8] & (1 << (i % 8))) {
                dep_indices[dep_count++] = smooth_idx[i]; /* Map back to original relation index */
            }
        }

        if (dep_count > 0) {
            factored = try_factor(dep_indices, dep_count, result);
            if (factored) {
                fprintf(stderr, "Found factor from dependency %d (size %d)!\n", try_idx - rank, dep_count);
            }
        }
    }

    /* Cleanup */
    for (int i = 0; i < batch_size; i++) mpz_clear(batch_vals[i]);
    free(batch_vals); free(batch_x);
    for (int i = 0; i < rows; i++) {
        free(matrix[i]);
        free(history[i]);
    }
    free(matrix); free(history); free(pivot_row); free(dep_indices); free(smooth_idx);

    return factored;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    mpz_t n, factor_result;
    mpz_init(n);
    mpz_init(N);
    mpz_init(factor_result);

    mpz_set_str(n, argv[1], 10);

    if (mpz_probab_prime_p(n, 25)) {
        gmp_printf("%Zd is prime\n", n);
        mpz_clear(n); mpz_clear(N); mpz_clear(factor_result);
        return 0;
    }

    if (factor(n, factor_result)) {
        mpz_t cofactor;
        mpz_init(cofactor);
        mpz_divexact(cofactor, n, factor_result);

        clock_gettime(CLOCK_MONOTONIC, &end);
        double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

        gmp_printf("%Zd = %Zd * %Zd\n", n, factor_result, cofactor);
        fprintf(stderr, "Time: %.3f seconds\n", elapsed);

        mpz_clear(cofactor);
    } else {
        clock_gettime(CLOCK_MONOTONIC, &end);
        double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
        fprintf(stderr, "FAILED to factor after %.3f seconds\n", elapsed);
        gmp_printf("FAIL %Zd\n", n);
    }

    mpz_clear(n); mpz_clear(N); mpz_clear(factor_result);
    if (fb) free(fb);
    if (primes) free(primes);

    return 0;
}
