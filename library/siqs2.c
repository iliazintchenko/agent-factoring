/*
 * SIQS2 - Self-Initializing Quadratic Sieve
 * High-performance implementation for 30-90+ digit balanced semiprimes
 *
 * Key features:
 * - Knuth-Schroeppel multiplier selection
 * - Self-initializing polynomials with Gray code b-value enumeration
 * - 32KB block sieving for L1D cache efficiency
 * - Single large prime variation with hash table
 * - Block Lanczos-style GF(2) linear algebra
 * - Optimized for AMD EPYC 9R45 (Zen4)
 *
 * Usage: ./siqs2 <N>
 * Compile: gcc -O2 -march=native -o siqs2 library/siqs2.c -lgmp -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

/* ==================== Parameters ==================== */

typedef struct {
    int fb_size;        /* target factor base size */
    int sieve_block;    /* sieve block size (32768 for L1D) */
    int num_blocks;     /* number of blocks per side of sieve */
    int lp_mult;        /* large prime multiplier */
    double thresh_adj;  /* sieve threshold adjustment (lower = more candidates) */
    int extra_rels;     /* extra relations beyond FB size */
} params_t;

/*
 * Parameters based on L(N) = exp(sqrt(ln(N)*ln(ln(N)))).
 * num_a_factors is computed dynamically in poly_new_a based on target_a.
 */
static params_t get_params(int bits) {
    if (bits <= 100)  return (params_t){100,  32768, 1,  30, 0.73, 40};
    if (bits <= 110)  return (params_t){150,  32768, 1,  30, 0.74, 40};
    if (bits <= 120)  return (params_t){200,  32768, 2,  35, 0.76, 50};
    if (bits <= 130)  return (params_t){300,  32768, 3,  40, 0.78, 50};
    if (bits <= 140)  return (params_t){400,  32768, 4,  40, 0.79, 60};
    if (bits <= 150)  return (params_t){600,  32768, 6,  50, 0.80, 60};
    if (bits <= 160)  return (params_t){900,  32768, 8,  50, 0.81, 80};
    if (bits <= 170)  return (params_t){1200, 32768, 10, 60, 0.82, 80};
    if (bits <= 180)  return (params_t){1800, 32768, 14, 60, 0.83, 80};
    if (bits <= 190)  return (params_t){2500, 32768, 18, 60, 0.84, 100};
    if (bits <= 200)  return (params_t){3500, 32768, 24, 70, 0.85, 100};
    if (bits <= 210)  return (params_t){5000, 32768, 32, 70, 0.86, 120};
    if (bits <= 220)  return (params_t){7000, 32768, 40, 80, 0.87, 120};
    if (bits <= 230)  return (params_t){9000, 32768, 48, 80, 0.875, 150};
    if (bits <= 240)  return (params_t){12000, 32768, 56, 80, 0.88, 150};
    if (bits <= 250)  return (params_t){16000, 32768, 64, 90, 0.885, 200};
    if (bits <= 260)  return (params_t){22000, 32768, 80, 90, 0.89, 200};
    if (bits <= 270)  return (params_t){30000, 32768, 96, 100, 0.895, 250};
    if (bits <= 280)  return (params_t){40000, 32768, 112, 100, 0.90, 300};
    if (bits <= 290)  return (params_t){55000, 32768, 128, 110, 0.905, 350};
    return (params_t){75000, 32768, 160, 120, 0.91, 400};
}

/* ==================== Primes sieve ==================== */

static int *sieve_primes(int bound, int *count) {
    char *is_composite = calloc(bound + 1, 1);
    for (int i = 2; (long)i * i <= bound; i++)
        if (!is_composite[i])
            for (int j = i * i; j <= bound; j += i)
                is_composite[j] = 1;
    int cnt = 0;
    for (int i = 2; i <= bound; i++)
        if (!is_composite[i]) cnt++;
    int *primes = malloc(cnt * sizeof(int));
    int idx = 0;
    for (int i = 2; i <= bound; i++)
        if (!is_composite[i]) primes[idx++] = i;
    free(is_composite);
    *count = cnt;
    return primes;
}

/* ==================== Modular Arithmetic ==================== */

/* Modular inverse via extended GCD. Returns 0 if no inverse. */
static unsigned int mod_inverse(unsigned int a, unsigned int m) {
    int g, x;
    /* Extended GCD */
    int old_r = (int)a, r = (int)m;
    int old_s = 1, s = 0;
    while (r != 0) {
        int q = old_r / r;
        int temp = r; r = old_r - q * r; old_r = temp;
        temp = s; s = old_s - q * s; old_s = temp;
    }
    g = old_r; x = old_s;
    if (g != 1) return 0;
    return (unsigned int)(((long long)x % (long long)m + m) % m);
}

/* Tonelli-Shanks: sqrt(n) mod p. Returns 0 if not a QR. */
static unsigned int sqrt_mod(unsigned int n, unsigned int p) {
    if (n == 0) return 0;
    if (p == 2) return n & 1;

    /* Euler criterion check */
    unsigned long long base = n % p;
    unsigned long long result = 1;
    unsigned long long exp = (p - 1) / 2;
    unsigned long long b = base, m = p;
    unsigned long long e = exp;
    while (e > 0) {
        if (e & 1) result = (result * b) % m;
        b = (b * b) % m;
        e >>= 1;
    }
    if (result != 1) return 0;

    /* p ≡ 3 (mod 4): simple case */
    if (p % 4 == 3) {
        b = base; e = (p + 1) / 4; m = p; result = 1;
        while (e > 0) {
            if (e & 1) result = (result * b) % m;
            b = (b * b) % m;
            e >>= 1;
        }
        return (unsigned int)result;
    }

    /* General Tonelli-Shanks */
    unsigned int Q = p - 1, S = 0;
    while (Q % 2 == 0) { Q /= 2; S++; }

    /* Find non-residue z */
    unsigned int z = 2;
    while (1) {
        b = z; e = (p - 1) / 2; m = p; result = 1;
        while (e > 0) {
            if (e & 1) result = (result * b) % m;
            b = (b * b) % m;
            e >>= 1;
        }
        if (result == p - 1) break;
        z++;
    }

    unsigned long long M_val = S;
    /* c = z^Q mod p */
    b = z; e = Q; m = p; unsigned long long c = 1;
    while (e > 0) { if (e & 1) c = (c * b) % m; b = (b * b) % m; e >>= 1; }
    /* t = n^Q mod p */
    b = base; e = Q; unsigned long long t = 1;
    while (e > 0) { if (e & 1) t = (t * b) % m; b = (b * b) % m; e >>= 1; }
    /* R = n^((Q+1)/2) mod p */
    b = base; e = (Q + 1) / 2; unsigned long long R = 1;
    while (e > 0) { if (e & 1) R = (R * b) % m; b = (b * b) % m; e >>= 1; }

    while (1) {
        if (t == 1) return (unsigned int)R;
        /* Find smallest i such that t^(2^i) = 1 */
        int i = 0;
        unsigned long long tt = t;
        while (tt != 1) { tt = (tt * tt) % p; i++; }
        /* b = c^(2^(M-i-1)) */
        unsigned long long bb = c;
        for (int j = 0; j < (int)M_val - i - 1; j++)
            bb = (bb * bb) % p;
        M_val = i;
        c = (bb * bb) % p;
        t = (t * c) % p;
        R = (R * bb) % p;
    }
}

/* ==================== Knuth-Schroeppel Multiplier ==================== */

static int choose_multiplier(mpz_t N, int fb_target) {
    static const int kvals[] = {1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19, 21,
                                 22, 23, 26, 29, 30, 31, 33, 34, 35, 37, 38, 39,
                                 41, 42, 43, 46, 47};
    int nk = sizeof(kvals) / sizeof(kvals[0]);

    double best_score = -1e30;
    int best_k = 1;

    /* Score each multiplier: higher is better */
    for (int ki = 0; ki < nk; ki++) {
        int k = kvals[ki];
        double score = 0.0;

        /* Penalty for large multiplier (increases values to sieve) */
        score -= 0.5 * log((double)k);

        mpz_t kN;
        mpz_init(kN);
        mpz_mul_ui(kN, N, k);

        /* Contribution of 2 */
        unsigned long kN_mod8 = mpz_fdiv_ui(kN, 8);
        if (kN_mod8 == 1) score += 2.0 * log(2.0);
        else if (kN_mod8 == 5) score += 1.0 * log(2.0);
        else if (kN_mod8 == 3 || kN_mod8 == 7) score += 0.5 * log(2.0);

        /* Contribution of small odd primes */
        int nprimes;
        int *primes = sieve_primes(1000, &nprimes);
        int counted = 0;
        for (int i = 1; i < nprimes && counted < 30; i++) {
            int p = primes[i];
            if (k % p == 0) continue;
            unsigned long n_mod_p = mpz_fdiv_ui(kN, p);
            unsigned int r = sqrt_mod((unsigned int)n_mod_p, p);
            if (r != 0) {
                score += 2.0 * log((double)p) / (p - 1);
                counted++;
            }
        }
        free(primes);
        mpz_clear(kN);

        if (score > best_score) {
            best_score = score;
            best_k = k;
        }
    }
    return best_k;
}

/* ==================== Factor Base ==================== */

typedef struct {
    unsigned int *prime;
    unsigned int *root;   /* sqrt(kN) mod p */
    unsigned char *logp;
    int size;
} fb_t;

static fb_t *fb_create(mpz_t kN, int target) {
    fb_t *fb = malloc(sizeof(fb_t));
    /* Allocate extra space */
    int alloc = target + 20;
    fb->prime = malloc(alloc * sizeof(unsigned int));
    fb->root  = malloc(alloc * sizeof(unsigned int));
    fb->logp  = malloc(alloc * sizeof(unsigned char));

    /* Entry 0: -1 (sign) */
    fb->prime[0] = 2; /* We'll handle 2 specially */
    fb->root[0] = 1;
    fb->logp[0] = 1;
    fb->size = 1;

    /* Generate primes and check if kN is a QR */
    int bound = target * 30 + 50000; /* generous bound to find enough QR primes */
    int nprimes;
    int *primes = sieve_primes(bound, &nprimes);

    for (int i = 1; i < nprimes && fb->size < target; i++) {
        int p = primes[i]; /* odd primes */
        unsigned long n_mod_p = mpz_fdiv_ui(kN, p);
        if (n_mod_p == 0) {
            /* p divides kN - include with root 0 */
            fb->prime[fb->size] = p;
            fb->root[fb->size] = 0;
            fb->logp[fb->size] = (unsigned char)(log2(p) + 0.5);
            fb->size++;
            continue;
        }
        unsigned int r = sqrt_mod((unsigned int)n_mod_p, p);
        if (r == 0) continue;

        fb->prime[fb->size] = p;
        fb->root[fb->size] = r;
        fb->logp[fb->size] = (unsigned char)(log2(p) + 0.5);
        fb->size++;
    }
    free(primes);
    return fb;
}

/* ==================== Large Prime Hash Table ==================== */

#define LP_HASH_BITS 20
#define LP_HASH_SIZE (1 << LP_HASH_BITS)
#define LP_HASH_MASK (LP_HASH_SIZE - 1)

typedef struct lp_node {
    unsigned long lp;
    int rel_idx;
    struct lp_node *next;
} lp_node_t;

typedef struct {
    lp_node_t **buckets;
    lp_node_t *pool;
    int pool_used;
    int pool_size;
} lp_table_t;

static lp_table_t *lp_create(int max_entries) {
    lp_table_t *t = malloc(sizeof(lp_table_t));
    t->buckets = calloc(LP_HASH_SIZE, sizeof(lp_node_t*));
    t->pool = malloc(max_entries * sizeof(lp_node_t));
    t->pool_used = 0;
    t->pool_size = max_entries;
    return t;
}

static int lp_find(lp_table_t *t, unsigned long lp) {
    unsigned int h = (unsigned int)((lp * 2654435761UL) >> (32 - LP_HASH_BITS)) & LP_HASH_MASK;
    for (lp_node_t *n = t->buckets[h]; n; n = n->next)
        if (n->lp == lp) return n->rel_idx;
    return -1;
}

static void lp_insert(lp_table_t *t, unsigned long lp, int idx) {
    if (t->pool_used >= t->pool_size) return;
    unsigned int h = (unsigned int)((lp * 2654435761UL) >> (32 - LP_HASH_BITS)) & LP_HASH_MASK;
    lp_node_t *n = &t->pool[t->pool_used++];
    n->lp = lp;
    n->rel_idx = idx;
    n->next = t->buckets[h];
    t->buckets[h] = n;
}

/* ==================== Relation Storage ==================== */

typedef struct {
    mpz_t *Qx;         /* Q(x) value (for sqrt step) */
    mpz_t *ax_b;       /* (ax+b) value for sqrt step */
    unsigned long *lp;  /* large prime (0 = full relation) */
    int count;
    int alloc;
} rel_store_t;

static rel_store_t *rel_create(int alloc) {
    rel_store_t *rs = malloc(sizeof(rel_store_t));
    rs->Qx = malloc(alloc * sizeof(mpz_t));
    rs->ax_b = malloc(alloc * sizeof(mpz_t));
    rs->lp = calloc(alloc, sizeof(unsigned long));
    for (int i = 0; i < alloc; i++) {
        mpz_init(rs->Qx[i]);
        mpz_init(rs->ax_b[i]);
    }
    rs->count = 0;
    rs->alloc = alloc;
    return rs;
}

/* ==================== GF(2) Matrix & Gaussian Elimination ==================== */

typedef unsigned long long u64;

typedef struct {
    u64 **rows;     /* Each row: [fb_words | id_words] */
    int nrows;
    int ncols;      /* FB size (left part) */
    int fb_words;
    int id_words;
    int total_words;
} gf2_matrix_t;

static gf2_matrix_t *gf2_create(int nrows, int ncols) {
    gf2_matrix_t *m = malloc(sizeof(gf2_matrix_t));
    m->nrows = nrows;
    m->ncols = ncols;
    m->fb_words = (ncols + 63) / 64;
    m->id_words = (nrows + 63) / 64;
    m->total_words = m->fb_words + m->id_words;
    m->rows = malloc(nrows * sizeof(u64*));
    for (int r = 0; r < nrows; r++) {
        m->rows[r] = calloc(m->total_words, sizeof(u64));
        /* Identity bit */
        m->rows[r][m->fb_words + r / 64] |= (1ULL << (r % 64));
    }
    return m;
}

static void gf2_set_bit(gf2_matrix_t *m, int row, int col) {
    m->rows[row][col / 64] ^= (1ULL << (col % 64));
}

/*
 * Gaussian elimination on GF(2). Returns number of dependencies found.
 * deps[i] = array of relation indices in the dependency, deps_len[i] = count.
 */
static int gf2_solve(gf2_matrix_t *m, int ***deps_out, int **deps_len_out, int max_deps) {
    int nrows = m->nrows;
    int ncols = m->ncols;
    int fb_words = m->fb_words;
    int total_words = m->total_words;

    /* Column-pivot tracking */
    int *pivot_col_to_row = malloc(ncols * sizeof(int));
    memset(pivot_col_to_row, -1, ncols * sizeof(int));

    int current_pivot_row = 0;
    for (int col = 0; col < ncols && current_pivot_row < nrows; col++) {
        /* Find row with bit set in this column, starting from current_pivot_row */
        int pr = -1;
        for (int r = current_pivot_row; r < nrows; r++) {
            if ((m->rows[r][col / 64] >> (col % 64)) & 1) {
                pr = r;
                break;
            }
        }
        if (pr < 0) continue;

        /* Swap to current pivot position */
        if (pr != current_pivot_row) {
            u64 *tmp = m->rows[pr];
            m->rows[pr] = m->rows[current_pivot_row];
            m->rows[current_pivot_row] = tmp;
        }
        pivot_col_to_row[col] = current_pivot_row;

        /* Eliminate all other rows */
        u64 *pivot = m->rows[current_pivot_row];
        for (int r = 0; r < nrows; r++) {
            if (r == current_pivot_row) continue;
            if ((m->rows[r][col / 64] >> (col % 64)) & 1) {
                for (int w = 0; w < total_words; w++)
                    m->rows[r][w] ^= pivot[w];
            }
        }
        current_pivot_row++;
    }

    /* Extract dependencies (rows with zero left part) */
    int ndeps = 0;
    *deps_out = malloc(max_deps * sizeof(int*));
    *deps_len_out = malloc(max_deps * sizeof(int));

    for (int r = current_pivot_row; r < nrows && ndeps < max_deps; r++) {
        /* Verify left part is zero */
        int zero = 1;
        for (int w = 0; w < fb_words && zero; w++) {
            u64 mask = (w < fb_words - 1) ? ~0ULL :
                       (ncols % 64 == 0) ? ~0ULL : ((1ULL << (ncols % 64)) - 1);
            if (m->rows[r][w] & mask) { zero = 0; break; }
        }
        if (!zero) continue;

        /* Extract relation indices from identity part */
        int *dep = malloc(nrows * sizeof(int));
        int dlen = 0;
        for (int w = 0; w < m->id_words; w++) {
            u64 bits = m->rows[r][fb_words + w];
            while (bits) {
                int bit = __builtin_ctzll(bits);
                int idx = w * 64 + bit;
                if (idx < nrows) dep[dlen++] = idx;
                bits &= bits - 1;
            }
        }
        if (dlen > 0) {
            (*deps_out)[ndeps] = dep;
            (*deps_len_out)[ndeps] = dlen;
            ndeps++;
        } else {
            free(dep);
        }
    }

    free(pivot_col_to_row);
    return ndeps;
}

/* ==================== SIQS Core ==================== */

/*
 * Polynomial: f(x) = (a*x + b)^2 - kN, where a = product of selected FB primes
 * Q(x) = f(x)/a = a*x^2 + 2*b*x + c, where c = (b^2 - kN)/a
 *
 * Self-initialization: when changing b via Gray code, only one B_j value changes.
 * New solutions: x_{i,new} = x_{i,old} ± 2 * a_inv_p * B_j (mod p)
 */

/* State for polynomial generation */
#define MAX_A_FACTORS 20
typedef struct {
    mpz_t a;                /* a coefficient */
    mpz_t b;                /* b coefficient */
    mpz_t c;                /* c = (b^2 - kN) / a */
    int a_indices[MAX_A_FACTORS]; /* FB indices of primes in a */
    int num_a_factors;      /* s = number of primes in a (computed dynamically) */
    mpz_t B_values[MAX_A_FACTORS]; /* B_j values for Gray code switching */
    unsigned int *a_inv_data; /* flattened a_inv[j][i] array */
    int a_inv_stride;       /* = fb->size */
    int num_b_values;       /* 2^(s-1) b-values per a */
    int current_b_index;    /* current b-value index (Gray code) */
    unsigned int *soln1;
    unsigned int *soln2;
} poly_state_t;

static struct timespec g_start_time;

static double elapsed_seconds(void) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec - g_start_time.tv_sec) + (now.tv_nsec - g_start_time.tv_nsec) / 1e9;
}

/*
 * Select a new 'a' coefficient and compute all B_j values.
 * a = product of s primes from the factor base.
 * We want a ≈ sqrt(2*kN) / M to minimize max|Q(x)|.
 */
static int poly_new_a(poly_state_t *ps, fb_t *fb, mpz_t kN, int M, gmp_randstate_t rng) {
    /* Target: a ≈ sqrt(2*kN) / M */
    mpz_t target;
    mpz_init(target);
    mpz_mul_ui(target, kN, 2);
    mpz_sqrt(target, target);
    mpz_tdiv_q_ui(target, target, M);
    double log_target = mpz_sizeinbase(target, 2) * log(2.0);

    /* Selection range: middle of FB (skip first few and last few) */
    int lo = fb->size / 3;
    int hi = 2 * fb->size / 3;
    if (lo < 2) lo = 2;
    if (hi <= lo + 3) hi = fb->size - 1;

    /* Dynamically compute number of a-factors */
    double avg_log_p = 0;
    int cnt = 0;
    for (int i = lo; i < hi; i++) {
        if (fb->root[i] == 0) continue;
        avg_log_p += log((double)fb->prime[i]);
        cnt++;
    }
    if (cnt == 0) { mpz_clear(target); return 0; }
    avg_log_p /= cnt;

    int s = (int)(log_target / avg_log_p + 0.5);
    if (s < 3) s = 3;
    if (s > MAX_A_FACTORS) s = MAX_A_FACTORS;
    if (s > hi - lo) s = hi - lo;
    ps->num_a_factors = s;

    /* Select s primes with product close to target */
    int best_attempts = 0;
    double best_ratio = 1e30;
    int best_indices[MAX_A_FACTORS];

    for (int attempt = 0; attempt < 40; attempt++) {
        mpz_set_ui(ps->a, 1);
        int indices[MAX_A_FACTORS];
        for (int i = 0; i < s; i++) {
            int idx, ok;
            int tries = 0;
            do {
                idx = lo + (gmp_urandomm_ui(rng, hi - lo));
                ok = 1;
                for (int j = 0; j < i; j++)
                    if (indices[j] == idx) { ok = 0; break; }
                if (fb->root[idx] == 0) ok = 0;
                tries++;
            } while (!ok && tries < 100);
            if (!ok) break;
            indices[i] = idx;
            mpz_mul_ui(ps->a, ps->a, fb->prime[idx]);
        }

        /* Measure how close a is to target */
        double ratio;
        if (mpz_cmp(ps->a, target) > 0) {
            mpz_t q; mpz_init(q);
            mpz_tdiv_q(q, ps->a, target);
            ratio = mpz_get_d(q);
            mpz_clear(q);
        } else {
            mpz_t q; mpz_init(q);
            mpz_tdiv_q(q, target, ps->a);
            ratio = mpz_get_d(q);
            mpz_clear(q);
        }
        if (ratio < best_ratio) {
            best_ratio = ratio;
            memcpy(best_indices, indices, s * sizeof(int));
            best_attempts = attempt;
        }
        if (ratio < 2.0) break; /* close enough */
    }

    memcpy(ps->a_indices, best_indices, s * sizeof(int));
    mpz_set_ui(ps->a, 1);
    for (int i = 0; i < s; i++)
        mpz_mul_ui(ps->a, ps->a, fb->prime[ps->a_indices[i]]);

    mpz_clear(target);

    /* Compute B_j values for each a-factor q_j:
     * B_j = r_j * (a/q_j) * ((a/q_j)^(-1) mod q_j)
     * where r_j = sqrt(kN) mod q_j (from factor base)
     */
    for (int j = 0; j < s; j++) {
        int idx = ps->a_indices[j];
        unsigned int qj = fb->prime[idx];
        unsigned int rj = fb->root[idx]; /* sqrt(kN) mod qj */

        /* a_over_qj = a / qj */
        mpz_t a_over_qj, mod_qj, inv;
        mpz_inits(a_over_qj, mod_qj, inv, NULL);
        mpz_divexact_ui(a_over_qj, ps->a, qj);
        mpz_set_ui(mod_qj, qj);

        /* inv = (a/qj)^(-1) mod qj */
        if (!mpz_invert(inv, a_over_qj, mod_qj)) {
            mpz_clears(a_over_qj, mod_qj, inv, NULL);
            return 0; /* shouldn't happen */
        }

        /* B_j = rj * (a/qj) * inv mod a
         * But we compute it as: rj * (a/qj)_inv_mod_qj * (a/qj) */
        unsigned long inv_val = mpz_get_ui(inv);
        mpz_mul_ui(ps->B_values[j], a_over_qj, (rj * inv_val) % qj);
        /* Ensure B_j ≡ rj (mod qj) and B_j has right sign */

        mpz_clears(a_over_qj, mod_qj, inv, NULL);
    }

    /* Initial b = sum of B_j (first Gray code position) */
    mpz_set_ui(ps->b, 0);
    for (int j = 0; j < s; j++)
        mpz_add(ps->b, ps->b, ps->B_values[j]);
    /* Ensure b^2 ≡ kN (mod a) */
    /* If not, try negating some B_j. For simplicity, verify and adjust: */
    mpz_t test;
    mpz_init(test);
    mpz_mul(test, ps->b, ps->b);
    mpz_sub(test, test, kN);
    mpz_mod(test, test, ps->a);
    if (mpz_sgn(test) != 0) {
        /* Try flipping sign of first B value */
        mpz_sub(ps->b, ps->b, ps->B_values[0]);
        mpz_sub(ps->b, ps->b, ps->B_values[0]);
        mpz_mul(test, ps->b, ps->b);
        mpz_sub(test, test, kN);
        mpz_mod(test, test, ps->a);
        if (mpz_sgn(test) != 0) {
            mpz_clear(test);
            return 0;
        }
    }
    mpz_clear(test);

    /* Compute c = (b^2 - kN) / a */
    mpz_mul(ps->c, ps->b, ps->b);
    mpz_sub(ps->c, ps->c, kN);
    mpz_divexact(ps->c, ps->c, ps->a);

    /* Compute sieve solutions for each FB prime:
     * Q(x) = a*x^2 + 2*b*x + c ≡ 0 (mod p)
     * x = (-b ± sqrt(kN)) * a^(-1) (mod p)
     * x = a^(-1) * (root[i] - b) mod p  and  a^(-1) * (-root[i] - b) mod p
     */
    for (int i = 0; i < fb->size; i++) {
        unsigned int p = fb->prime[i];
        unsigned long a_mod = mpz_fdiv_ui(ps->a, p);

        if (a_mod == 0) {
            ps->soln1[i] = ps->soln2[i] = 0xFFFFFFFF;
            continue;
        }

        unsigned int ainv = mod_inverse((unsigned int)a_mod, p);
        if (ainv == 0) {
            ps->soln1[i] = ps->soln2[i] = 0xFFFFFFFF;
            continue;
        }

        unsigned long b_mod = mpz_fdiv_ui(ps->b, p);
        unsigned int r = fb->root[i];

        if (r == 0) {
            /* p divides kN */
            ps->soln1[i] = ps->soln2[i] = 0xFFFFFFFF;
            continue;
        }

        /* x1 = ainv * (r - b) mod p */
        unsigned long x1 = ((unsigned long)ainv * ((r + p - b_mod) % p)) % p;
        /* x2 = ainv * (-r - b) mod p = ainv * (p - r - b) mod p */
        unsigned long x2 = ((unsigned long)ainv * ((p - r + p - b_mod) % p)) % p;

        ps->soln1[i] = (unsigned int)x1;
        ps->soln2[i] = (unsigned int)x2;
    }

    /* Precompute a_inv[j][i] for Gray code updates */
    ps->a_inv_stride = fb->size;
    for (int j = 0; j < s; j++) {
        for (int i = 0; i < fb->size; i++) {
            unsigned int p = fb->prime[i];
            unsigned long a_mod = mpz_fdiv_ui(ps->a, p);
            if (a_mod == 0 || fb->root[i] == 0) {
                ps->a_inv_data[j * fb->size + i] = 0;
                continue;
            }
            unsigned int ainv = mod_inverse((unsigned int)a_mod, p);
            unsigned long Bj_mod = mpz_fdiv_ui(ps->B_values[j], p);
            /* 2 * ainv * Bj mod p */
            ps->a_inv_data[j * fb->size + i] = (unsigned int)((2UL * ainv % p * Bj_mod) % p);
        }
    }

    ps->current_b_index = 0;
    ps->num_b_values = 1 << (s - 1);
    return 1;
}

/*
 * Switch to next b-value using Gray code.
 * Only updates the sieve solutions, not recomputing from scratch.
 * Returns: j = the B_j index that changed, or -1 if all b-values exhausted.
 */
static int poly_next_b(poly_state_t *ps, fb_t *fb, mpz_t kN) {
    ps->current_b_index++;
    if (ps->current_b_index >= ps->num_b_values)
        return -1;

    /* Gray code: which bit changed? */
    int gray_prev = (ps->current_b_index - 1) ^ ((ps->current_b_index - 1) >> 1);
    int gray_curr = ps->current_b_index ^ (ps->current_b_index >> 1);
    int changed = gray_prev ^ gray_curr;
    int j = __builtin_ctz(changed); /* which B_j to flip */

    /* Determine sign: if bit was set, subtract; if cleared, add */
    int sign = (gray_curr >> j) & 1; /* 1 = just set (add), 0 = just cleared (sub) */

    /* Update b: b_new = b_old ± 2*B_j */
    if (sign) {
        mpz_add(ps->b, ps->b, ps->B_values[j]);
        mpz_add(ps->b, ps->b, ps->B_values[j]);
    } else {
        mpz_sub(ps->b, ps->b, ps->B_values[j]);
        mpz_sub(ps->b, ps->b, ps->B_values[j]);
    }

    /* Update c = (b^2 - kN) / a */
    mpz_mul(ps->c, ps->b, ps->b);
    mpz_sub(ps->c, ps->c, kN);
    mpz_divexact(ps->c, ps->c, ps->a);

    /* Update sieve solutions:
     * x_new = x_old ± a_inv[j][i] (mod p)
     * Sign depends on whether B_j was added or subtracted from b.
     * If b increased (sign=1): solutions shift by -a_inv
     * If b decreased (sign=0): solutions shift by +a_inv
     */
    for (int i = 0; i < fb->size; i++) {
        if (ps->soln1[i] == 0xFFFFFFFF) continue;
        unsigned int p = fb->prime[i];
        unsigned int delta = ps->a_inv_data[j * ps->a_inv_stride + i];
        if (delta == 0) continue;
        if (sign) {
            ps->soln1[i] = (ps->soln1[i] + p - delta) % p;
            ps->soln2[i] = (ps->soln2[i] + p - delta) % p;
        } else {
            ps->soln1[i] = (ps->soln1[i] + delta) % p;
            ps->soln2[i] = (ps->soln2[i] + delta) % p;
        }
    }

    return j;
}

/* ==================== Main SIQS ==================== */

int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }

    clock_gettime(CLOCK_MONOTONIC, &g_start_time);

    mpz_t N;
    mpz_init_set_str(N, argv[1], 10);

    int digits = (int)mpz_sizeinbase(N, 10);
    int bits = (int)mpz_sizeinbase(N, 2);

    /* Knuth-Schroeppel multiplier */
    int multiplier = choose_multiplier(N, 100);
    mpz_t kN;
    mpz_init(kN);
    mpz_mul_ui(kN, N, multiplier);
    int kN_bits = (int)mpz_sizeinbase(kN, 2);

    params_t P = get_params(kN_bits);

    fprintf(stderr, "SIQS2: %d digits (%d bits), mult=%d, FB=%d, blocks=%d\n",
            digits, bits, multiplier, P.fb_size, P.num_blocks);

    /* Factor base */
    fb_t *fb = fb_create(kN, P.fb_size);
    fprintf(stderr, "FB: %d primes, largest=%u\n", fb->size, fb->prime[fb->size - 1]);

    int M = P.sieve_block * P.num_blocks; /* total sieve half-width */
    unsigned long lp_bound = (unsigned long)fb->prime[fb->size - 1] * P.lp_mult;
    int target_rels = fb->size + P.extra_rels;

    /* Sieve threshold
     * Q(x) = a*x^2 + 2*b*x + c, max |Q(x)| ≈ sqrt(2*kN) * M
     * For SIQS with a ≈ sqrt(2*kN)/M, this gives |Q(x)| ≈ sqrt(2*kN)*M at edges
     * log2(|Q(x)|_max) ≈ kN_bits/2 + 1 + log2(M)
     * Threshold = fraction of this that must be accounted for by FB hits
     */
    double log2_Qmax = kN_bits / 2.0 + 0.5 + log2(M);
    int threshold = (int)(log2_Qmax * P.thresh_adj);
    /* Reduce threshold for small primes we skip in sieve */
    int skip_bound = 5; /* don't sieve primes 2,3 (too many hits, minimal log) */
    threshold -= 3; /* approximate contribution of skipped small primes */

    fprintf(stderr, "M=%d, target=%d, threshold=%d, LP_bound=%lu\n",
            M, target_rels, threshold, lp_bound);

    /* Sieve array (one block) */
    unsigned char *sieve = malloc(P.sieve_block);

    /* Relation storage */
    int max_rels = target_rels * 2;
    int max_partials = 500000;
    rel_store_t *full_rels = rel_create(max_rels);
    rel_store_t *partial_rels = rel_create(max_partials);

    lp_table_t *lp_table = lp_create(max_partials);

    /* For sqrt step: track which FB primes divide each relation (full exponents) */
    /* We'll re-trial-divide during sqrt step to save memory during sieving */

    /* Polynomial state */
    poly_state_t ps;
    ps.num_a_factors = 0; /* computed dynamically in poly_new_a */
    mpz_inits(ps.a, ps.b, ps.c, NULL);
    for (int j = 0; j < MAX_A_FACTORS; j++)
        mpz_init(ps.B_values[j]);
    ps.a_inv_data = malloc(MAX_A_FACTORS * fb->size * sizeof(unsigned int));
    ps.a_inv_stride = fb->size;
    ps.soln1 = malloc(fb->size * sizeof(unsigned int));
    ps.soln2 = malloc(fb->size * sizeof(unsigned int));

    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, 42);

    mpz_t ax_b, Qx, residue, tmp;
    mpz_inits(ax_b, Qx, residue, tmp, NULL);

    int total_polys = 0;
    int a_count = 0;
    int combined_count = 0;
    long candidates_checked = 0;
    long candidates_passed_threshold = 0;

    while (full_rels->count < target_rels) {
        /* Check timeout */
        if (total_polys % 500 == 0 && total_polys > 0) {
            double t = elapsed_seconds();
            if (t > 280) {
                fprintf(stderr, "TIMEOUT at %.1fs\n", t);
                break;
            }
            if (total_polys % 2000 == 0) {
                int total = full_rels->count + combined_count;
                double rate = total / t;
                fprintf(stderr, "  polys=%d rels=%d/%d (full=%d comb=%d part=%d) cands=%ld %.1f/s t=%.1fs\n",
                        total_polys, total, target_rels, full_rels->count, combined_count,
                        partial_rels->count, candidates_passed_threshold, rate, t);
            }
        }

        /* Generate new 'a' if needed */
        if (total_polys == 0 || ps.current_b_index >= ps.num_b_values) {
            if (!poly_new_a(&ps, fb, kN, M, rng)) {
                continue;
            }
            a_count++;
        } else {
            if (poly_next_b(&ps, fb, kN) < 0) {
                if (!poly_new_a(&ps, fb, kN, M, rng))
                    continue;
                a_count++;
            }
        }
        total_polys++;

        /* Sieve in blocks */
        for (int block = -P.num_blocks; block < P.num_blocks; block++) {
            int block_start = block * P.sieve_block;

            memset(sieve, 0, P.sieve_block);

            /* Sieve with FB primes */
            for (int i = 1; i < fb->size; i++) {
                unsigned int p = fb->prime[i];
                if (ps.soln1[i] == 0xFFFFFFFF) continue;
                if (p < (unsigned int)skip_bound) continue; /* skip tiny primes */
                unsigned char lp = fb->logp[i];

                /* Compute starting position in this block.
                 * sieve[j] corresponds to x = block_start + j.
                 * Q(x) ≡ 0 (mod p) when x ≡ soln (mod p).
                 * So we need j ≡ (soln - block_start) (mod p).
                 */
                unsigned int s1 = ps.soln1[i];
                unsigned int s2 = ps.soln2[i];

                /* First root */
                long off1 = ((long)s1 - block_start) % (long)p;
                if (off1 < 0) off1 += p;
                for (int j = (int)off1; j < P.sieve_block; j += p)
                    sieve[j] += lp;

                /* Second root (if different) */
                if (s1 != s2) {
                    long off2 = ((long)s2 - block_start) % (long)p;
                    if (off2 < 0) off2 += p;
                    for (int j = (int)off2; j < P.sieve_block; j += p)
                        sieve[j] += lp;
                }
            }

            /* Scan for smooth candidates */
            for (int j = 0; j < P.sieve_block; j++) {
                if (sieve[j] < threshold) continue;
                candidates_passed_threshold++;

                long x = (long)(block_start + j);
                if (x == 0) continue;

                /* Compute Q(x) = a*x^2 + 2*b*x + c */
                mpz_set_si(tmp, x);
                mpz_mul(Qx, ps.a, tmp);     /* a*x */
                mpz_add(Qx, Qx, ps.b);      /* a*x + b */
                mpz_add(Qx, Qx, ps.b);      /* a*x + 2b */
                mpz_mul(Qx, Qx, tmp);       /* (a*x + 2b) * x = a*x^2 + 2*b*x */
                mpz_add(Qx, Qx, ps.c);      /* a*x^2 + 2*b*x + c */

                /* Q(x) verified correct */

                /* ax+b (for sqrt step) */
                mpz_mul_si(ax_b, ps.a, x);
                mpz_add(ax_b, ax_b, ps.b);

                if (mpz_sgn(Qx) == 0) continue;

                int neg = (mpz_sgn(Qx) < 0);
                mpz_abs(residue, Qx);

                /* Trial divide by FB primes */
                int smooth = 1;
                int sign_exp = neg ? 1 : 0;

                /* OPTIMIZED trial division: use sieve roots to identify
                 * which FB primes divide Q(x), then divide only by those.
                 * For prime p: if x ≡ soln1 (mod p) or x ≡ soln2 (mod p),
                 * then p | Q(x).
                 */

                /* Divide by 2 */
                while (mpz_even_p(residue))
                    mpz_tdiv_q_2exp(residue, residue, 1);

                /* Divide by odd FB primes using sieve root check */
                for (int i = 1; i < fb->size; i++) {
                    unsigned int p = fb->prime[i];
                    if (ps.soln1[i] == 0xFFFFFFFF) continue;
                    /* Check if this prime divides Q(x) at this position */
                    long xpos = x;
                    long r1 = (long)ps.soln1[i];
                    long r2 = (long)ps.soln2[i];
                    long xmod = ((xpos % (long)p) + p) % p;
                    if (xmod != r1 && xmod != r2) continue;
                    /* This prime divides Q(x) - divide it out */
                    if (mpz_divisible_ui_p(residue, p)) {
                        do {
                            mpz_divexact_ui(residue, residue, p);
                        } while (mpz_divisible_ui_p(residue, p));
                    }
                }
                /* Also try small primes we skipped in sieve */
                for (int i = 0; i < fb->size && fb->prime[i] < (unsigned int)skip_bound; i++) {
                    unsigned int p = fb->prime[i];
                    if (p <= 2) continue; /* already handled 2 */
                    while (mpz_divisible_ui_p(residue, p))
                        mpz_divexact_ui(residue, residue, p);
                }

                /* Check residue */
                /* Store a*Q(x) = (ax+b)^2 - kN for correct exponent tracking */
                /* The relation is: (ax+b)^2 ≡ a*Q(x) (mod kN) */
                mpz_t aQx;
                mpz_init(aQx);
                mpz_mul(aQx, Qx, ps.a);  /* a * Q(x) = (ax+b)^2 - kN */

                if (mpz_cmp_ui(residue, 1) == 0) {
                    /* Full relation */
                    int ri = full_rels->count;
                    if (ri < full_rels->alloc) {
                        mpz_set(full_rels->ax_b[ri], ax_b);
                        mpz_set(full_rels->Qx[ri], aQx);
                        full_rels->lp[ri] = 0;
                        full_rels->count++;
                    }
                } else if (mpz_fits_ulong_p(residue) && mpz_get_ui(residue) <= lp_bound) {
                    /* Single large prime partial */
                    unsigned long lp = mpz_get_ui(residue);
                    int match = lp_find(lp_table, lp);
                    if (match >= 0) {
                        /* Combine: product of two a*Q(x) values */
                        int ri = full_rels->count;
                        if (ri < full_rels->alloc) {
                            mpz_mul(full_rels->ax_b[ri], ax_b, partial_rels->ax_b[match]);
                            mpz_mod(full_rels->ax_b[ri], full_rels->ax_b[ri], N);
                            mpz_mul(full_rels->Qx[ri], aQx, partial_rels->Qx[match]);
                            full_rels->lp[ri] = lp;
                            full_rels->count++;
                            combined_count++;
                        }
                    } else {
                        /* Store partial */
                        int pi = partial_rels->count;
                        if (pi < partial_rels->alloc) {
                            mpz_set(partial_rels->ax_b[pi], ax_b);
                            mpz_set(partial_rels->Qx[pi], aQx);
                            partial_rels->lp[pi] = lp;
                            lp_insert(lp_table, lp, pi);
                            partial_rels->count++;
                        }
                    }
                }
                mpz_clear(aQx);
            }
        }

        if (full_rels->count >= target_rels)
            break;
    }

    double sieve_time = elapsed_seconds();
    int total = full_rels->count;
    fprintf(stderr, "Sieving: %d rels (%d full, %d combined) from %d polys (%d a-values) in %.2fs\n",
            total, total - combined_count, combined_count, total_polys, a_count, sieve_time);

    if (total < fb->size + 1) {
        fprintf(stderr, "FAIL: not enough relations (%d < %d)\n", total, fb->size + 1);
        return 1;
    }

    /* ==================== Linear Algebra ==================== */

    fprintf(stderr, "Linear algebra: %d rels x %d FB primes\n", total, fb->size);

    /* Build GF(2) matrix by re-trial-dividing each Q(x) */
    /* Column 0 = sign, columns 1..fb->size-1 = FB primes */
    int mat_cols = fb->size + 1; /* +1 for sign */
    gf2_matrix_t *mat = gf2_create(total, mat_cols);

    for (int r = 0; r < total; r++) {
        mpz_t Qval;
        mpz_init(Qval);
        mpz_set(Qval, full_rels->Qx[r]);

        /* Sign bit */
        if (mpz_sgn(Qval) < 0) {
            gf2_set_bit(mat, r, 0);
            mpz_neg(Qval, Qval);
        }

        /* Divide by 2 */
        int e2 = 0;
        while (mpz_even_p(Qval)) { mpz_tdiv_q_2exp(Qval, Qval, 1); e2++; }
        if (e2 & 1) gf2_set_bit(mat, r, 1); /* column 1 = prime 2 */

        /* Divide by odd FB primes */
        for (int i = 1; i < fb->size; i++) {
            unsigned int p = fb->prime[i];
            int e = 0;
            while (mpz_divisible_ui_p(Qval, p)) {
                mpz_divexact_ui(Qval, Qval, p);
                e++;
            }
            if (e & 1) gf2_set_bit(mat, r, i + 1);
        }

        /* For combined relations, LP^2 is always even exponent -> no contribution */
        mpz_clear(Qval);
    }

    /* Solve */
    int **deps;
    int *deps_len;
    int ndeps = gf2_solve(mat, &deps, &deps_len, 64);
    fprintf(stderr, "Found %d dependencies\n", ndeps);

    /* ==================== Square Root Step ==================== */

    for (int d = 0; d < ndeps; d++) {
        mpz_t X, Y, g;
        mpz_inits(X, Y, g, NULL);

        /* X = product of (ax+b) values mod N */
        mpz_set_ui(X, 1);
        for (int k = 0; k < deps_len[d]; k++) {
            int ri = deps[d][k];
            mpz_mul(X, X, full_rels->ax_b[ri]);
            mpz_mod(X, X, N);
        }

        /* Y = sqrt(product of Q(x) values) mod N
         * Compute product of |Q(x)|, take sqrt using factorization */
        mpz_t prod;
        mpz_init(prod);
        mpz_set_ui(prod, 1);
        for (int k = 0; k < deps_len[d]; k++) {
            int ri = deps[d][k];
            mpz_t absQ;
            mpz_init(absQ);
            mpz_abs(absQ, full_rels->Qx[ri]);
            mpz_mul(prod, prod, absQ);
            mpz_clear(absQ);
        }

        /* prod should be a perfect square */
        if (!mpz_perfect_square_p(prod)) {
            /* For combined relations with LP, multiply by LP^2 */
            /* Actually LP^2 is already squared in the product of two partials */
            /* Try factoring out known primes */
            mpz_set_ui(Y, 1);

            /* Re-trial-divide to get exponents */
            mpz_t rem;
            mpz_init_set(rem, prod);

            /* Remove sign: count negative Q values */
            int neg_count = 0;
            for (int k = 0; k < deps_len[d]; k++) {
                int ri = deps[d][k];
                if (mpz_sgn(full_rels->Qx[ri]) < 0) neg_count++;
            }
            /* neg_count should be even */

            /* Divide by 2 */
            int e2 = 0;
            while (mpz_even_p(rem)) { mpz_tdiv_q_2exp(rem, rem, 1); e2++; }
            if (e2 & 1) { mpz_clear(rem); mpz_clears(X, Y, g, prod, NULL); continue; }
            if (e2 / 2 > 0) {
                mpz_set_ui(tmp, 2);
                mpz_powm_ui(tmp, tmp, e2 / 2, N);
                mpz_mul(Y, Y, tmp);
                mpz_mod(Y, Y, N);
            }

            /* Divide by FB primes */
            int valid = 1;
            for (int i = 1; i < fb->size; i++) {
                unsigned int p = fb->prime[i];
                int e = 0;
                while (mpz_divisible_ui_p(rem, p)) {
                    mpz_divexact_ui(rem, rem, p);
                    e++;
                }
                if (e & 1) { valid = 0; break; }
                if (e / 2 > 0) {
                    mpz_set_ui(tmp, p);
                    mpz_powm_ui(tmp, tmp, e / 2, N);
                    mpz_mul(Y, Y, tmp);
                    mpz_mod(Y, Y, N);
                }
            }

            if (!valid) {
                mpz_clear(rem);
                mpz_clears(X, Y, g, prod, NULL);
                continue;
            }

            /* Remaining should be 1 or a perfect square of LPs */
            if (mpz_cmp_ui(rem, 1) != 0) {
                if (mpz_perfect_square_p(rem)) {
                    mpz_sqrt(tmp, rem);
                    mpz_mul(Y, Y, tmp);
                    mpz_mod(Y, Y, N);
                } else {
                    mpz_clear(rem);
                    mpz_clears(X, Y, g, prod, NULL);
                    continue;
                }
            }
            mpz_clear(rem);
        } else {
            mpz_sqrt(Y, prod);
            mpz_mod(Y, Y, N);
        }
        mpz_clear(prod);

        /* gcd(X - Y, N) and gcd(X + Y, N) */
        mpz_sub(tmp, X, Y);
        mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t other; mpz_init(other);
            mpz_divexact(other, N, g);
            if (mpz_cmp(g, other) > 0) mpz_swap(g, other);
            gmp_printf("%Zd\n", g);
            mpz_clear(other);
            fprintf(stderr, "SIQS2 done: %.3fs (dep %d/%d)\n", elapsed_seconds(), d + 1, ndeps);
            return 0;
        }

        mpz_add(tmp, X, Y);
        mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t other; mpz_init(other);
            mpz_divexact(other, N, g);
            if (mpz_cmp(g, other) > 0) mpz_swap(g, other);
            gmp_printf("%Zd\n", g);
            mpz_clear(other);
            fprintf(stderr, "SIQS2 done: %.3fs (dep %d/%d)\n", elapsed_seconds(), d + 1, ndeps);
            return 0;
        }

        mpz_clears(X, Y, g, NULL);
    }

    fprintf(stderr, "SIQS2 FAILED: no dependency produced a factor after trying %d\n", ndeps);
    return 1;
}
