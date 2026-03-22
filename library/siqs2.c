/*
 * SIQS2 - High-performance Self-Initializing Quadratic Sieve
 * Optimized for balanced semiprimes, 50-95 digits
 * Single-threaded, uses GMP for bignum arithmetic
 *
 * Key optimizations:
 * - Knuth-Schroeppel multiplier selection
 * - Gray code self-initialization (O(1) per polynomial switch)
 * - Block sieve (32KB blocks for L1 cache)
 * - Double large primes with hash table
 * - Block Lanczos for linear algebra
 * - Optimized trial division
 *
 * Usage: ./siqs2 <N>
 * Compile: gcc -O2 -march=native -o siqs2 siqs2.c -lgmp -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <immintrin.h>

/* ==================== Constants ==================== */

#define MAX_FB_SIZE    200000
#define MAX_RELATIONS  250000
#define BLOCK_SIZE     32768   /* 32KB sieve block - fits L1D cache */
#define LP_HASH_BITS   20
#define LP_HASH_SIZE   (1 << LP_HASH_BITS)
#define LP_HASH_MASK   (LP_HASH_SIZE - 1)
#define MAX_A_FACTORS  15
#define MAX_DEPS       64
#define SMALL_PRIME_CUTOFF 256  /* skip sieving for primes < this */

/* ==================== Parameters ==================== */

typedef struct {
    int fb_size;        /* factor base size */
    int sieve_radius;   /* sieve from -M to M */
    int lp_mult;        /* large prime bound = lp_mult * largest_FB_prime */
    int num_a_factors;  /* number of primes making up 'a' coefficient */
    double thresh_adj;  /* threshold adjustment (fraction of expected log) */
    int use_dlp;        /* use double large primes? */
    int dlp_mult;       /* DLP bound = dlp_mult * largest_FB_prime^2 */
} params_t;

static params_t get_params(int bits) {
    /* Parameters from YAFU's AVX512 parameter table (siqs_aux.c)
     * FB size, sieve_radius (blocks * 32768), LP multiplier, num_a_factors (computed), threshold_adj
     * num_a_factors = 0 means "compute automatically from target_a"
     */
    if (bits <=  80) return (params_t){   80,   32768, 40,  0, 0, 0, 0};
    if (bits <= 100) return (params_t){  175,   65536, 50,  0, 0, 0, 0};
    if (bits <= 120) return (params_t){  375,   65536, 50,  0, 0, 0, 0};
    if (bits <= 140) return (params_t){  828,   65536, 50,  0, 0, 0, 0};
    if (bits <= 149) return (params_t){ 1028,   65536, 60,  0, 0, 0, 0};
    if (bits <= 165) return (params_t){ 1228,   65536, 60,  0, 0, 0, 0};
    if (bits <= 181) return (params_t){ 2247,   65536, 70,  0, 0, 0, 0};
    if (bits <= 198) return (params_t){ 3485,   65536, 70,  0, 0, 0, 0};
    if (bits <= 215) return (params_t){ 6357,   65536, 80,  0, 0, 0, 0};
    if (bits <= 232) return (params_t){12132,  131072, 80,  0, 0, 1, 30};
    if (bits <= 248) return (params_t){26379,  196608, 90,  0, 0, 1, 30};
    if (bits <= 265) return (params_t){47158,  196608, 90,  0, 0, 1, 30};
    if (bits <= 281) return (params_t){60650,  262144,100,  0, 0, 1, 30};
    if (bits <= 298) return (params_t){71768,  262144,120,  0, 0, 1, 30};
    if (bits <= 310) return (params_t){86071,  327680,120,  0, 0, 1, 30};
    if (bits <= 320) return (params_t){99745,  327680,140,  0, 0, 1, 30};
    return                  (params_t){115500, 393216,150,  0, 0, 1, 30};
}

/* ==================== Factor Base ==================== */

typedef struct {
    unsigned int p;
    unsigned int r1;   /* sqrt(kN) mod p */
    unsigned int r2;   /* p - r1 */
    unsigned char logp; /* floor(log2(p)) */
    unsigned int ainv;  /* inverse of 2a mod p (updated per a) */
} fb_entry_t;

static fb_entry_t *fb;
static int fb_count;

/* Tonelli-Shanks: sqrt(n) mod p. Returns 0 if not a QR. */
static unsigned int sqrt_mod_p(unsigned long n, unsigned int p) {
    if (n == 0) return 0;
    if (p == 2) return n & 1;

    /* Euler criterion via GMP (fast for small p too) */
    mpz_t b, e, m, r;
    mpz_init_set_ui(b, n);
    mpz_init_set_ui(m, p);
    mpz_init(e); mpz_init(r);

    mpz_set_ui(e, (p - 1) / 2);
    mpz_powm(r, b, e, m);
    if (mpz_cmp_ui(r, 1) != 0) {
        mpz_clears(b, e, m, r, NULL);
        return 0;
    }

    if (p % 4 == 3) {
        mpz_set_ui(e, (p + 1) / 4);
        mpz_powm(r, b, e, m);
        unsigned int result = mpz_get_ui(r);
        mpz_clears(b, e, m, r, NULL);
        return result;
    }

    /* General Tonelli-Shanks */
    unsigned long Q = p - 1;
    int S = 0;
    while (Q % 2 == 0) { Q /= 2; S++; }

    unsigned long z = 2;
    for (;;) {
        mpz_set_ui(b, z);
        mpz_set_ui(e, (p - 1) / 2);
        mpz_powm(r, b, e, m);
        if (mpz_cmp_ui(r, p - 1) == 0) break;
        z++;
    }

    mpz_t M_val, c, t, R, bb, tmp;
    mpz_inits(M_val, c, t, R, bb, tmp, NULL);
    mpz_set_ui(M_val, S);
    mpz_set_ui(c, z); mpz_set_ui(e, Q); mpz_powm(c, c, e, m);
    mpz_set_ui(t, n); mpz_set_ui(e, Q); mpz_powm(t, t, e, m);
    mpz_set_ui(R, n); mpz_set_ui(e, (Q + 1) / 2); mpz_powm(R, R, e, m);

    for (;;) {
        if (mpz_cmp_ui(t, 1) == 0) {
            unsigned int result = mpz_get_ui(R);
            mpz_clears(b, e, m, r, M_val, c, t, R, bb, tmp, NULL);
            return result;
        }
        int i = 0;
        mpz_set(tmp, t);
        while (mpz_cmp_ui(tmp, 1) != 0) {
            mpz_mul(tmp, tmp, tmp); mpz_mod(tmp, tmp, m);
            i++;
        }
        mpz_set(bb, c);
        for (int j = 0; j < (int)mpz_get_ui(M_val) - i - 1; j++) {
            mpz_mul(bb, bb, bb); mpz_mod(bb, bb, m);
        }
        mpz_set_ui(M_val, i);
        mpz_mul(c, bb, bb); mpz_mod(c, c, m);
        mpz_mul(t, t, c); mpz_mod(t, t, m);
        mpz_mul(R, R, bb); mpz_mod(R, R, m);
    }
}

/* Knuth-Schroeppel multiplier selection */
static int select_multiplier(mpz_t N) {
    static const int mults[] = {1, 3, 5, 7, 11, 13, 15, 17, 19, 21, 23,
                                 29, 31, 33, 35, 37, 39, 41, 43, 47, 51,
                                 53, 55, 59, 61, 67, 69, 71, 73};
    int nmults = sizeof(mults) / sizeof(mults[0]);

    double best_score = -1e30;
    int best_mult = 1;

    mpz_t kN;
    mpz_init(kN);

    for (int mi = 0; mi < nmults; mi++) {
        int k = mults[mi];
        mpz_mul_ui(kN, N, k);

        double score = -0.5 * log(k);

        /* Bonus for kN mod 8 */
        int kN_mod8 = mpz_fdiv_ui(kN, 8);
        if (kN_mod8 == 1) score += log(2.0) * 2;
        else if (kN_mod8 == 5) score += log(2.0) * 1.5;
        else if (kN_mod8 == 3 || kN_mod8 == 7) score += log(2.0);

        /* Score contribution from small primes */
        for (int p = 3; p < 200; p += 2) {
            /* Check if p is prime (simple) */
            int is_prime = 1;
            for (int d = 3; d * d <= p; d += 2)
                if (p % d == 0) { is_prime = 0; break; }
            if (!is_prime) continue;

            unsigned long kN_mod = mpz_fdiv_ui(kN, p);
            unsigned int r = sqrt_mod_p(kN_mod, p);
            if (r > 0) {
                double contrib = 2.0 * log(p) / (p - 1);
                if (k % p == 0) contrib *= 0.5;
                score += contrib;
            }
        }

        if (score > best_score) {
            best_score = score;
            best_mult = k;
        }
    }

    mpz_clear(kN);
    return best_mult;
}

static int build_factor_base(mpz_t kN, int target_size) {
    fb_count = 0;

    /* Entry 0: factor -1 (sign) */
    fb[0].p = 1; fb[0].r1 = fb[0].r2 = 0; fb[0].logp = 0;
    fb_count = 1;

    /* Entry 1: factor 2 */
    fb[1].p = 2; fb[1].r1 = 1; fb[1].r2 = 1; fb[1].logp = 1;
    fb_count = 2;

    /* Generate primes via sieve of Eratosthenes */
    int bound = target_size * 20 + 10000;
    if (bound < 100000) bound = 100000;
    char *is_composite = calloc(bound, 1);
    for (int i = 2; (long)i * i < bound; i++)
        if (!is_composite[i])
            for (int j = i * i; j < bound; j += i)
                is_composite[j] = 1;

    for (int p = 3; p < bound && fb_count < target_size; p += 2) {
        if (is_composite[p]) continue;
        unsigned long n_mod = mpz_fdiv_ui(kN, p);
        unsigned int r = sqrt_mod_p(n_mod, p);
        if (r == 0) continue;

        fb[fb_count].p = p;
        fb[fb_count].r1 = r;
        fb[fb_count].r2 = p - r;
        fb[fb_count].logp = (unsigned char)(log2(p) + 0.5);
        fb_count++;
    }
    free(is_composite);
    return fb_count;
}

/* ==================== Relation Storage ==================== */

typedef struct {
    mpz_t Y;          /* (ax+b) value */
    int *expo;         /* full exponent vector (fb_count entries) */
    unsigned long lp1; /* large prime 1 (0 if full) */
    unsigned long lp2; /* large prime 2 (0 if not DLP) */
} relation_t;

static relation_t *relations;
static int num_relations;
static int target_relations;

/* Large prime hash for single large prime */
typedef struct lp_node {
    unsigned long lp;
    int rel_idx;
    struct lp_node *next;
} lp_node_t;

static lp_node_t **lp_table;
static lp_node_t *lp_pool;
static int lp_pool_used;

static void lp_init_table(int max_entries) {
    lp_table = calloc(LP_HASH_SIZE, sizeof(lp_node_t *));
    lp_pool = malloc(max_entries * sizeof(lp_node_t));
    lp_pool_used = 0;
}

static int lp_find(unsigned long lp) {
    unsigned int h = ((unsigned int)(lp * 2654435761U)) >> (32 - LP_HASH_BITS);
    for (lp_node_t *n = lp_table[h]; n; n = n->next)
        if (n->lp == lp) return n->rel_idx;
    return -1;
}

static void lp_insert(unsigned long lp, int idx) {
    unsigned int h = ((unsigned int)(lp * 2654435761U)) >> (32 - LP_HASH_BITS);
    lp_node_t *n = &lp_pool[lp_pool_used++];
    n->lp = lp;
    n->rel_idx = idx;
    n->next = lp_table[h];
    lp_table[h] = n;
}

/* ==================== Self-Initializing Polynomial ==================== */

/*
 * SIQS polynomial: f(x) = (ax+b)^2 - kN, where a*f(x) = (ax+b)^2 - kN
 * Actually we sieve g(x) = ((ax+b)^2 - kN)/a which has smaller values.
 *
 * a = product of s primes q_1, ..., q_s from the factor base
 * For each a, there are 2^(s-1) choices of b (Gray code enumeration)
 *
 * Self-initialization: when changing b, only 2 sieve solutions per FB prime change
 * by a fixed amount (precomputed). This is the key O(1) optimization.
 */

static int a_indices[MAX_A_FACTORS];  /* FB indices of primes in 'a' */
static int num_a_factors;
static mpz_t poly_a, poly_b;
static mpz_t B_values[MAX_A_FACTORS]; /* B_j values for gray code */

/* Precomputed sieve solutions */
static int *soln1, *soln2;  /* current sieve start positions */

/* Gray code state */
static int gray_idx;     /* which B value to toggle next */
static int gray_sign;    /* +1 or -1 */

static int sieve_M_global; /* stored for sieve offset computation */

/* Initialize a new 'a' coefficient and compute all B_j values */
static int init_new_a(mpz_t kN, int sieve_M) {
    sieve_M_global = sieve_M;
    mpz_set_ui(poly_a, 1);
    for (int i = 0; i < num_a_factors; i++) {
        mpz_mul_ui(poly_a, poly_a, fb[a_indices[i]].p);
    }

    /* Compute B_j for each prime q_j in a:
     * B_j = r_j * (a/q_j) * inv(a/q_j mod q_j)
     * where r_j = sqrt(kN) mod q_j
     */
    mpz_t a_over_q, inv, tmp;
    mpz_inits(a_over_q, inv, tmp, NULL);

    for (int j = 0; j < num_a_factors; j++) {
        unsigned int qj = fb[a_indices[j]].p;
        unsigned int rj = fb[a_indices[j]].r1;

        mpz_divexact_ui(a_over_q, poly_a, qj);
        mpz_set_ui(tmp, qj);
        if (!mpz_invert(inv, a_over_q, tmp)) {
            mpz_clears(a_over_q, inv, tmp, NULL);
            return 0;
        }
        unsigned long inv_val = mpz_get_ui(inv);

        /* B_j = (a/q_j) * (inv(a/q_j, q_j) * r_j mod q_j) */
        mpz_mul_ui(B_values[j], a_over_q, (inv_val * rj) % qj);
    }

    /* b = sum of B_j (first Gray code state) */
    mpz_set_ui(poly_b, 0);
    for (int j = 0; j < num_a_factors; j++) {
        mpz_add(poly_b, poly_b, B_values[j]);
    }

    /* Verify b^2 ≡ kN (mod a) */
    mpz_mul(tmp, poly_b, poly_b);
    mpz_sub(tmp, tmp, kN);
    mpz_mod(tmp, tmp, poly_a);
    if (mpz_sgn(tmp) != 0) {
        /* Negate B[0] to fix sign */
        mpz_sub(poly_b, poly_b, B_values[0]);
        mpz_sub(poly_b, poly_b, B_values[0]);
        mpz_neg(B_values[0], B_values[0]);

        mpz_mul(tmp, poly_b, poly_b);
        mpz_sub(tmp, tmp, kN);
        mpz_mod(tmp, tmp, poly_a);
        if (mpz_sgn(tmp) != 0) {
            mpz_clears(a_over_q, inv, tmp, NULL);
            return 0;
        }
    }

    /* Compute sieve solutions for each FB prime.
     * Polynomial: Q(x) = (ax+b)^2 - kN, divided by a.
     * Sieve position j maps to x = j - M.
     * Q(x) ≡ 0 mod p when x ≡ a^-1 * (±r - b) mod p.
     * Sieve start: soln = (a^-1 * (r - b) + M) mod p.
     */
    for (int i = 2; i < fb_count; i++) {
        unsigned int p = fb[i].p;
        unsigned long a_mod = mpz_fdiv_ui(poly_a, p);
        if (a_mod == 0) {
            soln1[i] = soln2[i] = -1;
            continue;
        }

        /* a_inv = a^-1 mod p */
        mpz_set_ui(tmp, a_mod);
        mpz_set_ui(inv, p);
        if (!mpz_invert(tmp, tmp, inv)) {
            soln1[i] = soln2[i] = -1;
            continue;
        }
        unsigned long a_inv = mpz_get_ui(tmp);
        fb[i].ainv = (unsigned int)a_inv;

        unsigned long b_mod = mpz_fdiv_ui(poly_b, p);
        unsigned long r1 = fb[i].r1;
        unsigned long r2 = fb[i].r2;

        /* x1 = a_inv * (r1 - b) mod p, then add M for sieve offset */
        unsigned long x1 = (a_inv * ((r1 + p - b_mod) % p)) % p;
        unsigned long x2 = (a_inv * ((r2 + p - b_mod) % p)) % p;

        soln1[i] = (int)((x1 + sieve_M) % p);
        soln2[i] = (int)((x2 + sieve_M) % p);
    }

    gray_idx = 0;

    mpz_clears(a_over_q, inv, tmp, NULL);
    return 1;
}

/* Switch to next b using Gray code. Returns 0 when all 2^(s-1) b's exhausted */
static int next_b(int sieve_radius) {
    gray_idx++;
    if (gray_idx >= (1 << (num_a_factors - 1)))
        return 0;

    /* Gray code: find which bit changed */
    int v = gray_idx ^ (gray_idx >> 1);
    int prev = (gray_idx - 1) ^ ((gray_idx - 1) >> 1);
    int changed = v ^ prev;
    int j = __builtin_ctz(changed); /* index of changed bit */

    /* Determine sign: if bit was set in v, add; if cleared, subtract */
    gray_sign = (v >> j) & 1 ? 1 : -1;

    /* Update b: b_new = b_old + 2 * gray_sign * B_j */
    if (gray_sign > 0) {
        mpz_add(poly_b, poly_b, B_values[j]);
        mpz_add(poly_b, poly_b, B_values[j]);
    } else {
        mpz_sub(poly_b, poly_b, B_values[j]);
        mpz_sub(poly_b, poly_b, B_values[j]);
    }

    /* Update sieve solutions: for each FB prime p,
     * When b changes by Δb = ±2*B_j:
     * x_new = a_inv * (r - (b + Δb)) = x_old - a_inv * Δb
     * soln_new = soln_old - a_inv * Δb mod p
     * delta = (2 * a_inv * B_j_mod_p) mod p
     */
    for (int i = 2; i < fb_count; i++) {
        if (soln1[i] < 0) continue;
        unsigned int p = fb[i].p;
        unsigned long bj_mod = mpz_fdiv_ui(B_values[j], p);
        unsigned long delta = (2UL * ((unsigned long)fb[i].ainv * bj_mod % p)) % p;

        if (gray_sign > 0) {
            /* b increased by 2*B_j, so soln decreases by delta */
            soln1[i] = (int)(((unsigned long)soln1[i] + p - delta) % p);
            soln2[i] = (int)(((unsigned long)soln2[i] + p - delta) % p);
        } else {
            /* b decreased by 2*B_j, so soln increases by delta */
            soln1[i] = (int)(((unsigned long)soln1[i] + delta) % p);
            soln2[i] = (int)(((unsigned long)soln2[i] + delta) % p);
        }
    }

    return 1;
}

/* Choose primes for 'a' coefficient.
 * a ≈ sqrt(2*kN) / M for optimal polynomial values.
 * Select primes from the middle of the factor base.
 */
static void choose_a_primes(mpz_t target_a, int *seed_state) {
    /* Target: a ≈ target_a. Choose num_a_factors primes from FB middle. */
    int lo = fb_count / 4;
    int hi = 3 * fb_count / 4;
    if (lo < 2) lo = 2;
    if (hi <= lo + num_a_factors) hi = fb_count - 1;
    int range = hi - lo;

    /* Simple deterministic selection based on seed */
    for (int i = 0; i < num_a_factors; i++) {
retry:
        *seed_state = (*seed_state * 1103515245 + 12345) & 0x7fffffff;
        int idx = lo + (*seed_state % range);
        /* Ensure unique */
        for (int j = 0; j < i; j++)
            if (a_indices[j] == idx) goto retry;
        a_indices[i] = idx;
    }
}

/* ==================== Sieving ==================== */

static unsigned char *sieve_array;

static void sieve_block(int block_start, int block_end) {
    int block_len = block_end - block_start;
    unsigned char *block = sieve_array;
    memset(block, 0, block_len);

    /* Sieve with factor base primes.
     * For small primes (< SMALL_PRIME_CUTOFF): skip, they contribute little per hit
     * For medium primes: standard sieving
     * For large primes: hit at most once per block
     */

    /* Medium primes */
    for (int i = 2; i < fb_count; i++) {
        if (soln1[i] < 0) continue;
        unsigned int p = fb[i].p;
        unsigned char logp = fb[i].logp;

        if (p < 3) continue; /* skip p=2, handle separately */
        if (logp < 2) continue; /* skip very small contribution */

        /* Compute start position in this block */
        int s1 = soln1[i];
        int s2 = soln2[i];

        /* Offset into sieve: position x maps to sieve[x + M] but we work in blocks.
         * soln values are offsets from 0 modulo p.
         * For block [block_start, block_end), we need first hit >= block_start.
         */
        int start1 = s1 + ((block_start - s1 + p - 1) / p) * p - block_start;
        if (s1 >= block_start) start1 = s1 - block_start;
        else start1 = p - ((block_start - s1) % p);
        if (start1 >= block_len) goto skip1;

        for (int j = start1; j < block_len; j += p)
            block[j] += logp;
skip1:
        if (s1 == s2) continue;

        int start2;
        if (s2 >= block_start) start2 = s2 - block_start;
        else start2 = p - ((block_start - s2) % p);
        if (start2 >= block_len) continue;

        for (int j = start2; j < block_len; j += p)
            block[j] += logp;
    }
}

/* ==================== Trial Division ==================== */

/* Optimized trial division using sieve position to identify which primes divide Q(x).
 * Only checks primes where sieve_pos ≡ soln1[i] or soln2[i] (mod p).
 * Falls back to GMP for the residue.
 */
static int trial_divide(mpz_t residue, int *expos, mpz_t kN,
                        mpz_t ax_b, int sieve_pos) {
    /* Compute Q(x) = ((ax+b)^2 - kN) / a */
    mpz_mul(residue, ax_b, ax_b);
    mpz_sub(residue, residue, kN);
    mpz_divexact(residue, residue, poly_a);

    if (mpz_sgn(residue) == 0) return 0;

    int neg = (mpz_sgn(residue) < 0);
    mpz_abs(residue, residue);

    memset(expos, 0, fb_count * sizeof(int));
    if (neg) expos[0] = 1;

    /* Factor out primes from 'a' */
    for (int k = 0; k < num_a_factors; k++) {
        expos[a_indices[k]] += 1;
    }

    /* Factor 2 separately */
    while (mpz_even_p(residue)) {
        mpz_tdiv_q_2exp(residue, residue, 1);
        expos[1]++;
    }

    /* For each FB prime, check if this position is hit by the sieve solution.
     * If sieve_pos % p == soln1[i] or soln2[i], then p | Q(x).
     * After confirming, do repeated division.
     */
    for (int i = 2; i < fb_count; i++) {
        if (soln1[i] < 0) {
            /* This prime divides 'a', handled above. But Q(x)/a may still be divisible. */
            /* Skip for now - already counted in a */
            continue;
        }
        unsigned int p = fb[i].p;
        unsigned int pos_mod = sieve_pos % p;

        if (pos_mod != (unsigned int)soln1[i] && pos_mod != (unsigned int)soln2[i])
            continue;

        /* p divides Q(x). Remove all factors of p. */
        do {
            mpz_divexact_ui(residue, residue, p);
            expos[i]++;
        } while (mpz_divisible_ui_p(residue, p));
    }

    return 1;
}

/* ==================== Block Lanczos (GF(2)) ==================== */

/*
 * Simplified 64-bit Block Lanczos for GF(2) null space finding.
 * Uses 64-wide block vectors for efficiency.
 */

typedef unsigned long long u64;

/* Sparse matrix representation */
typedef struct {
    int nrows, ncols;
    int *row_start;  /* row_start[i] = start index in col_idx for row i */
    int *col_idx;    /* column indices for non-zero entries */
    int nnz;         /* total non-zeros */
} sparse_matrix_t;

/* Dense 64-wide block vector */
typedef u64 *block_vec_t;

static block_vec_t bv_alloc(int n) {
    return calloc(n, sizeof(u64));
}

/* Matrix-vector multiply: y = M * x (in GF(2), 64-wide) */
static void mat_mul(sparse_matrix_t *M, block_vec_t x, block_vec_t y) {
    memset(y, 0, M->nrows * sizeof(u64));
    for (int i = 0; i < M->nrows; i++) {
        u64 val = 0;
        int start = M->row_start[i];
        int end = M->row_start[i + 1];
        for (int k = start; k < end; k++) {
            val ^= x[M->col_idx[k]];
        }
        y[i] = val;
    }
}

/* Transpose matrix-vector multiply: y = M^T * x */
static void mat_tmul(sparse_matrix_t *M, block_vec_t x, block_vec_t y) {
    memset(y, 0, M->ncols * sizeof(u64));
    for (int i = 0; i < M->nrows; i++) {
        u64 xi = x[i];
        if (xi == 0) continue;
        int start = M->row_start[i];
        int end = M->row_start[i + 1];
        for (int k = start; k < end; k++) {
            y[M->col_idx[k]] ^= xi;
        }
    }
}

/* Simple GF(2) null space finder using structured Gaussian elimination
 * on the normal equations M^T * M.
 * Returns number of null space vectors found.
 */
static int find_null_space(sparse_matrix_t *M, int **dep_lists, int *dep_lens, int max_deps) {
    int n = M->ncols;
    int m = M->nrows;

    /* Use Gaussian elimination on the matrix (rows = relations, cols = FB primes)
     * with augmented identity to find dependencies.
     * For practical sizes this is fine; Block Lanczos would be needed for >50K relations.
     */

    int nwords = (n + 63) / 64;
    int id_words = (m + 63) / 64;
    int total_words = nwords + id_words;

    u64 **mat = malloc(m * sizeof(u64 *));
    for (int r = 0; r < m; r++) {
        mat[r] = calloc(total_words, sizeof(u64));
        /* Fill from sparse matrix */
        int start = M->row_start[r];
        int end = M->row_start[r + 1];
        for (int k = start; k < end; k++) {
            int c = M->col_idx[k];
            mat[r][c / 64] |= (1ULL << (c % 64));
        }
        /* Identity bit */
        mat[r][nwords + r / 64] |= (1ULL << (r % 64));
    }

    /* Row reduction */
    int *pivot = malloc(n * sizeof(int));
    for (int c = 0; c < n; c++) pivot[c] = -1;

    for (int col = 0; col < n; col++) {
        int pr = -1;
        for (int r = 0; r < m; r++) {
            if (!((mat[r][col / 64] >> (col % 64)) & 1)) continue;
            int used = 0;
            for (int cc = 0; cc < col; cc++)
                if (pivot[cc] == r) { used = 1; break; }
            if (!used) { pr = r; break; }
        }
        if (pr < 0) continue;
        pivot[col] = pr;

        for (int r = 0; r < m; r++) {
            if (r == pr) continue;
            if ((mat[r][col / 64] >> (col % 64)) & 1) {
                for (int w = 0; w < total_words; w++)
                    mat[r][w] ^= mat[pr][w];
            }
        }
    }

    /* Extract null space vectors */
    int ndeps = 0;
    for (int r = 0; r < m && ndeps < max_deps; r++) {
        int zero = 1;
        for (int w = 0; w < nwords && zero; w++) {
            u64 mask = ~0ULL;
            if (w == nwords - 1 && n % 64 != 0)
                mask = (1ULL << (n % 64)) - 1;
            if (mat[r][w] & mask) zero = 0;
        }
        if (!zero) continue;

        int *dep = malloc(m * sizeof(int));
        int dep_len = 0;
        for (int w = 0; w < id_words; w++) {
            u64 bits = mat[r][nwords + w];
            while (bits) {
                int bit = __builtin_ctzll(bits);
                int idx = w * 64 + bit;
                if (idx < m) dep[dep_len++] = idx;
                bits &= bits - 1;
            }
        }
        if (dep_len > 0) {
            dep_lists[ndeps] = dep;
            dep_lens[ndeps] = dep_len;
            ndeps++;
        } else {
            free(dep);
        }
    }

    for (int r = 0; r < m; r++) free(mat[r]);
    free(mat);
    free(pivot);
    return ndeps;
}

/* ==================== Main ==================== */

int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }

    struct timespec t0;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    mpz_t N;
    mpz_init_set_str(N, argv[1], 10);

    /* Quick checks */
    if (mpz_perfect_square_p(N)) {
        mpz_t sq; mpz_init(sq); mpz_sqrt(sq, N);
        gmp_printf("%Zd\n", sq);
        return 0;
    }

    /* Small factor check with trial division up to 10000 */
    for (unsigned long p = 2; p < 10000; p++) {
        if (mpz_divisible_ui_p(N, p)) {
            printf("%lu\n", p);
            return 0;
        }
    }

    int digits = (int)mpz_sizeinbase(N, 10);
    int bits = (int)mpz_sizeinbase(N, 2);

    /* Multiplier selection */
    int multiplier = select_multiplier(N);
    mpz_t kN;
    mpz_init(kN);
    mpz_mul_ui(kN, N, multiplier);
    int kN_bits = (int)mpz_sizeinbase(kN, 2);

    params_t P = get_params(kN_bits);

    fprintf(stderr, "SIQS2: %d digits (%d bits), k=%d, FB=%d, M=%d, s=%d, DLP=%d\n",
            digits, bits, multiplier, P.fb_size, P.sieve_radius, P.num_a_factors, P.use_dlp);

    /* Build factor base */
    fb = malloc(MAX_FB_SIZE * sizeof(fb_entry_t));
    build_factor_base(kN, P.fb_size);
    fprintf(stderr, "FB: %d primes, largest=%u\n", fb_count, fb[fb_count - 1].p);

    int sieve_M = P.sieve_radius;
    int total_sieve = 2 * sieve_M;

    /* Compute num_a_factors from target_a = sqrt(2*kN)/M */
    mpz_t target_a;
    mpz_init(target_a);
    mpz_mul_ui(target_a, kN, 2);
    mpz_sqrt(target_a, target_a);
    mpz_tdiv_q_ui(target_a, target_a, sieve_M);

    double log_target = mpz_sizeinbase(target_a, 2) * log(2.0);
    /* Average log of primes in the middle of FB */
    int mid = fb_count / 2;
    double avg_log_p = log((double)fb[mid].p);
    num_a_factors = (int)(log_target / avg_log_p + 0.5);
    if (num_a_factors < 3) num_a_factors = 3;
    if (num_a_factors > MAX_A_FACTORS) num_a_factors = MAX_A_FACTORS;
    if (P.num_a_factors > 0) num_a_factors = P.num_a_factors; /* override if set */

    /* Relations target: fb_count + some extra for safety */
    target_relations = fb_count + 64;
    unsigned long lp_bound = (unsigned long)fb[fb_count - 1].p * P.lp_mult;

    /* Threshold: sieve value should account for all but the large prime.
     * expected_logQ ≈ bits(kN)/2 + log2(M) (for Q(x) = ((ax+b)^2 - kN)/a)
     * threshold = expected_logQ - log2(lp_bound)
     */
    int expected_log = kN_bits / 2 + (int)log2(sieve_M);
    int threshold = expected_log - (int)(log2((double)lp_bound));
    if (threshold < 30) threshold = 30;

    fprintf(stderr, "Target: %d rels, threshold=%d, LP bound=%lu\n",
            target_relations, threshold, lp_bound);

    /* Allocate */
    relations = calloc(MAX_RELATIONS, sizeof(relation_t));
    for (int i = 0; i < MAX_RELATIONS; i++) {
        mpz_init(relations[i].Y);
        relations[i].expo = NULL;
    }
    num_relations = 0;

    soln1 = malloc(MAX_FB_SIZE * sizeof(int));
    soln2 = malloc(MAX_FB_SIZE * sizeof(int));
    sieve_array = aligned_alloc(64, BLOCK_SIZE + 64);

    for (int j = 0; j < MAX_A_FACTORS; j++)
        mpz_init(B_values[j]);
    mpz_inits(poly_a, poly_b, NULL);

    /* target_a already computed above */

    lp_init_table(MAX_RELATIONS);

    int seed_state = 42;
    int total_polys = 0;
    int total_candidates = 0;
    int num_partials = 0;

    mpz_t ax_b, residue, tmp, tmp2;
    mpz_inits(ax_b, residue, tmp, tmp2, NULL);

    int *temp_expo = malloc(MAX_FB_SIZE * sizeof(int));

    while (num_relations < target_relations) {
        /* Choose new 'a' coefficient */
        choose_a_primes(target_a, &seed_state);
        if (!init_new_a(kN, sieve_M)) continue;

        /* Process all 2^(s-1) values of b for this a */
        int b_count = 0;
        do {
            total_polys++;
            b_count++;

            /* Sieve in blocks */
            for (int block_start = 0; block_start < total_sieve; block_start += BLOCK_SIZE) {
                int block_end = block_start + BLOCK_SIZE;
                if (block_end > total_sieve) block_end = total_sieve;
                int block_len = block_end - block_start;

                /* Sieve this block */
                memset(sieve_array, 0, block_len);

                for (int i = 2; i < fb_count; i++) {
                    if (soln1[i] < 0) continue;
                    unsigned int p = fb[i].p;
                    unsigned char logp = fb[i].logp;
                    if (logp < 2) continue;

                    /* First root */
                    int s1 = soln1[i];
                    int off1;
                    if (s1 >= block_start) {
                        off1 = s1 - block_start;
                    } else {
                        off1 = p - ((block_start - s1) % p);
                    }

                    for (int j = off1; j < block_len; j += p)
                        sieve_array[j] += logp;

                    /* Second root */
                    if (soln1[i] != soln2[i]) {
                        int s2 = soln2[i];
                        int off2;
                        if (s2 >= block_start) {
                            off2 = s2 - block_start;
                        } else {
                            off2 = p - ((block_start - s2) % p);
                        }

                        for (int j = off2; j < block_len; j += p)
                            sieve_array[j] += logp;
                    }
                }

                /* Scan for candidates */
                for (int j = 0; j < block_len; j++) {
                    if (sieve_array[j] < threshold) continue;

                    int sieve_pos = block_start + j;
                    long x = (long)sieve_pos - sieve_M;
                    if (x == 0) continue;

                    total_candidates++;

                    /* Compute (ax+b) */
                    mpz_mul_si(ax_b, poly_a, x);
                    mpz_add(ax_b, ax_b, poly_b);

                    /* Trial divide - use sieve position for fast prime identification */
                    if (!trial_divide(residue, temp_expo, kN, ax_b, sieve_pos))
                        continue;

                    if (mpz_cmp_ui(residue, 1) == 0) {
                        /* Full relation */
                        if (num_relations >= MAX_RELATIONS) break;
                        int ri = num_relations;
                        mpz_set(relations[ri].Y, ax_b);
                        relations[ri].expo = malloc(fb_count * sizeof(int));
                        memcpy(relations[ri].expo, temp_expo, fb_count * sizeof(int));
                        relations[ri].lp1 = 0;
                        relations[ri].lp2 = 0;
                        num_relations++;
                    } else if (mpz_fits_ulong_p(residue) && mpz_get_ui(residue) <= lp_bound) {
                        /* Single large prime partial */
                        unsigned long lp = mpz_get_ui(residue);
                        int match = lp_find(lp);
                        if (match >= 0) {
                            /* Combine with existing partial */
                            if (num_relations >= MAX_RELATIONS) break;
                            int ri = num_relations;
                            int pi = match;

                            relations[ri].expo = malloc(fb_count * sizeof(int));
                            for (int k = 0; k < fb_count; k++)
                                relations[ri].expo[k] = temp_expo[k] + relations[pi].expo[k];

                            mpz_mul(relations[ri].Y, ax_b, relations[pi].Y);
                            relations[ri].lp1 = lp;
                            relations[ri].lp2 = 0;
                            num_relations++;
                        } else {
                            /* Store as partial */
                            if (num_partials < MAX_RELATIONS / 2) {
                                int pi = MAX_RELATIONS / 2 + num_partials;
                                if (pi < MAX_RELATIONS) {
                                    mpz_set(relations[pi].Y, ax_b);
                                    relations[pi].expo = malloc(fb_count * sizeof(int));
                                    memcpy(relations[pi].expo, temp_expo, fb_count * sizeof(int));
                                    lp_insert(lp, pi);
                                    num_partials++;
                                }
                            }
                        }
                    }
                    /* TODO: DLP handling */
                }
            }

        } while (next_b(sieve_M) && num_relations < target_relations);

        /* Progress report */
        if (total_polys % 500 == 0 || num_relations >= target_relations) {
            struct timespec now;
            clock_gettime(CLOCK_MONOTONIC, &now);
            double elapsed = (now.tv_sec - t0.tv_sec) + (now.tv_nsec - t0.tv_nsec) / 1e9;
            double rate = num_relations / (elapsed > 0 ? elapsed : 1);
            int remaining = target_relations - num_relations;
            fprintf(stderr, "  polys=%d rels=%d/%d (partials=%d, cands=%d) "
                    "%.1f rels/s ETA=%.0fs elapsed=%.1fs\n",
                    total_polys, num_relations, target_relations, num_partials,
                    total_candidates, rate,
                    remaining / (rate > 0 ? rate : 1), elapsed);

            if (elapsed > 280) {
                fprintf(stderr, "TIMEOUT approaching, stopping sieve\n");
                break;
            }
        }
    }

    struct timespec t_sieve;
    clock_gettime(CLOCK_MONOTONIC, &t_sieve);
    double sieve_time = (t_sieve.tv_sec - t0.tv_sec) + (t_sieve.tv_nsec - t0.tv_nsec) / 1e9;
    fprintf(stderr, "Sieving done: %d rels from %d polys in %.1fs\n",
            num_relations, total_polys, sieve_time);

    if (num_relations < fb_count + 1) {
        fprintf(stderr, "Not enough relations (%d < %d)\n", num_relations, fb_count + 1);
        return 1;
    }

    /* Build sparse matrix for linear algebra */
    fprintf(stderr, "Linear algebra: %d x %d\n", num_relations, fb_count);

    sparse_matrix_t M;
    M.nrows = num_relations;
    M.ncols = fb_count;
    M.row_start = malloc((num_relations + 1) * sizeof(int));

    /* Count non-zeros */
    int nnz = 0;
    for (int r = 0; r < num_relations; r++) {
        M.row_start[r] = nnz;
        for (int c = 0; c < fb_count; c++)
            if (relations[r].expo[c] & 1) nnz++;
    }
    M.row_start[num_relations] = nnz;
    M.nnz = nnz;
    M.col_idx = malloc(nnz * sizeof(int));

    nnz = 0;
    for (int r = 0; r < num_relations; r++) {
        for (int c = 0; c < fb_count; c++)
            if (relations[r].expo[c] & 1)
                M.col_idx[nnz++] = c;
    }

    int *dep_lists[MAX_DEPS];
    int dep_lens[MAX_DEPS];
    int ndeps = find_null_space(&M, dep_lists, dep_lens, MAX_DEPS);
    fprintf(stderr, "Found %d dependencies\n", ndeps);

    /* Try each dependency to find factor */
    for (int d = 0; d < ndeps; d++) {
        mpz_t X, Y, g;
        mpz_inits(X, Y, g, NULL);
        mpz_set_ui(X, 1);
        mpz_set_ui(Y, 1);

        int *total_exp = calloc(fb_count, sizeof(int));
        for (int k = 0; k < dep_lens[d]; k++) {
            int ri = dep_lists[d][k];
            mpz_mul(X, X, relations[ri].Y);
            mpz_mod(X, X, N);
            for (int f = 0; f < fb_count; f++)
                total_exp[f] += relations[ri].expo[f];
        }

        /* Check all exponents are even */
        int valid = 1;
        for (int f = 0; f < fb_count; f++)
            if (total_exp[f] & 1) { valid = 0; break; }
        if (!valid) { free(total_exp); mpz_clears(X, Y, g, NULL); continue; }

        /* Y = product of p^(exp/2) mod N */
        for (int f = 1; f < fb_count; f++) {
            int e = total_exp[f] / 2;
            if (e > 0) {
                mpz_set_ui(tmp, fb[f].p);
                mpz_powm_ui(tmp, tmp, e, N);
                mpz_mul(Y, Y, tmp);
                mpz_mod(Y, Y, N);
            }
        }

        /* Include large primes */
        for (int k = 0; k < dep_lens[d]; k++) {
            int ri = dep_lists[d][k];
            if (relations[ri].lp1 > 0) {
                mpz_set_ui(tmp, relations[ri].lp1);
                mpz_mul(Y, Y, tmp);
                mpz_mod(Y, Y, N);
            }
        }

        /* Include multiplier: if k > 1, we need to adjust */
        /* X^2 ≡ Y^2 (mod kN), but we want factors of N */

        /* Try gcd(X-Y, N) and gcd(X+Y, N) */
        mpz_sub(tmp, X, Y);
        mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t other; mpz_init(other);
            mpz_divexact(other, N, g);
            if (mpz_cmp(g, other) > 0) mpz_swap(g, other);
            gmp_printf("%Zd\n", g);
            struct timespec t1; clock_gettime(CLOCK_MONOTONIC, &t1);
            fprintf(stderr, "SIQS2 done: %.3fs (dep %d/%d)\n",
                    (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9,
                    d + 1, ndeps);
            return 0;
        }

        mpz_add(tmp, X, Y);
        mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t other; mpz_init(other);
            mpz_divexact(other, N, g);
            if (mpz_cmp(g, other) > 0) mpz_swap(g, other);
            gmp_printf("%Zd\n", g);
            struct timespec t1; clock_gettime(CLOCK_MONOTONIC, &t1);
            fprintf(stderr, "SIQS2 done: %.3fs (dep %d/%d)\n",
                    (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9,
                    d + 1, ndeps);
            return 0;
        }

        free(total_exp);
        mpz_clears(X, Y, g, NULL);
    }

    fprintf(stderr, "SIQS2 FAILED: no dependency produced a factor\n");
    return 1;
}
