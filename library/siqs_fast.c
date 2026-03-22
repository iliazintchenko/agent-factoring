/*
 * High-performance SIQS (Self-Initializing Quadratic Sieve)
 * For balanced semiprimes, 30-80+ digits.
 *
 * Key optimizations:
 * - AVX2/AVX512 sieve scanning (64 bytes at a time)
 * - 32KB sieve blocks (fits L1 cache)
 * - Gray-code self-initialization for O(1) poly switching
 * - Single Large Prime (SLP) + Double Large Prime (DLP)
 * - Block Lanczos or GF(2) Gaussian elimination
 * - Knuth-Schroeppel multiplier selection
 *
 * Compile: gcc -O3 -march=native -mavx512bw -o siqs_fast library/siqs_fast.c -lgmp -lm
 * Usage:   ./siqs_fast <N>
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <immintrin.h>
#include <gmp.h>

/* ==================== Configuration ==================== */
#define SIEVE_SIZE       32768   /* 32KB = L1 cache */
#define MAX_FB           200000
#define MAX_AFACTORS     16
#define MAX_RELATIONS    600000
#define LP_HASH_BITS     21
#define LP_HASH_SIZE     (1 << LP_HASH_BITS)
#define LP_HASH_MASK     (LP_HASH_SIZE - 1)
#define DLP_HASH_BITS    20
#define DLP_HASH_SIZE    (1 << DLP_HASH_BITS)
#define DLP_HASH_MASK    (DLP_HASH_SIZE - 1)

/* ==================== Parameter Table ==================== */
typedef struct {
    int fb_size;        /* Factor base size */
    int num_blocks;     /* Number of 32KB sieve blocks per side */
    int lp_mult;        /* Large prime multiplier (lp_bound = lp_mult * largest_fb_prime) */
    int num_a_factors;  /* Number of primes in A polynomial */
    int extra_rels;     /* Extra relations over FB size */
} params_t;

static params_t get_params(int digits) {
    /* Tuned for balanced semiprimes */
    if (digits <= 30) return (params_t){150,   1, 40, 3, 60};
    if (digits <= 35) return (params_t){300,   1, 50, 3, 80};
    if (digits <= 40) return (params_t){600,   1, 60, 4, 100};
    if (digits <= 45) return (params_t){1200,  2, 80, 5, 100};
    if (digits <= 50) return (params_t){2200,  3, 100, 6, 120};
    if (digits <= 55) return (params_t){4000,  4, 120, 7, 120};
    if (digits <= 60) return (params_t){7000,  5, 150, 7, 140};
    if (digits <= 65) return (params_t){12000, 8, 200, 8, 150};
    if (digits <= 70) return (params_t){22000, 10, 250, 9, 180};
    if (digits <= 75) return (params_t){38000, 14, 300, 10, 200};
    if (digits <= 80) return (params_t){60000, 18, 350, 11, 220};
    if (digits <= 85) return (params_t){72000, 22, 400, 11, 200};
    if (digits <= 90) return (params_t){100000, 28, 500, 12, 250};
    return (params_t){150000, 36, 600, 13, 300};
}

/* ==================== Global State ==================== */
static mpz_t N_orig, kN;
static gmp_randstate_t g_rng;
static int g_multiplier;
static uint64_t g_lp_bound;

/* Factor base */
static uint32_t *fb_prime;    /* primes */
static uint32_t *fb_sqrtn;   /* sqrt(kN) mod p */
static uint8_t  *fb_logp;    /* ceil(log2(p)) */
static int fb_count;
static int fb_sieve_start;   /* skip tiny primes from sieve (< ~7) */
static int fb_med_start;     /* primes above this use bucket sieve */

/* Polynomial state */
static mpz_t poly_a, poly_b;
static int a_indices[MAX_AFACTORS];
static int num_a_factors;
static mpz_t B_values[MAX_AFACTORS];
static int gray_count, gray_max;

/* Per-prime solution offsets */
static int32_t *soln1, *soln2;  /* roots for current poly */
static uint32_t *fb_ainv;       /* modular inverse of a mod p */
static uint32_t **bi_delta;     /* precomputed 2*ainv*B[j] mod p */

/* Sieve array */
static uint8_t sieve[SIEVE_SIZE] __attribute__((aligned(64)));

/* Bucket sieve for large primes */
typedef struct { uint16_t offset; uint16_t fb_idx; } bucket_entry_t;
static bucket_entry_t **buckets;
static int *bucket_count;
static int *bucket_alloc;
static int num_buckets;

/* ==================== Relations ==================== */
typedef struct {
    mpz_t Y;            /* (ax+b) */
    int *exponent;      /* FB exponent vector (full counts) */
    mpz_t lp_prod;      /* product of large primes (for sqrt step) */
} relation_t;

static relation_t *relations;
static int num_relations;
static int target_relations;

/* SLP hash table */
typedef struct slp_entry {
    uint64_t lp;
    mpz_t Y;
    int *exponent;
    struct slp_entry *next;
} slp_entry_t;
static slp_entry_t **slp_ht;
static int slp_count, slp_combined;

/* DLP hash table */
typedef struct dlp_entry {
    uint64_t lp1, lp2;  /* two large primes */
    mpz_t Y;
    int *exponent;
    int used;
    struct dlp_entry *next;
} dlp_entry_t;
static dlp_entry_t **dlp_ht;
static int dlp_stored, dlp_combined;

/* Pool allocator for entries */
#define POOL_CHUNK 16384
typedef struct pool_chunk {
    void *data;
    struct pool_chunk *next;
} pool_chunk_t;

/* ==================== Math Utilities ==================== */

static uint32_t mod_inverse(uint32_t a, uint32_t m) {
    int64_t g0 = m, g1 = a, u0 = 0, u1 = 1;
    while (g1 != 0) {
        int64_t q = g0 / g1;
        int64_t t = g0 - q * g1; g0 = g1; g1 = t;
        t = u0 - q * u1; u0 = u1; u1 = t;
    }
    if (u0 < 0) u0 += m;
    return (uint32_t)u0;
}

/* Tonelli-Shanks: find sqrt(n) mod p */
static uint32_t tonelli_shanks(uint32_t n, uint32_t p) {
    if (p == 2) return n & 1;
    if (n == 0) return 0;

    /* Find Q, S such that p-1 = Q * 2^S */
    uint32_t Q = p - 1, S = 0;
    while ((Q & 1) == 0) { Q >>= 1; S++; }

    if (S == 1) {
        /* p ≡ 3 mod 4: sqrt = n^((p+1)/4) */
        uint64_t r = 1, base = n, exp = (p + 1) / 4;
        while (exp > 0) {
            if (exp & 1) r = (r * base) % p;
            base = (base * base) % p;
            exp >>= 1;
        }
        return (uint32_t)r;
    }

    /* Find a non-residue z */
    uint32_t z = 2;
    while (1) {
        uint64_t r = 1, base = z, exp = (p - 1) / 2;
        while (exp > 0) {
            if (exp & 1) r = (r * base) % p;
            base = (base * base) % p;
            exp >>= 1;
        }
        if (r == p - 1) break;
        z++;
    }

    uint32_t M = S;
    uint64_t c = 1, base_c = z, exp_c = Q;
    while (exp_c > 0) {
        if (exp_c & 1) c = (c * base_c) % p;
        base_c = (base_c * base_c) % p;
        exp_c >>= 1;
    }

    uint64_t t = 1, base_t = n, exp_t = Q;
    while (exp_t > 0) {
        if (exp_t & 1) t = (t * base_t) % p;
        base_t = (base_t * base_t) % p;
        exp_t >>= 1;
    }

    uint64_t R = 1, base_R = n, exp_R = (Q + 1) / 2;
    while (exp_R > 0) {
        if (exp_R & 1) R = (R * base_R) % p;
        base_R = (base_R * base_R) % p;
        exp_R >>= 1;
    }

    while (1) {
        if (t == 1) return (uint32_t)R;

        uint32_t i = 0;
        uint64_t tmp = t;
        while (tmp != 1) {
            tmp = (tmp * tmp) % p;
            i++;
        }

        uint64_t b = c;
        for (uint32_t j = 0; j < M - i - 1; j++)
            b = (b * b) % p;

        M = i;
        c = (b * b) % p;
        t = (t * c) % p;
        R = (R * b) % p;
    }
}

/* Knuth-Schroeppel multiplier selection */
static int select_multiplier(mpz_t n) {
    static const int mult_candidates[] = {1, 3, 5, 7, 11, 13, 15, 17, 19, 21,
        23, 29, 31, 33, 35, 37, 39, 41, 43, 47, 0};

    double best_score = -1e30;
    int best_mult = 1;

    for (int mi = 0; mult_candidates[mi] != 0; mi++) {
        int k = mult_candidates[mi];
        mpz_t kn;
        mpz_init(kn);
        mpz_mul_ui(kn, n, k);

        double score = -0.5 * log((double)k);
        uint32_t kn_mod8 = mpz_fdiv_ui(kn, 8);

        /* Bonus for kN ≡ 1 mod 8 (2 is in FB twice) */
        if (kn_mod8 == 1) score += 2.0 * log(2.0);
        else if (kn_mod8 == 5) score += log(2.0);

        /* Score based on small primes in factor base */
        static const uint32_t small_primes[] = {3,5,7,11,13,17,19,23,29,31,37,41,43,47,53};
        for (int i = 0; i < 15; i++) {
            uint32_t p = small_primes[i];
            if (k % p == 0) {
                score += log((double)p);
            } else {
                uint32_t kn_mod_p = mpz_fdiv_ui(kn, p);
                /* Euler criterion: kn^((p-1)/2) mod p */
                uint64_t r = 1, b = kn_mod_p, e = (p - 1) / 2;
                while (e > 0) { if (e & 1) r = r * b % p; b = b * b % p; e >>= 1; }
                if (r == 1) score += 2.0 * log((double)p) / (p - 1);
            }
        }

        if (score > best_score) {
            best_score = score;
            best_mult = k;
        }
        mpz_clear(kn);
    }
    return best_mult;
}

/* ==================== Factor Base ==================== */
static int build_factor_base(int target_size) {
    /* Sieve of Eratosthenes up to a reasonable bound */
    int bound = target_size * 20 + 10000;
    if (bound < 100000) bound = 100000;
    uint8_t *is_prime = calloc(bound + 1, 1);
    for (int i = 2; i <= bound; i++) is_prime[i] = 1;
    for (int i = 2; i * i <= bound; i++)
        if (is_prime[i])
            for (int j = i * i; j <= bound; j += i)
                is_prime[j] = 0;

    fb_count = 0;
    /* Always include -1 and 2 */
    fb_prime[fb_count] = 1; /* represents -1 */
    fb_sqrtn[fb_count] = 0;
    fb_logp[fb_count] = 0;
    fb_count++;

    fb_prime[fb_count] = 2;
    fb_sqrtn[fb_count] = mpz_fdiv_ui(kN, 2); /* kN mod 2 */
    fb_logp[fb_count] = 1;
    fb_count++;

    uint32_t kn_mod;
    for (int p = 3; p <= bound && fb_count < target_size; p += 2) {
        if (!is_prime[p]) continue;

        kn_mod = mpz_fdiv_ui(kN, p);
        if (kn_mod == 0) {
            /* p divides kN — include it but it's special */
            fb_prime[fb_count] = p;
            fb_sqrtn[fb_count] = 0;
            fb_logp[fb_count] = (uint8_t)(log2((double)p) + 0.5);
            fb_count++;
            continue;
        }

        /* Check if kN is a QR mod p using Euler criterion */
        uint64_t r = 1, base = kn_mod, exp = (p - 1) / 2;
        while (exp > 0) {
            if (exp & 1) r = r * base % p;
            base = base * base % p;
            exp >>= 1;
        }
        if (r != 1) continue; /* kN is not a QR mod p */

        fb_prime[fb_count] = p;
        fb_sqrtn[fb_count] = tonelli_shanks(kn_mod, p);
        fb_logp[fb_count] = (uint8_t)(log2((double)p) + 0.5);
        fb_count++;
    }

    free(is_prime);

    /* Set sieve start (skip primes < 7 from sieve, account for them in threshold) */
    fb_sieve_start = 2;
    while (fb_sieve_start < fb_count && fb_prime[fb_sieve_start] < 7)
        fb_sieve_start++;

    /* Medium prime boundary: primes above this go to bucket sieve */
    fb_med_start = fb_count;
    for (int i = fb_sieve_start; i < fb_count; i++) {
        if (fb_prime[i] > SIEVE_SIZE) {
            fb_med_start = i;
            break;
        }
    }

    return fb_count;
}

/* ==================== Polynomial Management ==================== */

static void compute_poly_b(void) {
    mpz_set_ui(poly_b, 0);
    for (int j = 0; j < num_a_factors; j++) {
        mpz_init_set_ui(B_values[j], 0);
        uint32_t qj = fb_prime[a_indices[j]];
        uint32_t tj = fb_sqrtn[a_indices[j]];

        mpz_t Aj, inv_Aj;
        mpz_inits(Aj, inv_Aj, NULL);
        mpz_divexact_ui(Aj, poly_a, qj);

        mpz_t mod_q;
        mpz_init_set_ui(mod_q, qj);
        if (!mpz_invert(inv_Aj, Aj, mod_q)) {
            mpz_clears(Aj, inv_Aj, mod_q, NULL);
            continue;
        }

        uint64_t gamma = ((uint64_t)tj * mpz_get_ui(inv_Aj)) % qj;
        if (gamma > qj / 2) gamma = qj - gamma;

        mpz_mul_ui(B_values[j], Aj, (uint32_t)gamma);
        mpz_add(poly_b, poly_b, B_values[j]);
        mpz_clears(Aj, inv_Aj, mod_q, NULL);
    }

    /* Verify b^2 ≡ kN mod a */
    mpz_t test;
    mpz_init(test);
    mpz_mul(test, poly_b, poly_b);
    mpz_sub(test, test, kN);
    mpz_mod(test, test, poly_a);
    if (mpz_sgn(test) != 0) {
        mpz_sub(poly_b, poly_a, poly_b);
    }
    mpz_clear(test);

    gray_count = 0;
    gray_max = 1 << (num_a_factors - 1);
}

static void compute_ainv(void) {
    for (int i = 2; i < fb_count; i++) {
        uint32_t p = fb_prime[i];
        uint32_t am = mpz_fdiv_ui(poly_a, p);
        fb_ainv[i] = (am == 0) ? 0 : mod_inverse(am, p);
    }
}

static void compute_roots(int M) {
    for (int i = 2; i < fb_count; i++) {
        uint32_t p = fb_prime[i];
        uint32_t ainv = fb_ainv[i];
        if (ainv == 0) { soln1[i] = soln2[i] = -1; continue; }

        uint32_t bm = mpz_fdiv_ui(poly_b, p);
        uint32_t r1 = fb_sqrtn[i];
        uint32_t r2 = p - r1;

        int64_t x1 = ((int64_t)ainv * ((int64_t)((r1 + p - bm) % p))) % p;
        int64_t x2 = ((int64_t)ainv * ((int64_t)((r2 + p - bm) % p))) % p;

        soln1[i] = (int32_t)(((x1 % p) + M) % p);
        soln2[i] = (int32_t)(((x2 % p) + M) % p);
    }
}

static void precompute_bi_deltas(void) {
    for (int j = 0; j < num_a_factors; j++) {
        for (int i = 2; i < fb_count; i++) {
            uint32_t p = fb_prime[i];
            uint32_t ainv = fb_ainv[i];
            if (ainv == 0) { bi_delta[j][i] = 0; continue; }
            uint32_t bm = mpz_fdiv_ui(B_values[j], p);
            bi_delta[j][i] = (uint32_t)((2ULL * ainv * (uint64_t)bm) % p);
        }
    }
}

static int new_poly_a(params_t *par) {
    num_a_factors = par->num_a_factors;
    int M = par->num_blocks * SIEVE_SIZE;

    mpz_t target;
    mpz_init(target);
    mpz_mul_ui(target, kN, 2);
    mpz_sqrt(target, target);
    mpz_tdiv_q_ui(target, target, M);

    double td = mpz_get_d(target);
    double ideal_prime = pow(td, 1.0 / num_a_factors);

    int center = 2;
    for (int i = 2; i < fb_count; i++) {
        if (fb_prime[i] > ideal_prime) { center = i; break; }
    }

    int spread = num_a_factors * 6;
    if (spread < 50) spread = 50;
    int lo = center - spread, hi = center + spread;
    if (lo < 2) lo = 2;
    if (hi >= fb_count) hi = fb_count - 1;
    if (hi <= lo + num_a_factors) hi = lo + num_a_factors + 30;

    /* Try random A selections until we get one close to target */
    for (int att = 0; att < 500; att++) {
        mpz_set_ui(poly_a, 1);
        for (int i = 0; i < num_a_factors; i++) {
            int idx, dup;
            do {
                idx = lo + (int)gmp_urandomm_ui(g_rng, hi - lo);
                dup = 0;
                for (int j = 0; j < i; j++)
                    if (a_indices[j] == idx) { dup = 1; break; }
            } while (dup);
            a_indices[i] = idx;
            mpz_mul_ui(poly_a, poly_a, fb_prime[idx]);
        }
        double ratio = mpz_get_d(poly_a) / td;
        if (ratio > 0.7 && ratio < 1.5) {
            mpz_clear(target);
            return 1;
        }
    }
    mpz_clear(target);
    return 1; /* best effort */
}

static int next_poly_b(int M) {
    gray_count++;
    if (gray_count >= gray_max) return 0;

    int bit = __builtin_ctz(gray_count);
    int neg = (gray_count >> (bit + 1)) & 1;

    if (neg)
        mpz_submul_ui(poly_b, B_values[bit], 2);
    else
        mpz_addmul_ui(poly_b, B_values[bit], 2);

    /* Update roots using precomputed deltas */
    for (int i = 2; i < fb_count; i++) {
        if (soln1[i] < 0) continue;
        uint32_t p = fb_prime[i];
        uint32_t d = bi_delta[bit][i];
        if (d == 0) continue;

        uint32_t s1 = (uint32_t)soln1[i];
        uint32_t s2 = (uint32_t)soln2[i];

        if (neg) {
            s1 += d; if (s1 >= p) s1 -= p;
            s2 += d; if (s2 >= p) s2 -= p;
        } else {
            s1 = (s1 >= d) ? s1 - d : s1 + p - d;
            s2 = (s2 >= d) ? s2 - d : s2 + p - d;
        }
        soln1[i] = (int32_t)s1;
        soln2[i] = (int32_t)s2;
    }
    return 1;
}

/* ==================== Bucket Sieve ==================== */

static void init_buckets(int M) {
    num_buckets = (2 * M + SIEVE_SIZE - 1) / SIEVE_SIZE;
    buckets = malloc(num_buckets * sizeof(bucket_entry_t *));
    bucket_count = calloc(num_buckets, sizeof(int));
    bucket_alloc = malloc(num_buckets * sizeof(int));
    for (int i = 0; i < num_buckets; i++) {
        bucket_alloc[i] = 4096;
        buckets[i] = malloc(4096 * sizeof(bucket_entry_t));
    }
}

static inline void bucket_add(int bn, uint16_t off, uint16_t fi) {
    if (bucket_count[bn] >= bucket_alloc[bn]) {
        bucket_alloc[bn] *= 2;
        buckets[bn] = realloc(buckets[bn], bucket_alloc[bn] * sizeof(bucket_entry_t));
    }
    buckets[bn][bucket_count[bn]++] = (bucket_entry_t){off, fi};
}

static void fill_buckets(int M) {
    for (int i = 0; i < num_buckets; i++) bucket_count[i] = 0;

    uint32_t total = (uint32_t)(2 * M);
    for (int i = fb_med_start; i < fb_count; i++) {
        uint32_t p = fb_prime[i];
        if (soln1[i] < 0) continue;
        uint16_t fi = (uint16_t)(i < 65536 ? i : 65535);

        for (uint32_t pos = (uint32_t)soln1[i]; pos < total; pos += p)
            bucket_add(pos / SIEVE_SIZE, (uint16_t)(pos % SIEVE_SIZE), fi);
        if (soln1[i] != soln2[i])
            for (uint32_t pos = (uint32_t)soln2[i]; pos < total; pos += p)
                bucket_add(pos / SIEVE_SIZE, (uint16_t)(pos % SIEVE_SIZE), fi);
    }
}

/* ==================== SIEVE KERNEL ==================== */

/* Compute threshold correction for unsieved tiny primes */
static int tiny_prime_correction(void) {
    double c = 0;
    for (int i = 2; i < fb_sieve_start; i++)
        c += fb_logp[i] * 2.0 / fb_prime[i];
    return (int)(c + 0.5);
}

/* Sieve one 32KB block with medium primes */
static void sieve_block_medium(int block_start, int block_end) {
    int blen = block_end - block_start;
    memset(sieve, 0, SIEVE_SIZE);

    for (int i = fb_sieve_start; i < fb_med_start; i++) {
        uint32_t p = fb_prime[i];
        uint8_t lp = fb_logp[i];
        if (soln1[i] < 0) continue;

        /* Compute starting positions within this block */
        int32_t s1 = soln1[i] - block_start;
        int32_t s2 = soln2[i] - block_start;

        /* Adjust to be within block */
        if (s1 < 0) { int k = (-s1 + p - 1) / p; s1 += k * p; }
        if (s2 < 0) { int k = (-s2 + p - 1) / p; s2 += k * p; }

        uint32_t r1 = (uint32_t)s1, r2 = (uint32_t)s2;
        uint32_t end = (uint32_t)blen;

        /* Interleave two roots for small primes (improves ILP) */
        if (p < 256) {
            if (r1 > r2) { uint32_t t = r1; r1 = r2; r2 = t; }
            while (r2 < end) {
                sieve[r1] += lp;
                sieve[r2] += lp;
                r1 += p;
                r2 += p;
            }
            if (r1 < end) sieve[r1] += lp;
        } else {
            for (uint32_t j = r1; j < end; j += p) sieve[j] += lp;
            if (r1 != r2)
                for (uint32_t j = r2; j < end; j += p) sieve[j] += lp;
        }
    }
}

/* Apply bucket sieve hits for large primes */
static void apply_bucket_hits(int blk) {
    int cnt = bucket_count[blk];
    bucket_entry_t *be = buckets[blk];
    for (int i = 0; i < cnt; i++)
        sieve[be[i].offset] += fb_logp[be[i].fb_idx];
}

/* Scan sieve for candidates using AVX512 (or AVX2 fallback) */
static int scan_sieve(uint8_t threshold, int *candidates) {
    int nc = 0;

#ifdef __AVX512BW__
    __m512i vthresh = _mm512_set1_epi8((char)(threshold - 1));
    for (int i = 0; i < SIEVE_SIZE; i += 64) {
        __m512i v = _mm512_load_si512((__m512i *)(sieve + i));
        uint64_t mask = _mm512_cmpgt_epu8_mask(v, vthresh);
        while (mask) {
            int bit = __builtin_ctzll(mask);
            candidates[nc++] = i + bit;
            mask &= mask - 1;
        }
    }
#elif defined(__AVX2__)
    __m256i vthresh = _mm256_set1_epi8((char)threshold);
    for (int i = 0; i < SIEVE_SIZE; i += 32) {
        __m256i v = _mm256_load_si256((__m256i *)(sieve + i));
        /* Compare >= threshold: subtract threshold, check high bit */
        __m256i diff = _mm256_subs_epu8(v, vthresh);
        int mask = _mm256_movemask_epi8(_mm256_cmpgt_epi8(diff, _mm256_setzero_si256()));
        /* Also include exact matches */
        mask |= _mm256_movemask_epi8(_mm256_cmpeq_epi8(v, vthresh));
        while (mask) {
            int bit = __builtin_ctz(mask);
            candidates[nc++] = i + bit;
            mask &= mask - 1;
        }
    }
#else
    for (int i = 0; i < SIEVE_SIZE; i++)
        if (sieve[i] >= threshold) candidates[nc++] = i;
#endif

    return nc;
}

/* ==================== Trial Division ==================== */

/* Try to factor Q(x) over the factor base.
 * Returns 1 if smooth (or has acceptable large prime residue). */
static int trial_divide(mpz_t Qval, int *exponents, uint64_t *lp1_out, uint64_t *lp2_out) {
    mpz_t q;
    mpz_init_set(q, Qval);

    memset(exponents, 0, fb_count * sizeof(int));

    /* Handle sign */
    if (mpz_sgn(q) < 0) {
        exponents[0] = 1;  /* -1 factor */
        mpz_neg(q, q);
    }

    /* Divide by factor base primes */
    for (int i = 1; i < fb_count; i++) {
        uint32_t p = fb_prime[i];
        while (mpz_divisible_ui_p(q, p)) {
            mpz_divexact_ui(q, q, p);
            exponents[i]++;
        }
    }

    *lp1_out = 0;
    *lp2_out = 0;

    if (mpz_cmp_ui(q, 1) == 0) {
        mpz_clear(q);
        return 1;  /* fully smooth */
    }

    uint64_t residue = 0;
    if (mpz_fits_ulong_p(q)) {
        residue = mpz_get_ui(q);
    } else if (mpz_sizeinbase(q, 2) <= 64) {
        residue = mpz_get_ui(q);
    }

    uint64_t lp_bound = g_lp_bound;

    if (residue > 0 && residue <= lp_bound) {
        /* Single large prime */
        if (mpz_probab_prime_p(q, 1)) {
            *lp1_out = residue;
            mpz_clear(q);
            return 2;
        }
    }

    /* Check for DLP: residue = lp1 * lp2 where both <= lp_bound */
    if (mpz_sizeinbase(q, 2) <= 52 && !mpz_probab_prime_p(q, 1)) {
        uint64_t qval = 0;
        if (mpz_fits_ulong_p(q)) {
            qval = mpz_get_ui(q);
            /* Quick trial division first */
            uint64_t d;
            for (d = 2; d * d <= qval && d < 50000; d++) {
                if (qval % d == 0) {
                    uint64_t cofactor = qval / d;
                    if (d <= lp_bound && cofactor <= lp_bound) {
                        *lp1_out = (d < cofactor) ? d : cofactor;
                        *lp2_out = (d < cofactor) ? cofactor : d;
                        mpz_clear(q);
                        return 3;  /* DLP */
                    }
                    break;
                }
            }
            /* If trial division didn't find it, try simple Rho */
            if (d * d > qval) {
                /* qval is prime, can't split */
            } else if (qval > 0) {
                /* Rho on qval */
                uint64_t x = 2, y = 2, c = 1, g = 1;
                for (int iter = 0; iter < 100000 && g == 1; iter++) {
                    x = ((__int128)x * x + c) % qval;
                    y = ((__int128)y * y + c) % qval;
                    y = ((__int128)y * y + c) % qval;
                    uint64_t diff = (x > y) ? x - y : y - x;
                    /* Batch GCD every 32 iters */
                    mpz_t gv;
                    mpz_init(gv);
                    mpz_set_ui(gv, diff);
                    mpz_t qv;
                    mpz_init_set_ui(qv, qval);
                    mpz_gcd(gv, gv, qv);
                    g = mpz_get_ui(gv);
                    mpz_clears(gv, qv, NULL);
                }
                if (g > 1 && g < qval) {
                    uint64_t cofactor = qval / g;
                    if (g <= lp_bound && cofactor <= lp_bound) {
                        *lp1_out = (g < cofactor) ? g : cofactor;
                        *lp2_out = (g < cofactor) ? cofactor : g;
                        mpz_clear(q);
                        return 3;
                    }
                }
            }
        }
    }

    mpz_clear(q);
    return 0;  /* not smooth */
}

/* ==================== Relation Management ==================== */

static void add_full_relation(mpz_t Y, int *exp, mpz_t lp_prod) {
    if (num_relations >= MAX_RELATIONS) return;
    relation_t *r = &relations[num_relations];
    mpz_init_set(r->Y, Y);
    r->exponent = malloc(fb_count * sizeof(int));
    memcpy(r->exponent, exp, fb_count * sizeof(int));
    mpz_init_set(r->lp_prod, lp_prod);
    num_relations++;
}

static slp_entry_t *find_slp(uint64_t lp) {
    uint32_t h = (uint32_t)((lp * 2654435761ULL) >> 12) & LP_HASH_MASK;
    for (slp_entry_t *e = slp_ht[h]; e; e = e->next)
        if (e->lp == lp) return e;
    return NULL;
}

static void process_slp(mpz_t Y, int *exp, uint64_t lp) {
    slp_entry_t *existing = find_slp(lp);
    if (existing) {
        /* Combine two SLPs into a full relation */
        int *combined = malloc(fb_count * sizeof(int));
        for (int i = 0; i < fb_count; i++)
            combined[i] = exp[i] + existing->exponent[i];
        /* lp^2 is already a perfect square, so the combined exponents
         * have even parity for all FB primes (plus the lp^2) */

        mpz_t combined_Y, lp_val;
        mpz_inits(combined_Y, lp_val, NULL);
        mpz_mul(combined_Y, Y, existing->Y);
        mpz_set_ui(lp_val, lp);
        /* Check that combined relation is non-trivial (has odd exponents) */
        int has_odd = 0;
        for (int k = 0; k < fb_count; k++) {
            if (combined[k] % 2 != 0) { has_odd = 1; break; }
        }
        if (!has_odd) {
            /* Trivial combination (same Q value) — try gcd directly instead */
            mpz_t gval;
            mpz_init(gval);
            mpz_sub(gval, Y, existing->Y);
            mpz_gcd(gval, gval, N_orig);
            if (mpz_cmp_ui(gval, 1) > 0 && mpz_cmp(gval, N_orig) < 0) {
                mpz_t cof;
                mpz_init(cof);
                mpz_divexact(cof, N_orig, gval);
                if (mpz_cmp(gval, cof) > 0) mpz_swap(gval, cof);
                gmp_printf("%Zd %Zd\n", gval, cof);
                mpz_clears(cof, gval, combined_Y, lp_val, NULL);
                free(combined);
                exit(0);
            }
            mpz_clear(gval);
            /* Don't add trivial relation */
            mpz_clears(combined_Y, lp_val, NULL);
            free(combined);
            slp_combined++; /* count but don't add */
            return;
        }
        add_full_relation(combined_Y, combined, lp_val);
        mpz_clears(combined_Y, lp_val, NULL);
        free(combined);
        slp_combined++;
    } else {
        /* Debug: check stored SLP exponents */
        if (0) {
            int nz = 0;
            for (int k = 0; k < fb_count; k++) if (exp[k]) nz++;
            fprintf(stderr, "  SLP store: %d nonzero exponents, lp=%lu\n", nz, lp);
        }
        /* Store new SLP */
        slp_entry_t *e = malloc(sizeof(slp_entry_t));
        e->lp = lp;
        mpz_init_set(e->Y, Y);
        e->exponent = malloc(fb_count * sizeof(int));
        memcpy(e->exponent, exp, fb_count * sizeof(int));
        uint32_t h = (uint32_t)((lp * 2654435761ULL) >> 12) & LP_HASH_MASK;
        e->next = slp_ht[h];
        slp_ht[h] = e;
        slp_count++;
    }
}

static void process_dlp(mpz_t Y, int *exp, uint64_t lp1, uint64_t lp2) {
    slp_entry_t *s1 = find_slp(lp1);
    slp_entry_t *s2 = find_slp(lp2);

    if (s1 && s2) {
        /* Triangle: combine all three */
        int *combined = malloc(fb_count * sizeof(int));
        for (int i = 0; i < fb_count; i++)
            combined[i] = exp[i] + s1->exponent[i] + s2->exponent[i];

        mpz_t combined_Y, lp_val;
        mpz_inits(combined_Y, lp_val, NULL);
        mpz_mul(combined_Y, Y, s1->Y);
        mpz_mul(combined_Y, combined_Y, s2->Y);
        mpz_set_ui(lp_val, lp1);
        mpz_mul_ui(lp_val, lp_val, lp2);
        add_full_relation(combined_Y, combined, lp_val);
        mpz_clears(combined_Y, lp_val, NULL);
        free(combined);
        dlp_combined++;
        return;
    }

    /* Store for later */
    uint32_t h1 = (uint32_t)((lp1 * 2654435761ULL) >> 12) & DLP_HASH_MASK;
    dlp_entry_t *e = malloc(sizeof(dlp_entry_t));
    e->lp1 = lp1;
    e->lp2 = lp2;
    mpz_init_set(e->Y, Y);
    e->exponent = malloc(fb_count * sizeof(int));
    memcpy(e->exponent, exp, fb_count * sizeof(int));
    e->used = 0;
    e->next = dlp_ht[h1];
    dlp_ht[h1] = e;
    dlp_stored++;
}

/* ==================== Linear Algebra (GF(2) Gaussian Elimination) ==================== */

/* For small matrices, use simple GF(2) Gaussian elimination with 64-bit words */
static int solve_matrix(int *dep_rows, int max_deps) {
    int nrows = num_relations;
    int ncols = fb_count;
    int nwords = (ncols + 63) / 64;

    /* Build matrix: each row has fb_count columns + identity columns for tracking */
    int total_words = nwords + (nrows + 63) / 64;
    uint64_t **matrix = malloc(nrows * sizeof(uint64_t *));
    for (int i = 0; i < nrows; i++) {
        matrix[i] = calloc(total_words, sizeof(uint64_t));
        /* Set exponent columns (parity of full exponent counts) */
        for (int j = 0; j < ncols; j++) {
            if (relations[i].exponent[j] % 2 != 0)
                matrix[i][j / 64] |= (1ULL << (j % 64));
        }
        /* Set identity column for this row */
        matrix[i][nwords + i / 64] |= (1ULL << (i % 64));
    }

    /* Debug: count nonzero bits per row */
    {
        int zero_rows = 0;
        for (int r = 0; r < nrows; r++) {
            int bits = 0;
            for (int w = 0; w < nwords; w++)
                bits += __builtin_popcountll(matrix[r][w]);
            if (bits == 0) zero_rows++;
        }
        fprintf(stderr, "LA debug: %d/%d rows have all-zero exponent parity before GE\n", zero_rows, nrows);
    }

    /* Forward elimination with proper pivot tracking */
    int *is_pivot = calloc(nrows, sizeof(int));  /* 1 if row is used as pivot */
    int cur_pivot = 0;

    for (int col = 0; col < ncols; col++) {
        /* Find a non-pivot row with a 1 in this column */
        int found = -1;
        for (int row = cur_pivot; row < nrows; row++) {
            if (matrix[row][col / 64] & (1ULL << (col % 64))) {
                found = row;
                break;
            }
        }
        if (found < 0) continue;

        /* Swap found row to cur_pivot position */
        if (found != cur_pivot) {
            uint64_t *tmp = matrix[cur_pivot];
            matrix[cur_pivot] = matrix[found];
            matrix[found] = tmp;
        }
        is_pivot[cur_pivot] = 1;

        /* Eliminate this column from all other rows */
        for (int row = 0; row < nrows; row++) {
            if (row == cur_pivot) continue;
            if (matrix[row][col / 64] & (1ULL << (col % 64))) {
                for (int w = 0; w < total_words; w++)
                    matrix[row][w] ^= matrix[cur_pivot][w];
            }
        }
        cur_pivot++;
    }

    /* Non-pivot rows (from cur_pivot to nrows-1) are in the null space */
    int ndeps = 0;
    for (int row = cur_pivot; row < nrows && ndeps < max_deps; row++) {
        /* Verify exponent part is all zero */
        int all_zero = 1;
        for (int w = 0; w < nwords; w++) {
            if (matrix[row][w]) { all_zero = 0; break; }
        }
        if (!all_zero) continue;

        /* Extract dependency: which original rows are combined */
        int cnt = 0;
        for (int i = 0; i < nrows; i++) {
            if (matrix[row][nwords + i / 64] & (1ULL << (i % 64))) {
                dep_rows[ndeps * (nrows + 1) + cnt] = i;
                cnt++;
            }
        }
        dep_rows[ndeps * (nrows + 1) + cnt] = -1;
        ndeps++;
    }
    free(is_pivot);

    for (int i = 0; i < nrows; i++) free(matrix[i]);
    free(matrix);
    /* pivot_row removed — now using is_pivot + cur_pivot */

    return ndeps;
}

/* ==================== Square Root Step ==================== */

static int extract_factor(int *dep_indices, mpz_t factor) {
    mpz_t X, Y;
    mpz_inits(X, Y, NULL);
    mpz_set_ui(X, 1);
    mpz_set_ui(Y, 1);

    /* Sum full exponent vectors across all relations in dependency */
    int *total_exp = calloc(fb_count, sizeof(int));

    for (int i = 0; dep_indices[i] >= 0; i++) {
        int ri = dep_indices[i];
        /* X = product of Y values mod N */
        mpz_mul(X, X, relations[ri].Y);
        mpz_mod(X, X, N_orig);

        /* Sum exponents */
        for (int j = 0; j < fb_count; j++)
            total_exp[j] += relations[ri].exponent[j];

        /* Include large prime products in Y */
        mpz_mul(Y, Y, relations[ri].lp_prod);
        mpz_mod(Y, Y, N_orig);
    }

    /* All exponents should be even */
    int odd_count = 0;
    for (int j = 0; j < fb_count; j++) {
        if (total_exp[j] % 2 != 0) {
            if (odd_count < 5)
                fprintf(stderr, "  DEBUG: exp[%d] (p=%u) = %d (odd!)\n", j, fb_prime[j], total_exp[j]);
            odd_count++;
        }
    }
    if (odd_count > 0) {
        fprintf(stderr, "  DEBUG: %d odd exponents, skipping dependency\n", odd_count);
        free(total_exp);
        mpz_clears(X, Y, NULL);
        return 0;
    }

    /* Y = product of p^(total_exp[j]/2) * product of large_prime products (already in Y) */
    /* Handle -1: if total_exp[0] is odd, X and Y differ by sign — but we checked it's even */
    for (int j = 1; j < fb_count; j++) {
        int half_exp = total_exp[j] / 2;
        if (half_exp == 0) continue;
        mpz_t pk;
        mpz_init(pk);
        mpz_ui_pow_ui(pk, fb_prime[j], half_exp);
        mpz_mul(Y, Y, pk);
        mpz_mod(Y, Y, N_orig);
        mpz_clear(pk);
    }

    /* Verify: X^2 ≡ Y^2 mod N */
    {
        mpz_t x2, y2;
        mpz_inits(x2, y2, NULL);
        mpz_mul(x2, X, X);
        mpz_mod(x2, x2, N_orig);
        mpz_mul(y2, Y, Y);
        mpz_mod(y2, y2, N_orig);
        if (mpz_cmp(x2, y2) != 0) {
            fprintf(stderr, "  VERIFY FAIL: X^2 != Y^2 mod N\n");
            mpz_clears(x2, y2, NULL);
            free(total_exp);
            mpz_clears(X, Y, NULL);
            return 0;
        }
        mpz_clears(x2, y2, NULL);
    }

    /* Debug: print X and Y for first few deps */
    static int dep_dbg = 0;
    if (0) {
        gmp_fprintf(stderr, "  X = %Zd\n  Y = %Zd\n  N = %Zd\n", X, Y, N_orig);
        mpz_t xmy;
        mpz_init(xmy);
        mpz_sub(xmy, X, Y);
        mpz_mod(xmy, xmy, N_orig);
        gmp_fprintf(stderr, "  X-Y mod N = %Zd\n", xmy);
        mpz_add(xmy, X, Y);
        mpz_mod(xmy, xmy, N_orig);
        gmp_fprintf(stderr, "  X+Y mod N = %Zd\n", xmy);
        mpz_clear(xmy);
        dep_dbg++;
    }

    /* Now X^2 ≡ Y^2 mod N. Try gcd(X-Y, N) and gcd(X+Y, N). */
    mpz_t diff;
    mpz_init(diff);

    mpz_sub(diff, X, Y);
    mpz_gcd(factor, diff, N_orig);
    if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N_orig) < 0) {
        mpz_clear(diff);
        free(total_exp);
        mpz_clears(X, Y, NULL);
        return 1;
    }

    mpz_add(diff, X, Y);
    mpz_gcd(factor, diff, N_orig);
    int success = (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N_orig) < 0);

    mpz_clear(diff);
    free(total_exp);
    mpz_clears(X, Y, NULL);
    return success;
}

/* ==================== Main SIQS ==================== */

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    struct timespec t_start, t_now;
    clock_gettime(CLOCK_MONOTONIC, &t_start);

    mpz_init(N_orig);
    mpz_set_str(N_orig, argv[1], 10);

    int digits = mpz_sizeinbase(N_orig, 10);
    int bits = mpz_sizeinbase(N_orig, 2);

    /* Init RNG with seed 42 */
    gmp_randinit_default(g_rng);
    gmp_randseed_ui(g_rng, 42);

    /* Multiplier disabled for now — causes trivial congruences in sqrt step */
    g_multiplier = 1;
    mpz_init(kN);
    mpz_mul_ui(kN, N_orig, g_multiplier);

    /* Get parameters */
    params_t par = get_params(digits);

    /* Allocate arrays */
    fb_prime = calloc(MAX_FB, sizeof(uint32_t));
    fb_sqrtn = calloc(MAX_FB, sizeof(uint32_t));
    fb_logp = calloc(MAX_FB, sizeof(uint8_t));
    soln1 = calloc(MAX_FB, sizeof(int32_t));
    soln2 = calloc(MAX_FB, sizeof(int32_t));
    fb_ainv = calloc(MAX_FB, sizeof(uint32_t));
    bi_delta = malloc(MAX_AFACTORS * sizeof(uint32_t *));
    for (int j = 0; j < MAX_AFACTORS; j++)
        bi_delta[j] = calloc(MAX_FB, sizeof(uint32_t));

    mpz_init(poly_a);
    mpz_init(poly_b);
    for (int j = 0; j < MAX_AFACTORS; j++)
        mpz_init(B_values[j]);

    /* Build factor base */
    build_factor_base(par.fb_size);

    fprintf(stderr, "SIQS: %dd (%d bits), mult=%d, FB=%d, blocks=%d, s=%d\n",
            digits, bits, g_multiplier, fb_count, par.num_blocks, par.num_a_factors);
    fprintf(stderr, "  FB: %d primes, largest=%u, med_start=%d (prime=%u)\n",
            fb_count, fb_prime[fb_count - 1], fb_med_start,
            fb_med_start < fb_count ? fb_prime[fb_med_start] : 0);

    /* Allocate relations */
    relations = calloc(MAX_RELATIONS, sizeof(relation_t));
    num_relations = 0;
    target_relations = fb_count + par.extra_rels;
    g_lp_bound = (uint64_t)fb_prime[fb_count - 1] * par.lp_mult;

    /* Allocate hash tables */
    slp_ht = calloc(LP_HASH_SIZE, sizeof(slp_entry_t *));
    dlp_ht = calloc(DLP_HASH_SIZE, sizeof(dlp_entry_t *));

    /* Init buckets */
    int M = par.num_blocks * SIEVE_SIZE;
    init_buckets(M);

    uint8_t threshold;
    {
        double log_M = log2((double)M);
        double log_N = bits * 0.5;  /* log2(sqrt(kN)) approximately */
        double log_target = log_N + log_M;
        int correction = tiny_prime_correction();
        threshold = (uint8_t)(log_target * 0.72 - correction);
        if (threshold < 30) threshold = 30;
        if (threshold > 130) threshold = 130;
    }

    fprintf(stderr, "  target=%d rels, threshold=%d\n", target_relations, threshold);

    /* Main sieve loop */
    int total_polys = 0;
    int total_smooth = 0;
    int total_candidates = 0;
    int *candidates = malloc(SIEVE_SIZE * sizeof(int));

    while (num_relations < target_relations) {
        /* Generate new A polynomial */
        new_poly_a(&par);
        compute_ainv();
        compute_poly_b();
        compute_roots(M);
        precompute_bi_deltas();

        /* Loop over B polynomials */
        do {
            total_polys++;
            fill_buckets(M);

            /* Sieve each block */
            for (int blk = 0; blk < num_buckets; blk++) {
                int block_start = blk * SIEVE_SIZE;
                int block_end = block_start + SIEVE_SIZE;
                if (block_end > 2 * M) block_end = 2 * M;

                sieve_block_medium(block_start, block_end);
                apply_bucket_hits(blk);

                int nc = scan_sieve(threshold, candidates);
                total_candidates += nc;

                /* Trial divide each candidate */
                for (int ci = 0; ci < nc; ci++) {
                    int x_offset = block_start + candidates[ci] - M;

                    /* Compute Q(x) = (ax + b)^2 - kN */
                    mpz_t ax_b, Qval;
                    mpz_inits(ax_b, Qval, NULL);

                    mpz_mul_si(ax_b, poly_a, x_offset);
                    mpz_add(ax_b, ax_b, poly_b);

                    /* Skip mirror positions: only process Y > 0
                     * (mirror at -Y gives same Q but trivial combination) */
                    if (mpz_sgn(ax_b) <= 0) {
                        mpz_clears(ax_b, Qval, NULL);
                        continue;
                    }

                    mpz_mul(Qval, ax_b, ax_b);
                    mpz_sub(Qval, Qval, kN);

                    /* If using multiplier, need to handle: Q = (ax+b)^2 - kN
                     * The actual value we factor is Q/k or something...
                     * Actually for SIQS with multiplier k, Q(x) = (ax+b)^2 - kN
                     * and we need Q(x) to be smooth. The relation is
                     * (ax+b)^2 ≡ Q(x) mod N (not kN).
                     * So Y = (ax+b), and Y^2 - kN = Q(x).
                     * When we combine: Y1*Y2... ≡ sqrt(Q1*Q2...) mod kN
                     * Then gcd(Y1*Y2... - sqrt(Q1*Q2...), N) might give factor.
                     * But we need to be careful with the multiplier in sqrt step. */

                    int *exp = calloc(fb_count, sizeof(int));
                    uint64_t lp1 = 0, lp2 = 0;
                    int result = trial_divide(Qval, exp, &lp1, &lp2);

                    if (result == 1) {
                        /* Fully smooth — verify */
                        if (0) {
                            fprintf(stderr, "  VERIFY sm: Q has %d bits, exps:", (int)mpz_sizeinbase(Qval, 2));
                            int nz = 0;
                            for (int k = 0; k < fb_count && nz < 8; k++)
                                if (exp[k]) { fprintf(stderr, " [%d]=%d", k, exp[k]); nz++; }
                            fprintf(stderr, "\n");
                            /* Verify: product of p^exp should equal |Q| */
                            mpz_t check;
                            mpz_init_set_ui(check, 1);
                            for (int k = 1; k < fb_count; k++)
                                for (int e = 0; e < exp[k]; e++)
                                    mpz_mul_ui(check, check, fb_prime[k]);
                            mpz_t absQ;
                            mpz_init(absQ);
                            mpz_abs(absQ, Qval);
                            if (mpz_cmp(check, absQ) != 0) {
                                gmp_fprintf(stderr, "  MISMATCH: check=%Zd, |Q|=%Zd\n", check, absQ);
                            } else {
                                fprintf(stderr, "  MATCH OK\n");
                            }
                            mpz_clears(check, absQ, NULL);
                        }
                        mpz_t one;
                        mpz_init_set_ui(one, 1);
                        add_full_relation(ax_b, exp, one);
                        mpz_clear(one);
                        total_smooth++;
                    } else if (result == 2) {
                        /* Single large prime */
                        process_slp(ax_b, exp, lp1);
                    } else if (result == 3) {
                        /* Double large prime */
                        process_dlp(ax_b, exp, lp1, lp2);
                    }

                    free(exp);
                    mpz_clears(ax_b, Qval, NULL);
                }
            }

            /* Progress reporting */
            if (total_polys % 200 == 0) {
                clock_gettime(CLOCK_MONOTONIC, &t_now);
                double elapsed = (t_now.tv_sec - t_start.tv_sec) +
                                 (t_now.tv_nsec - t_start.tv_nsec) * 1e-9;
                fprintf(stderr, "  poly=%d rels=%d/%d (sm=%d slp=%d/%d dlp=%d/%d) %.1fs\n",
                        total_polys, num_relations, target_relations,
                        total_smooth, slp_combined, slp_count,
                        dlp_combined, dlp_stored, elapsed);
            }

        } while (next_poly_b(M) && num_relations < target_relations);
    }

    clock_gettime(CLOCK_MONOTONIC, &t_now);
    double sieve_time = (t_now.tv_sec - t_start.tv_sec) +
                        (t_now.tv_nsec - t_start.tv_nsec) * 1e-9;
    fprintf(stderr, "Sieve done: %d rels (%d smooth, %d SLP, %d DLP) in %.1fs\n",
            num_relations, total_smooth, slp_combined, dlp_combined, sieve_time);

    /* Linear algebra */
    fprintf(stderr, "LA: %d x %d matrix\n", num_relations, fb_count);

    int max_deps = 64;
    int *dep_rows = malloc(max_deps * (num_relations + 1) * sizeof(int));
    int ndeps = solve_matrix(dep_rows, max_deps);

    fprintf(stderr, "Found %d dependencies\n", ndeps);

    /* Try each dependency */
    mpz_t factor, cofactor;
    mpz_inits(factor, cofactor, NULL);
    int factored = 0;

    for (int d = 0; d < ndeps && !factored; d++) {
        int *dep = &dep_rows[d * (num_relations + 1)];
        int dep_size = 0;
        while (dep[dep_size] >= 0) dep_size++;
        if (d < 3) fprintf(stderr, "  Dep %d: %d rows\n", d, dep_size);
        if (dep_size == 0) continue;
        if (extract_factor(dep, factor)) {
            mpz_divexact(cofactor, N_orig, factor);
            if (mpz_cmp(factor, cofactor) > 0) mpz_swap(factor, cofactor);
            gmp_printf("%Zd %Zd\n", factor, cofactor);
            factored = 1;
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &t_now);
    double total_time = (t_now.tv_sec - t_start.tv_sec) +
                        (t_now.tv_nsec - t_start.tv_nsec) * 1e-9;

    if (!factored) {
        fprintf(stderr, "FAIL: no factor found from %d dependencies in %.1fs\n", ndeps, total_time);
        return 1;
    }

    fprintf(stderr, "Done: %.3fs (sieve=%.3fs, LA=%.3fs)\n",
            total_time, sieve_time, total_time - sieve_time);

    /* Cleanup */
    free(dep_rows);
    free(candidates);
    mpz_clears(factor, cofactor, N_orig, kN, poly_a, poly_b, NULL);
    gmp_randclear(g_rng);

    return 0;
}
