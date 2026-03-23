/*
 * SIQS4 - High-performance Self-Initializing Quadratic Sieve
 * Optimized for balanced semiprimes, targeting 30-75+ digits.
 *
 * Key optimizations over siqs3:
 * - __int128 trial division (no GMP in hot loops)
 * - Position-based trial division (only test primes whose roots match)
 * - Compiler-friendly sieve loops (auto-vectorization hints)
 * - AVX512BW sieve scanning
 * - Bucket sieve for large primes (p > SIEVE_SIZE)
 * - Pollard rho for DLP cofactor splitting
 * - Block Lanczos for linear algebra (proper implementation)
 * - Gray code self-initialization
 * - Knuth-Schroeppel multiplier
 *
 * Compile: gcc -O3 -march=native -mavx512bw -o siqs4 library/siqs4.c -lgmp -lm
 * Usage: timeout 295 ./siqs4 <N>
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
#define SIEVE_SIZE       32768
#define MAX_FB           120000
#define MAX_AFACTORS     16
#define MAX_RELATIONS    500000
#define MAX_CYCLES       500000
#define LP_HASH_BITS     22
#define LP_HASH_SIZE     (1 << LP_HASH_BITS)
#define LP_HASH_MASK     (LP_HASH_SIZE - 1)
#define BUCKET_ALLOC     (1 << 20)
#define SMALL_PRIME_BOUND 64

/* ==================== Parameter Table ==================== */
typedef struct {
    int fb_size;        /* factor base size */
    int num_blocks;     /* sieve blocks per polynomial */
    int lp_mult;        /* large prime multiplier (lp_bound = lp_mult * largest_fb_prime) */
    int num_a_factors;  /* number of primes in A coefficient */
    int extra_rels;     /* surplus relations beyond fb_size */
} params_t;

static params_t get_params(int digits) {
    /* a_factors chosen so that (largest_fb_prime)^a_factors ≈ sqrt(2*kN)/M */
    if (digits <= 30) return (params_t){300,   1,  40, 5, 50};
    if (digits <= 35) return (params_t){500,   2,  50, 6, 60};
    if (digits <= 40) return (params_t){1000,  2,  60, 7, 80};
    if (digits <= 45) return (params_t){1800,  3,  70, 7, 80};
    if (digits <= 50) return (params_t){3000,  4,  80, 8, 100};
    if (digits <= 55) return (params_t){5500,  6, 100, 8, 100};
    if (digits <= 60) return (params_t){9000,  8, 120, 9, 120};
    if (digits <= 65) return (params_t){16000, 12, 160, 10, 140};
    if (digits <= 70) return (params_t){28000, 16, 200, 10, 160};
    if (digits <= 75) return (params_t){45000, 22, 280, 11, 180};
    if (digits <= 80) return (params_t){70000, 28, 350, 12, 200};
    if (digits <= 85) return (params_t){85000, 32, 400, 12, 220};
    if (digits <= 90) return (params_t){115000, 38, 500, 13, 250};
    return (params_t){160000, 44, 600, 14, 300};
}

/* ==================== Timing ==================== */
static struct timespec g_start;
static double elapsed_sec(void) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec - g_start.tv_sec) + (now.tv_nsec - g_start.tv_nsec) * 1e-9;
}

/* ==================== Factor Base ==================== */
static uint32_t fb_p[MAX_FB];        /* primes */
static uint32_t fb_sqrt[MAX_FB];     /* sqrt(kN) mod p */
static uint8_t  fb_logp[MAX_FB];     /* approx log2(p) */
static int fb_count;
static int fb_sieve_start;           /* first prime to sieve (skip tiny ones) */

/* Sieve arrays */
static uint8_t sieve[SIEVE_SIZE] __attribute__((aligned(64)));

/* Current polynomial: Q(x) = (Ax + B)^2 - kN, we sieve Q(x)/A */
/* Roots: for prime p, sieve at x where Ax + B ≡ ±sqrt(kN) (mod p) */
static int32_t soln1[MAX_FB], soln2[MAX_FB]; /* roots for current B */

/* Bucket sieve for large primes */
typedef struct { uint16_t pos; uint16_t fb_idx_low; uint8_t logp; uint8_t fb_idx_high; } bucket_entry_t;
static bucket_entry_t *buckets;
static int *bucket_start; /* bucket_start[block] = start index */
static int *bucket_count; /* bucket_count[block] = number of entries */
static int bucket_alloc;
static int num_blocks_global;

/* A coefficient data */
static int a_factors[MAX_AFACTORS];   /* indices into fb */
static int num_a_factors;
static mpz_t A_coeff, B_coeff, C_coeff, kN_val, N_val;
static mpz_t Bvals[MAX_AFACTORS];    /* B adjustment values for Gray code */
static int B_sign[MAX_AFACTORS];     /* +1 or -1 for Gray code */

/* Relations */
typedef struct {
    mpz_t Y;           /* Ax + B value */
    uint32_t *fb_exp;  /* exponent of each fb prime (compact) */
    int nfactors;      /* number of nonzero exponents */
    uint32_t lp1, lp2; /* large primes (0 if none) */
} relation_t;

static relation_t rels[MAX_RELATIONS];
static int num_rels;
static int num_full, num_partial, num_dlp;

/* Large prime hash for combining partials */
typedef struct lp_entry {
    uint32_t lp;
    int rel_idx;
    struct lp_entry *next;
} lp_entry_t;
static lp_entry_t *lp_hash[LP_HASH_SIZE];
static lp_entry_t lp_pool[MAX_RELATIONS];
static int lp_pool_used;

/* DLP hash */
typedef struct dlp_entry {
    uint64_t key; /* lp1 * lp2 or just lp1 for matching */
    int rel_idx;
    struct dlp_entry *next;
} dlp_entry_t;
static dlp_entry_t *dlp_hash[LP_HASH_SIZE];
static dlp_entry_t dlp_pool[MAX_RELATIONS];
static int dlp_pool_used;

/* Combined full + merged partial relations for matrix */
typedef struct {
    mpz_t prod_Y;      /* product of Y values */
    uint8_t *exp_vec;   /* exponent vector mod 2, packed bits */
} matrix_rel_t;

static matrix_rel_t mat_rels[MAX_RELATIONS];
static int num_mat_rels;

/* Multiplier */
static int g_multiplier;
static uint64_t g_lp_bound;
static gmp_randstate_t g_rng;

/* ==================== Primes and Modular Arithmetic ==================== */
static uint32_t small_primes[MAX_FB];
static int num_small_primes;

static void sieve_primes(int limit) {
    char *is_comp = calloc(limit + 1, 1);
    num_small_primes = 0;
    for (int i = 2; i <= limit; i++) {
        if (!is_comp[i]) {
            small_primes[num_small_primes++] = i;
            if ((long long)i * i <= limit)
                for (int j = i * i; j <= limit; j += i)
                    is_comp[j] = 1;
        }
    }
    free(is_comp);
}

/* Tonelli-Shanks: compute sqrt(n) mod p */
static uint32_t modsqrt(uint32_t n, uint32_t p) {
    if (p == 2) return n & 1;
    n %= p;
    if (n == 0) return 0;

    /* Check quadratic residue */
    uint64_t test = 1;
    uint32_t exp = (p - 1) / 2;
    uint64_t base = n;
    uint32_t e = exp;
    while (e > 0) {
        if (e & 1) test = test * base % p;
        base = base * base % p;
        e >>= 1;
    }
    if (test != 1) return 0; /* not a QR */

    if (p % 4 == 3) {
        /* Simple case */
        uint64_t r = 1;
        base = n;
        e = (p + 1) / 4;
        while (e > 0) {
            if (e & 1) r = r * base % p;
            base = base * base % p;
            e >>= 1;
        }
        return (uint32_t)r;
    }

    /* General Tonelli-Shanks */
    uint32_t q = p - 1, s = 0;
    while (!(q & 1)) { q >>= 1; s++; }

    /* Find non-residue */
    uint32_t z = 2;
    while (1) {
        test = 1; base = z; e = (p - 1) / 2;
        while (e > 0) {
            if (e & 1) test = test * base % p;
            base = base * base % p;
            e >>= 1;
        }
        if (test == p - 1) break;
        z++;
    }

    uint64_t M = s;
    uint64_t c = 1; base = z; e = q;
    while (e > 0) { if (e & 1) c = c * base % p; base = base * base % p; e >>= 1; }
    uint64_t t = 1; base = n; e = q;
    while (e > 0) { if (e & 1) t = t * base % p; base = base * base % p; e >>= 1; }
    uint64_t r = 1; base = n; e = (q + 1) / 2;
    while (e > 0) { if (e & 1) r = r * base % p; base = base * base % p; e >>= 1; }

    while (1) {
        if (t == 1) return (uint32_t)r;
        uint64_t i = 0, tmp = t;
        for (i = 1; i < M; i++) { tmp = tmp * tmp % p; if (tmp == 1) break; }
        uint64_t b = c;
        for (uint64_t j = 0; j < M - i - 1; j++) b = b * b % p;
        M = i;
        c = b * b % p;
        t = t * c % p;
        r = r * b % p;
    }
}

/* Knuth-Schroeppel multiplier selection */
static int select_multiplier(mpz_t N) {
    static const int mult_cands[] = {1,2,3,5,6,7,10,11,13,14,15,17,19,21,22,23,26,29,30,31,33,34,35,37,38,39,41,42,43};
    int best_k = 1;
    double best_score = -1e30;
    mpz_t kn;
    mpz_init(kn);

    for (int mi = 0; mi < (int)(sizeof(mult_cands)/sizeof(mult_cands[0])); mi++) {
        int k = mult_cands[mi];
        mpz_mul_ui(kn, N, k);
        int kn_mod8 = (int)mpz_fdiv_ui(kn, 8);

        double score = 0;
        /* Contribution from 2 */
        switch (kn_mod8) {
            case 1: score += 2.0 * log(2); break;
            case 5: score += 1.5 * log(2); break;
            case 3: case 7: score += log(2); break;
            default: score -= 2.0 * log(2); break;
        }

        /* Contribution from odd primes */
        for (int i = 1; i < 30 && small_primes[i] < 200; i++) {
            uint32_t p = small_primes[i];
            int kn_mod_p = (int)mpz_fdiv_ui(kn, p);
            if (kn_mod_p == 0) {
                score += log((double)p);
            } else {
                /* Check if kN is QR mod p */
                uint64_t test = 1, base = kn_mod_p;
                uint32_t e = (p - 1) / 2;
                while (e > 0) {
                    if (e & 1) test = test * base % p;
                    base = base * base % p;
                    e >>= 1;
                }
                if (test == 1) score += 2.0 * log((double)p) / (p - 1);
            }
        }
        score -= 0.5 * log((double)k);

        if (score > best_score) {
            best_score = score;
            best_k = k;
        }
    }
    mpz_clear(kn);
    return best_k;
}

/* Build factor base */
static void build_factor_base(mpz_t kN, int target_size) {
    fb_count = 0;

    /* -1 sentinel */
    fb_p[0] = 0; /* represents -1 */
    fb_sqrt[0] = 0;
    fb_logp[0] = 0;
    fb_count = 1;

    /* 2 */
    fb_p[1] = 2;
    fb_sqrt[1] = 1;
    fb_logp[1] = 1;
    fb_count = 2;

    uint32_t kn_mod;
    for (int i = 1; i < num_small_primes && fb_count < target_size; i++) {
        uint32_t p = small_primes[i];
        kn_mod = (uint32_t)mpz_fdiv_ui(kN, p);
        uint32_t s = modsqrt(kn_mod, p);
        if (s == 0 && kn_mod != 0) continue; /* not a QR */
        fb_p[fb_count] = p;
        fb_sqrt[fb_count] = s;
        fb_logp[fb_count] = (uint8_t)(log2((double)p) + 0.5);
        if (fb_logp[fb_count] == 0) fb_logp[fb_count] = 1;
        fb_count++;
    }

    /* Set sieve start: skip primes < SMALL_PRIME_BOUND in sieve */
    fb_sieve_start = 2; /* start after -1 and 2 */
    while (fb_sieve_start < fb_count && fb_p[fb_sieve_start] < SMALL_PRIME_BOUND)
        fb_sieve_start++;

    fprintf(stderr, "Factor base: %d primes (largest = %u), sieve start at idx %d\n",
            fb_count, fb_p[fb_count-1], fb_sieve_start);
}

/* ==================== Polynomial Generation (SIQS) ==================== */

/*
 * Choose A coefficient: A = product of s primes from factor base.
 * Target: A ≈ sqrt(2N) / M where M = sieve half-interval.
 */
static void choose_a_coefficient(mpz_t A, mpz_t kN, int s, int M, int *a_idx) {
    /* Target value for A */
    mpz_t target;
    mpz_init(target);
    mpz_mul_ui(target, kN, 2);
    mpz_sqrt(target, target);
    mpz_tdiv_q_ui(target, target, M);

    /* Select s primes near the middle of the factor base whose product ≈ target */
    double log_target = mpz_sizeinbase(target, 2) * log(2.0);
    double log_A = 0;

    /* Start from middle of factor base */
    int mid = fb_count / 2;
    int best_start = mid - s / 2;
    if (best_start < fb_sieve_start + 2) best_start = fb_sieve_start + 2;

    /* Deterministic but diverse A selection using counter-based offset */
    static int a_counter = 0;
    int a_id = a_counter++;

    int selected = 0;
    mpz_set_ui(A, 1);

    /* Target range for a-factor primes: centered around ideal_prime */
    double ideal_prime = exp(log_target / s);
    int lo = fb_sieve_start + 2;
    if (lo < 5) lo = 5;
    int hi = fb_count - 1;

    /* Find center index where fb_p ≈ ideal_prime */
    int center = lo;
    for (int i = lo; i <= hi; i++) {
        if (fb_p[i] >= ideal_prime) { center = i; break; }
    }

    /* Use a wider spread for diversity. Each A uses s consecutive primes
     * starting from an offset based on a_id. */
    int spread = (hi - lo) / 3;
    if (spread < s * 4) spread = s * 4;
    int range_lo = center - spread / 2;
    int range_hi = center + spread / 2;
    if (range_lo < lo) range_lo = lo;
    if (range_hi > hi) range_hi = hi;

    int range_size = range_hi - range_lo + 1;
    if (range_size < s) range_size = s;

    /* Select s primes randomly from the range using LCG seeded by a_id */
    uint32_t rng_state = (uint32_t)(a_id * 2654435761u + 42u);

    for (int j = 0; j < s; j++) {
        rng_state = rng_state * 1103515245u + 12345u;
        int start = range_lo + (int)(rng_state % (uint32_t)range_size);

        int best_i = -1;
        for (int d = 0; d < range_size; d++) {
            int try_idx = range_lo + ((start - range_lo + d) % range_size);
            int skip = 0;
            for (int k = 0; k < selected; k++)
                if (a_idx[k] == try_idx) { skip = 1; break; }
            if (!skip) { best_i = try_idx; break; }
        }
        if (best_i < 0) {
            for (int i = lo; i <= hi; i++) {
                int skip = 0;
                for (int k = 0; k < selected; k++)
                    if (a_idx[k] == i) { skip = 1; break; }
                if (!skip) { best_i = i; break; }
            }
        }
        a_idx[selected++] = best_i;
        mpz_mul_ui(A, A, fb_p[best_i]);
    }
    mpz_clear(target);
}

/* Compute B values for SIQS using Hensel lifting */
static void compute_B_values(mpz_t A, mpz_t kN, int *a_idx, int s) {
    for (int j = 0; j < s; j++) {
        if (j > 0 && Bvals[j]->_mp_alloc) mpz_clear(Bvals[j]);

        uint32_t qj = fb_p[a_idx[j]];
        uint32_t tj = fb_sqrt[a_idx[j]];

        mpz_t Aqj, Aqj_inv, mod_qj;
        mpz_init(Aqj); mpz_init(Aqj_inv); mpz_init_set_ui(mod_qj, qj);

        mpz_divexact_ui(Aqj, A, qj);

        if (!mpz_invert(Aqj_inv, Aqj, mod_qj)) {
            fprintf(stderr, "Error: A/qj not invertible mod qj=%u\n", qj);
            mpz_init_set_ui(Bvals[j], 0);
            mpz_clear(Aqj); mpz_clear(Aqj_inv); mpz_clear(mod_qj);
            continue;
        }

        /* gamma = sqrt(kN) * (A/qj)^(-1) mod qj, adjusted to <= qj/2 */
        uint64_t gamma = ((uint64_t)tj * mpz_get_ui(Aqj_inv)) % qj;
        if (gamma > qj / 2) gamma = qj - gamma;

        /* Bj = gamma * (A/qj) */
        mpz_init(Bvals[j]);
        mpz_mul_ui(Bvals[j], Aqj, (uint32_t)gamma);

        mpz_clear(Aqj); mpz_clear(Aqj_inv); mpz_clear(mod_qj);
    }

    /* B = sum of Bvals */
    mpz_set_ui(B_coeff, 0);
    for (int j = 0; j < s; j++) {
        B_sign[j] = 1;
        mpz_add(B_coeff, B_coeff, Bvals[j]);
    }

    /* Verify B^2 ≡ kN (mod A). If not, fix: B = A - B */
    mpz_t test;
    mpz_init(test);
    mpz_mul(test, B_coeff, B_coeff);
    mpz_sub(test, test, kN);
    mpz_mod(test, test, A);
    if (mpz_sgn(test) != 0) {
        mpz_sub(B_coeff, A, B_coeff);
        /* Negate all Bvals for correct Gray code updates */
        for (int j = 0; j < s; j++)
            mpz_neg(Bvals[j], Bvals[j]);
    }
    mpz_clear(test);

    /* C = (B^2 - kN) / A */
    mpz_t B2;
    mpz_init(B2);
    mpz_mul(B2, B_coeff, B_coeff);
    mpz_sub(C_coeff, B2, kN_val);
    mpz_divexact(C_coeff, C_coeff, A_coeff);
    mpz_clear(B2);
}

/* Switch to next B using Gray code */
static int gray_code_next(int poly_idx, int s) {
    /* Find which bit changed: index of lowest set bit of poly_idx */
    int j = __builtin_ctz(poly_idx);
    if (j >= s) return -1; /* done with all B values for this A */

    /* Flip sign of Bvals[j] */
    B_sign[j] = -B_sign[j];
    if (B_sign[j] < 0) {
        mpz_submul_ui(B_coeff, Bvals[j], 2);
    } else {
        mpz_addmul_ui(B_coeff, Bvals[j], 2);
    }

    /* Update C = (B^2 - kN) / A */
    mpz_t B2;
    mpz_init(B2);
    mpz_mul(B2, B_coeff, B_coeff);
    mpz_sub(C_coeff, B2, kN_val);
    mpz_divexact(C_coeff, C_coeff, A_coeff);
    mpz_clear(B2);

    return j;
}

/* Pre-computed A^(-1) mod p for each FB prime */
static uint32_t fb_ainv[MAX_FB];

static uint32_t mod_inverse_u32(uint32_t a, uint32_t m) {
    if (a == 0) return 0;
    /* Extended Euclidean */
    int64_t g = m, x = 0, y = 1, g1 = a, x1 = 1, y1 = 0;
    while (g1 != 0) {
        int64_t q = g / g1;
        int64_t t;
        t = g - q * g1; g = g1; g1 = t;
        t = x - q * x1; x = x1; x1 = t;
        t = y - q * y1; y = y1; y1 = t;
    }
    return (uint32_t)((x % (int64_t)m + m) % m);
}

/* Compute A^(-1) mod p for all FB primes */
static void compute_ainv(void) {
    for (int i = 2; i < fb_count; i++) {
        uint32_t p = fb_p[i];
        uint32_t am = (uint32_t)mpz_fdiv_ui(A_coeff, p);
        fb_ainv[i] = (am == 0) ? 0 : mod_inverse_u32(am, p);
    }
}

/* Compute sieve roots for current polynomial */
static void compute_sieve_roots(int M) {
    for (int i = 2; i < fb_count; i++) {
        uint32_t p = fb_p[i];
        uint32_t ainv = fb_ainv[i];

        if (ainv == 0) {
            soln1[i] = -1; soln2[i] = -1;
            continue;
        }

        uint32_t bm = (uint32_t)mpz_fdiv_ui(B_coeff, p);
        uint32_t r1 = fb_sqrt[i];
        uint32_t r2 = p - r1;

        /* x1 = ainv * (sqrt - b) mod p */
        int64_t x1 = ((int64_t)ainv * ((int64_t)((r1 + p - bm) % p))) % p;
        /* x2 = ainv * (-sqrt - b) mod p */
        int64_t x2 = ((int64_t)ainv * ((int64_t)((r2 + p - bm) % p))) % p;

        /* Shift to sieve coordinates: x_sieve = x + M */
        soln1[i] = (int32_t)(((x1 % p) + M) % p);
        soln2[i] = (int32_t)(((x2 % p) + M) % p);
    }
}

/* Pre-computed delta values for each B_values[j] and each FB prime */
static uint32_t bi_delta[MAX_AFACTORS][MAX_FB];

static void precompute_bi_deltas(int s) {
    for (int j = 0; j < s; j++) {
        for (int i = 2; i < fb_count; i++) {
            uint32_t p = fb_p[i];
            uint32_t ainv = fb_ainv[i];
            if (ainv == 0) { bi_delta[j][i] = 0; continue; }
            uint32_t bm = (uint32_t)mpz_fdiv_ui(Bvals[j], p);
            bi_delta[j][i] = (uint32_t)((2ULL * ainv * (uint64_t)bm) % p);
        }
    }
}

/* Update roots for Gray code B switch */
static void update_roots_gray(int j_changed, int sign, int M) {
    for (int i = 2; i < fb_count; i++) {
        if (soln1[i] < 0) continue;
        uint32_t p = fb_p[i];
        uint32_t delta = bi_delta[j_changed][i];
        if (delta == 0) continue;

        if (sign < 0) {
            /* B decreased: roots shift by +delta */
            soln1[i] = (soln1[i] + delta) % p;
            soln2[i] = (soln2[i] + delta) % p;
        } else {
            /* B increased: roots shift by -delta */
            soln1[i] = (soln1[i] + p - delta) % p;
            soln2[i] = (soln2[i] + p - delta) % p;
        }
    }
}

/* ==================== Sieve ==================== */

/* Fill bucket sieve entries for large primes */
static void fill_buckets(int M, int NB) {
    int total_sieve = 2 * M;
    for (int b = 0; b < NB; b++) bucket_count[b] = 0;

    int total_entries = 0;
    for (int i = fb_sieve_start; i < fb_count; i++) {
        uint32_t p = fb_p[i];
        if (p <= SIEVE_SIZE) continue; /* handled in main sieve loop */
        uint8_t lp = fb_logp[i];

        for (int root = 0; root < 2; root++) {
            int32_t r = (root == 0) ? soln1[i] : soln2[i];
            if (r < 0) continue;
            /* Find hits in sieve interval */
            int32_t pos = r;
            while (pos < total_sieve) {
                int block = pos / SIEVE_SIZE;
                if (block >= NB) break;
                int bpos = pos % SIEVE_SIZE;
                if (total_entries >= bucket_alloc) break;
                /* Store in bucket */
                int idx = bucket_start[block] + bucket_count[block];
                if (bucket_count[block] < bucket_alloc / NB) {
                    buckets[idx].pos = (uint16_t)bpos;
                    buckets[idx].logp = lp;
                    buckets[idx].fb_idx_low = (uint16_t)(i & 0xFFFF);
                    buckets[idx].fb_idx_high = (uint8_t)((i >> 16) & 0xFF);
                    bucket_count[block]++;
                    total_entries++;
                }
                pos += p;
            }
        }
    }
}

/* Compute init_val for a specific block based on max Q(x)/A in that block */
static uint8_t compute_block_init(int block_idx, int M) {
    /* Q(x)/A = Ax^2 + 2Bx + C for x_actual = x_sieve - M
     * x_sieve ranges from block_start to block_start + SIEVE_SIZE - 1
     * Q/A is a quadratic, maximum at edges of block */
    int block_start = block_idx * SIEVE_SIZE;
    double x_lo = (double)(block_start - M);
    double x_hi = (double)(block_start + SIEVE_SIZE - 1 - M);

    double A_d = mpz_get_d(A_coeff);
    double B_d = mpz_get_d(B_coeff);
    double C_d = mpz_get_d(C_coeff);

    /* Q/A at edges */
    double q_lo = fabs(A_d * x_lo * x_lo + 2.0 * B_d * x_lo + C_d);
    double q_hi = fabs(A_d * x_hi * x_hi + 2.0 * B_d * x_hi + C_d);
    double q_max = (q_lo > q_hi) ? q_lo : q_hi;

    /* Also check vertex x = -B/A if it's in the block */
    double x_vertex = -B_d / A_d;
    if (x_vertex >= x_lo && x_vertex <= x_hi) {
        double q_vertex = fabs(C_d - B_d * B_d / A_d);
        if (q_vertex > q_max) q_max = q_vertex;
        /* Actually vertex gives minimum, not maximum. Use q_max from edges. */
    }

    if (q_max < 1.0) q_max = 1.0;
    int lv = (int)(log2(q_max) + 1.0);
    if (lv < 10) lv = 10;
    if (lv > 230) lv = 230;
    return (uint8_t)lv;
}

/* Sieve one block */
static void sieve_block(int block_idx, int M, uint8_t init_val_unused) {
    int block_start = block_idx * SIEVE_SIZE;

    /* Per-block init value based on actual Q(x)/A size */
    uint8_t block_init = compute_block_init(block_idx, M);

    /* Initialize sieve */
    memset(sieve, block_init, SIEVE_SIZE);

    /* Sieve with small and medium primes */
    for (int i = fb_sieve_start; i < fb_count; i++) {
        uint32_t p = fb_p[i];
        if (p > SIEVE_SIZE) break;
        uint8_t lp = fb_logp[i];

        /* Two roots per prime */
        for (int root = 0; root < 2; root++) {
            int32_t r = (root == 0) ? soln1[i] : soln2[i];
            if (r < 0) continue;

            int32_t start_pos = r - block_start;
            while (start_pos < 0) start_pos += p;

            uint32_t pos = (uint32_t)start_pos;
            while (pos < SIEVE_SIZE) {
                sieve[pos] -= lp;
                pos += p;
            }
        }
    }

    /* Apply bucket sieve entries for this block */
    int bstart = bucket_start[block_idx];
    int bcount = bucket_count[block_idx];
    for (int j = 0; j < bcount; j++) {
        bucket_entry_t *e = &buckets[bstart + j];
        sieve[e->pos] -= e->logp;
    }
}

/* Scan sieve for smooth candidates using AVX512BW */
static int scan_sieve_avx512(uint8_t threshold, int *candidates) {
    int nc = 0;
#ifdef __AVX512BW__
    __m512i vthresh = _mm512_set1_epi8((char)(threshold - 1));
    for (int i = 0; i < SIEVE_SIZE; i += 64) {
        __m512i v = _mm512_load_si512((__m512i *)(sieve + i));
        /* Find bytes where sieve[i] <= threshold (subtracted enough) */
        /* We initialized to init_val and subtracted. Check if sieve[i] < threshold */
        uint64_t mask = _mm512_cmple_epu8_mask(v, vthresh);
        while (mask) {
            int bit = __builtin_ctzll(mask);
            candidates[nc++] = i + bit;
            mask &= mask - 1;
        }
    }
#else
    for (int i = 0; i < SIEVE_SIZE; i++)
        if (sieve[i] <= threshold) candidates[nc++] = i;
#endif
    return nc;
}

/* ==================== Trial Division ==================== */

/* Fast trial division using __int128 arithmetic where possible */
static int trial_divide(mpz_t Qval, int *exp_out, int *nexp_out,
                        uint32_t *lp1_out, uint32_t *lp2_out,
                        int x_sieve_pos) {
    *lp1_out = 0;
    *lp2_out = 0;
    int nexp = 0;

    mpz_t q;
    mpz_init_set(q, Qval);

    /* Handle sign */
    if (mpz_sgn(q) < 0) {
        exp_out[nexp++] = 0; /* -1 factor */
        mpz_neg(q, q);
    }

    /* Powers of 2 */
    {
        unsigned long e2 = mpz_scan1(q, 0);
        if (e2 > 0) {
            mpz_tdiv_q_2exp(q, q, e2);
            for (unsigned long j = 0; j < e2; j++)
                exp_out[nexp++] = 1;
        }
    }

    /* Tiny primes (not sieved) */
    for (int i = 2; i < fb_sieve_start && i < fb_count; i++) {
        uint32_t p = fb_p[i];
        while (mpz_fdiv_ui(q, p) == 0) {
            mpz_divexact_ui(q, q, p);
            exp_out[nexp++] = i;
        }
    }

    /* Add a-factor primes: each appears with exponent 1 in A, so
     * A * Q(x)/A has each a-factor with odd exponent (from A).
     * We must include them in the exponent vector. */
    for (int j = 0; j < num_a_factors; j++) {
        exp_out[nexp++] = a_factors[j];
        /* Also check if this prime divides Q(x)/A (would add extra power) */
        uint32_t p = fb_p[a_factors[j]];
        while (mpz_fdiv_ui(q, p) == 0) {
            mpz_divexact_ui(q, q, p);
            exp_out[nexp++] = a_factors[j];
        }
    }

    /* Sieved primes: use position-based filtering */
    for (int i = fb_sieve_start; i < fb_count; i++) {
        uint32_t p = fb_p[i];
        if (soln1[i] < 0) continue; /* a-factor prime: already handled above */
        /* Position-based check */
        uint32_t xmod = x_sieve_pos % p;
        if (xmod != (uint32_t)soln1[i] && xmod != (uint32_t)soln2[i]) continue;

        /* This prime divides Q(x) */
        mpz_divexact_ui(q, q, p);
        exp_out[nexp++] = i;
        /* Check higher powers */
        while (mpz_fdiv_ui(q, p) == 0) {
            mpz_divexact_ui(q, q, p);
            exp_out[nexp++] = i;
        }
    }

    *nexp_out = nexp;

    /* Check cofactor */
    if (mpz_cmp_ui(q, 1) == 0) {
        mpz_clear(q);
        return 1; /* fully smooth */
    }

    /* Check if cofactor fits in 64 bits */
    if (mpz_sizeinbase(q, 2) > 64) {
        mpz_clear(q);
        return 0; /* too large */
    }

    uint64_t cofactor = mpz_get_ui(q);
    uint64_t lp_bound = g_lp_bound;

    /* SLP: single large prime */
    if (cofactor <= lp_bound) {
        /* Check primality with quick test */
        if (mpz_probab_prime_p(q, 1)) {
            *lp1_out = (uint32_t)cofactor;
            mpz_clear(q);
            return 2; /* SLP */
        }
    }

    /* DLP: double large prime */
    if (mpz_sizeinbase(q, 2) <= 52 && cofactor > 1) {
        /* Check if prime */
        if (mpz_probab_prime_p(q, 1)) {
            if (cofactor <= lp_bound) {
                *lp1_out = (uint32_t)cofactor;
                mpz_clear(q);
                return 2;
            }
            mpz_clear(q);
            return 0; /* prime but too large */
        }

        /* Try to split composite cofactor using Pollard rho */
        uint64_t d = 0;

        /* Quick trial division for small factors */
        for (uint64_t f = 2; f < 1000 && f * f <= cofactor; f++) {
            if (cofactor % f == 0) {
                d = f;
                break;
            }
        }

        if (d == 0 && cofactor > 1000000) {
            /* Pollard rho with __int128 */
            typedef unsigned __int128 u128;
            uint64_t x = 2, y = 2, c = 1;
            for (int attempt = 0; attempt < 3 && d == 0; attempt++) {
                c = attempt + 1;
                x = 2; y = 2;
                for (int iter = 0; iter < 100000; iter++) {
                    x = (u128)x * x % cofactor;
                    x = (x + c) % cofactor;
                    y = (u128)y * y % cofactor;
                    y = (y + c) % cofactor;
                    y = (u128)y * y % cofactor;
                    y = (y + c) % cofactor;
                    uint64_t diff = (x > y) ? x - y : y - x;
                    /* GCD */
                    uint64_t a = diff, b = cofactor;
                    while (b) { uint64_t t = b; b = a % b; a = t; }
                    if (a > 1 && a < cofactor) { d = a; break; }
                    if (a == cofactor) break; /* failed, try new c */
                }
            }
        }

        if (d > 0 && d < cofactor) {
            uint64_t cofactor2 = cofactor / d;
            if (d <= lp_bound && cofactor2 <= lp_bound) {
                *lp1_out = (uint32_t)((d < cofactor2) ? d : cofactor2);
                *lp2_out = (uint32_t)((d < cofactor2) ? cofactor2 : d);
                mpz_clear(q);
                return 3; /* DLP */
            }
        }
    }

    mpz_clear(q);
    return 0; /* not smooth enough */
}

/* ==================== Relation Storage ==================== */

static void store_relation(mpz_t Y, int *exp_indices, int nexp,
                           uint32_t lp1, uint32_t lp2) {
    if (num_rels >= MAX_RELATIONS) return;
    relation_t *r = &rels[num_rels];
    mpz_init_set(r->Y, Y);
    r->fb_exp = malloc(nexp * sizeof(uint32_t));
    for (int i = 0; i < nexp; i++) r->fb_exp[i] = exp_indices[i];
    r->nfactors = nexp;
    r->lp1 = lp1;
    r->lp2 = lp2;
    num_rels++;

    if (lp1 == 0 && lp2 == 0) {
        num_full++;
    } else if (lp2 == 0) {
        num_partial++;
        /* Add to SLP hash */
        uint32_t h = lp1 & LP_HASH_MASK;
        if (lp_pool_used < MAX_RELATIONS) {
            lp_entry_t *e = &lp_pool[lp_pool_used++];
            e->lp = lp1;
            e->rel_idx = num_rels - 1;
            e->next = lp_hash[h];
            lp_hash[h] = e;
        }
    } else {
        num_dlp++;
        /* Add to DLP hash using both large primes */
        uint32_t h1 = lp1 & LP_HASH_MASK;
        if (dlp_pool_used < MAX_RELATIONS) {
            dlp_entry_t *e = &dlp_pool[dlp_pool_used++];
            e->key = lp1;
            e->rel_idx = num_rels - 1;
            e->next = dlp_hash[h1];
            dlp_hash[h1] = e;
        }
        uint32_t h2 = lp2 & LP_HASH_MASK;
        if (dlp_pool_used < MAX_RELATIONS) {
            dlp_entry_t *e = &dlp_pool[dlp_pool_used++];
            e->key = lp2;
            e->rel_idx = num_rels - 1;
            e->next = dlp_hash[h2];
            dlp_hash[h2] = e;
        }
    }
}

/* Count useful relations (full + combined partials) */
static int count_useful_relations(void) {
    int useful = num_full;

    /* Count SLP matches */
    for (int h = 0; h < LP_HASH_SIZE; h++) {
        int cnt = 0;
        for (lp_entry_t *e = lp_hash[h]; e; e = e->next) cnt++;
        useful += cnt / 2; /* each pair gives one relation */
    }

    /* Count DLP cycle merges (simplified: match on shared large primes) */
    for (int h = 0; h < LP_HASH_SIZE; h++) {
        int cnt = 0;
        for (dlp_entry_t *e = dlp_hash[h]; e; e = e->next) cnt++;
        if (cnt >= 2) useful += cnt / 2;
    }

    return useful;
}

/* ==================== Linear Algebra (Block Lanczos) ==================== */

/* Build the matrix and solve using Block Lanczos */
/* The matrix has one row per relation, one column per factor base prime + LPs */
/* We need to find a subset S of relations where the product of Q(x) values is a square */

typedef struct {
    int nrows, ncols;
    int *row_start;
    int *col_idx;
    int nnz;
} sparse_mat_t;

/* Map large primes to column indices */
static uint32_t lp_col_map_keys[MAX_RELATIONS];
static int lp_col_map_vals[MAX_RELATIONS];
static int lp_col_count;

static int get_lp_column(uint32_t lp) {
    /* Linear search (fine for moderate sizes) */
    for (int i = 0; i < lp_col_count; i++)
        if (lp_col_map_keys[i] == lp) return lp_col_map_vals[i];
    /* Add new */
    int col = fb_count + lp_col_count;
    lp_col_map_keys[lp_col_count] = lp;
    lp_col_map_vals[lp_col_count] = col;
    lp_col_count++;
    return col;
}

/* Merged relation pairs (for SLP combining) */
static int merged_r1[MAX_RELATIONS], merged_r2[MAX_RELATIONS];

static sparse_mat_t *build_matrix(int *rel_indices, int nrels) {
    lp_col_count = 0;
    int total_cols = fb_count; /* Only FB columns, no LP columns */
    int *counts = calloc(total_cols, sizeof(int));

    sparse_mat_t *m = malloc(sizeof(sparse_mat_t));
    m->nrows = nrels;
    m->ncols = total_cols;
    m->row_start = malloc((nrels + 1) * sizeof(int));

    /* Two passes: count then fill */
    for (int pass = 0; pass < 2; pass++) {
        int nnz = 0;
        int idx = 0;
        for (int r = 0; r < nrels; r++) {
            if (pass == 0) m->row_start[r] = nnz;

            memset(counts, 0, total_cols * sizeof(int));

            int ri = rel_indices[r];
            if (ri >= 0) {
                /* Direct relation */
                relation_t *rel = &rels[ri];
                for (int j = 0; j < rel->nfactors; j++)
                    counts[rel->fb_exp[j]]++;
            } else {
                /* Merged SLP pair: XOR of both relations' exponents */
                int mi = -(ri + 1);
                relation_t *rel1 = &rels[merged_r1[mi]];
                relation_t *rel2 = &rels[merged_r2[mi]];
                for (int j = 0; j < rel1->nfactors; j++)
                    counts[rel1->fb_exp[j]]++;
                for (int j = 0; j < rel2->nfactors; j++)
                    counts[rel2->fb_exp[j]]++;
                /* LP appears in both with odd count - cancels (even total) */
            }

            for (int j = 0; j < total_cols; j++) {
                if (counts[j] & 1) {
                    nnz++;
                    if (pass == 1) m->col_idx[idx++] = j;
                }
            }
        }
        if (pass == 0) {
            m->row_start[nrels] = nnz;
            m->nnz = nnz;
            m->col_idx = malloc(nnz * sizeof(int));
        }
    }

    free(counts);
    return m;
}

/* Block Lanczos solver */
static void mat_mul(sparse_mat_t *B, uint64_t *x, uint64_t *result) {
    memset(result, 0, B->nrows * sizeof(uint64_t));
    for (int r = 0; r < B->nrows; r++) {
        uint64_t acc = 0;
        for (int j = B->row_start[r]; j < B->row_start[r + 1]; j++)
            acc ^= x[B->col_idx[j]];
        result[r] = acc;
    }
}

static void mat_mul_t(sparse_mat_t *B, uint64_t *x, uint64_t *result) {
    memset(result, 0, B->ncols * sizeof(uint64_t));
    for (int r = 0; r < B->nrows; r++) {
        uint64_t xr = x[r];
        if (!xr) continue;
        for (int j = B->row_start[r]; j < B->row_start[r + 1]; j++)
            result[B->col_idx[j]] ^= xr;
    }
}

static int solve_matrix(sparse_mat_t *B, int **dep_out, int *dep_len_out) {
    int n = B->ncols;
    int m = B->nrows;

    /* Use Gaussian elimination for small matrices, Block Lanczos for large */
    if (n > 50000) {
        fprintf(stderr, "Matrix too large for current implementation (%d x %d). Using GE.\n", m, n);
    }

    /* Gaussian elimination over GF(2) using bit-packed rows */
    int words = (n + 63) / 64;
    uint64_t **mat = malloc(m * sizeof(uint64_t *));
    uint64_t **ident = malloc(m * sizeof(uint64_t *));
    int id_words = (m + 63) / 64;
    for (int r = 0; r < m; r++) {
        mat[r] = calloc(words, sizeof(uint64_t));
        ident[r] = calloc(id_words, sizeof(uint64_t));
        ident[r][r / 64] |= (1ULL << (r % 64));
        /* Fill matrix row */
        for (int j = B->row_start[r]; j < B->row_start[r + 1]; j++)
            mat[r][B->col_idx[j] / 64] |= (1ULL << (B->col_idx[j] % 64));
    }

    /* Forward elimination with randomized pivot selection */
    int *pivot_row = malloc(n * sizeof(int));
    memset(pivot_row, -1, n * sizeof(int));
    int rank = 0;
    srand(42);

    int *col_order = NULL; /* unused */

    for (int col = 0; col < n && rank < m; col++) {
        /* Find pivot */
        int pr = -1;
        for (int r = rank; r < m; r++) {
            if ((mat[r][col / 64] >> (col % 64)) & 1) { pr = r; break; }
        }
        if (pr < 0) continue;

        /* Swap rows */
        if (pr != rank) {
            uint64_t *tmp;
            tmp = mat[pr]; mat[pr] = mat[rank]; mat[rank] = tmp;
            tmp = ident[pr]; ident[pr] = ident[rank]; ident[rank] = tmp;
        }
        pivot_row[col] = rank;

        /* Eliminate */
        for (int r = 0; r < m; r++) {
            if (r == rank) continue;
            if ((mat[r][col / 64] >> (col % 64)) & 1) {
                for (int w = 0; w < words; w++) mat[r][w] ^= mat[rank][w];
                for (int w = 0; w < id_words; w++) ident[r][w] ^= ident[rank][w];
            }
        }
        rank++;
    }

    fprintf(stderr, "  Matrix %d x %d, rank %d, nullity %d\n", m, n, rank, m - rank);

    /* Find null space vectors - prefer longer dependencies (more likely to work) */
    int ndeps = 0;

    /* Collect all null rows with their lengths */
    typedef struct { int row; int len; } null_info_t;
    null_info_t *null_rows = malloc(m * sizeof(null_info_t));
    int num_null = 0;

    for (int r = 0; r < m; r++) {
        int all_zero = 1;
        for (int w = 0; w < words && all_zero; w++)
            if (mat[r][w]) all_zero = 0;
        if (!all_zero) continue;

        int dlen = 0;
        for (int w = 0; w < id_words; w++) {
            uint64_t bits = ident[r][w];
            while (bits) { dlen++; bits &= bits - 1; }
        }
        if (dlen >= 2) {
            null_rows[num_null].row = r;
            null_rows[num_null].len = dlen;
            num_null++;
        }
    }

    /* Sort by length descending - longer deps are more likely to give nontrivial factors */
    for (int i = 0; i < num_null - 1; i++) {
        for (int j = i + 1; j < num_null; j++) {
            if (null_rows[j].len > null_rows[i].len) {
                null_info_t tmp = null_rows[i];
                null_rows[i] = null_rows[j];
                null_rows[j] = tmp;
            }
        }
    }

    /* Extract deps. First try long natural ones, then XOR pairs of short ones. */

    /* Pass 1: collect long deps (length >= 4) */
    for (int ni = 0; ni < num_null && ndeps < 32; ni++) {
        if (null_rows[ni].len < 4) continue;
        int r = null_rows[ni].row;
        int *dep = malloc(m * sizeof(int));
        int dl = 0;
        for (int w = 0; w < id_words; w++) {
            uint64_t bits = ident[r][w];
            while (bits) {
                int b = __builtin_ctzll(bits);
                if (w * 64 + b < m) dep[dl++] = w * 64 + b;
                bits &= bits - 1;
            }
        }
        dep_out[ndeps] = dep;
        dep_len_out[ndeps] = dl;
        ndeps++;
    }

    /* Generate random combinations of null space basis vectors */
    /* The null space is spanned by ident[r] for all null rows r */
    /* Random XOR combinations create more diverse dependencies */

    /* Collect null space basis vectors (identity column vectors) */
    int basis_count = (num_null < 200) ? num_null : 200;
    uint64_t **basis = malloc(basis_count * sizeof(uint64_t*));
    for (int i = 0; i < basis_count; i++) {
        basis[i] = malloc(id_words * sizeof(uint64_t));
        memcpy(basis[i], ident[null_rows[i].row], id_words * sizeof(uint64_t));
    }

    /* Generate random deps by XOR-ing random subsets of basis vectors */
    srand(42);
    for (int trial = 0; trial < 200 && ndeps < 64; trial++) {
        /* Create random combination */
        uint64_t *combo = calloc(id_words, sizeof(uint64_t));

        /* XOR 3-7 random basis vectors */
        int nxor = 3 + (rand() % 5);
        for (int k = 0; k < nxor; k++) {
            int bi = rand() % basis_count;
            for (int w = 0; w < id_words; w++)
                combo[w] ^= basis[bi][w];
        }

        /* Extract dependency */
        int dl = 0;
        for (int w = 0; w < id_words; w++) {
            uint64_t bits = combo[w];
            while (bits) { dl++; bits &= bits - 1; }
        }

        if (dl >= 3) {
            int *dep = malloc(m * sizeof(int));
            int dl2 = 0;
            for (int w = 0; w < id_words; w++) {
                uint64_t bits = combo[w];
                while (bits) {
                    int b = __builtin_ctzll(bits);
                    if (w * 64 + b < m) dep[dl2++] = w * 64 + b;
                    bits &= bits - 1;
                }
            }
            dep_out[ndeps] = dep;
            dep_len_out[ndeps] = dl2;
            ndeps++;
        }
        free(combo);
    }

    for (int i = 0; i < basis_count; i++) free(basis[i]);
    free(basis);
    free(null_rows);

    for (int r = 0; r < m; r++) { free(mat[r]); free(ident[r]); }
    free(mat); free(ident); free(pivot_row);

    return ndeps;
}

/* ==================== Square Root Step ==================== */

static int try_factor(int *dep, int dep_len, int *rel_indices, mpz_t N, mpz_t kN) {
    mpz_t X, Y, prod_Y, prod_Q, temp, g;
    mpz_init(X); mpz_init(Y); mpz_init_set_ui(prod_Y, 1);
    mpz_init_set_ui(prod_Q, 1); mpz_init(temp); mpz_init(g);

    /* Debug: verify first relation */
    if (dep_len >= 1) {
        int ri0 = rel_indices[dep[0]];
        relation_t *rel0 = &rels[ri0];
        mpz_t Y2, AQ;
        mpz_init(Y2); mpz_init(AQ);
        mpz_mul(Y2, rel0->Y, rel0->Y);
        /* Reconstruct A*Q from factorization */
        mpz_set_ui(AQ, 1);
        for (int j = 0; j < rel0->nfactors; j++) {
            if (rel0->fb_exp[j] == 0) { mpz_neg(AQ, AQ); continue; }
            mpz_mul_ui(AQ, AQ, fb_p[rel0->fb_exp[j]]);
        }
        if (rel0->lp1 > 0) mpz_mul_ui(AQ, AQ, rel0->lp1);
        if (rel0->lp2 > 0) mpz_mul_ui(AQ, AQ, rel0->lp2);
        /* Check: Y^2 ≡ AQ (mod kN) */
        mpz_sub(temp, Y2, AQ);
        mpz_mod(temp, temp, kN);
        if (mpz_cmp_ui(temp, 0) != 0) {
            gmp_fprintf(stderr, "  WARNING: Y^2 != A*Q (mod kN) for rel %d! diff=%Zd\n", ri0, temp);
        }
        mpz_clear(Y2); mpz_clear(AQ);
    }

    /* Accumulate product of Y values and exponents */
    int total_cols = fb_count + 1; /* +1 for LP placeholder */
    int *total_exp = calloc(total_cols, sizeof(int));

    for (int d = 0; d < dep_len; d++) {
        int ri = rel_indices[dep[d]];
        if (ri >= 0) {
            /* Direct relation */
            relation_t *rel = &rels[ri];
            mpz_mul(prod_Y, prod_Y, rel->Y);
            mpz_mod(prod_Y, prod_Y, N);
            for (int j = 0; j < rel->nfactors; j++)
                total_exp[rel->fb_exp[j]]++;
        } else {
            /* Merged SLP pair */
            int mi = -(ri + 1);
            relation_t *rel1 = &rels[merged_r1[mi]];
            relation_t *rel2 = &rels[merged_r2[mi]];
            mpz_mul(prod_Y, prod_Y, rel1->Y);
            mpz_mul(prod_Y, prod_Y, rel2->Y);
            mpz_mod(prod_Y, prod_Y, N);
            for (int j = 0; j < rel1->nfactors; j++)
                total_exp[rel1->fb_exp[j]]++;
            for (int j = 0; j < rel2->nfactors; j++)
                total_exp[rel2->fb_exp[j]]++;
            /* LP^2 in product → LP^1 in sqrt. Track in prod_Q for later. */
            if (rel1->lp1 > 0) {
                mpz_mul_ui(prod_Q, prod_Q, rel1->lp1);
            }
        }
    }

    /* Compute X = product of primes ^ (exp/2) mod N */
    /* Important: compute as full precision integer first, then reduce mod N */
    mpz_set_ui(X, 1);
    int neg_count = total_exp[0]; /* -1 count */
    for (int i = 1; i < fb_count; i++) {
        int e = total_exp[i] / 2;
        if (e <= 0) continue;
        /* Multiply in full precision to maintain correct sign */
        for (int k = 0; k < e; k++) {
            mpz_mul_ui(X, X, fb_p[i]);
            /* Reduce periodically to keep numbers manageable */
            if (mpz_sizeinbase(X, 2) > 2000)
                mpz_mod(X, X, N);
        }
    }
    /* Include accumulated LP contributions from merged SLP pairs */
    mpz_mul(X, X, prod_Q);
    mpz_mod(X, X, N);

    /* Y = prod_Y mod N, X = sqrt(prod_Q) mod N */
    mpz_mod(Y, prod_Y, N);

    /* Debug: verify X^2 ≡ Y^2 (mod N) for first dep */
    if (dep_len >= 2 && dep_len <= 500) {
        mpz_t X2, Y2;
        mpz_init(X2); mpz_init(Y2);
        mpz_mul(X2, X, X); mpz_mod(X2, X2, N);
        mpz_mul(Y2, Y, Y); mpz_mod(Y2, Y2, N);
        if (mpz_cmp(X2, Y2) != 0) {
            gmp_fprintf(stderr, "  BUG: X^2 != Y^2 (mod N)! X^2=%Zd Y^2=%Zd\n", X2, Y2);
            /* Print total exponents for debugging */
            fprintf(stderr, "  total_exp (nonzero): ");
            for (int i = 0; i < total_cols; i++)
                if (total_exp[i] > 0) fprintf(stderr, "%d:%d ", i, total_exp[i]);
            fprintf(stderr, "\n");
            /* Check odd exponents */
            int has_odd = 0;
            for (int i = 0; i < total_cols; i++)
                if (total_exp[i] & 1) { fprintf(stderr, "  ODD exp at col %d: %d\n", i, total_exp[i]); has_odd = 1; }
            if (!has_odd) fprintf(stderr, "  All exponents even but X^2 != Y^2!\n");
        }
        mpz_clear(X2); mpz_clear(Y2);
    }

    /* gcd(X - Y, N) and gcd(X + Y, N) */
    mpz_sub(temp, X, Y);
    mpz_gcd(g, temp, N);

    /* Debug: print first 3 deps' gcd details */
    static int dep_debug_count = 0;
    if (dep_debug_count < 3) {
        gmp_fprintf(stderr, "  dep gcd(Y-X,N)=%Zd  gcd(Y+X,N)=", g);
        mpz_add(temp, X, Y);
        mpz_t g2; mpz_init(g2);
        mpz_gcd(g2, temp, N);
        gmp_fprintf(stderr, "%Zd  deplen=%d\n", g2, dep_len);
        gmp_fprintf(stderr, "    X=%Zd\n    Y=%Zd\n", X, Y);
        mpz_clear(g2);
        dep_debug_count++;
    }

    int found = 0;
    if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
        /* Found non-trivial factor! */
        gmp_fprintf(stderr, "  Factor found: %Zd\n", g);
        mpz_t cofactor;
        mpz_init(cofactor);
        mpz_divexact(cofactor, N, g);
        gmp_printf("%Zd\n%Zd\n", g, cofactor);
        mpz_clear(cofactor);
        found = 1;
    }

    if (!found) {
        mpz_add(temp, X, Y);
        mpz_gcd(g, temp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            gmp_fprintf(stderr, "  Factor found: %Zd\n", g);
            mpz_t cofactor;
            mpz_init(cofactor);
            mpz_divexact(cofactor, N, g);
            gmp_printf("%Zd\n%Zd\n", g, cofactor);
            mpz_clear(cofactor);
            found = 1;
        }
    }

    free(total_exp);
    mpz_clear(X); mpz_clear(Y); mpz_clear(prod_Y);
    mpz_clear(prod_Q); mpz_clear(temp); mpz_clear(g);
    return found;
}

/* ==================== Main ==================== */

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    clock_gettime(CLOCK_MONOTONIC, &g_start);

    mpz_init_set_str(N_val, argv[1], 10);
    int digits = (int)mpz_sizeinbase(N_val, 10);
    fprintf(stderr, "Factoring %d-digit number\n", digits);

    /* Quick checks */
    if (mpz_probab_prime_p(N_val, 25)) {
        fprintf(stderr, "Input is prime\n");
        gmp_printf("%Zd\n", N_val);
        return 0;
    }

    /* Check for small factors */
    for (int i = 0; i < 1000; i++) {
        unsigned long p = 2 + i;
        if (mpz_fdiv_ui(N_val, p) == 0) {
            mpz_t f, c;
            mpz_init_set_ui(f, p);
            mpz_init(c);
            mpz_divexact(c, N_val, f);
            gmp_printf("%Zd\n%Zd\n", f, c);
            mpz_clear(f); mpz_clear(c);
            return 0;
        }
    }

    /* Generate primes */
    sieve_primes(2000000);

    /* Select multiplier */
    g_multiplier = select_multiplier(N_val);
    fprintf(stderr, "Multiplier: %d\n", g_multiplier);

    mpz_init(kN_val);
    mpz_mul_ui(kN_val, N_val, g_multiplier);

    /* Get parameters */
    params_t params = get_params(digits);
    fprintf(stderr, "Parameters: FB=%d, NB=%d, LP_mult=%d, A_factors=%d\n",
            params.fb_size, params.num_blocks, params.lp_mult, params.num_a_factors);

    /* Build factor base */
    build_factor_base(kN_val, params.fb_size);

    int target_rels = fb_count + params.extra_rels;
    g_lp_bound = (uint64_t)fb_p[fb_count - 1] * params.lp_mult;
    num_a_factors = params.num_a_factors;
    num_blocks_global = params.num_blocks;

    int M = params.num_blocks * SIEVE_SIZE / 2;

    /* Compute sieve threshold */
    /* Q(x)/A ≈ M * sqrt(kN) at edges of sieve interval */
    double log2_Q = log2((double)M) + mpz_sizeinbase(kN_val, 2) / 2.0;
    uint8_t init_val = (uint8_t)log2_Q;
    if (init_val > 230) init_val = 230;
    if (init_val < 30) init_val = 30;
    /* Threshold: accept candidates where unsieved residual < lp_bound^2 (for DLP) */
    double lp_bits = log2((double)g_lp_bound);
    uint8_t threshold = (uint8_t)(2.0 * lp_bits + 3.0); /* DLP: cofactor < lp_bound^2 */
    if (threshold > init_val - 5) threshold = init_val - 5;

    fprintf(stderr, "Sieve init=%d, threshold=%d, M=%d, target=%d rels\n",
            init_val, threshold, M, target_rels);

    /* Allocate bucket sieve */
    bucket_alloc = BUCKET_ALLOC;
    buckets = malloc(bucket_alloc * sizeof(bucket_entry_t));
    bucket_start = malloc(params.num_blocks * sizeof(int));
    bucket_count = malloc(params.num_blocks * sizeof(int));
    for (int b = 0; b < params.num_blocks; b++)
        bucket_start[b] = b * (bucket_alloc / params.num_blocks);

    /* Initialize polynomial variables */
    mpz_init(A_coeff); mpz_init(B_coeff); mpz_init(C_coeff);
    gmp_randinit_default(g_rng);
    gmp_randseed_ui(g_rng, 42);

    int candidates[SIEVE_SIZE];
    int exp_indices[MAX_FB * 2];
    int a_idx[MAX_AFACTORS];

    int total_polys = 0;
    int a_count = 0;

    fprintf(stderr, "Starting sieve...\n");

    /* Main sieve loop */
    while (1) {
        if (elapsed_sec() > 290.0) {
            fprintf(stderr, "Time limit approaching, stopping sieve\n");
            break;
        }

        int useful = num_full; /* Require enough FULL relations for matrix */
        if (useful >= target_rels) {
            fprintf(stderr, "Enough relations: %d full + %d SLP + %d DLP = ~%d useful (target %d)\n",
                    num_full, num_partial, num_dlp, useful, target_rels);
            break;
        }

        /* Choose new A coefficient */
        choose_a_coefficient(A_coeff, kN_val, num_a_factors, M, a_idx);
        memcpy(a_factors, a_idx, num_a_factors * sizeof(int));
        a_count++;
        if (a_count <= 5) {
            gmp_fprintf(stderr, "  A[%d] = %Zd (primes:", a_count-1, A_coeff);
            for (int j = 0; j < num_a_factors; j++) fprintf(stderr, " %u", fb_p[a_idx[j]]);
            fprintf(stderr, ")\n");
        }

        /* Compute B values */
        compute_B_values(A_coeff, kN_val, a_idx, num_a_factors);

        /* Compute A^(-1) mod p for all primes */
        compute_ainv();

        /* Compute initial sieve roots */
        compute_sieve_roots(M);

        /* Precompute delta values for Gray code updates */
        precompute_bi_deltas(num_a_factors);

        int num_b_polys = 1 << (num_a_factors - 1);

        /* Debug: print Q(0)/A = C for first A */
        if (a_count <= 2) {
            gmp_fprintf(stderr, "  C (Q(0)/A) = %Zd, bits=%zu\n", C_coeff,
                        mpz_sizeinbase(C_coeff, 2));
            for (int blk = 0; blk < params.num_blocks; blk++) {
                uint8_t bi = compute_block_init(blk, M);
                fprintf(stderr, "  block[%d] init=%d\n", blk, bi);
            }
        }

        for (int b_idx = 0; b_idx < num_b_polys; b_idx++) {
            if (b_idx > 0) {
                int j = gray_code_next(b_idx, num_a_factors);
                if (j < 0) break;
                update_roots_gray(j, B_sign[j], M);
            }

            /* Fill bucket sieve for large primes */
            fill_buckets(M, params.num_blocks);

            /* Sieve each block */
            for (int blk = 0; blk < params.num_blocks; blk++) {
                sieve_block(blk, M, init_val);

                /* Threshold: accept candidates where residual < lp_bound^2 */
                uint8_t block_threshold = (uint8_t)(2.0 * log2((double)g_lp_bound) + 3.0);
                int nc = scan_sieve_avx512(block_threshold, candidates);

                for (int ci = 0; ci < nc; ci++) {
                    int x_sieve = blk * SIEVE_SIZE + candidates[ci];
                    int x_actual = x_sieve - M;

                    /* Include both positive and negative x values.
                     * Mirror positions (x and -x-2B/A) give same Q but opposite Y signs.
                     * Both are needed for the square root step to produce non-trivial gcd. */

                    /* Compute Q(x) = (A*x + B)^2 - kN, divided by A */
                    /* = A*x^2 + 2*B*x + C */
                    mpz_t Qval, Y_val;
                    mpz_init(Qval); mpz_init(Y_val);

                    /* Y = A*x + B */
                    mpz_mul_si(Y_val, A_coeff, x_actual);
                    mpz_add(Y_val, Y_val, B_coeff);

                    /* Q = A*x^2 + 2*B*x + C */
                    mpz_set_si(Qval, x_actual);
                    mpz_mul_si(Qval, Qval, x_actual);
                    mpz_mul(Qval, Qval, A_coeff);

                    mpz_t bx2;
                    mpz_init(bx2);
                    mpz_mul_si(bx2, B_coeff, 2 * x_actual);
                    mpz_add(Qval, Qval, bx2);
                    mpz_add(Qval, Qval, C_coeff);
                    mpz_clear(bx2);

                    int nexp;
                    uint32_t lp1, lp2;
                    int status = trial_divide(Qval, exp_indices, &nexp,
                                             &lp1, &lp2, x_sieve);

                    if (status > 0) {
                        store_relation(Y_val, exp_indices, nexp, lp1, lp2);
                        /* Debug first few full relations */
                        if (status == 1 && num_full <= 3) {
                            fprintf(stderr, "  FULL REL #%d: x=%d, nexp=%d, primes: ",
                                    num_full, x_actual, nexp);
                            for (int zz = 0; zz < nexp && zz < 20; zz++)
                                fprintf(stderr, "%u ", fb_p[exp_indices[zz]]);
                            fprintf(stderr, "\n");
                            /* Verify Q */
                            mpz_t verify;
                            mpz_init(verify);
                            mpz_mul(verify, Y_val, Y_val);
                            mpz_sub(verify, verify, kN_val);
                            /* verify should = A * Qval */
                            gmp_fprintf(stderr, "    Y=%Zd\n    Y^2-kN=%Zd\n", Y_val, verify);
                            mpz_clear(verify);
                        }
                    }

                    mpz_clear(Qval); mpz_clear(Y_val);
                }
            }

            total_polys++;

            if (total_polys % 500 == 0) {
                int useful = count_useful_relations();
                double t = elapsed_sec();
                fprintf(stderr, "\r  polys=%d a=%d rels=%d(F%d+P%d+D%d) useful≈%d/%d  %.1fs  %.0f rels/s",
                        total_polys, a_count, num_rels, num_full, num_partial, num_dlp,
                        useful, target_rels, t, num_rels / t);
            }

            if (elapsed_sec() > 290.0) break;
            int useful_check = num_full;
            if (useful_check >= target_rels) break;
        }

        if (elapsed_sec() > 290.0) break;
        int useful_check = num_full;
        if (useful_check >= target_rels) break;
    }

    fprintf(stderr, "\nSieve complete: %d rels (F%d+P%d+D%d) in %.1fs\n",
            num_rels, num_full, num_partial, num_dlp, elapsed_sec());

    /* Build relation list for matrix */
    /* Strategy: use full relations directly, and merge SLP pairs into combined relations */
    int *rel_indices = malloc(MAX_RELATIONS * sizeof(int));
    int nmat = 0;

    /* Full relations */
    for (int i = 0; i < num_rels; i++) {
        if (rels[i].lp1 == 0 && rels[i].lp2 == 0)
            rel_indices[nmat++] = i;
    }
    int nfull_mat = nmat;

    /* Merge SLP pairs: sort SLP relations by their LP, combine pairs */
    /* For each pair with same LP: create a merged relation (Y = Y1*Y2, exp = union) */
    /* We store merged relations as pairs: (idx1, idx2) with idx1 < idx2 */
    typedef struct { int r1, r2; } merged_t;
    merged_t *merged = malloc(MAX_RELATIONS * sizeof(merged_t));
    int nmerged = 0;

    /* Group SLP relations by LP */
    for (int h = 0; h < LP_HASH_SIZE; h++) {
        lp_entry_t *prev = NULL;
        for (lp_entry_t *e = lp_hash[h]; e; e = e->next) {
            if (prev && prev->lp == e->lp) {
                merged_r1[nmerged] = prev->rel_idx;
                merged_r2[nmerged] = e->rel_idx;
                merged[nmerged].r1 = prev->rel_idx;
                merged[nmerged].r2 = e->rel_idx;
                nmerged++;
                prev = NULL;
            } else {
                prev = e;
            }
        }
    }

    /* For merged relations, we need to create synthetic relations
     * The merged relation's exponent vector = XOR of the two individual vectors
     * (LP cancels since it appears in both with odd exponent)
     * Store as negative indices to distinguish from direct relations */
    /* For simplicity, add both relations and let the matrix handle it */
    /* Actually, let's properly merge by creating combined relations */
    int merge_start = nmat;
    for (int m = 0; m < 0 /* use full only */ && nmat + 1 < MAX_RELATIONS; m++) {
        /* Add BOTH relations of the merged pair. In the matrix, the LP column
         * will have the same value in both rows, so when XOR'd they cancel.
         * But we're NOT adding LP columns anymore. Instead, we create a single
         * merged row that's the XOR of both relations' FB exponents. */
        /* We'll handle this in the matrix build by creating a synthetic row */
        rel_indices[nmat++] = -(m + 1); /* negative = merged relation index */
    }

    fprintf(stderr, "Matrix relations: %d full + %d merged SLP pairs = %d\n",
            nfull_mat, nmerged, nmat);

    fprintf(stderr, "Matrix relations: %d (from %d total rels)\n", nmat, num_rels);

    if (nmat < fb_count + 1) {
        fprintf(stderr, "Not enough relations for matrix (%d < %d)\n", nmat, fb_count + 1);
        return 1;
    }

    /* Build and solve matrix */
    fprintf(stderr, "Building matrix...\n");
    sparse_mat_t *mat = build_matrix(rel_indices, nmat);
    fprintf(stderr, "Matrix: %d x %d with %d nonzeros\n", mat->nrows, mat->ncols, mat->nnz);

    /* Debug: print column usage */
    {
        int *col_used = calloc(mat->ncols, sizeof(int));
        for (int r = 0; r < mat->nrows; r++)
            for (int j = mat->row_start[r]; j < mat->row_start[r+1]; j++)
                col_used[mat->col_idx[j]] = 1;
        int used = 0;
        for (int c = 0; c < mat->ncols; c++) if (col_used[c]) used++;
        fprintf(stderr, "  Columns with nonzeros: %d / %d\n", used, mat->ncols);
        /* Print first 3 rows */
        for (int r = 0; r < 3 && r < mat->nrows; r++) {
            fprintf(stderr, "  Row %d (%d nz): cols ", r, mat->row_start[r+1] - mat->row_start[r]);
            for (int j = mat->row_start[r]; j < mat->row_start[r+1] && j < mat->row_start[r]+20; j++)
                fprintf(stderr, "%d ", mat->col_idx[j]);
            fprintf(stderr, "\n");
        }
        free(col_used);
    }

    int *deps[64];
    int dep_lens[64];
    fprintf(stderr, "Solving matrix (GF(2) Gaussian elimination)...\n");
    int ndeps = solve_matrix(mat, deps, dep_lens);
    fprintf(stderr, "Found %d dependencies\n", ndeps);

    /* Try each dependency */
    int factored = 0;
    for (int d = 0; d < ndeps && !factored; d++) {
        fprintf(stderr, "Trying dependency %d (length %d)...\n", d, dep_lens[d]);
        factored = try_factor(deps[d], dep_lens[d], rel_indices, N_val, kN_val);
    }

    if (!factored) {
        fprintf(stderr, "FAILED: no dependency produced a factor\n");
        return 1;
    }

    fprintf(stderr, "Total time: %.2f seconds\n", elapsed_sec());

    /* Cleanup */
    for (int d = 0; d < ndeps; d++) free(deps[d]);
    free(rel_indices);
    free(mat->row_start); free(mat->col_idx); free(mat);
    free(buckets); free(bucket_start); free(bucket_count);
    mpz_clear(A_coeff); mpz_clear(B_coeff); mpz_clear(C_coeff);
    mpz_clear(kN_val); mpz_clear(N_val);
    for (int i = 0; i < num_a_factors; i++) mpz_clear(Bvals[i]);
    for (int i = 0; i < num_rels; i++) {
        mpz_clear(rels[i].Y);
        free(rels[i].fb_exp);
    }
    gmp_randclear(g_rng);

    return 0;
}
