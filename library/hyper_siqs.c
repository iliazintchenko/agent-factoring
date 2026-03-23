/*
 * Hyper SIQS - Batch polynomial SIQS with all known optimizations
 *
 * Key features:
 * 1. Batch polynomial sieving: BATCH=4 polynomials per sieve block pass
 *    - Allocate BATCH sieve arrays (each BLOCK_SIZE bytes)
 *    - FB prime loop iterates ONCE but updates ALL BATCH sieve arrays
 *    - Amortizes outer FB prime loop overhead
 * 2. 48KB sieve blocks (L1d cache = 48KB per core on AMD EPYC 9R45)
 * 3. Gray code self-initialization for polynomial enumeration
 * 4. Bucket sieving for large primes (p > BLOCK_SIZE)
 * 5. SLP (single large prime) matching with hash table
 * 6. DLP (double large prime) with cofactor splitting via Pollard rho
 * 7. Sieve-informed trial division (only check primes whose roots match x)
 * 8. Gaussian elimination for GF(2) linear algebra
 * 9. Standard SIQS square root step
 *
 * Compile: gcc -O3 -march=native -o hyper_siqs library/hyper_siqs.c -lgmp -lm
 * Usage: ./hyper_siqs <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <gmp.h>

/* ==================== Constants ==================== */
#define SEED          42
#define BLOCK_SIZE    49152       /* 48KB = L1d cache */
#define BATCH         4           /* polynomials per sieve pass */
#define MAX_FB        50000
#define MAX_A_FACTORS 25
#define MAX_FULL_RELS 500000
#define MAX_PARTIAL   4000000
#define BUCKET_ALLOC  65536

/* ==================== Timing ==================== */
static struct timespec g_start;
static double elapsed(void) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec - g_start.tv_sec) + (now.tv_nsec - g_start.tv_nsec) / 1e9;
}

/* ==================== Modular Arithmetic ==================== */
static inline uint32_t mod_inv32(uint32_t a, uint32_t m) {
    int64_t old_r = a, r = m, old_s = 1, s = 0;
    while (r) {
        int64_t q = old_r / r, t;
        t = r; r = old_r - q * r; old_r = t;
        t = s; s = old_s - q * s; old_s = t;
    }
    return old_r == 1 ? (uint32_t)(((old_s % (int64_t)m) + m) % m) : 0;
}

static uint32_t sqrt_mod_p(uint32_t n, uint32_t p) {
    if (n == 0) return 0;
    if (p == 2) return n & 1;
    uint64_t mod = p;
    uint64_t b = n % p, e = (p - 1) / 2, r = 1;
    while (e) { if (e & 1) r = r * b % mod; b = b * b % mod; e >>= 1; }
    if (r != 1) return 0;
    if ((p & 3) == 3) {
        b = n % p; e = (p + 1) / 4; r = 1;
        while (e) { if (e & 1) r = r * b % mod; b = b * b % mod; e >>= 1; }
        return (uint32_t)r;
    }
    uint32_t Q = p - 1, S = 0;
    while (!(Q & 1)) { Q >>= 1; S++; }
    uint32_t z = 2;
    while (1) {
        b = z; e = (p - 1) / 2; r = 1;
        while (e) { if (e & 1) r = r * b % mod; b = b * b % mod; e >>= 1; }
        if (r == (uint64_t)(p - 1)) break;
        z++;
    }
    uint64_t M = S;
    b = z; e = Q; uint64_t c = 1;
    while (e) { if (e & 1) c = c * b % mod; b = b * b % mod; e >>= 1; }
    b = n % p; e = Q; uint64_t t = 1;
    while (e) { if (e & 1) t = t * b % mod; b = b * b % mod; e >>= 1; }
    b = n % p; e = (Q + 1) / 2; uint64_t R = 1;
    while (e) { if (e & 1) R = R * b % mod; b = b * b % mod; e >>= 1; }
    while (1) {
        if (t == 1) return (uint32_t)R;
        uint64_t i = 0, tt = t;
        while (tt != 1) { tt = tt * tt % p; i++; }
        uint64_t bb = c;
        for (uint64_t j = 0; j < M - i - 1; j++) bb = bb * bb % p;
        M = i; c = bb * bb % p; t = t * c % p; R = R * bb % p;
    }
}

/* ==================== Factor Base ==================== */
typedef struct {
    uint32_t *prime, *root;
    uint8_t  *logp;
    int size;
} fb_t;

static fb_t *fb_create(mpz_t kN, int target) {
    fb_t *fb = calloc(1, sizeof(fb_t));
    fb->prime = malloc((target + 100) * sizeof(uint32_t));
    fb->root  = malloc((target + 100) * sizeof(uint32_t));
    fb->logp  = malloc((target + 100) * sizeof(uint8_t));
    fb->prime[0] = 2; fb->root[0] = 1; fb->logp[0] = 1; fb->size = 1;
    int bound = target * 30 + 100000;
    char *sv = calloc(bound + 1, 1);
    for (int i = 2; (long)i * i <= bound; i++)
        if (!sv[i]) for (int j = i * i; j <= bound; j += i) sv[j] = 1;
    for (int i = 3; i <= bound && fb->size < target; i += 2) {
        if (sv[i]) continue;
        uint32_t p = (uint32_t)i;
        unsigned long nm = mpz_fdiv_ui(kN, p);
        if (nm == 0) {
            fb->prime[fb->size] = p; fb->root[fb->size] = 0;
            fb->logp[fb->size] = (uint8_t)(log2(p) + 0.5); fb->size++;
            continue;
        }
        uint32_t r = sqrt_mod_p((uint32_t)nm, p);
        if (!r) continue;
        fb->prime[fb->size] = p; fb->root[fb->size] = r;
        fb->logp[fb->size] = (uint8_t)(log2(p) + 0.5); fb->size++;
    }
    free(sv);
    return fb;
}

/* ==================== Multiplier ==================== */
static int choose_multiplier(mpz_t N) {
    static const int ks[] = {1,2,3,5,6,7,10,11,13,14,15,17,19,21,23,29,31,37,41,43,0};
    double best = -1e30; int best_k = 1;
    for (int ki = 0; ks[ki]; ki++) {
        int k = ks[ki]; mpz_t kN; mpz_init(kN); mpz_mul_ui(kN, N, k);
        double s = -0.5 * log((double)k);
        unsigned long m8 = mpz_fdiv_ui(kN, 8);
        if (m8 == 1) s += 2*log(2.0); else if (m8 == 5) s += log(2.0);
        else if (m8 == 3 || m8 == 7) s += 0.5*log(2.0);
        int ps[] = {3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67};
        for (int i = 0; i < 18; i++) {
            if (k % ps[i] == 0) { s += log((double)ps[i])/ps[i]; continue; }
            if (sqrt_mod_p(mpz_fdiv_ui(kN, ps[i]), ps[i])) s += 2.0*log(ps[i])/(ps[i]-1);
        }
        if (s > best) { best = s; best_k = k; }
        mpz_clear(kN);
    }
    return best_k;
}

/* ==================== Bucket Sieve ==================== */
typedef struct { uint32_t pos; uint8_t logp; } bucket_entry_t;
typedef struct { bucket_entry_t *entries; int count, alloc; } bucket_t;

static inline void bucket_add(bucket_t *b, uint32_t pos, uint8_t logp) {
    if (b->count >= b->alloc) {
        b->alloc *= 2;
        b->entries = realloc(b->entries, b->alloc * sizeof(bucket_entry_t));
    }
    b->entries[b->count++] = (bucket_entry_t){pos, logp};
}

/* ==================== LP Hash for SLP matching ==================== */
#define LP_HASH_BITS 22
#define LP_HASH_SIZE (1 << LP_HASH_BITS)
typedef struct lp_e { uint64_t lp; int idx; struct lp_e *next; } lp_e_t;
typedef struct { lp_e_t **b; lp_e_t *pool; int used, max; } lp_hash_t;

static lp_hash_t *lp_create(int m) {
    lp_hash_t *t = calloc(1, sizeof(lp_hash_t));
    t->b = calloc(LP_HASH_SIZE, sizeof(lp_e_t*));
    t->pool = calloc(m, sizeof(lp_e_t));
    t->max = m;
    return t;
}
static int lp_find(lp_hash_t *t, uint64_t lp) {
    uint32_t h = (uint32_t)((lp * 0x9E3779B97F4A7C15ULL) >> (64 - LP_HASH_BITS));
    for (lp_e_t *e = t->b[h]; e; e = e->next) if (e->lp == lp) return e->idx;
    return -1;
}
static void lp_insert(lp_hash_t *t, uint64_t lp, int idx) {
    if (t->used >= t->max) return;
    uint32_t h = (uint32_t)((lp * 0x9E3779B97F4A7C15ULL) >> (64 - LP_HASH_BITS));
    lp_e_t *e = &t->pool[t->used++];
    e->lp = lp; e->idx = idx; e->next = t->b[h]; t->b[h] = e;
}

/* ==================== DLP hash for pairwise matching ==================== */
#define DLP_HASH_BITS 22
#define DLP_HASH_SIZE (1 << DLP_HASH_BITS)
typedef struct dp_e { uint64_t lp; int rel_idx; struct dp_e *next; } dp_e_t;

/* ==================== Relation Storage ==================== */
typedef struct {
    mpz_t *ax_b, *Qx;
    uint64_t *lp1, *lp2;
    int count, alloc;
} rels_t;

static rels_t *rels_create(int n) {
    rels_t *r = calloc(1, sizeof(rels_t));
    r->alloc = n;
    r->ax_b = malloc(n * sizeof(mpz_t));
    r->Qx   = malloc(n * sizeof(mpz_t));
    r->lp1  = calloc(n, sizeof(uint64_t));
    r->lp2  = calloc(n, sizeof(uint64_t));
    for (int i = 0; i < n; i++) { mpz_init(r->ax_b[i]); mpz_init(r->Qx[i]); }
    return r;
}

static int rels_add(rels_t *r, mpz_t ax_b, mpz_t Qx, uint64_t lp1, uint64_t lp2) {
    if (r->count >= r->alloc) return -1;
    int i = r->count++;
    mpz_set(r->ax_b[i], ax_b); mpz_set(r->Qx[i], Qx);
    r->lp1[i] = lp1; r->lp2[i] = lp2;
    return i;
}

/* ==================== Pollard Rho (64-bit cofactor splitting) ==================== */
static int is_prime64(uint64_t n) {
    if (n < 2) return 0; if (n < 4) return 1;
    if (n % 2 == 0 || n % 3 == 0) return 0;
    uint64_t d = n - 1; int r = 0;
    while (!(d & 1)) { d >>= 1; r++; }
    uint64_t witnesses[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (int wi = 0; wi < 12; wi++) {
        uint64_t a = witnesses[wi];
        if (a >= n) continue;
        __uint128_t x = 1, base = a; uint64_t exp = d;
        while (exp) { if (exp & 1) x = x * base % n; base = base * base % n; exp >>= 1; }
        uint64_t xv = (uint64_t)x;
        if (xv == 1 || xv == n - 1) continue;
        int found = 0;
        for (int i = 0; i < r - 1; i++) {
            xv = (uint64_t)((__uint128_t)xv * xv % n);
            if (xv == n - 1) { found = 1; break; }
        }
        if (!found) return 0;
    }
    return 1;
}

static int split64(uint64_t n, uint64_t *f1, uint64_t *f2) {
    if (n < 4) return 0;
    if (n % 2 == 0) { *f1 = 2; *f2 = n / 2; return 1; }
    if (is_prime64(n)) return 0;
    for (uint64_t c = 1; c < 200; c++) {
        uint64_t x = 2, y = 2, p = 1, ys = 2, q = 1;
        int m = 128;
        do {
            x = y;
            for (int i = 0; i < m; i++) y = (uint64_t)((__uint128_t)y * y % n) + c;
            int k = 0;
            do {
                ys = y;
                int lim = (m - k < 32) ? m - k : 32;
                for (int i = 0; i < lim; i++) {
                    y = (uint64_t)((__uint128_t)y * y % n) + c;
                    uint64_t diff = x > y ? x - y : y - x;
                    q = (uint64_t)((__uint128_t)q * diff % n);
                }
                { uint64_t a = q, b = n; while (b) { uint64_t t = b; b = a % b; a = t; } p = a; }
                k += 32;
            } while (p == 1 && k < m);
            m *= 2;
        } while (p == 1);
        if (p == n) {
            do {
                ys = (uint64_t)((__uint128_t)ys * ys % n) + c;
                uint64_t diff = x > ys ? x - ys : ys - x;
                uint64_t a = diff, b = n;
                while (b) { uint64_t t = b; b = a % b; a = t; }
                p = a;
            } while (p == 1);
        }
        if (p != n && p != 1) { *f1 = p; *f2 = n / p; return 1; }
    }
    return 0;
}

/* ==================== GF(2) Gaussian Elimination ==================== */
typedef uint64_t gf2w;

static int gf2_solve(int nrels, int ncols, gf2w **rows, int fb_words, int id_words,
                     int ***deps_out, int **dlen_out, int max_deps) {
    int total_words = fb_words + id_words;
    int pivot = 0;
    for (int c = 0; c < ncols && pivot < nrels; c++) {
        int pr = -1;
        for (int r = pivot; r < nrels; r++)
            if ((rows[r][c/64] >> (c % 64)) & 1) { pr = r; break; }
        if (pr < 0) continue;
        if (pr != pivot) { gf2w *tmp = rows[pr]; rows[pr] = rows[pivot]; rows[pivot] = tmp; }
        for (int r = 0; r < nrels; r++) {
            if (r == pivot) continue;
            if ((rows[r][c/64] >> (c % 64)) & 1)
                for (int w = 0; w < total_words; w++) rows[r][w] ^= rows[pivot][w];
        }
        pivot++;
    }
    int nd = 0;
    *deps_out = malloc(max_deps * sizeof(int*));
    *dlen_out = malloc(max_deps * sizeof(int));
    for (int r = pivot; r < nrels && nd < max_deps; r++) {
        int zero = 1;
        for (int w = 0; w < fb_words && zero; w++) {
            gf2w mask = (w < fb_words - 1) ? ~0ULL
                      : (ncols % 64 == 0 ? ~0ULL : (1ULL << (ncols % 64)) - 1);
            if (rows[r][w] & mask) zero = 0;
        }
        if (!zero) continue;
        int *d = malloc(nrels * sizeof(int)); int dl = 0;
        for (int w = 0; w < id_words; w++) {
            gf2w bits = rows[r][fb_words + w];
            while (bits) {
                int bit = __builtin_ctzll(bits);
                int idx = w * 64 + bit;
                if (idx < nrels) d[dl++] = idx;
                bits &= bits - 1;
            }
        }
        if (dl > 0) { (*deps_out)[nd] = d; (*dlen_out)[nd] = dl; nd++; } else free(d);
    }
    return nd;
}

/* ==================== Parameters ==================== */
typedef struct {
    int fb_size, nblocks, lp_mult, extra;
    double thresh;
} params_t;

static params_t get_params(int bits) {
    if (bits <= 100) return (params_t){120,   1,  40,  30, 0.73};
    if (bits <= 110) return (params_t){170,   1,  50,  35, 0.74};
    if (bits <= 120) return (params_t){240,   2,  60,  45, 0.76};
    if (bits <= 130) return (params_t){340,   3,  70,  50, 0.77};
    if (bits <= 140) return (params_t){480,   4,  80,  60, 0.78};
    if (bits <= 150) return (params_t){650,   5,  90,  70, 0.79};
    if (bits <= 160) return (params_t){900,   7, 100,  80, 0.80};
    if (bits <= 170) return (params_t){1200, 10, 110,  90, 0.81};
    if (bits <= 180) return (params_t){1700, 14, 120, 100, 0.82};
    if (bits <= 190) return (params_t){2300, 18, 140, 120, 0.825};
    if (bits <= 200) return (params_t){3000, 24, 160, 140, 0.83};
    if (bits <= 210) return (params_t){4000, 30, 180, 160, 0.835};
    if (bits <= 220) return (params_t){5200, 38, 200, 180, 0.84};
    if (bits <= 230) return (params_t){6500, 46, 220, 200, 0.845};
    if (bits <= 240) return (params_t){7500, 56, 250, 220, 0.85};
    if (bits <= 250) return (params_t){9000, 66, 280, 250, 0.855};
    if (bits <= 260) return (params_t){11000,78, 320, 280, 0.86};
    if (bits <= 270) return (params_t){13000,92, 360, 320, 0.865};
    return              (params_t){16000,108, 400, 360, 0.87};
}

/* ==================== MAIN ==================== */
int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }
    clock_gettime(CLOCK_MONOTONIC, &g_start);

    mpz_t N, kN;
    mpz_inits(N, kN, NULL);
    mpz_set_str(N, argv[1], 10);

    int digits = (int)mpz_sizeinbase(N, 10);
    int bits   = (int)mpz_sizeinbase(N, 2);

    /* Trial division */
    for (unsigned long p = 2; p < 100000; p++) {
        if (p > 2 && p % 2 == 0) continue;
        if (mpz_divisible_ui_p(N, p)) { printf("%lu\n", p); return 0; }
    }
    { mpz_t sq; mpz_init(sq);
      if (mpz_perfect_square_p(N)) { mpz_sqrt(sq, N); gmp_printf("%Zd\n", sq); return 0; }
      mpz_clear(sq); }

    int mult = choose_multiplier(N);
    mpz_mul_ui(kN, N, mult);
    int kN_bits = (int)mpz_sizeinbase(kN, 2);
    params_t P = get_params(kN_bits);

    fb_t *fb = fb_create(kN, P.fb_size);
    int M = BLOCK_SIZE * P.nblocks;            /* half-sieve interval */
    int total_blocks = 2 * P.nblocks;          /* sieve from -M..+M */
    uint64_t lp_bound  = (uint64_t)fb->prime[fb->size - 1] * P.lp_mult;
    uint64_t dlp_max   = lp_bound * lp_bound;
    int target = fb->size + P.extra;

    double log2_Qmax = kN_bits / 2.0 + 0.5 + log2((double)M);
    int threshold = (int)(log2_Qmax * P.thresh) - 3;
    if (threshold < 20) threshold = 20;

    /* Bucket sieve threshold: primes > BLOCK_SIZE */
    int bucket_thresh = fb->size;
    for (int i = 0; i < fb->size; i++) {
        if (fb->prime[i] > (uint32_t)BLOCK_SIZE) { bucket_thresh = i; break; }
    }

    /*
     * Per-polynomial, per-block buckets for large primes.
     * We use BATCH separate sets of buckets (one per polynomial in the batch).
     */
    bucket_t *buckets[BATCH];
    for (int bi = 0; bi < BATCH; bi++) {
        buckets[bi] = calloc(total_blocks, sizeof(bucket_t));
        for (int i = 0; i < total_blocks; i++) {
            buckets[bi][i].alloc = BUCKET_ALLOC;
            buckets[bi][i].entries = malloc(BUCKET_ALLOC * sizeof(bucket_entry_t));
        }
    }

    fprintf(stderr,
            "HyperSIQS: %dd (%db), k=%d, FB=%d, M=%d, thresh=%d, LP=%lu [DLP], target=%d\n",
            digits, bits, mult, fb->size, M, threshold, lp_bound, target);

    /* BATCH sieve arrays (each BLOCK_SIZE bytes) */
    uint8_t *sieves[BATCH];
    for (int bi = 0; bi < BATCH; bi++) sieves[bi] = malloc(BLOCK_SIZE);

    rels_t *full = rels_create(MAX_FULL_RELS);
    rels_t *part = rels_create(MAX_PARTIAL);
    lp_hash_t *slp_ht = lp_create(MAX_PARTIAL);

    /* DLP pairwise matching hash */
    dp_e_t **dp_hash = calloc(DLP_HASH_SIZE, sizeof(dp_e_t*));
    dp_e_t *dp_pool  = calloc(MAX_PARTIAL * 2, sizeof(dp_e_t));
    int dp_pool_used = 0;

    /* Polynomial state */
    mpz_t a, b_val, c_val, B_vals[MAX_A_FACTORS];
    mpz_inits(a, b_val, c_val, NULL);
    for (int j = 0; j < MAX_A_FACTORS; j++) mpz_init(B_vals[j]);

    /* Per-batch polynomial b/c values */
    mpz_t batch_b[BATCH], batch_c[BATCH];
    for (int bi = 0; bi < BATCH; bi++) { mpz_init(batch_b[bi]); mpz_init(batch_c[bi]); }

    /*
     * Global sieve roots (updated via Gray code, tracks the "current" polynomial).
     * We also keep per-batch copies for the BATCH polynomials we're processing.
     */
    uint32_t *soln1 = malloc(fb->size * sizeof(uint32_t));
    uint32_t *soln2 = malloc(fb->size * sizeof(uint32_t));
    uint32_t *bs_soln1[BATCH], *bs_soln2[BATCH];
    for (int bi = 0; bi < BATCH; bi++) {
        bs_soln1[bi] = malloc(fb->size * sizeof(uint32_t));
        bs_soln2[bi] = malloc(fb->size * sizeof(uint32_t));
    }

    /* ainv[j][i] = 2 * (a/qj)^{-1} * B_j (mod p_i), for Gray code updates */
    uint32_t *ainv_data = malloc((size_t)MAX_A_FACTORS * fb->size * sizeof(uint32_t));
#define AINV(j,i) ainv_data[(j)*fb->size+(i)]

    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, SEED);

    mpz_t ax_b, Qx, aQx, residue, tmp, tmp2;
    mpz_inits(ax_b, Qx, aQx, residue, tmp, tmp2, NULL);

    int total_polys = 0, combined_slp = 0, combined_dlp = 0, dlp_found = 0;
    int a_idx[MAX_A_FACTORS], num_a_factors = 0;
    /* Track whether soln1/soln2 are initialized for this 'a' */
    int soln_initialized = 0;
    int prev_gray_global = 0; /* tracks Gray code index of last processed polynomial */

    /* ==================== Main Sieving Loop ==================== */
    while (full->count < target) {
        double t0 = elapsed();
        if (t0 > 275.0) {
            fprintf(stderr, "TIMEOUT at %.1fs with %d/%d rels\n", t0, full->count, target);
            break;
        }

        /* ========== Generate new 'a' ========== */
        {
            mpz_t tgt; mpz_init(tgt);
            mpz_mul_ui(tgt, kN, 2); mpz_sqrt(tgt, tgt); mpz_tdiv_q_ui(tgt, tgt, M);
            double log_tgt = mpz_sizeinbase(tgt, 2) * log(2.0);

            int lo = fb->size / 4, hi = 3 * fb->size / 4;
            if (lo < 2) lo = 2;
            if (hi <= lo + 5) hi = fb->size - 1;

            double avg_logp = 0; int cnt = 0;
            for (int i = lo; i < hi; i++) {
                if (fb->root[i] == 0) continue;
                avg_logp += log(fb->prime[i]); cnt++;
            }
            if (cnt == 0) { mpz_clear(tgt); break; }
            avg_logp /= cnt;

            int s = (int)(log_tgt / avg_logp + 0.5);
            if (s < 3) s = 3;
            if (s > MAX_A_FACTORS) s = MAX_A_FACTORS;
            if (s > hi - lo) s = hi - lo;
            num_a_factors = s;

            double best_ratio = 1e30;
            int best_idx[MAX_A_FACTORS];

            for (int att = 0; att < 50; att++) {
                mpz_set_ui(a, 1);
                int idx[MAX_A_FACTORS]; int ok = 1;
                for (int i = 0; i < s && ok; i++) {
                    int tries = 0, good;
                    do {
                        idx[i] = lo + (int)gmp_urandomm_ui(rng, hi - lo);
                        good = 1;
                        for (int j = 0; j < i; j++) if (idx[j] == idx[i]) { good = 0; break; }
                        if (fb->root[idx[i]] == 0) good = 0;
                        tries++;
                    } while (!good && tries < 100);
                    if (!good) { ok = 0; break; }
                    mpz_mul_ui(a, a, fb->prime[idx[i]]);
                }
                if (!ok) continue;
                double ratio;
                if (mpz_cmp(a, tgt) > 0) { mpz_tdiv_q(tmp, a, tgt); ratio = mpz_get_d(tmp); }
                else { mpz_tdiv_q(tmp, tgt, a); ratio = mpz_get_d(tmp); }
                if (ratio < best_ratio) { best_ratio = ratio; memcpy(best_idx, idx, s * sizeof(int)); }
                if (ratio < 1.5) break;
            }
            memcpy(a_idx, best_idx, s * sizeof(int));
            mpz_set_ui(a, 1);
            for (int i = 0; i < s; i++) mpz_mul_ui(a, a, fb->prime[a_idx[i]]);
            mpz_clear(tgt);
        }

        /* ========== Compute B_j values ========== */
        for (int j = 0; j < num_a_factors; j++) {
            int idx = a_idx[j];
            uint32_t qj = fb->prime[idx], rj = fb->root[idx];
            mpz_t a_q, mod_q, inv; mpz_inits(a_q, mod_q, inv, NULL);
            mpz_divexact_ui(a_q, a, qj);
            mpz_set_ui(mod_q, qj);
            mpz_invert(inv, a_q, mod_q);
            uint64_t iv = mpz_get_ui(inv);
            mpz_mul_ui(B_vals[j], a_q, (uint32_t)(((uint64_t)rj * iv) % qj));
            mpz_clears(a_q, mod_q, inv, NULL);
        }

        /* ========== Precompute ainv for Gray code updates ========== */
        for (int j = 0; j < num_a_factors; j++) {
            for (int i = 0; i < fb->size; i++) {
                uint32_t p = fb->prime[i];
                unsigned long am = mpz_fdiv_ui(a, p);
                if (am == 0 || fb->root[i] == 0) { AINV(j, i) = 0; continue; }
                uint32_t ai = mod_inv32((uint32_t)am, p);
                unsigned long Bm = mpz_fdiv_ui(B_vals[j], p);
                AINV(j, i) = (uint32_t)((uint64_t)2 * ai % p * Bm % p);
            }
        }

        /* ========== Gray code enumeration of b-values ========== */
        int num_b = 1 << (num_a_factors - 1);
        soln_initialized = 0;
        prev_gray_global = 0;

        /* ========== Process b-values in batches of BATCH ========== */
        for (int b_start = 0; b_start < num_b && full->count < target; b_start += BATCH) {
            int batch = BATCH;
            if (b_start + batch > num_b) batch = num_b - b_start;

            /*
             * For each polynomial in this batch, compute b value and sieve roots.
             *
             * Strategy:
             *  - For b_start==0 bi==0: compute from scratch (b = -sum B_j).
             *  - Otherwise: apply Gray code update from the previous polynomial.
             *
             * We maintain soln1/soln2 as a "cursor" that always reflects
             * the Gray code polynomial just before the current batch.
             */
            for (int bi = 0; bi < batch; bi++) {
                int b_idx = b_start + bi;
                int gray  = b_idx ^ (b_idx >> 1);

                if (!soln_initialized) {
                    /* Very first polynomial for this 'a': compute from scratch */
                    /* gray = 0, b = -sum(B_j) */
                    mpz_set_ui(b_val, 0);
                    for (int j = 0; j < num_a_factors; j++) mpz_sub(b_val, b_val, B_vals[j]);
                    if (mpz_sgn(b_val) <= 0) mpz_add(b_val, b_val, a);

                    for (int i = 0; i < fb->size; i++) {
                        uint32_t p = fb->prime[i];
                        unsigned long am = mpz_fdiv_ui(a, p);
                        if (am == 0 || fb->root[i] == 0) { soln1[i] = soln2[i] = 0xFFFFFFFF; continue; }
                        uint32_t ai = mod_inv32((uint32_t)am, p);
                        if (ai == 0) { soln1[i] = soln2[i] = 0xFFFFFFFF; continue; }
                        unsigned long bm = mpz_fdiv_ui(b_val, p);
                        uint32_t r = fb->root[i];
                        soln1[i] = (uint32_t)((uint64_t)ai * ((r + p - (uint32_t)bm) % p) % p);
                        soln2[i] = (uint32_t)((uint64_t)ai * ((p - r + p - (uint32_t)bm) % p) % p);
                    }
                    soln_initialized = 1;
                    prev_gray_global = 0;
                    /* Copy to batch slot 0 */
                    memcpy(bs_soln1[0], soln1, fb->size * sizeof(uint32_t));
                    memcpy(bs_soln2[0], soln2, fb->size * sizeof(uint32_t));
                    mpz_set(batch_b[0], b_val);
                } else if (bi == 0) {
                    /*
                     * First poly in a new batch (not the very first).
                     * Apply Gray code update from prev_gray_global to gray.
                     * soln1/soln2/b_val currently reflect prev_gray_global.
                     */
                    int changed = gray ^ prev_gray_global;
                    int j = __builtin_ctz(changed);
                    int sign = (gray >> j) & 1;

                    if (sign) mpz_addmul_ui(b_val, B_vals[j], 2);
                    else      mpz_submul_ui(b_val, B_vals[j], 2);

                    for (int i = 0; i < fb->size; i++) {
                        if (soln1[i] == 0xFFFFFFFF) continue;
                        uint32_t p = fb->prime[i];
                        uint32_t delta = AINV(j, i);
                        if (delta == 0) continue;
                        if (sign) {
                            soln1[i] = (soln1[i] >= delta) ? soln1[i] - delta : soln1[i] + p - delta;
                            soln2[i] = (soln2[i] >= delta) ? soln2[i] - delta : soln2[i] + p - delta;
                        } else {
                            soln1[i] += delta; if (soln1[i] >= p) soln1[i] -= p;
                            soln2[i] += delta; if (soln2[i] >= p) soln2[i] -= p;
                        }
                    }
                    prev_gray_global = gray;
                    memcpy(bs_soln1[0], soln1, fb->size * sizeof(uint32_t));
                    memcpy(bs_soln2[0], soln2, fb->size * sizeof(uint32_t));
                    mpz_set(batch_b[0], b_val);
                } else {
                    /*
                     * Subsequent poly within the same batch.
                     * Apply Gray code update from (b_idx-1) to b_idx,
                     * based on the previous batch slot.
                     */
                    int prev_gray = (b_idx - 1) ^ ((b_idx - 1) >> 1);
                    int changed = gray ^ prev_gray;
                    int j = __builtin_ctz(changed);
                    int sign = (gray >> j) & 1;

                    memcpy(bs_soln1[bi], bs_soln1[bi-1], fb->size * sizeof(uint32_t));
                    memcpy(bs_soln2[bi], bs_soln2[bi-1], fb->size * sizeof(uint32_t));
                    mpz_set(batch_b[bi], batch_b[bi-1]);

                    if (sign) mpz_addmul_ui(batch_b[bi], B_vals[j], 2);
                    else      mpz_submul_ui(batch_b[bi], B_vals[j], 2);

                    uint32_t *s1 = bs_soln1[bi];
                    uint32_t *s2 = bs_soln2[bi];
                    for (int i = 0; i < fb->size; i++) {
                        if (s1[i] == 0xFFFFFFFF) continue;
                        uint32_t p = fb->prime[i];
                        uint32_t delta = AINV(j, i);
                        if (delta == 0) continue;
                        if (sign) {
                            s1[i] = (s1[i] >= delta) ? s1[i] - delta : s1[i] + p - delta;
                            s2[i] = (s2[i] >= delta) ? s2[i] - delta : s2[i] + p - delta;
                        } else {
                            s1[i] += delta; if (s1[i] >= p) s1[i] -= p;
                            s2[i] += delta; if (s2[i] >= p) s2[i] -= p;
                        }
                    }
                }

                /* Compute c = (b^2 - kN) / a */
                mpz_mul(batch_c[bi], batch_b[bi], batch_b[bi]);
                mpz_sub(batch_c[bi], batch_c[bi], kN);
                mpz_divexact(batch_c[bi], batch_c[bi], a);
            }

            /* After setting up this batch, advance soln1/soln2/b_val/prev_gray_global
             * to reflect the last polynomial in the batch, ready for next b_start. */
            {
                int last_b_idx = b_start + batch - 1;
                int last_gray  = last_b_idx ^ (last_b_idx >> 1);
                if (last_gray != prev_gray_global) {
                    /* Apply remaining Gray code steps to bring cursor up to date */
                    /* Actually: the cursor should already be at b_start's gray after bi==0.
                     * For bi>0 updates we kept them in batch slots only.
                     * We need to advance the cursor from b_start gray to last_b_idx gray.
                     */
                    /* Re-sync: cursor is at b_start's gray (set during bi==0 above).
                     * We need to walk from b_start to b_start+batch-1. */
                    for (int bi = 1; bi < batch; bi++) {
                        int cur_b_idx  = b_start + bi;
                        int cur_gray   = cur_b_idx ^ (cur_b_idx >> 1);
                        int prev_g     = (cur_b_idx-1) ^ ((cur_b_idx-1) >> 1);
                        int changed    = cur_gray ^ prev_g;
                        int j          = __builtin_ctz(changed);
                        int sign       = (cur_gray >> j) & 1;

                        if (sign) mpz_addmul_ui(b_val, B_vals[j], 2);
                        else      mpz_submul_ui(b_val, B_vals[j], 2);

                        for (int i = 0; i < fb->size; i++) {
                            if (soln1[i] == 0xFFFFFFFF) continue;
                            uint32_t p = fb->prime[i];
                            uint32_t delta = AINV(j, i);
                            if (delta == 0) continue;
                            if (sign) {
                                soln1[i] = (soln1[i] >= delta) ? soln1[i] - delta : soln1[i] + p - delta;
                                soln2[i] = (soln2[i] >= delta) ? soln2[i] - delta : soln2[i] + p - delta;
                            } else {
                                soln1[i] += delta; if (soln1[i] >= p) soln1[i] -= p;
                                soln2[i] += delta; if (soln2[i] >= p) soln2[i] -= p;
                            }
                        }
                    }
                    prev_gray_global = last_gray;
                }
            }

            total_polys += batch;

            /* ========== Fill buckets for large primes (per polynomial) ========== */
            for (int bi = 0; bi < batch; bi++) {
                for (int k = 0; k < total_blocks; k++) buckets[bi][k].count = 0;
                uint32_t *s1 = bs_soln1[bi], *s2 = bs_soln2[bi];
                for (int i = bucket_thresh; i < fb->size; i++) {
                    if (s1[i] == 0xFFFFFFFF) continue;
                    uint32_t p = fb->prime[i]; uint8_t lp = fb->logp[i];
                    for (int root = 0; root < 2; root++) {
                        uint32_t sv = (root == 0) ? s1[i] : s2[i];
                        if (root == 1 && s1[i] == s2[i]) continue;
                        int64_t pos = ((int64_t)sv - (-(int64_t)M)) % (int64_t)p;
                        if (pos < 0) pos += p;
                        int64_t x = -(int64_t)M + pos;
                        while (x < (int64_t)M) {
                            int blk = (int)((x + M) / BLOCK_SIZE);
                            if (blk >= 0 && blk < total_blocks) {
                                uint32_t bp = (uint32_t)((x + M) - (int64_t)blk * BLOCK_SIZE);
                                if (bp < (uint32_t)BLOCK_SIZE) bucket_add(&buckets[bi][blk], bp, lp);
                            }
                            x += p;
                        }
                    }
                }
            }

            /* ========== KEY: Batch sieve over each block ========== */
            /*
             * For each block: iterate through FB primes ONCE, update ALL BATCH sieves.
             * This amortizes the outer FB prime loop overhead by factor BATCH.
             */
            for (int blk = 0; blk < total_blocks; blk++) {
                int64_t bbase = (int64_t)blk * BLOCK_SIZE - (int64_t)M;

                /* Zero all BATCH sieve arrays */
                for (int bi = 0; bi < batch; bi++)
                    memset(sieves[bi], 0, BLOCK_SIZE);

                /*
                 * Core batch sieve loop:
                 * For each FB prime p, compute offsets for ALL batch polynomials
                 * and update their respective sieve arrays.
                 * The outer prime loop runs ONCE; per-prime work is BATCH times.
                 */
                for (int i = 1; i < bucket_thresh; i++) {
                    uint32_t p = fb->prime[i];
                    if (p < 4) continue;
                    uint8_t lp = fb->logp[i];

                    /* Unrolled for BATCH=4 */
                    uint32_t *s10 = bs_soln1[0], *s20 = bs_soln2[0];
                    if (batch > 0 && s10[i] != 0xFFFFFFFF) {
                        int64_t off1 = ((int64_t)s10[i] - bbase) % (int64_t)p;
                        if (off1 < 0) off1 += p;
                        uint8_t *sv = sieves[0];
                        for (int64_t j = off1; j < BLOCK_SIZE; j += p) sv[j] += lp;
                        if (s10[i] != s20[i]) {
                            int64_t off2 = ((int64_t)s20[i] - bbase) % (int64_t)p;
                            if (off2 < 0) off2 += p;
                            for (int64_t j = off2; j < BLOCK_SIZE; j += p) sv[j] += lp;
                        }
                    }
                    if (batch > 1) {
                        uint32_t *s1 = bs_soln1[1], *s2 = bs_soln2[1];
                        if (s1[i] != 0xFFFFFFFF) {
                            int64_t off1 = ((int64_t)s1[i] - bbase) % (int64_t)p;
                            if (off1 < 0) off1 += p;
                            uint8_t *sv = sieves[1];
                            for (int64_t j = off1; j < BLOCK_SIZE; j += p) sv[j] += lp;
                            if (s1[i] != s2[i]) {
                                int64_t off2 = ((int64_t)s2[i] - bbase) % (int64_t)p;
                                if (off2 < 0) off2 += p;
                                for (int64_t j = off2; j < BLOCK_SIZE; j += p) sv[j] += lp;
                            }
                        }
                    }
                    if (batch > 2) {
                        uint32_t *s1 = bs_soln1[2], *s2 = bs_soln2[2];
                        if (s1[i] != 0xFFFFFFFF) {
                            int64_t off1 = ((int64_t)s1[i] - bbase) % (int64_t)p;
                            if (off1 < 0) off1 += p;
                            uint8_t *sv = sieves[2];
                            for (int64_t j = off1; j < BLOCK_SIZE; j += p) sv[j] += lp;
                            if (s1[i] != s2[i]) {
                                int64_t off2 = ((int64_t)s2[i] - bbase) % (int64_t)p;
                                if (off2 < 0) off2 += p;
                                for (int64_t j = off2; j < BLOCK_SIZE; j += p) sv[j] += lp;
                            }
                        }
                    }
                    if (batch > 3) {
                        uint32_t *s1 = bs_soln1[3], *s2 = bs_soln2[3];
                        if (s1[i] != 0xFFFFFFFF) {
                            int64_t off1 = ((int64_t)s1[i] - bbase) % (int64_t)p;
                            if (off1 < 0) off1 += p;
                            uint8_t *sv = sieves[3];
                            for (int64_t j = off1; j < BLOCK_SIZE; j += p) sv[j] += lp;
                            if (s1[i] != s2[i]) {
                                int64_t off2 = ((int64_t)s2[i] - bbase) % (int64_t)p;
                                if (off2 < 0) off2 += p;
                                for (int64_t j = off2; j < BLOCK_SIZE; j += p) sv[j] += lp;
                            }
                        }
                    }
                }

                /* Apply bucket entries for each polynomial */
                for (int bi = 0; bi < batch; bi++) {
                    bucket_t *bkt = &buckets[bi][blk];
                    uint8_t *sv = sieves[bi];
                    int cnt = bkt->count;
                    bucket_entry_t *ent = bkt->entries;
                    for (int e = 0; e < cnt; e++) sv[ent[e].pos] += ent[e].logp;
                }

                /* ========== Scan ALL BATCH sieves for candidates ========== */
                for (int bi = 0; bi < batch; bi++) {
                    uint8_t *sv = sieves[bi];
                    uint32_t *s1 = bs_soln1[bi];
                    uint32_t *s2 = bs_soln2[bi];

                    for (int j = 0; j < BLOCK_SIZE; j++) {
                        if (sv[j] < (uint8_t)threshold) continue;
                        int64_t x = bbase + j;
                        if (x == 0) continue;

                        /* Compute Q(x) = ax^2 + 2bx + c and ax+b */
                        mpz_set_si(tmp, (long)x);
                        mpz_mul_si(ax_b, a, (long)x);
                        mpz_add(ax_b, ax_b, batch_b[bi]);

                        mpz_mul(Qx, tmp, tmp); mpz_mul(Qx, Qx, a);
                        mpz_mul(tmp2, batch_b[bi], tmp); mpz_addmul_ui(Qx, tmp2, 2);
                        mpz_add(Qx, Qx, batch_c[bi]);

                        if (mpz_sgn(Qx) == 0) continue;
                        mpz_abs(residue, Qx);

                        /* Divide out 2 */
                        while (mpz_even_p(residue)) mpz_tdiv_q_2exp(residue, residue, 1);

                        /* Sieve-informed trial division:
                         * only test primes p where x ≡ soln1 or soln2 (mod p) */
                        for (int i = 1; i < fb->size; i++) {
                            if (s1[i] == 0xFFFFFFFF) continue;
                            uint32_t p = fb->prime[i];
                            int64_t xmod = ((x % (int64_t)p) + p) % p;
                            if (xmod != (int64_t)s1[i] && xmod != (int64_t)s2[i]) continue;
                            while (mpz_divisible_ui_p(residue, p))
                                mpz_divexact_ui(residue, residue, p);
                        }
                        /* a-factor primes (marked 0xFFFF) - must divide Q(x) */
                        for (int i = 0; i < num_a_factors; i++) {
                            uint32_t p = fb->prime[a_idx[i]];
                            while (mpz_divisible_ui_p(residue, p))
                                mpz_divexact_ui(residue, residue, p);
                        }

                        mpz_mul(aQx, Qx, a);

                        if (mpz_cmp_ui(residue, 1) == 0) {
                            /* Full relation */
                            rels_add(full, ax_b, aQx, 0, 0);
                        } else if (mpz_fits_ulong_p(residue)) {
                            uint64_t cof = mpz_get_ui(residue);

                            if (cof <= lp_bound) {
                                /* SLP: single large prime */
                                int match = lp_find(slp_ht, cof);
                                if (match >= 0) {
                                    mpz_mul(tmp, ax_b, part->ax_b[match]);
                                    mpz_mod(tmp, tmp, N);
                                    mpz_mul(tmp2, aQx, part->Qx[match]);
                                    rels_add(full, tmp, tmp2, 0, 0);
                                    combined_slp++;
                                } else {
                                    int pi = rels_add(part, ax_b, aQx, cof, 0);
                                    if (pi >= 0) lp_insert(slp_ht, cof, pi);
                                }
                            } else if (cof <= dlp_max) {
                                /* DLP: split cofactor via Pollard rho */
                                uint64_t f1, f2;
                                if (split64(cof, &f1, &f2)) {
                                    if (f1 > f2) { uint64_t t = f1; f1 = f2; f2 = t; }
                                    if (f1 <= lp_bound && f2 <= lp_bound) {
                                        int pi = rels_add(part, ax_b, aQx, f1, f2);
                                        dlp_found++;
                                        if (pi >= 0) {
                                            /* Look for exact match: another partial with same {f1,f2} */
                                            uint32_t h1 = (uint32_t)((f1 * 0x9E3779B97F4A7C15ULL) >> (64 - DLP_HASH_BITS));
                                            uint32_t h2 = (uint32_t)((f2 * 0x9E3779B97F4A7C15ULL) >> (64 - DLP_HASH_BITS));

                                            int matched = 0;
                                            for (dp_e_t *e = dp_hash[h1]; e; e = e->next) {
                                                if (e->lp == f1) {
                                                    int oi = e->rel_idx;
                                                    if ((part->lp1[oi] == f1 || part->lp2[oi] == f1) &&
                                                        (part->lp1[oi] == f2 || part->lp2[oi] == f2)) {
                                                        mpz_mul(tmp, ax_b, part->ax_b[oi]);
                                                        mpz_mod(tmp, tmp, N);
                                                        mpz_mul(tmp2, aQx, part->Qx[oi]);
                                                        rels_add(full, tmp, tmp2, 0, 0);
                                                        combined_dlp++;
                                                        matched = 1;
                                                        break;
                                                    }
                                                }
                                            }
                                            if (!matched) {
                                                for (dp_e_t *e = dp_hash[h2]; e; e = e->next) {
                                                    if (e->lp == f2) {
                                                        int oi = e->rel_idx;
                                                        if ((part->lp1[oi] == f1 || part->lp2[oi] == f1) &&
                                                            (part->lp1[oi] == f2 || part->lp2[oi] == f2)) {
                                                            mpz_mul(tmp, ax_b, part->ax_b[oi]);
                                                            mpz_mod(tmp, tmp, N);
                                                            mpz_mul(tmp2, aQx, part->Qx[oi]);
                                                            rels_add(full, tmp, tmp2, 0, 0);
                                                            combined_dlp++;
                                                            break;
                                                        }
                                                    }
                                                }
                                            }

                                            /* Insert this DLP partial into hash by both LPs */
                                            if (dp_pool_used + 2 <= (int)(MAX_PARTIAL * 2)) {
                                                dp_e_t *e1 = &dp_pool[dp_pool_used++];
                                                e1->lp = f1; e1->rel_idx = pi;
                                                e1->next = dp_hash[h1]; dp_hash[h1] = e1;
                                                dp_e_t *e2 = &dp_pool[dp_pool_used++];
                                                e2->lp = f2; e2->rel_idx = pi;
                                                e2->next = dp_hash[h2]; dp_hash[h2] = e2;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    } /* j (sieve position) loop */
                } /* bi (polynomial in batch) loop */
            } /* blk loop */

            /* Progress reporting */
            if (total_polys % 500 == 0) {
                double t = elapsed();
                if (t > 275.0) break;
                if (total_polys % 2000 == 0)
                    fprintf(stderr,
                            "  poly=%d rels=%d/%d (full=%d slp=%d dlp=%d) part=%d t=%.1fs\n",
                            total_polys, full->count, target,
                            full->count - combined_slp - combined_dlp,
                            combined_slp, combined_dlp, part->count, t);
            }
        } /* b_start loop */
    } /* while full->count < target */

    double sieve_time = elapsed();
    fprintf(stderr,
            "Sieve done: %d rels in %.2fs (%d polys, SLP=%d, DLP=%d/%d found)\n",
            full->count, sieve_time, total_polys, combined_slp, combined_dlp, dlp_found);

    if (full->count < fb->size + 2) {
        fprintf(stderr, "FAIL: not enough relations (%d < %d)\n", full->count, fb->size + 2);
        return 1;
    }

    /* ==================== Linear Algebra ==================== */
    int nrels = full->count;
    /*
     * Column layout:
     *   col 0       = sign bit (negative Q)
     *   col 1       = prime 2
     *   col i+1     = fb->prime[i]  for i = 1..fb->size-1
     * Total: fb->size + 1 columns
     */
    int ncols    = fb->size + 1;
    int fb_words = (ncols + 63) / 64;
    int id_words = (nrels + 63) / 64;
    int total_words = fb_words + id_words;

    fprintf(stderr, "Building matrix: %d x %d\n", nrels, ncols);

    gf2w **rows = malloc(nrels * sizeof(gf2w*));
    for (int r = 0; r < nrels; r++) {
        rows[r] = calloc(total_words, sizeof(gf2w));
        /* Identity part: bit r marks this relation */
        rows[r][fb_words + r / 64] |= (1ULL << (r % 64));

        mpz_abs(residue, full->Qx[r]);
        /* Sign bit */
        if (mpz_sgn(full->Qx[r]) < 0) rows[r][0] |= 1ULL;

        /* Prime 2 parity */
        int exp2 = 0;
        while (mpz_even_p(residue)) { mpz_tdiv_q_2exp(residue, residue, 1); exp2++; }
        if (exp2 & 1) rows[r][0] |= (1ULL << 1);

        /* Odd FB primes */
        for (int i = 1; i < fb->size; i++) {
            uint32_t p = fb->prime[i]; int exp = 0;
            while (mpz_divisible_ui_p(residue, p)) { mpz_divexact_ui(residue, residue, p); exp++; }
            if (exp & 1) rows[r][(i + 1) / 64] |= (1ULL << ((i + 1) % 64));
        }
    }

    int **deps; int *dlen;
    int ndeps = gf2_solve(nrels, ncols, rows, fb_words, id_words, &deps, &dlen, 128);
    fprintf(stderr, "LA: %d deps from %dx%d in %.2fs\n", ndeps, nrels, ncols, elapsed());

    if (ndeps == 0) { fprintf(stderr, "FAIL: no dependencies\n"); return 1; }

    /* ==================== Square Root ==================== */
    mpz_t X, Y, g;
    mpz_inits(X, Y, g, NULL);
    int factored = 0;

    for (int d = 0; d < ndeps && !factored; d++) {
        /*
         * X = product of (ax+b) mod N
         * Y = sqrt(product of |Q(x)|) mod N
         *
         * We factor out the product of |Q| step by step:
         *  - extract power of 2
         *  - extract each FB prime
         *  - verify remainder is 1 (should be, by construction of GF2 null space)
         */
        mpz_set_ui(X, 1);
        for (int k = 0; k < dlen[d]; k++) {
            int ri = deps[d][k];
            mpz_mul(X, X, full->ax_b[ri]);
            mpz_mod(X, X, N);
        }

        mpz_set_ui(Y, 1);
        mpz_t prod; mpz_init(prod); mpz_set_ui(prod, 1);
        for (int k = 0; k < dlen[d]; k++) {
            mpz_abs(tmp, full->Qx[deps[d][k]]);
            mpz_mul(prod, prod, tmp);
        }

        /* Power of 2 */
        int e2 = 0;
        while (mpz_even_p(prod)) { mpz_tdiv_q_2exp(prod, prod, 1); e2++; }
        if (e2 & 1) { mpz_clear(prod); continue; }
        if (e2 / 2 > 0) {
            mpz_set_ui(tmp, 2); mpz_powm_ui(tmp, tmp, e2 / 2, N);
            mpz_mul(Y, Y, tmp); mpz_mod(Y, Y, N);
        }

        /* FB odd primes */
        int valid = 1;
        for (int i = 1; i < fb->size && valid; i++) {
            uint32_t p = fb->prime[i]; int exp = 0;
            while (mpz_divisible_ui_p(prod, p)) { mpz_divexact_ui(prod, prod, p); exp++; }
            if (exp & 1) { valid = 0; break; }
            if (exp / 2 > 0) {
                mpz_set_ui(tmp, p); mpz_powm_ui(tmp, tmp, exp / 2, N);
                mpz_mul(Y, Y, tmp); mpz_mod(Y, Y, N);
            }
        }
        if (!valid) { mpz_clear(prod); continue; }

        /* Remainder should be 1 */
        if (mpz_cmp_ui(prod, 1) != 0) {
            if (mpz_perfect_square_p(prod)) {
                mpz_sqrt(tmp, prod); mpz_mod(tmp, tmp, N);
                mpz_mul(Y, Y, tmp); mpz_mod(Y, Y, N);
            } else { mpz_clear(prod); continue; }
        }
        mpz_clear(prod);

        mpz_sub(tmp, X, Y); mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            gmp_printf("%Zd\n", g);
            fprintf(stderr, "HyperSIQS: factored in %.3fs\n", elapsed());
            factored = 1;
        } else {
            mpz_add(tmp, X, Y); mpz_gcd(g, tmp, N);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                gmp_printf("%Zd\n", g);
                fprintf(stderr, "HyperSIQS: factored in %.3fs\n", elapsed());
                factored = 1;
            }
        }
    }

    if (!factored) {
        fprintf(stderr, "FAIL: no factor found from %d dependencies\n", ndeps);
        return 1;
    }

    fprintf(stderr, "Total: %.3fs\n", elapsed());

    /* Cleanup */
    for (int r = 0; r < nrels; r++) free(rows[r]);
    free(rows);
    mpz_clears(N, kN, a, b_val, c_val, ax_b, Qx, aQx, residue, tmp, tmp2, X, Y, g, NULL);
    for (int j = 0; j < MAX_A_FACTORS; j++) mpz_clear(B_vals[j]);
    for (int bi = 0; bi < BATCH; bi++) { mpz_clear(batch_b[bi]); mpz_clear(batch_c[bi]); }
    gmp_randclear(rng);

    return 0;
}
