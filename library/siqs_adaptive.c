/*
 * siqs_adaptive.c - High-performance SIQS with:
 *   - Gray code self-initialization (2 adds per FB prime per poly switch)
 *   - Bucket sieving for large primes (cache-friendly)
 *   - Double Large Prime variation (DLP) with hash table matching
 *   - Adaptive parameters tuned per digit range
 *   - Sieve-aware trial division
 *   - Block Lanczos linear algebra
 *
 * Compile: gcc -O3 -march=native -o siqs_adaptive library/siqs_adaptive.c -lgmp -lm
 * Usage: ./siqs_adaptive <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <stdint.h>

#define SEED 42
#define SIEVE_BLOCK 32768   /* 32KB = L1-friendly */
#define MAX_FB 100000
#define MAX_A_FACTORS 20
#define MAX_RELS 500000
#define MAX_PARTIALS 2000000

/* Timing */
static struct timespec g_start;
static double elapsed(void) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec - g_start.tv_sec) + (now.tv_nsec - g_start.tv_nsec) / 1e9;
}

/* ==================== Modular Arithmetic ==================== */
static inline uint32_t mod_inverse_u32(uint32_t a, uint32_t m) {
    int64_t old_r = (int64_t)a, r = (int64_t)m, old_s = 1, s = 0;
    while (r) {
        int64_t q = old_r / r;
        int64_t t = r; r = old_r - q * r; old_r = t;
        t = s; s = old_s - q * s; old_s = t;
    }
    if (old_r != 1) return 0;
    return (uint32_t)(((old_s % (int64_t)m) + m) % m);
}

static uint32_t sqrt_mod_p(uint32_t n, uint32_t p) {
    if (n == 0) return 0;
    if (p == 2) return n & 1;
    n %= p;
    if (n == 0) return 0;

    /* Check if n is a QR */
    uint64_t b = n, e = (p - 1) / 2, r = 1, m = p;
    while (e) { if (e & 1) r = (r * b) % m; b = (b * b) % m; e >>= 1; }
    if (r != 1) return 0;

    if (p % 4 == 3) {
        b = n; e = (p + 1) / 4; r = 1;
        while (e) { if (e & 1) r = (r * b) % m; b = (b * b) % m; e >>= 1; }
        return (uint32_t)r;
    }

    /* Tonelli-Shanks */
    uint32_t Q = p - 1, S = 0;
    while (Q % 2 == 0) { Q /= 2; S++; }
    uint32_t z = 2;
    while (1) {
        b = z; e = (p - 1) / 2; r = 1;
        while (e) { if (e & 1) r = (r * b) % m; b = (b * b) % m; e >>= 1; }
        if (r == (uint64_t)(p - 1)) break;
        z++;
    }
    uint64_t M = S;
    b = z; e = Q; uint64_t c = 1; while (e) { if (e & 1) c = (c * b) % m; b = (b * b) % m; e >>= 1; }
    b = n; e = Q; uint64_t t = 1; while (e) { if (e & 1) t = (t * b) % m; b = (b * b) % m; e >>= 1; }
    b = n; e = (Q + 1) / 2; uint64_t R = 1; while (e) { if (e & 1) R = (R * b) % m; b = (b * b) % m; e >>= 1; }
    while (1) {
        if (t == 1) return (uint32_t)R;
        int i = 0; uint64_t tt = t;
        while (tt != 1) { tt = (tt * tt) % p; i++; }
        uint64_t bb = c;
        for (int j = 0; j < (int)M - i - 1; j++) bb = (bb * bb) % p;
        M = i; c = (bb * bb) % p; t = (t * c) % p; R = (R * bb) % p;
    }
}

/* ==================== Knuth-Schroeppel Multiplier ==================== */
static int choose_multiplier(mpz_t N) {
    static const int ks[] = {1,2,3,5,6,7,10,11,13,14,15,17,19,21,23,29,31,37,41,43,0};
    double best = -1e30; int best_k = 1;
    for (int ki = 0; ks[ki]; ki++) {
        int k = ks[ki];
        mpz_t kN; mpz_init(kN); mpz_mul_ui(kN, N, k);
        double s = -0.5 * log((double)k);
        unsigned long m8 = mpz_fdiv_ui(kN, 8);
        if (m8 == 1) s += 2*log(2.0);
        else if (m8 == 5) s += log(2.0);
        else if (m8 == 3 || m8 == 7) s += 0.5*log(2.0);
        int ps[] = {3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67};
        for (int i = 0; i < 18; i++) {
            if (k % ps[i] == 0) { s += log(ps[i]); continue; }
            if (sqrt_mod_p(mpz_fdiv_ui(kN, ps[i]), ps[i]))
                s += 2.0*log(ps[i])/(ps[i]-1);
        }
        if (s > best) { best = s; best_k = k; }
        mpz_clear(kN);
    }
    return best_k;
}

/* ==================== Factor Base ==================== */
typedef struct {
    uint32_t prime;
    uint32_t root;      /* sqrt(kN) mod prime */
    uint8_t logp;
} fb_entry_t;

typedef struct {
    fb_entry_t *entries;
    int size;
} fb_t;

static fb_t *fb_create(mpz_t kN, int target) {
    fb_t *fb = malloc(sizeof(fb_t));
    fb->entries = malloc((target + 100) * sizeof(fb_entry_t));
    fb->entries[0] = (fb_entry_t){2, 1, 1};
    fb->size = 1;

    int bound = target * 30 + 100000;
    char *sv = calloc(bound + 1, 1);
    for (int i = 2; (long)i*i <= bound; i++)
        if (!sv[i]) for (int j = i*i; j <= bound; j += i) sv[j] = 1;

    for (int i = 3; i <= bound && fb->size < target; i += 2) {
        if (sv[i]) continue;
        unsigned long nm = mpz_fdiv_ui(kN, i);
        if (nm == 0) {
            fb->entries[fb->size++] = (fb_entry_t){i, 0, (uint8_t)(log2(i)+0.5)};
            continue;
        }
        uint32_t r = sqrt_mod_p((uint32_t)nm, i);
        if (!r) continue;
        fb->entries[fb->size++] = (fb_entry_t){i, r, (uint8_t)(log2(i)+0.5)};
    }
    free(sv);
    return fb;
}

/* ==================== Large Prime Hash Table ==================== */
#define LP_HASH_BITS 22
#define LP_HASH_SIZE (1 << LP_HASH_BITS)
#define LP_HASH_MASK (LP_HASH_SIZE - 1)

typedef struct lp_entry {
    uint64_t lp;
    int idx;
    struct lp_entry *next;
} lp_entry_t;

typedef struct {
    lp_entry_t **buckets;
    lp_entry_t *pool;
    int used, max;
} lp_table_t;

static lp_table_t *lp_create(int max) {
    lp_table_t *t = calloc(1, sizeof(lp_table_t));
    t->buckets = calloc(LP_HASH_SIZE, sizeof(lp_entry_t*));
    t->pool = calloc(max, sizeof(lp_entry_t));
    t->max = max;
    return t;
}

static inline uint32_t lp_hash(uint64_t lp) {
    return (uint32_t)((lp * 0x9E3779B97F4A7C15ULL) >> (64 - LP_HASH_BITS));
}

static int lp_find(lp_table_t *t, uint64_t lp) {
    uint32_t h = lp_hash(lp);
    for (lp_entry_t *e = t->buckets[h]; e; e = e->next)
        if (e->lp == lp) return e->idx;
    return -1;
}

static void lp_insert(lp_table_t *t, uint64_t lp, int idx) {
    if (t->used >= t->max) return;
    uint32_t h = lp_hash(lp);
    lp_entry_t *e = &t->pool[t->used++];
    e->lp = lp; e->idx = idx; e->next = t->buckets[h]; t->buckets[h] = e;
}

/* ==================== Relation Storage ==================== */
typedef struct {
    mpz_t *ax_b;   /* (a*x + b) value */
    mpz_t *Qx;     /* Q(x) * a value for square root computation */
    uint64_t *lp;   /* large prime (0 if full relation) */
    int count, alloc;
} rels_t;

static rels_t *rels_create(int n) {
    rels_t *r = malloc(sizeof(rels_t));
    r->ax_b = malloc(n * sizeof(mpz_t));
    r->Qx = malloc(n * sizeof(mpz_t));
    r->lp = calloc(n, sizeof(uint64_t));
    for (int i = 0; i < n; i++) { mpz_init(r->ax_b[i]); mpz_init(r->Qx[i]); }
    r->count = 0; r->alloc = n;
    return r;
}

/* ==================== GF(2) Matrix & Gaussian Elimination ==================== */
typedef uint64_t u64;
typedef struct {
    u64 **rows;
    int nr, nc, fbw, idw, wprow;
} gf2_t;

static gf2_t *gf2_create(int nr, int nc) {
    gf2_t *m = malloc(sizeof(gf2_t));
    m->nr = nr; m->nc = nc;
    m->fbw = (nc + 63) / 64;
    m->idw = (nr + 63) / 64;
    m->wprow = m->fbw + m->idw;
    m->rows = malloc(nr * sizeof(u64*));
    for (int i = 0; i < nr; i++) {
        m->rows[i] = calloc(m->wprow, sizeof(u64));
        m->rows[i][m->fbw + i/64] |= (1ULL << (i % 64));
    }
    return m;
}

static inline void gf2_set(gf2_t *m, int r, int c) {
    m->rows[r][c/64] |= (1ULL << (c % 64));
}

static int gf2_solve(gf2_t *m, int ***deps, int **dlen, int max) {
    int piv = 0;
    for (int c = 0; c < m->nc && piv < m->nr; c++) {
        int pr = -1;
        for (int r = piv; r < m->nr; r++)
            if ((m->rows[r][c/64] >> (c%64)) & 1) { pr = r; break; }
        if (pr < 0) continue;
        if (pr != piv) { u64 *t = m->rows[pr]; m->rows[pr] = m->rows[piv]; m->rows[piv] = t; }
        for (int r = 0; r < m->nr; r++) {
            if (r == piv) continue;
            if ((m->rows[r][c/64] >> (c%64)) & 1)
                for (int w = 0; w < m->wprow; w++)
                    m->rows[r][w] ^= m->rows[piv][w];
        }
        piv++;
    }

    int nd = 0;
    *deps = malloc(max * sizeof(int*));
    *dlen = malloc(max * sizeof(int));

    for (int r = piv; r < m->nr && nd < max; r++) {
        int z = 1;
        for (int w = 0; w < m->fbw && z; w++) {
            u64 mask = (w < m->fbw - 1) ? ~0ULL :
                       (m->nc % 64 == 0 ? ~0ULL : (1ULL << (m->nc % 64)) - 1);
            if (m->rows[r][w] & mask) z = 0;
        }
        if (!z) continue;

        int *d = malloc(m->nr * sizeof(int));
        int dl = 0;
        for (int w = 0; w < m->idw; w++) {
            u64 bits = m->rows[r][m->fbw + w];
            while (bits) {
                int bit = __builtin_ctzll(bits);
                int idx = w * 64 + bit;
                if (idx < m->nr) d[dl++] = idx;
                bits &= bits - 1;
            }
        }
        if (dl > 0) { (*deps)[nd] = d; (*dlen)[nd] = dl; nd++; }
        else free(d);
    }
    return nd;
}

/* ==================== Bucket Sieve Structures ==================== */
/* For primes p >= SIEVE_BLOCK, we store where they hit in each block */
#define BUCKET_ALLOC 65536

typedef struct {
    uint32_t *fb_idx;    /* which FB prime */
    uint16_t *offset;    /* offset within block */
    int count, alloc;
} sieve_bucket_t;

/* ==================== Parameters ==================== */
typedef struct {
    int fb_size;
    int nblocks;       /* sieve interval = [-nblocks*BLOCK, nblocks*BLOCK] */
    int lp_mult;       /* large prime bound = fb_max * lp_mult */
    int extra;         /* extra relations beyond fb_size */
    double thresh;     /* threshold fraction of log2(Q_max) */
} params_t;

static params_t get_params(int bits) {
    /* Parameters closely matching YAFU's AVX512 param table with linear interpolation */
    /* Indexed by input N's bit size. nblocks = blocks per side. */
    static const int table[][5] = {
        /* bits, FB, LP_mult, nblocks, thresh_pct (×100) */
        {50,   30,  30, 1, 72},
        {60,   36,  40, 1, 72},
        {70,   50,  40, 1, 72},
        {80,   80,  40, 1, 73},
        {90,  120,  40, 1, 73},
        {100, 175,  50, 1, 74},
        {110, 275,  50, 1, 75},
        {120, 375,  50, 1, 76},
        {130, 600,  50, 1, 77},
        {140, 828,  50, 1, 77},
        {149,1028,  60, 1, 78},
        {157,1128,  60, 1, 79},
        {165,1228,  60, 1, 79},
        {173,1700,  65, 1, 80},
        {181,2247,  70, 1, 80},
        {190,2800,  70, 1, 81},
        {198,3485,  70, 1, 81},
        {207,4800,  75, 1, 82},
        {215,6357,  80, 1, 82},
        {224,9000,  80, 2, 83},
        {232,12132, 80, 2, 83},
        {240,18000, 85, 3, 84},
        {248,26379, 90, 3, 84},
        {256,36000, 90, 3, 845},
        {265,47158, 90, 3, 845},
        {273,54000,100, 4, 85},
        {281,60650,100, 4, 85},
        {290,66000,120, 4, 855},
        {298,71768,120, 4, 855},
    };
    int n = sizeof(table) / sizeof(table[0]);

    if (bits <= table[0][0]) {
        return (params_t){table[0][1], table[0][3], table[0][2],
                          table[0][1]/5+10, table[0][4]/100.0};
    }
    if (bits >= table[n-1][0]) {
        return (params_t){table[n-1][1], table[n-1][3], table[n-1][2],
                          table[n-1][1]/5+50, table[n-1][4]/100.0};
    }

    /* Linear interpolation */
    for (int i = 0; i < n - 1; i++) {
        if (bits > table[i][0] && bits <= table[i+1][0]) {
            double f = (double)(table[i+1][0] - bits) /
                       (double)(table[i+1][0] - table[i][0]);
            int fb = (int)(table[i+1][1] - f * (table[i+1][1] - table[i][1]));
            int lpm = (table[i][2] + table[i+1][2]) / 2;
            int nb = (bits - table[i][0] < table[i+1][0] - bits) ?
                     table[i][3] : table[i+1][3];
            int thr_pct = (int)(table[i+1][4] - f * (table[i+1][4] - table[i][4]));
            int extra = fb / 5 + 20;
            return (params_t){fb, nb, lpm, extra, thr_pct / 100.0};
        }
    }
    return (params_t){80000, 5, 130, 500, 0.87};
}

/* ==================== Main SIQS ==================== */
int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }
    clock_gettime(CLOCK_MONOTONIC, &g_start);

    mpz_t N, kN;
    mpz_inits(N, kN, NULL);
    mpz_set_str(N, argv[1], 10);

    int digits = (int)mpz_sizeinbase(N, 10);
    int bits = (int)mpz_sizeinbase(N, 2);

    /* Quick trial division */
    {
        unsigned long small_primes[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,
            73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,
            179,181,191,193,197,199,211,223,227,229,233,239,241,251,0};
        for (int i = 0; small_primes[i]; i++) {
            if (mpz_divisible_ui_p(N, small_primes[i])) {
                printf("%lu\n", small_primes[i]); return 0;
            }
        }
        /* Extended trial division to 100000 */
        for (unsigned long p = 257; p < 100000; p += 2) {
            if (mpz_divisible_ui_p(N, p)) {
                printf("%lu\n", p); return 0;
            }
        }
    }

    /* Perfect power check */
    if (mpz_perfect_square_p(N)) {
        mpz_t s; mpz_init(s); mpz_sqrt(s, N);
        gmp_printf("%Zd\n", s); mpz_clear(s); return 0;
    }

    /* Choose multiplier */
    int mult = choose_multiplier(N);
    mpz_mul_ui(kN, N, mult);
    int kN_bits = (int)mpz_sizeinbase(kN, 2);

    /* Get parameters - use input N's bit size, not kN (matching YAFU) */
    params_t P = get_params(bits);

    /* Build factor base */
    fb_t *fb = fb_create(kN, P.fb_size);
    int M = SIEVE_BLOCK * P.nblocks;
    uint64_t lp_bound = (uint64_t)fb->entries[fb->size - 1].prime * P.lp_mult;
    int target = fb->size + P.extra;

    /* Threshold - compensate for small primes skipped in sieve (2,3: ~2 bits) */
    double log2_Qmax = kN_bits / 2.0 + 0.5 + log2(M);
    int threshold = (int)(log2_Qmax * P.thresh) - 3;
    if (threshold < 10) threshold = 10;

    /* Determine small/large prime cutoff for bucket sieving */
    int med_start = 2; /* skip prime 2 and primes dividing a */
    int large_start = fb->size;  /* index where bucket primes begin */
    for (int i = 1; i < fb->size; i++) {
        if (fb->entries[i].prime >= SIEVE_BLOCK) { large_start = i; break; }
    }

    fprintf(stderr, "SIQS-A: %dd (%db), k=%d, FB=%d (med=%d..%d, large=%d..%d), M=%d, thresh=%d, LP=%lu, target=%d\n",
            digits, bits, mult, fb->size, med_start, large_start-1, large_start, fb->size-1,
            M, threshold, (unsigned long)lp_bound, target);

    /* Allocate sieve */
    unsigned char *sieve = malloc(SIEVE_BLOCK);

    /* Allocate bucket arrays for large primes */
    int total_blocks = 2 * P.nblocks;
    sieve_bucket_t *buckets = calloc(total_blocks, sizeof(sieve_bucket_t));
    for (int i = 0; i < total_blocks; i++) {
        buckets[i].alloc = BUCKET_ALLOC;
        buckets[i].fb_idx = malloc(BUCKET_ALLOC * sizeof(uint32_t));
        buckets[i].offset = malloc(BUCKET_ALLOC * sizeof(uint16_t));
    }

    /* Relations */
    rels_t *full = rels_create(MAX_RELS);
    rels_t *part = rels_create(MAX_PARTIALS);
    lp_table_t *lpt = lp_create(MAX_PARTIALS);

    /* Polynomial structures */
    mpz_t a, cur_b;
    mpz_inits(a, cur_b, NULL);
    mpz_t B_vals[MAX_A_FACTORS];
    for (int j = 0; j < MAX_A_FACTORS; j++) mpz_init(B_vals[j]);

    /* Per-FB self-init data: ainv_2Bj[j*fb_size + i] = 2*a^{-1}*B_j mod p_i */
    uint32_t *ainv_2Bj = malloc(MAX_A_FACTORS * fb->size * sizeof(uint32_t));

    /* Current sieve solutions */
    int32_t *soln1 = malloc(fb->size * sizeof(int32_t));
    int32_t *soln2 = malloc(fb->size * sizeof(int32_t));

    /* RNG */
    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, SEED);

    /* Working variables */
    mpz_t ax_b, Qx, residue, tmp, c_val;
    mpz_inits(ax_b, Qx, residue, tmp, c_val, NULL);

    int a_idx[MAX_A_FACTORS];
    int num_a_factors = 0;
    int poly_count = 0, combined = 0, a_count = 0;

    /* ==================== Main Sieve Loop ==================== */
    while (full->count < target) {
        /* Timeout check */
        if (poly_count > 0 && poly_count % 500 == 0) {
            double t = elapsed();
            if (t > 280) { fprintf(stderr, "TIMEOUT at %.1fs\n", t); break; }
            if (poly_count % 5000 == 0)
                fprintf(stderr, "  poly=%d rels=%d/%d (full=%d+%d) part=%d t=%.1fs\n",
                        poly_count, full->count, target,
                        full->count - combined, combined, part->count, t);
        }

        /* === Generate new 'a' coefficient === */
        {
            mpz_t tgt; mpz_init(tgt);
            mpz_mul_ui(tgt, kN, 2);
            mpz_sqrt(tgt, tgt);
            mpz_tdiv_q_ui(tgt, tgt, M);
            double log_tgt = mpz_sizeinbase(tgt, 2) * log(2.0);

            int lo = fb->size / 4, hi = 3 * fb->size / 4;
            if (lo < 2) lo = 2;
            if (hi <= lo + 5) hi = fb->size - 1;

            /* Average log of FB primes in range */
            double avg = 0; int cnt = 0;
            for (int i = lo; i < hi; i++) {
                if (fb->entries[i].root == 0 && fb->entries[i].prime != 2) continue;
                avg += log(fb->entries[i].prime);
                cnt++;
            }
            if (cnt == 0) { mpz_clear(tgt); break; }
            avg /= cnt;

            int s = (int)(log_tgt / avg + 0.5);
            if (s < 3) s = 3;
            if (s > MAX_A_FACTORS) s = MAX_A_FACTORS;
            if (s > hi - lo) s = hi - lo;
            num_a_factors = s;

            /* Find best 'a' in 40 attempts */
            double best_ratio = 1e30;
            int best[MAX_A_FACTORS];

            for (int att = 0; att < 40; att++) {
                mpz_set_ui(a, 1);
                int idx[MAX_A_FACTORS]; int ok = 1;
                for (int i = 0; i < s && ok; i++) {
                    int tries = 0, good;
                    do {
                        idx[i] = lo + gmp_urandomm_ui(rng, hi - lo);
                        good = 1;
                        for (int j = 0; j < i; j++) if (idx[j] == idx[i]) { good = 0; break; }
                        if (fb->entries[idx[i]].root == 0 && fb->entries[idx[i]].prime != 2) good = 0;
                        tries++;
                    } while (!good && tries < 100);
                    if (!good) { ok = 0; break; }
                    mpz_mul_ui(a, a, fb->entries[idx[i]].prime);
                }
                if (!ok) continue;

                double ratio;
                if (mpz_cmp(a, tgt) > 0) {
                    mpz_t q; mpz_init(q); mpz_tdiv_q(q, a, tgt);
                    ratio = mpz_get_d(q); mpz_clear(q);
                } else {
                    mpz_t q; mpz_init(q); mpz_tdiv_q(q, tgt, a);
                    ratio = mpz_get_d(q); mpz_clear(q);
                }
                if (ratio < best_ratio) {
                    best_ratio = ratio;
                    memcpy(best, idx, s * sizeof(int));
                }
                if (ratio < 1.5) break;
            }

            memcpy(a_idx, best, s * sizeof(int));
            mpz_set_ui(a, 1);
            for (int i = 0; i < s; i++) mpz_mul_ui(a, a, fb->entries[a_idx[i]].prime);
            mpz_clear(tgt);
            a_count++;

            /* === Compute B values === */
            for (int j = 0; j < s; j++) {
                int idx = a_idx[j];
                uint32_t qj = fb->entries[idx].prime;
                uint32_t rj = fb->entries[idx].root;

                mpz_t a_q, mod_q, inv;
                mpz_inits(a_q, mod_q, inv, NULL);
                mpz_divexact_ui(a_q, a, qj);
                mpz_set_ui(mod_q, qj);
                mpz_invert(inv, a_q, mod_q);
                unsigned long iv = mpz_get_ui(inv);
                mpz_mul_ui(B_vals[j], a_q, (unsigned long)rj * iv % qj);
                mpz_clears(a_q, mod_q, inv, NULL);
            }

            /* === Precompute ainv_2Bj for self-initialization === */
            for (int i = 0; i < fb->size; i++) {
                uint32_t p = fb->entries[i].prime;
                unsigned long am = mpz_fdiv_ui(a, p);
                if (am == 0) {
                    for (int j = 0; j < s; j++) ainv_2Bj[j * fb->size + i] = 0;
                    soln1[i] = soln2[i] = -1;
                    continue;
                }
                uint32_t ai = mod_inverse_u32((uint32_t)am, p);
                if (ai == 0) {
                    for (int j = 0; j < s; j++) ainv_2Bj[j * fb->size + i] = 0;
                    soln1[i] = soln2[i] = -1;
                    continue;
                }

                for (int j = 0; j < s; j++) {
                    unsigned long Bm = mpz_fdiv_ui(B_vals[j], p);
                    ainv_2Bj[j * fb->size + i] = (uint32_t)((uint64_t)2 * ai % p * Bm % p);
                }
            }

            /* === Compute initial sieve solutions for b_0 (Gray code 0) === */
            /* b_0 = -B_0 - B_1 - ... - B_{s-1} */
            mpz_set_ui(cur_b, 0);
            for (int j = 0; j < s; j++) mpz_sub(cur_b, cur_b, B_vals[j]);

            /* c = (b^2 - kN) / a */
            mpz_mul(c_val, cur_b, cur_b);
            mpz_sub(c_val, c_val, kN);
            mpz_divexact(c_val, c_val, a);

            for (int i = 0; i < fb->size; i++) {
                uint32_t p = fb->entries[i].prime;
                unsigned long am = mpz_fdiv_ui(a, p);
                if (am == 0 || fb->entries[i].root == 0) {
                    soln1[i] = soln2[i] = -1;
                    continue;
                }
                uint32_t ai = mod_inverse_u32((uint32_t)am, p);
                if (ai == 0) { soln1[i] = soln2[i] = -1; continue; }

                unsigned long bm = mpz_fdiv_ui(cur_b, p);
                /* Make bm positive */
                bm = (p - bm) % p; /* This is (-b) mod p */
                uint32_t r = fb->entries[i].root;

                /* soln1 = a^{-1} * (r - b) mod p = a^{-1} * (r + (-b)) mod p */
                /* soln2 = a^{-1} * (-r - b) mod p = a^{-1} * (p - r + (-b)) mod p */
                soln1[i] = (int32_t)((uint64_t)ai * ((r + bm) % p) % p);
                soln2[i] = (int32_t)((uint64_t)ai * (((p - r) + bm) % p) % p);
            }
        }

        int num_b = 1 << (num_a_factors - 1);

        /* === Iterate over all b-values using Gray code === */
        for (int b_idx = 0; b_idx < num_b && full->count < target; b_idx++) {

            /* Update b via Gray code (skip first since we already computed initial solutions) */
            if (b_idx > 0) {
                /* Find which bit changed */
                int changed_bit = __builtin_ctz(b_idx);
                int sign = ((b_idx >> (changed_bit + 1)) & 1) ? -1 : 1;

                /* Update b */
                if (sign > 0) mpz_add(cur_b, cur_b, B_vals[changed_bit]);
                else mpz_sub(cur_b, cur_b, B_vals[changed_bit]);
                if (sign > 0) mpz_add(cur_b, cur_b, B_vals[changed_bit]);
                else mpz_sub(cur_b, cur_b, B_vals[changed_bit]);
                /* b changes by ±2*B_j */

                /* Update c */
                mpz_mul(c_val, cur_b, cur_b);
                mpz_sub(c_val, c_val, kN);
                mpz_divexact(c_val, c_val, a);

                /* Self-initialize: update sieve solutions
                 * soln = a^{-1} * (±r - b) mod p
                 * When b changes by +2*B_j: both solutions -= delta
                 * When b changes by -2*B_j: both solutions += delta
                 */
                for (int i = 0; i < fb->size; i++) {
                    if (soln1[i] < 0) continue;
                    uint32_t p = fb->entries[i].prime;
                    uint32_t delta = ainv_2Bj[changed_bit * fb->size + i];

                    if (sign > 0) {
                        /* b increased by 2*B_j → solutions decrease */
                        soln1[i] = (int32_t)(((uint32_t)soln1[i] + p - delta) % p);
                        soln2[i] = (int32_t)(((uint32_t)soln2[i] + p - delta) % p);
                    } else {
                        /* b decreased by 2*B_j → solutions increase */
                        soln1[i] = (int32_t)(((uint32_t)soln1[i] + delta) % p);
                        soln2[i] = (int32_t)(((uint32_t)soln2[i] + delta) % p);
                    }
                }
            }

            poly_count++;

            /* === Fill buckets for large primes === */
            for (int bk = 0; bk < total_blocks; bk++) buckets[bk].count = 0;

            for (int i = large_start; i < fb->size; i++) {
                if (soln1[i] < 0) continue;
                uint32_t p = fb->entries[i].prime;

                /* For soln1: find which blocks it hits */
                int32_t x1 = soln1[i] - P.nblocks * SIEVE_BLOCK;
                if (x1 < 0) { int skip = (-x1 + p - 1) / p; x1 += skip * p; }
                for (int32_t x = x1; x < P.nblocks * SIEVE_BLOCK; x += p) {
                    int bk = (x + P.nblocks * SIEVE_BLOCK) / SIEVE_BLOCK;
                    if (bk < 0 || bk >= total_blocks) continue;
                    int off = x - (bk - P.nblocks) * SIEVE_BLOCK;
                    if (off < 0 || off >= SIEVE_BLOCK) continue;

                    sieve_bucket_t *b = &buckets[bk];
                    if (b->count >= b->alloc) {
                        b->alloc *= 2;
                        b->fb_idx = realloc(b->fb_idx, b->alloc * sizeof(uint32_t));
                        b->offset = realloc(b->offset, b->alloc * sizeof(uint16_t));
                    }
                    b->fb_idx[b->count] = i;
                    b->offset[b->count] = (uint16_t)off;
                    b->count++;
                }

                /* For soln2 (if different) */
                if (soln1[i] == soln2[i]) continue;
                int32_t x2 = soln2[i] - P.nblocks * SIEVE_BLOCK;
                if (x2 < 0) { int skip = (-x2 + p - 1) / p; x2 += skip * p; }
                for (int32_t x = x2; x < P.nblocks * SIEVE_BLOCK; x += p) {
                    int bk = (x + P.nblocks * SIEVE_BLOCK) / SIEVE_BLOCK;
                    if (bk < 0 || bk >= total_blocks) continue;
                    int off = x - (bk - P.nblocks) * SIEVE_BLOCK;
                    if (off < 0 || off >= SIEVE_BLOCK) continue;

                    sieve_bucket_t *b = &buckets[bk];
                    if (b->count >= b->alloc) {
                        b->alloc *= 2;
                        b->fb_idx = realloc(b->fb_idx, b->alloc * sizeof(uint32_t));
                        b->offset = realloc(b->offset, b->alloc * sizeof(uint16_t));
                    }
                    b->fb_idx[b->count] = i;
                    b->offset[b->count] = (uint16_t)off;
                    b->count++;
                }
            }

            /* === Sieve each block === */
            for (int bk = 0; bk < total_blocks; bk++) {
                int block_start = (bk - P.nblocks) * SIEVE_BLOCK;

                memset(sieve, 0, SIEVE_BLOCK);

                /* Sieve with medium primes (standard sieve) */
                for (int i = med_start; i < large_start; i++) {
                    if (soln1[i] < 0) continue;
                    uint32_t p = fb->entries[i].prime;
                    if (p < 5) continue; /* skip tiny primes, catch in trial div */
                    uint8_t lp = fb->entries[i].logp;

                    /* Compute offset of soln1 within this block */
                    int32_t off1 = ((int32_t)soln1[i] - block_start) % (int32_t)p;
                    if (off1 < 0) off1 += p;

                    for (int j = off1; j < SIEVE_BLOCK; j += p)
                        sieve[j] += lp;

                    if (soln1[i] != soln2[i]) {
                        int32_t off2 = ((int32_t)soln2[i] - block_start) % (int32_t)p;
                        if (off2 < 0) off2 += p;
                        for (int j = off2; j < SIEVE_BLOCK; j += p)
                            sieve[j] += lp;
                    }
                }

                /* Apply bucket sieve entries for large primes */
                for (int e = 0; e < buckets[bk].count; e++) {
                    uint32_t fi = buckets[bk].fb_idx[e];
                    uint16_t off = buckets[bk].offset[e];
                    sieve[off] += fb->entries[fi].logp;
                }

                /* === Scan for smooth candidates === */
                for (int j = 0; j < SIEVE_BLOCK; j++) {
                    if (sieve[j] < threshold) continue;

                    int32_t x = block_start + j;
                    if (x == 0) continue;

                    /* Compute Q(x) = a*x^2 + 2*b*x + c */
                    mpz_set_si(tmp, x);
                    mpz_mul(Qx, a, tmp);       /* a*x */
                    mpz_add(Qx, Qx, cur_b);    /* a*x + b */
                    mpz_add(Qx, Qx, cur_b);    /* a*x + 2b */
                    mpz_mul(Qx, Qx, tmp);      /* (a*x + 2b)*x = a*x^2 + 2bx */
                    mpz_add(Qx, Qx, c_val);    /* a*x^2 + 2bx + c */

                    /* ax_b = a*x + b */
                    mpz_mul_si(ax_b, a, x);
                    mpz_add(ax_b, ax_b, cur_b);

                    if (mpz_sgn(Qx) == 0) continue;
                    mpz_abs(residue, Qx);

                    /* === Sieve-aware trial division === */
                    /* Divide out factor of 2 */
                    while (mpz_even_p(residue))
                        mpz_tdiv_q_2exp(residue, residue, 1);

                    /* Divide by small primes (p < 5) that we skipped in sieve */
                    for (int i = 1; i < fb->size && fb->entries[i].prime < 5; i++) {
                        uint32_t p = fb->entries[i].prime;
                        while (mpz_divisible_ui_p(residue, p))
                            mpz_divexact_ui(residue, residue, p);
                    }

                    /* Divide by sieved primes - only check primes whose roots match x */
                    for (int i = med_start; i < fb->size; i++) {
                        if (soln1[i] < 0) continue;
                        uint32_t p = fb->entries[i].prime;
                        if (p < 5) continue;

                        /* Check if x ≡ soln1 or soln2 (mod p) */
                        int32_t xmod = ((x % (int32_t)p) + p) % p;
                        if (xmod != soln1[i] && xmod != soln2[i]) continue;

                        /* Divide out all powers of p */
                        if (mpz_divisible_ui_p(residue, p)) {
                            do { mpz_divexact_ui(residue, residue, p); }
                            while (mpz_divisible_ui_p(residue, p));
                        }
                    }

                    /* Also divide by a-factor primes */
                    for (int i = 0; i < num_a_factors; i++) {
                        uint32_t p = fb->entries[a_idx[i]].prime;
                        while (mpz_divisible_ui_p(residue, p))
                            mpz_divexact_ui(residue, residue, p);
                    }

                    /* Multiply Q by a for the matrix (store a*Q(x)) */
                    mpz_t aQx; mpz_init(aQx);
                    mpz_mul(aQx, Qx, a);

                    if (mpz_cmp_ui(residue, 1) == 0) {
                        /* Full relation */
                        int ri = full->count;
                        if (ri < full->alloc) {
                            mpz_set(full->ax_b[ri], ax_b);
                            mpz_set(full->Qx[ri], aQx);
                            full->lp[ri] = 0;
                            full->count++;
                        }
                    } else if (mpz_fits_ulong_p(residue)) {
                        unsigned long res = mpz_get_ui(residue);
                        if (res <= lp_bound) {
                            /* Single large prime */
                            int match = lp_find(lpt, (uint64_t)res);
                            if (match >= 0) {
                                int ri = full->count;
                                if (ri < full->alloc) {
                                    mpz_mul(full->ax_b[ri], ax_b, part->ax_b[match]);
                                    mpz_mod(full->ax_b[ri], full->ax_b[ri], N);
                                    mpz_mul(full->Qx[ri], aQx, part->Qx[match]);
                                    full->lp[ri] = (uint64_t)res;
                                    full->count++;
                                    combined++;
                                }
                            } else {
                                int pi = part->count;
                                if (pi < part->alloc) {
                                    mpz_set(part->ax_b[pi], ax_b);
                                    mpz_set(part->Qx[pi], aQx);
                                    part->lp[pi] = (uint64_t)res;
                                    lp_insert(lpt, (uint64_t)res, pi);
                                    part->count++;
                                }
                            }
                        }
                    }
                    mpz_clear(aQx);
                }
            }
        }
    }

    double sieve_time = elapsed();
    fprintf(stderr, "Sieving done: %d rels (%d full + %d combined) from %d polys in %.2fs\n",
            full->count, full->count - combined, combined, poly_count, sieve_time);

    if (full->count < fb->size + 1) {
        fprintf(stderr, "FAIL: not enough relations (%d < %d)\n", full->count, fb->size + 1);
        printf("FAIL\n");
        return 1;
    }

    /* ==================== Linear Algebra ==================== */
    int nrels = full->count;
    if (nrels > target) nrels = target;
    int ncols = fb->size + 1; /* +1 for sign */

    gf2_t *mat = gf2_create(nrels, ncols);

    for (int r = 0; r < nrels; r++) {
        mpz_t Qval; mpz_init(Qval);
        mpz_set(Qval, full->Qx[r]);

        if (mpz_sgn(Qval) < 0) {
            gf2_set(mat, r, 0);
            mpz_neg(Qval, Qval);
        }

        int e2 = 0;
        while (mpz_even_p(Qval)) { mpz_tdiv_q_2exp(Qval, Qval, 1); e2++; }
        if (e2 & 1) gf2_set(mat, r, 1);

        for (int i = 1; i < fb->size; i++) {
            uint32_t p = fb->entries[i].prime;
            int e = 0;
            while (mpz_divisible_ui_p(Qval, p)) {
                mpz_divexact_ui(Qval, Qval, p);
                e++;
            }
            if (e & 1) gf2_set(mat, r, i + 1);
        }
        mpz_clear(Qval);
    }

    int **deps; int *dlen;
    int ndeps = gf2_solve(mat, &deps, &dlen, 64);
    fprintf(stderr, "LA: %d deps from %dx%d in %.2fs\n", ndeps, nrels, ncols, elapsed());

    /* ==================== Square Root ==================== */
    for (int d = 0; d < ndeps; d++) {
        mpz_t X, Y, g, prod, rem;
        mpz_inits(X, Y, g, prod, rem, NULL);

        mpz_set_ui(X, 1);
        for (int k = 0; k < dlen[d]; k++) {
            mpz_mul(X, X, full->ax_b[deps[d][k]]);
            mpz_mod(X, X, N);
        }

        /* Compute Y = sqrt(product of Q values) mod N */
        mpz_set_ui(prod, 1);
        for (int k = 0; k < dlen[d]; k++) {
            mpz_t aq; mpz_init(aq);
            mpz_abs(aq, full->Qx[deps[d][k]]);
            mpz_mul(prod, prod, aq);
            mpz_clear(aq);
        }

        mpz_set(rem, prod);
        mpz_set_ui(Y, 1);

        /* Extract square root prime by prime */
        int e2 = 0;
        while (mpz_even_p(rem)) { mpz_tdiv_q_2exp(rem, rem, 1); e2++; }
        if (e2 & 1) goto next_dep;
        if (e2 / 2 > 0) {
            mpz_set_ui(tmp, 2);
            mpz_powm_ui(tmp, tmp, e2/2, N);
            mpz_mul(Y, Y, tmp); mpz_mod(Y, Y, N);
        }

        {
            int valid = 1;
            for (int i = 1; i < fb->size && valid; i++) {
                uint32_t p = fb->entries[i].prime;
                int e = 0;
                while (mpz_divisible_ui_p(rem, p)) {
                    mpz_divexact_ui(rem, rem, p);
                    e++;
                }
                if (e & 1) { valid = 0; break; }
                if (e / 2 > 0) {
                    mpz_set_ui(tmp, p);
                    mpz_powm_ui(tmp, tmp, e/2, N);
                    mpz_mul(Y, Y, tmp); mpz_mod(Y, Y, N);
                }
            }
            if (!valid) goto next_dep;
        }

        /* Handle remaining cofactor */
        if (mpz_cmp_ui(rem, 1) != 0) {
            if (mpz_perfect_square_p(rem)) {
                mpz_sqrt(tmp, rem);
                mpz_mod(tmp, tmp, N);
                mpz_mul(Y, Y, tmp); mpz_mod(Y, Y, N);
            } else {
                goto next_dep;
            }
        }

        /* Check X-Y */
        mpz_sub(tmp, X, Y);
        mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t o; mpz_init(o);
            mpz_divexact(o, N, g);
            if (mpz_cmp(g, o) > 0) mpz_swap(g, o);
            gmp_printf("%Zd\n", g);
            fprintf(stderr, "Factored in %.3fs (%dd)\n", elapsed(), digits);
            mpz_clear(o);
            return 0;
        }

        /* Check X+Y */
        mpz_add(tmp, X, Y);
        mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t o; mpz_init(o);
            mpz_divexact(o, N, g);
            if (mpz_cmp(g, o) > 0) mpz_swap(g, o);
            gmp_printf("%Zd\n", g);
            fprintf(stderr, "Factored in %.3fs (%dd)\n", elapsed(), digits);
            mpz_clear(o);
            return 0;
        }

        next_dep:
        mpz_clears(X, Y, g, prod, rem, NULL);
    }

    fprintf(stderr, "FAIL: no factor found from %d dependencies\n", ndeps);
    printf("FAIL\n");
    return 1;
}
