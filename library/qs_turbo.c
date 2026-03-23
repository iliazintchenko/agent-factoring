/*
 * QS-Turbo: High-performance SIQS with bucket sieving and DLP
 *
 * Key optimizations over SPQS:
 * 1. Bucket sieving: large FB primes stored in per-block buckets
 * 2. YAFU-calibrated parameters (larger FB, smaller sieve range)
 * 3. Double Large Prime with SQUFOF cofactorization
 * 4. Sieve-informed trial division (check only matching primes)
 * 5. Gray code SIQS self-initialization (O(1) poly switch)
 *
 * Compile: gcc -O3 -march=native -o qs_turbo library/qs_turbo.c -lgmp -lm
 * Usage: ./qs_turbo <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <stdint.h>

#define SEED 42
#define BLOCK_SIZE 32768          /* 32KB = L1 cache */
#define MAX_FB 120000
#define MAX_A_FACTORS 20
#define MAX_RELS 500000
#define MAX_PARTIALS 2000000
#define BUCKET_ALLOC 4096         /* entries per bucket slice */

/* ==================== Timing ==================== */
static struct timespec g_start;
static double elapsed(void) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec - g_start.tv_sec) + (now.tv_nsec - g_start.tv_nsec) / 1e9;
}

/* ==================== Modular Arithmetic ==================== */
static uint32_t mod_inverse_u32(uint32_t a, uint32_t m) {
    int64_t old_r = a, r = m, old_s = 1, s = 0;
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
    uint64_t b, e, m = p, r;
    /* Check QR */
    b = n % p; e = (p - 1) / 2; r = 1;
    { uint64_t bb = b, ee = e; while (ee) { if (ee & 1) r = (r * bb) % m; bb = (bb * bb) % m; ee >>= 1; } }
    if (r != 1) return 0;
    /* p ≡ 3 (mod 4): direct formula */
    if (p % 4 == 3) {
        b = n % p; e = (p + 1) / 4; r = 1;
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
    uint64_t M_val = S;
    b = z; e = Q; uint64_t c = 1;
    while (e) { if (e & 1) c = (c * b) % m; b = (b * b) % m; e >>= 1; }
    b = n % p; e = Q; uint64_t t = 1;
    while (e) { if (e & 1) t = (t * b) % m; b = (b * b) % m; e >>= 1; }
    b = n % p; e = (Q + 1) / 2; uint64_t R = 1;
    while (e) { if (e & 1) R = (R * b) % m; b = (b * b) % m; e >>= 1; }
    while (1) {
        if (t == 1) return (uint32_t)R;
        int i = 0; uint64_t tt = t;
        while (tt != 1) { tt = (tt * tt) % p; i++; }
        uint64_t bb = c;
        for (int j = 0; j < (int)M_val - i - 1; j++) bb = (bb * bb) % p;
        M_val = i; c = (bb * bb) % p; t = (t * c) % p; R = (R * bb) % p;
    }
}

/* ==================== Multiplier Selection ==================== */
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
                s += 2.0 * log(ps[i]) / (ps[i] - 1);
        }
        if (s > best) { best = s; best_k = k; }
        mpz_clear(kN);
    }
    return best_k;
}

/* ==================== Factor Base ==================== */
typedef struct {
    uint32_t *prime;
    uint32_t *root;       /* sqrt(kN) mod p */
    uint8_t  *logp;
    int size;
    int med_B;            /* index where primes exceed BLOCK_SIZE */
} fb_t;

static fb_t *fb_create(mpz_t kN, int target) {
    fb_t *fb = malloc(sizeof(fb_t));
    int alloc = target + 100;
    fb->prime = malloc(alloc * sizeof(uint32_t));
    fb->root  = malloc(alloc * sizeof(uint32_t));
    fb->logp  = malloc(alloc * sizeof(uint8_t));
    fb->prime[0] = 2; fb->root[0] = 1; fb->logp[0] = 1; fb->size = 1;
    fb->med_B = -1;

    /* Sieve of Eratosthenes for prime generation */
    int bound = (int)((double)target * 15.0) + 200000;
    char *sv = calloc(bound + 1, 1);
    for (int i = 2; (long)i*i <= bound; i++)
        if (!sv[i]) for (int j = i*i; j <= bound; j += i) sv[j] = 1;

    for (int i = 3; i <= bound && fb->size < target; i += 2) {
        if (sv[i]) continue;
        unsigned long nm = mpz_fdiv_ui(kN, i);
        if (nm == 0) {
            fb->prime[fb->size] = i; fb->root[fb->size] = 0;
            fb->logp[fb->size] = (uint8_t)(log2(i) + 0.5);
            fb->size++; continue;
        }
        uint32_t r = sqrt_mod_p((uint32_t)nm, i);
        if (!r) continue;
        fb->prime[fb->size] = i;
        fb->root[fb->size] = r;
        fb->logp[fb->size] = (uint8_t)(log2(i) + 0.5);
        if (fb->med_B < 0 && i > BLOCK_SIZE) fb->med_B = fb->size;
        fb->size++;
    }
    free(sv);
    if (fb->med_B < 0) fb->med_B = fb->size;
    return fb;
}

/* ==================== Large Prime Hash Table ==================== */
#define LP_HASH_BITS 22
#define LP_HASH_SIZE (1 << LP_HASH_BITS)
#define LP_HASH_MASK (LP_HASH_SIZE - 1)

typedef struct lp_entry {
    uint64_t lp;
    int rel_idx;
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
        if (e->lp == lp) return e->rel_idx;
    return -1;
}

static void lp_insert(lp_table_t *t, uint64_t lp, int idx) {
    if (t->used >= t->max) return;
    uint32_t h = lp_hash(lp);
    lp_entry_t *e = &t->pool[t->used++];
    e->lp = lp; e->rel_idx = idx; e->next = t->buckets[h]; t->buckets[h] = e;
}

/* ==================== DLP Hash Table ==================== */
/* For double large primes: key = min(lp1,lp2) * LARGE + max(lp1,lp2) */
/* Actually, we use a union-find on LP values. Each partial with LP=p is
 * an edge in the LP graph. When two partials share an LP, they can be
 * combined. Cycles give full relations. */

/* Simpler approach: for DLP, store (lp1, lp2, rel_idx). When we find
 * another relation with the same lp1 or lp2, we can combine. */

/* For simplicity, use SLP hash for each large prime in a DLP relation.
 * When both LPs of a DLP match different partials, combine all three. */

/* ==================== Relation Storage ==================== */
typedef struct {
    mpz_t *ax_b;      /* (ax+b) values */
    mpz_t *Qx;        /* Q(x) = (ax+b)^2 - kN (or a*Q) values */
    uint64_t *lp1;     /* first large prime (0 if full) */
    uint64_t *lp2;     /* second large prime (0 if SLP) */
    int count, alloc;
} rels_t;

static rels_t *rels_create(int n) {
    rels_t *r = malloc(sizeof(rels_t));
    r->ax_b = malloc(n * sizeof(mpz_t));
    r->Qx = malloc(n * sizeof(mpz_t));
    r->lp1 = calloc(n, sizeof(uint64_t));
    r->lp2 = calloc(n, sizeof(uint64_t));
    for (int i = 0; i < n; i++) { mpz_init(r->ax_b[i]); mpz_init(r->Qx[i]); }
    r->count = 0; r->alloc = n;
    return r;
}

/* ==================== GF(2) Matrix and Gaussian Elimination ==================== */
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

static void gf2_set(gf2_t *m, int r, int c) {
    m->rows[r][c/64] |= (1ULL << (c % 64));
}

static int gf2_solve(gf2_t *m, int ***deps, int **dlen, int max) {
    int piv = 0;
    for (int c = 0; c < m->nc && piv < m->nr; c++) {
        int pr = -1;
        for (int r = piv; r < m->nr; r++)
            if ((m->rows[r][c/64] >> (c % 64)) & 1) { pr = r; break; }
        if (pr < 0) continue;
        if (pr != piv) { u64 *t = m->rows[pr]; m->rows[pr] = m->rows[piv]; m->rows[piv] = t; }
        for (int r = 0; r < m->nr; r++) {
            if (r == piv) continue;
            if ((m->rows[r][c/64] >> (c % 64)) & 1)
                for (int w = 0; w < m->wprow; w++) m->rows[r][w] ^= m->rows[piv][w];
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

/* ==================== SQUFOF for DLP cofactorization ==================== */
/* Based on Jason Papadopoulos's implementation (public domain) */
/* Races multiple multipliers, uses saved-Q list to avoid trivial factors */

#define SQUFOF_MAX_MULT 16
#define SQUFOF_QSIZE 50
#define SQUFOF_MAX_CYCLES 100000
#define SQUFOF_ONE_CYCLE 300

static uint64_t gcd64(uint64_t a, uint64_t b) {
    while (b) { uint64_t t = b; b = a % b; a = t; }
    return a;
}

static uint64_t isqrt64(uint64_t n) {
    if (n == 0) return 0;
    uint64_t x = (uint64_t)sqrt((double)n);
    while (x > 0 && x * x > n) x--;
    while ((x + 1) * (x + 1) <= n) x++;
    return x;
}

static uint64_t squfof_factor(uint64_t n) {
    if (n <= 1) return n;
    if (n % 2 == 0) return 2;
    { uint64_t s = isqrt64(n); if (s * s == n) return s; }

    static const uint32_t mults[] = {1,3,5,7,11,15,21,33,35,55,77,105,165,231,385,1155};

    uint32_t sqrtn[SQUFOF_MAX_MULT], cutoff[SQUFOF_MAX_MULT];
    uint32_t q0[SQUFOF_MAX_MULT], p1[SQUFOF_MAX_MULT], q1[SQUFOF_MAX_MULT];
    uint32_t num_saved[SQUFOF_MAX_MULT];
    uint16_t saved[SQUFOF_MAX_MULT][SQUFOF_QSIZE];
    uint8_t failed[SQUFOF_MAX_MULT];
    int num_mult = 0;

    for (int i = 0; i < SQUFOF_MAX_MULT; i++) {
        __uint128_t sn = (__uint128_t)n * mults[i];
        if (sn >> 62) break; /* must fit in 62 bits */
        uint64_t scaled = (uint64_t)sn;
        uint64_t sq = isqrt64(scaled);

        sqrtn[i] = (uint32_t)sq;
        cutoff[i] = (uint32_t)sqrt(2.0 * (double)sq);
        q0[i] = 1;
        p1[i] = (uint32_t)sq;
        q1[i] = (uint32_t)(scaled - sq * sq);
        num_saved[i] = 0;
        failed[i] = 0;
        if (q1[i] == 0) { /* perfect square */
            uint64_t g = gcd64(n, sq);
            if (g > 1 && g < n) return g;
            failed[i] = 1;
        }
        num_mult = i + 1;
    }
    if (num_mult == 0) return 0;

    uint32_t total_iter = 0;
    while (total_iter < SQUFOF_MAX_CYCLES) {
        for (int mi = num_mult - 1; mi >= 0; mi--) {
            if (failed[mi]) continue;

            uint32_t sq_n = sqrtn[mi], cut = cutoff[mi];
            uint32_t multiplier = 2 * mults[mi];
            uint32_t coarse_cut = cut * multiplier;
            uint32_t Q0 = q0[mi], P1 = p1[mi], Q1 = q1[mi];
            uint32_t ns = num_saved[mi];
            uint32_t P0 = 0, sqrtq = 0;
            int found_sq = 0;

            for (uint32_t it = 0; it < SQUFOF_ONE_CYCLE; it++) {
                /* Even step */
                uint32_t tmp = sq_n + P1 - Q1;
                uint32_t q = 1;
                if (tmp >= Q1) q += tmp / Q1;
                P0 = q * Q1 - P1;
                Q0 = Q0 + (P1 - P0) * q;

                /* Save small Q1 values */
                if (Q1 < coarse_cut) {
                    uint32_t t = Q1 / gcd64(Q1, multiplier);
                    if (t < cut) {
                        if (ns >= SQUFOF_QSIZE) { failed[mi] = 1; break; }
                        saved[mi][ns++] = (uint16_t)t;
                    }
                }

                /* Check Q0 for perfect square */
                uint32_t bits = 0; tmp = Q0;
                while (tmp && !(tmp & 1)) { bits++; tmp >>= 1; }
                if (!(bits & 1) && ((tmp & 7) == 1)) {
                    sqrtq = (uint32_t)sqrt((double)Q0);
                    while (sqrtq > 0 && sqrtq * sqrtq > Q0) sqrtq--;
                    while ((sqrtq + 1) * (sqrtq + 1) <= Q0) sqrtq++;
                    if (sqrtq * sqrtq == Q0) {
                        uint32_t j;
                        for (j = 0; j < ns; j++)
                            if (saved[mi][j] == sqrtq) break;
                        if (j == ns) { found_sq = 1; break; }
                    }
                }

                /* Odd step */
                tmp = sq_n + P0 - Q0;
                q = 1;
                if (tmp >= Q0) q += tmp / Q0;
                P1 = q * Q0 - P0;
                Q1 = Q1 + (P0 - P1) * q;

                if (Q0 < coarse_cut) {
                    uint32_t t = Q0 / gcd64(Q0, multiplier);
                    if (t < cut) {
                        if (ns >= SQUFOF_QSIZE) { failed[mi] = 1; break; }
                        saved[mi][ns++] = (uint16_t)t;
                    }
                }
                total_iter++;
            }

            if (failed[mi]) continue;

            if (!found_sq) {
                q0[mi] = Q0; p1[mi] = P1; q1[mi] = Q1;
                num_saved[mi] = ns;
                continue;
            }

            if (sqrtq <= 1) { failed[mi] = 1; continue; }

            /* Inverse cycle */
            uint64_t scaledn = (uint64_t)n * mults[mi];
            Q0 = sqrtq;
            P1 = P0 + sqrtq * ((sq_n - P0) / sqrtq);
            Q1 = (uint32_t)((scaledn - (uint64_t)P1 * P1) / Q0);

            while (1) {
                uint32_t tmp = sq_n + P1 - Q1;
                uint32_t q = 1;
                if (tmp >= Q1) q += tmp / Q1;
                P0 = q * Q1 - P1;
                Q0 = Q0 + (P1 - P0) * q;
                if (P0 == P1) { Q0 = Q1; break; }

                tmp = sq_n + P0 - Q0;
                q = 1;
                if (tmp >= Q0) q += tmp / Q0;
                P1 = q * Q0 - P0;
                Q1 = Q1 + (P0 - P1) * q;
                if (P0 == P1) break;
            }

            Q0 = Q0 / (uint32_t)gcd64(Q0, multiplier);
            if (Q0 > 1 && Q0 < n && n % Q0 == 0) return Q0;
        }
    }
    return 0;
}

/* ==================== Parameters (YAFU-calibrated) ==================== */
typedef struct {
    int fb_size;
    int num_blocks;
    int lp_mult;       /* LP bound = FB_max * lp_mult */
    double dlp_exp;    /* DLP: cofactor < FB_max^dlp_exp */
    double thresh;     /* sieve threshold fraction */
    int extra_rels;
} params_t;

static params_t get_params(int bits) {
    /* Tuned for scalar sieve: smaller FB, more blocks, aggressive DLP */
    if (bits <= 60)  return (params_t){30,    1, 30,  0,    0.72, 20};
    if (bits <= 70)  return (params_t){40,    1, 40,  0,    0.72, 25};
    if (bits <= 80)  return (params_t){60,    1, 40,  0,    0.73, 30};
    if (bits <= 90)  return (params_t){90,    1, 40,  0,    0.74, 30};
    if (bits <= 100) return (params_t){130,   1, 50,  0,    0.75, 40};
    if (bits <= 110) return (params_t){200,   1, 50,  1.75, 0.76, 40};
    if (bits <= 120) return (params_t){280,   2, 50,  1.75, 0.77, 50};
    if (bits <= 130) return (params_t){400,   2, 50,  1.75, 0.78, 50};
    if (bits <= 140) return (params_t){550,   3, 60,  1.75, 0.79, 60};
    if (bits <= 150) return (params_t){750,   4, 60,  1.75, 0.80, 60};
    if (bits <= 160) return (params_t){1000,  5, 60,  1.80, 0.81, 80};
    if (bits <= 170) return (params_t){1400,  6, 70,  1.80, 0.82, 80};
    if (bits <= 180) return (params_t){2000,  8, 70,  1.80, 0.83, 100};
    if (bits <= 190) return (params_t){2800, 10, 80,  1.80, 0.84, 100};
    if (bits <= 200) return (params_t){4000, 12, 80,  1.80, 0.85, 120};
    if (bits <= 215) return (params_t){6000, 16, 90,  1.80, 0.86, 120};
    if (bits <= 232) return (params_t){10000, 20, 90, 1.80, 0.87, 150};
    if (bits <= 248) return (params_t){18000, 28, 100, 1.80, 0.88, 150};
    if (bits <= 265) return (params_t){30000, 36, 100, 1.85, 0.89, 200};
    if (bits <= 281) return (params_t){45000, 44, 110, 1.85, 0.90, 200};
    if (bits <= 298) return (params_t){65000, 52, 120, 1.85, 0.91, 250};
    return (params_t){90000, 64, 140, 1.85, 0.91, 300};
}

/* ==================== Bucket Sieve Data Structure ==================== */
typedef struct {
    uint32_t *hits;    /* packed: (position << 8) | logp */
    int *count;        /* hits per block */
    int num_blocks;
    int alloc_per_block;
} bucket_t;

static bucket_t *bucket_create(int nblocks, int alloc) {
    bucket_t *b = malloc(sizeof(bucket_t));
    b->num_blocks = nblocks * 2; /* positive + negative */
    b->alloc_per_block = alloc;
    b->hits = malloc((size_t)b->num_blocks * alloc * sizeof(uint32_t));
    b->count = calloc(b->num_blocks, sizeof(int));
    return b;
}

static void bucket_reset(bucket_t *b) {
    memset(b->count, 0, b->num_blocks * sizeof(int));
}

static inline void bucket_add(bucket_t *b, int block_idx, uint32_t pos, uint8_t logp) {
    if (block_idx < 0 || block_idx >= b->num_blocks) return;
    int c = b->count[block_idx];
    if (c >= b->alloc_per_block) return;
    b->hits[(size_t)block_idx * b->alloc_per_block + c] = (pos << 8) | logp;
    b->count[block_idx] = c + 1;
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
    for (unsigned long p = 2; p < 100000; p++) {
        if (p > 2 && p % 2 == 0) continue;
        if (p > 3 && p % 3 == 0) continue;
        if (mpz_divisible_ui_p(N, p)) {
            printf("%lu\n", p);
            fprintf(stderr, "Trial division: %lu in %.3fs\n", p, elapsed());
            return 0;
        }
    }

    int mult = choose_multiplier(N);
    mpz_mul_ui(kN, N, mult);
    int kN_bits = (int)mpz_sizeinbase(kN, 2);
    params_t P = get_params(kN_bits);

    fb_t *fb = fb_create(kN, P.fb_size);
    int M = BLOCK_SIZE * P.num_blocks;
    uint64_t lp_bound = (uint64_t)fb->prime[fb->size - 1] * P.lp_mult;
    uint64_t dlp_bound = 0;
    if (P.dlp_exp > 0) {
        double fb_max_d = (double)fb->prime[fb->size - 1];
        dlp_bound = (uint64_t)pow(fb_max_d, P.dlp_exp);
        if (dlp_bound < lp_bound) dlp_bound = lp_bound;
    }

    int target = fb->size + P.extra_rels;

    /* Compute sieve threshold */
    double log2_Qmax = kN_bits / 2.0 + 0.5 + log2((double)M);
    int threshold = (int)(log2_Qmax * P.thresh);
    if (threshold < 30) threshold = 30;

    fprintf(stderr, "QS-Turbo: %dd (%db), k=%d, FB=%d (med_B=%d), M=%d, "
            "thresh=%d, LP=%lu, DLP=%lu, target=%d\n",
            digits, bits, mult, fb->size, fb->med_B, M, threshold,
            (unsigned long)lp_bound, (unsigned long)dlp_bound, target);

    /* Allocate sieve and buckets */
    uint8_t *sieve = aligned_alloc(64, BLOCK_SIZE);
    int total_blocks = P.num_blocks * 2;
    bucket_t *buckets = bucket_create(P.num_blocks, BUCKET_ALLOC);

    /* Relation storage */
    rels_t *full_rels = rels_create(MAX_RELS);
    rels_t *part_rels = rels_create(MAX_PARTIALS);
    lp_table_t *slp_table = lp_create(MAX_PARTIALS);

    /* Polynomial state */
    mpz_t a, b_val, c_val;
    mpz_t B_vals[MAX_A_FACTORS];
    mpz_inits(a, b_val, c_val, NULL);
    for (int j = 0; j < MAX_A_FACTORS; j++) mpz_init(B_vals[j]);

    /* Sieve solution arrays */
    uint32_t *soln1 = malloc(fb->size * sizeof(uint32_t));
    uint32_t *soln2 = malloc(fb->size * sizeof(uint32_t));

    /* ainv data for Gray code update */
    uint32_t *ainv = malloc(MAX_A_FACTORS * fb->size * sizeof(uint32_t));

    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, SEED);

    mpz_t ax_b, Qx, residue, tmp, aQx;
    mpz_inits(ax_b, Qx, residue, tmp, aQx, NULL);

    int total_polys = 0, combined_slp = 0, combined_dlp = 0;
    int a_count = 0;
    int a_idx[MAX_A_FACTORS];
    int num_a_factors = 0;
    int total_sieve_hits = 0;

    /* Main sieving loop */
    while (full_rels->count < target) {
        double t = elapsed();
        if (t > 280) { fprintf(stderr, "TIMEOUT at %.1fs\n", t); break; }

        if (total_polys > 0 && total_polys % 500 == 0) {
            fprintf(stderr, "  poly=%d rels=%d/%d (full=%d slp=%d dlp=%d) part=%d hits=%d t=%.1fs\n",
                    total_polys, full_rels->count, target,
                    full_rels->count - combined_slp - combined_dlp,
                    combined_slp, combined_dlp, part_rels->count,
                    total_sieve_hits, t);
        }

        /* === Generate new 'a' value === */
        {
            mpz_t tgt; mpz_init(tgt);
            mpz_mul_ui(tgt, kN, 2);
            mpz_sqrt(tgt, tgt);
            mpz_tdiv_q_ui(tgt, tgt, M > 0 ? M : 1);
            double log_tgt = mpz_sizeinbase(tgt, 2) * log(2.0);

            /* Select FB primes for 'a' from middle third */
            int lo = fb->size / 4, hi = 3 * fb->size / 4;
            if (lo < 2) lo = 2;
            if (hi <= lo + 3) hi = fb->size - 1;

            double avg = 0; int cnt = 0;
            for (int i = lo; i < hi; i++) {
                if (fb->root[i] == 0) continue;
                avg += log(fb->prime[i]); cnt++;
            }
            if (cnt == 0) { mpz_clear(tgt); break; }
            avg /= cnt;

            int s = (int)(log_tgt / avg + 0.5);
            if (s < 3) s = 3;
            if (s > MAX_A_FACTORS) s = MAX_A_FACTORS;
            if (s > hi - lo) s = hi - lo;
            num_a_factors = s;

            /* Try multiple random selections, keep best */
            double best_ratio = 1e30;
            int best[MAX_A_FACTORS];

            for (int att = 0; att < 50; att++) {
                mpz_set_ui(a, 1);
                int idx[MAX_A_FACTORS]; int ok = 1;
                for (int i = 0; i < s && ok; i++) {
                    int tries = 0, good;
                    do {
                        idx[i] = lo + gmp_urandomm_ui(rng, hi - lo);
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
            for (int i = 0; i < s; i++)
                mpz_mul_ui(a, a, fb->prime[a_idx[i]]);
            mpz_clear(tgt);
            a_count++;

            /* Compute B values (Tonelli-Shanks for each a-factor) */
            for (int j = 0; j < s; j++) {
                int idx = a_idx[j];
                uint32_t qj = fb->prime[idx], rj = fb->root[idx];
                mpz_t a_q, mod_q, inv;
                mpz_inits(a_q, mod_q, inv, NULL);
                mpz_divexact_ui(a_q, a, qj);
                mpz_set_ui(mod_q, qj);
                mpz_invert(inv, a_q, mod_q);
                unsigned long iv = mpz_get_ui(inv);
                mpz_mul_ui(B_vals[j], a_q, (rj * iv) % qj);
                mpz_clears(a_q, mod_q, inv, NULL);
            }

            /* Precompute ainv for Gray code polynomial switching */
            for (int j = 0; j < s; j++) {
                for (int i = 0; i < fb->size; i++) {
                    uint32_t p = fb->prime[i];
                    unsigned long am = mpz_fdiv_ui(a, p);
                    if (am == 0 || fb->root[i] == 0) {
                        ainv[j * fb->size + i] = 0;
                        continue;
                    }
                    uint32_t ai = mod_inverse_u32((uint32_t)am, p);
                    unsigned long Bm = mpz_fdiv_ui(B_vals[j], p);
                    ainv[j * fb->size + i] = (uint32_t)((2ULL * ai % p * Bm) % p);
                }
            }
        }

        /* === Iterate over all b-values using Gray code === */
        int num_b = 1 << (num_a_factors - 1);

        /* Initialize first b-value */
        mpz_set_ui(b_val, 0);
        for (int j = 0; j < num_a_factors; j++)
            mpz_add(b_val, b_val, B_vals[j]);

        /* Compute initial sieve solutions */
        for (int i = 0; i < fb->size; i++) {
            uint32_t p = fb->prime[i];
            unsigned long am = mpz_fdiv_ui(a, p);
            if (am == 0 || fb->root[i] == 0) {
                soln1[i] = soln2[i] = 0xFFFFFFFF;
                continue;
            }
            uint32_t ai = mod_inverse_u32((uint32_t)am, p);
            if (ai == 0) { soln1[i] = soln2[i] = 0xFFFFFFFF; continue; }
            unsigned long bm = mpz_fdiv_ui(b_val, p);
            uint32_t r = fb->root[i];
            soln1[i] = (uint32_t)((uint64_t)ai * ((r + p - bm) % p) % p);
            soln2[i] = (uint32_t)((uint64_t)ai * ((p - r + p - bm) % p) % p);
        }

        /* c = (b^2 - kN) / a */
        mpz_mul(c_val, b_val, b_val);
        mpz_sub(c_val, c_val, kN);

        /* Verify b^2 ≡ kN (mod a) */
        mpz_mod(tmp, c_val, a);
        if (mpz_sgn(tmp) != 0) {
            /* Try negation */
            mpz_neg(b_val, b_val);
            mpz_mul(c_val, b_val, b_val);
            mpz_sub(c_val, c_val, kN);
            mpz_mod(tmp, c_val, a);
            if (mpz_sgn(tmp) != 0) continue;
            /* Recompute solutions */
            for (int i = 0; i < fb->size; i++) {
                uint32_t p = fb->prime[i];
                unsigned long am = mpz_fdiv_ui(a, p);
                if (am == 0 || fb->root[i] == 0) { soln1[i] = soln2[i] = 0xFFFFFFFF; continue; }
                uint32_t ai = mod_inverse_u32((uint32_t)am, p);
                if (ai == 0) { soln1[i] = soln2[i] = 0xFFFFFFFF; continue; }
                unsigned long bm = mpz_fdiv_ui(b_val, p);
                uint32_t r = fb->root[i];
                soln1[i] = (uint32_t)((uint64_t)ai * ((r + p - bm) % p) % p);
                soln2[i] = (uint32_t)((uint64_t)ai * ((p - r + p - bm) % p) % p);
            }
        }
        mpz_divexact(c_val, c_val, a);

        for (int b_idx = 0; b_idx < num_b && full_rels->count < target; b_idx++) {
            /* For b_idx > 0, update b via Gray code */
            if (b_idx > 0) {
                int gray_prev = (b_idx - 1) ^ ((b_idx - 1) >> 1);
                int gray_curr = b_idx ^ (b_idx >> 1);
                int changed_bit = gray_prev ^ gray_curr;
                int j = __builtin_ctz(changed_bit);
                int sign = (gray_curr >> j) & 1; /* 1 = add, 0 = subtract */

                if (sign) {
                    mpz_add(b_val, b_val, B_vals[j]);
                    mpz_add(b_val, b_val, B_vals[j]);
                } else {
                    mpz_sub(b_val, b_val, B_vals[j]);
                    mpz_sub(b_val, b_val, B_vals[j]);
                }

                /* Update sieve solutions: soln += ±ainv[j] */
                for (int i = 0; i < fb->size; i++) {
                    if (soln1[i] == 0xFFFFFFFF) continue;
                    uint32_t p = fb->prime[i];
                    uint32_t delta = ainv[j * fb->size + i];
                    if (delta == 0) continue;
                    if (sign) {
                        /* Adding: soln = (soln - delta) mod p */
                        soln1[i] = (soln1[i] + p - delta) % p;
                        soln2[i] = (soln2[i] + p - delta) % p;
                    } else {
                        /* Subtracting: soln = (soln + delta) mod p */
                        soln1[i] = (soln1[i] + delta) % p;
                        soln2[i] = (soln2[i] + delta) % p;
                    }
                }

                /* Update c */
                mpz_mul(c_val, b_val, b_val);
                mpz_sub(c_val, c_val, kN);
                mpz_divexact(c_val, c_val, a);
            }

            total_polys++;

            /* === Fill buckets for large primes === */
            bucket_reset(buckets);

            for (int i = fb->med_B; i < fb->size; i++) {
                if (soln1[i] == 0xFFFFFFFF) continue;
                uint32_t p = fb->prime[i];
                uint8_t lp = fb->logp[i];

                /* Positive side */
                for (int r = 0; r < 2; r++) {
                    uint32_t root = (r == 0) ? soln1[i] : soln2[i];
                    if (r == 1 && soln1[i] == soln2[i]) break;

                    /* Hits in positive blocks [0, M) */
                    uint32_t pos = root;
                    while (pos < (uint32_t)M) {
                        int blk = pos / BLOCK_SIZE;
                        uint32_t off = pos % BLOCK_SIZE;
                        bucket_add(buckets, P.num_blocks + blk, off, lp);
                        pos += p;
                    }
                    /* Hits in negative blocks [-M, 0) mapped to [0, M) */
                    /* For negative x: need (soln - x) ≡ 0 mod p where x in [-M,0)
                     * i.e., positions in the negative sieve block */
                    long neg_start = -((long)p - (long)root);
                    if (neg_start >= 0) neg_start -= p;
                    for (long nx = neg_start; nx > -(long)M; nx -= p) {
                        long abs_nx = -nx;
                        if (abs_nx <= 0 || abs_nx > M) continue;
                        int blk = (int)((abs_nx - 1) / BLOCK_SIZE);
                        uint32_t off = BLOCK_SIZE - 1 - (int)((abs_nx - 1) % BLOCK_SIZE);
                        if (blk < P.num_blocks)
                            bucket_add(buckets, blk, off, lp);
                    }
                }
            }

            /* === Sieve each block === */
            for (int side = 0; side < 2; side++) {
                for (int blk = 0; blk < P.num_blocks; blk++) {
                    int block_start;
                    if (side == 1) {
                        /* Positive side: x in [blk*BLOCK_SIZE, (blk+1)*BLOCK_SIZE) */
                        block_start = blk * BLOCK_SIZE;
                    } else {
                        /* Negative side: x in [-(blk+1)*BLOCK_SIZE, -blk*BLOCK_SIZE) */
                        block_start = -(blk + 1) * BLOCK_SIZE;
                    }

                    memset(sieve, 0, BLOCK_SIZE);

                    /* Small/medium primes: standard sieve */
                    for (int i = 1; i < fb->med_B; i++) {
                        uint32_t p = fb->prime[i];
                        if (p < 5) continue;
                        uint8_t lp_val = fb->logp[i];

                        for (int r = 0; r < 2; r++) {
                            uint32_t root = (r == 0) ? soln1[i] : soln2[i];
                            if (soln1[i] == 0xFFFFFFFF) break;
                            if (r == 1 && soln1[i] == soln2[i]) break;

                            long off = ((long)root - (long)block_start) % (long)p;
                            if (off < 0) off += p;

                            /* Unrolled inner sieve loop */
                            int j = (int)off;
                            int limit = BLOCK_SIZE - (int)p;
                            while (j <= limit) {
                                sieve[j] += lp_val;
                                j += p;
                            }
                            if (j < BLOCK_SIZE) sieve[j] += lp_val;
                        }
                    }

                    /* Large primes: apply bucket hits */
                    int bucket_idx = side ? (P.num_blocks + blk) : blk;
                    int nhits = buckets->count[bucket_idx];
                    uint32_t *hits = buckets->hits + (size_t)bucket_idx * buckets->alloc_per_block;
                    for (int h = 0; h < nhits; h++) {
                        uint32_t packed = hits[h];
                        uint32_t pos = packed >> 8;
                        uint8_t lp_val = packed & 0xFF;
                        if (pos < BLOCK_SIZE)
                            sieve[pos] += lp_val;
                    }

                    /* === Scan for smooth candidates === */
                    for (int j = 0; j < BLOCK_SIZE; j++) {
                        if (sieve[j] < threshold) continue;

                        long x = (long)block_start + j;
                        if (x == 0) continue;
                        total_sieve_hits++;

                        /* Compute Q(x) = a*x^2 + 2*b*x + c */
                        mpz_set_si(tmp, x);
                        mpz_mul(Qx, a, tmp);
                        mpz_add(Qx, Qx, b_val);
                        mpz_add(Qx, Qx, b_val);
                        mpz_mul(Qx, Qx, tmp);
                        mpz_add(Qx, Qx, c_val);

                        /* ax_b = a*x + b */
                        mpz_mul_si(ax_b, a, x);
                        mpz_add(ax_b, ax_b, b_val);

                        if (mpz_sgn(Qx) == 0) continue;
                        mpz_abs(residue, Qx);

                        /* === Sieve-informed trial division === */
                        /* Divide by 2 */
                        while (mpz_even_p(residue))
                            mpz_tdiv_q_2exp(residue, residue, 1);

                        /* Only trial-divide primes whose roots match x mod p */
                        for (int i = 1; i < fb->size; i++) {
                            if (soln1[i] == 0xFFFFFFFF) continue;
                            uint32_t p = fb->prime[i];

                            /* Check if this prime divides Q(x) via sieve root matching */
                            long xmod = ((x % (long)p) + p) % p;
                            if (xmod != (long)soln1[i] && xmod != (long)soln2[i]) continue;

                            /* Divide out all powers */
                            if (mpz_divisible_ui_p(residue, p))
                                do { mpz_divexact_ui(residue, residue, p); }
                                while (mpz_divisible_ui_p(residue, p));
                        }

                        /* Handle small primes not caught by roots */
                        for (int i = 1; i < fb->size && fb->prime[i] < 5; i++) {
                            uint32_t p = fb->prime[i];
                            while (mpz_divisible_ui_p(residue, p))
                                mpz_divexact_ui(residue, residue, p);
                        }

                        /* Compute a*Q(x) for the congruence (ax+b)^2 ≡ a*Q(x) (mod N) */
                        mpz_mul(aQx, Qx, a);

                        /* === Classify relation === */
                        if (mpz_cmp_ui(residue, 1) == 0) {
                            /* Full relation */
                            int ri = full_rels->count;
                            if (ri < full_rels->alloc) {
                                mpz_set(full_rels->ax_b[ri], ax_b);
                                mpz_set(full_rels->Qx[ri], aQx);
                                full_rels->lp1[ri] = 0;
                                full_rels->lp2[ri] = 0;
                                full_rels->count++;
                            }
                        } else if (mpz_fits_ulong_p(residue)) {
                            uint64_t cofactor = mpz_get_ui(residue);

                            if (cofactor <= lp_bound) {
                                /* Single Large Prime */
                                int match = lp_find(slp_table, cofactor);
                                if (match >= 0) {
                                    int ri = full_rels->count;
                                    if (ri < full_rels->alloc) {
                                        mpz_mul(full_rels->ax_b[ri], ax_b, part_rels->ax_b[match]);
                                        mpz_mod(full_rels->ax_b[ri], full_rels->ax_b[ri], N);
                                        mpz_mul(full_rels->Qx[ri], aQx, part_rels->Qx[match]);
                                        full_rels->lp1[ri] = cofactor;
                                        full_rels->count++;
                                        combined_slp++;
                                    }
                                } else {
                                    int pi = part_rels->count;
                                    if (pi < part_rels->alloc) {
                                        mpz_set(part_rels->ax_b[pi], ax_b);
                                        mpz_set(part_rels->Qx[pi], aQx);
                                        part_rels->lp1[pi] = cofactor;
                                        lp_insert(slp_table, cofactor, pi);
                                        part_rels->count++;
                                    }
                                }
                            } else if (dlp_bound > 0 && cofactor <= dlp_bound) {
                                /* Double Large Prime candidate - try to split */
                                uint64_t f1 = squfof_factor(cofactor);
                                if (f1 > 1 && f1 < cofactor) {
                                    uint64_t f2 = cofactor / f1;
                                    if (f1 > f2) { uint64_t t2 = f1; f1 = f2; f2 = t2; }
                                    if (f1 <= lp_bound && f2 <= lp_bound) {
                                        /* Both factors within LP bound - treat as SLP for each */
                                        /* Store with lp1 as the product key */
                                        int match1 = lp_find(slp_table, f1);
                                        int match2 = lp_find(slp_table, f2);

                                        if (match1 >= 0 && match2 >= 0) {
                                            /* Both LPs matched - combine three relations */
                                            int ri = full_rels->count;
                                            if (ri < full_rels->alloc) {
                                                mpz_mul(full_rels->ax_b[ri], ax_b, part_rels->ax_b[match1]);
                                                mpz_mul(full_rels->ax_b[ri], full_rels->ax_b[ri], part_rels->ax_b[match2]);
                                                mpz_mod(full_rels->ax_b[ri], full_rels->ax_b[ri], N);
                                                mpz_mul(full_rels->Qx[ri], aQx, part_rels->Qx[match1]);
                                                mpz_mul(full_rels->Qx[ri], full_rels->Qx[ri], part_rels->Qx[match2]);
                                                full_rels->lp1[ri] = f1;
                                                full_rels->lp2[ri] = f2;
                                                full_rels->count++;
                                                combined_dlp++;
                                            }
                                        } else if (match1 >= 0) {
                                            /* Only f1 matched - combine and store as partial with f2 */
                                            int pi = part_rels->count;
                                            if (pi < part_rels->alloc) {
                                                mpz_mul(part_rels->ax_b[pi], ax_b, part_rels->ax_b[match1]);
                                                mpz_mod(part_rels->ax_b[pi], part_rels->ax_b[pi], N);
                                                mpz_mul(part_rels->Qx[pi], aQx, part_rels->Qx[match1]);
                                                part_rels->lp1[pi] = f2;
                                                lp_insert(slp_table, f2, pi);
                                                part_rels->count++;
                                            }
                                        } else if (match2 >= 0) {
                                            int pi = part_rels->count;
                                            if (pi < part_rels->alloc) {
                                                mpz_mul(part_rels->ax_b[pi], ax_b, part_rels->ax_b[match2]);
                                                mpz_mod(part_rels->ax_b[pi], part_rels->ax_b[pi], N);
                                                mpz_mul(part_rels->Qx[pi], aQx, part_rels->Qx[match2]);
                                                part_rels->lp1[pi] = f1;
                                                lp_insert(slp_table, f1, pi);
                                                part_rels->count++;
                                            }
                                        } else {
                                            /* Neither matched - store both as partials */
                                            int pi = part_rels->count;
                                            if (pi < part_rels->alloc) {
                                                mpz_set(part_rels->ax_b[pi], ax_b);
                                                mpz_set(part_rels->Qx[pi], aQx);
                                                part_rels->lp1[pi] = f1;
                                                lp_insert(slp_table, f1, pi);
                                                part_rels->count++;
                                            }
                                            pi = part_rels->count;
                                            if (pi < part_rels->alloc) {
                                                mpz_set(part_rels->ax_b[pi], ax_b);
                                                mpz_set(part_rels->Qx[pi], aQx);
                                                part_rels->lp1[pi] = f2;
                                                lp_insert(slp_table, f2, pi);
                                                part_rels->count++;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    double sieve_time = elapsed();
    fprintf(stderr, "Sieving done: %d rels (%d full + %d SLP + %d DLP) in %.2fs, %d polys\n",
            full_rels->count,
            full_rels->count - combined_slp - combined_dlp,
            combined_slp, combined_dlp, sieve_time, total_polys);

    if (full_rels->count < fb->size + 1) {
        fprintf(stderr, "FAIL: not enough relations (%d < %d)\n",
                full_rels->count, fb->size + 1);
        printf("FAIL\n");
        return 1;
    }

    /* === Linear Algebra === */
    int nrels = full_rels->count;
    if (nrels > target) nrels = target;
    int ncols = fb->size + 1; /* +1 for sign */
    gf2_t *mat = gf2_create(nrels, ncols);

    for (int r = 0; r < nrels; r++) {
        mpz_t Qval;
        mpz_init(Qval);
        mpz_set(Qval, full_rels->Qx[r]);
        if (mpz_sgn(Qval) < 0) { gf2_set(mat, r, 0); mpz_neg(Qval, Qval); }

        int e2 = 0;
        while (mpz_even_p(Qval)) { mpz_tdiv_q_2exp(Qval, Qval, 1); e2++; }
        if (e2 & 1) gf2_set(mat, r, 1);

        for (int i = 1; i < fb->size; i++) {
            uint32_t p = fb->prime[i];
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
    fprintf(stderr, "LA: %d deps from %dx%d in %.2fs\n",
            ndeps, nrels, ncols, elapsed() - sieve_time);

    /* === Square Root === */
    for (int d = 0; d < ndeps; d++) {
        mpz_t X, Y, g, prod, rem;
        mpz_inits(X, Y, g, prod, rem, NULL);

        mpz_set_ui(X, 1);
        for (int k = 0; k < dlen[d]; k++) {
            mpz_mul(X, X, full_rels->ax_b[deps[d][k]]);
            mpz_mod(X, X, N);
        }

        mpz_set_ui(prod, 1);
        for (int k = 0; k < dlen[d]; k++) {
            mpz_t aq;
            mpz_init(aq);
            mpz_abs(aq, full_rels->Qx[deps[d][k]]);
            mpz_mul(prod, prod, aq);
            mpz_clear(aq);
        }

        mpz_set(rem, prod);
        int e2 = 0;
        while (mpz_even_p(rem)) { mpz_tdiv_q_2exp(rem, rem, 1); e2++; }
        if (e2 & 1) goto next_dep;

        mpz_set_ui(Y, 1);
        if (e2 / 2 > 0) {
            mpz_set_ui(tmp, 2);
            mpz_powm_ui(tmp, tmp, e2 / 2, N);
            mpz_mul(Y, Y, tmp);
            mpz_mod(Y, Y, N);
        }

        {
            int valid = 1;
            for (int i = 1; i < fb->size && valid; i++) {
                uint32_t p = fb->prime[i];
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
            if (!valid) goto next_dep;
        }

        /* Handle remaining cofactor */
        if (mpz_cmp_ui(rem, 1) != 0) {
            if (mpz_perfect_square_p(rem)) {
                mpz_sqrt(tmp, rem);
                mpz_mod(tmp, tmp, N);
                mpz_mul(Y, Y, tmp);
                mpz_mod(Y, Y, N);
            } else {
                goto next_dep;
            }
        }

        /* Check GCD */
        mpz_sub(tmp, X, Y); mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t o; mpz_init(o); mpz_divexact(o, N, g);
            if (mpz_cmp(g, o) > 0) mpz_swap(g, o);
            gmp_printf("%Zd\n", g);
            fprintf(stderr, "QS-Turbo: factored %dd in %.3fs (%d polys, %d rels)\n",
                    digits, elapsed(), total_polys, full_rels->count);
            mpz_clear(o);
            return 0;
        }
        mpz_add(tmp, X, Y); mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t o; mpz_init(o); mpz_divexact(o, N, g);
            if (mpz_cmp(g, o) > 0) mpz_swap(g, o);
            gmp_printf("%Zd\n", g);
            fprintf(stderr, "QS-Turbo: factored %dd in %.3fs (%d polys, %d rels)\n",
                    digits, elapsed(), total_polys, full_rels->count);
            mpz_clear(o);
            return 0;
        }

        next_dep:
        mpz_clears(X, Y, g, prod, rem, NULL);
    }

    fprintf(stderr, "QS-Turbo: FAILED after %.2fs\n", elapsed());
    printf("FAIL\n");
    return 1;
}
