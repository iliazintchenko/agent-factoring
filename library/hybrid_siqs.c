/*
 * hybrid_siqs.c - Hybrid SIQS combining multi-polynomial batching (SPQS) with bucket sieve
 *
 * Key ideas:
 * 1. SPQS multi-polynomial batching: for each 'a', generate BATCH_POLYS b-values and
 *    sieve them simultaneously, amortizing factor base iteration cost over BATCH_POLYS polys.
 * 2. Bucket sieve for large primes (p >= BLOCKSIZE): pre-scatter (block_idx, offset, logp)
 *    entries to per-block, per-polynomial buckets during polynomial setup, then apply them
 *    during block processing instead of inner-looping over large primes.
 *
 * This combines the strengths of both approaches:
 * - SPQS batching is best at 30-55d (small FB, loop amortization wins)
 * - Bucket sieve is best at 60-65d+ (large primes dominate, direct sieve is slow)
 *
 * Compile: gcc -O3 -march=native -o hybrid_siqs library/hybrid_siqs.c -lgmp -lm
 * Usage: ./hybrid_siqs <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

/* ==================== Configuration ==================== */
#define SEED        42
#define BLOCKSIZE   32768
#define BLOCKBITS   15
#define BLOCKMASK   (BLOCKSIZE - 1)
#define BATCH_POLYS 4   /* polynomials sieved simultaneously per block pass */

#define MAX_FB        100000
#define MAX_RELS      500000
#define MAX_PARTIALS  2000000
#define MAX_A_FACTORS 25
#define MAX_DEPS      64

/* Bucket sieve */
#define BUCKET_INIT_ALLOC 8192

/* ==================== Timing ==================== */
static struct timespec g_start;
static double elapsed(void) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec - g_start.tv_sec) + (now.tv_nsec - g_start.tv_nsec) * 1e-9;
}

/* ==================== Parameters ==================== */
typedef struct {
    int fb_size;
    int num_blocks;    /* blocks per side */
    int lp_mult;       /* large prime bound = fb_max * lp_mult */
    int extra_rels;    /* relations beyond FB size to collect */
    double thresh_adj; /* sieve threshold multiplier */
} params_t;

static params_t get_params(int bits) {
    if (bits <= 80)  return (params_t){60,    1,  30,  30,  0.68};
    if (bits <= 90)  return (params_t){80,    1,  30,  30,  0.70};
    if (bits <= 100) return (params_t){120,   1,  40,  35,  0.72};
    if (bits <= 110) return (params_t){180,   1,  40,  35,  0.73};
    if (bits <= 120) return (params_t){230,   1,  40,  40,  0.74};
    if (bits <= 130) return (params_t){300,   2,  40,  40,  0.75};
    if (bits <= 140) return (params_t){400,   2,  40,  50,  0.76};
    if (bits <= 150) return (params_t){500,   2,  40,  50,  0.77};
    if (bits <= 160) return (params_t){650,   3,  40,  60,  0.78};
    if (bits <= 170) return (params_t){900,   3,  40,  60,  0.79};
    if (bits <= 180) return (params_t){1200,  4,  50,  70,  0.80};
    if (bits <= 190) return (params_t){1700,  5,  50,  80,  0.81};
    if (bits <= 200) return (params_t){2200,  6,  50,  80,  0.82};
    if (bits <= 210) return (params_t){3000,  8,  50,  90,  0.83};
    if (bits <= 220) return (params_t){4000,  10, 60,  100, 0.84};
    if (bits <= 230) return (params_t){5000,  12, 60,  100, 0.84};
    if (bits <= 240) return (params_t){6500,  16, 60,  120, 0.85};
    if (bits <= 250) return (params_t){9000,  20, 70,  120, 0.86};
    if (bits <= 260) return (params_t){12000, 26, 70,  150, 0.86};
    if (bits <= 270) return (params_t){16000, 32, 80,  150, 0.87};
    if (bits <= 280) return (params_t){22000, 40, 80,  200, 0.87};
    if (bits <= 290) return (params_t){30000, 48, 90,  200, 0.88};
    if (bits <= 300) return (params_t){40000, 56, 90,  250, 0.88};
    return                  (params_t){55000, 64, 100, 300, 0.89};
}

/* ==================== Modular Arithmetic ==================== */
static inline unsigned int mod_inverse(unsigned int a, unsigned int m) {
    int old_r = (int)a, r = (int)m, old_s = 1, s = 0;
    while (r) {
        int q = old_r / r;
        int t = r; r = old_r - q * r; old_r = t;
        t = s; s = old_s - q * s; old_s = t;
    }
    if (old_r != 1) return 0;
    return (unsigned int)(((long long)old_s % m + m) % m);
}

/* Tonelli-Shanks modular square root */
static unsigned int sqrt_mod(unsigned int n, unsigned int p) {
    if (n == 0) return 0;
    if (p == 2) return n & 1;
    unsigned long long nn = n % p, m = p;

    unsigned long long r = 1, b = nn, e = (p - 1) / 2;
    while (e) { if (e & 1) r = r * b % m; b = b * b % m; e >>= 1; }
    if (r != 1) return 0;

    if (p % 4 == 3) {
        r = 1; b = nn; e = (p + 1) / 4;
        while (e) { if (e & 1) r = r * b % m; b = b * b % m; e >>= 1; }
        return (unsigned int)r;
    }

    unsigned int Q = p - 1, S = 0;
    while (Q % 2 == 0) { Q /= 2; S++; }

    unsigned int z = 2;
    for (;;) {
        r = 1; b = z; e = (p - 1) / 2;
        while (e) { if (e & 1) r = r * b % m; b = b * b % m; e >>= 1; }
        if (r == m - 1) break;
        z++;
    }

    unsigned long long M_val = S;
    r = 1; b = z; e = Q;
    while (e) { if (e & 1) r = r * b % m; b = b * b % m; e >>= 1; }
    unsigned long long c = r;

    r = 1; b = nn; e = Q;
    while (e) { if (e & 1) r = r * b % m; b = b * b % m; e >>= 1; }
    unsigned long long t = r;

    r = 1; b = nn; e = (Q + 1) / 2;
    while (e) { if (e & 1) r = r * b % m; b = b * b % m; e >>= 1; }
    unsigned long long R = r;

    for (;;) {
        if (t == 1) return (unsigned int)R;
        int i = 0;
        unsigned long long tt = t;
        while (tt != 1) { tt = tt * tt % p; i++; }
        unsigned long long bb = c;
        for (int j = 0; j < (int)M_val - i - 1; j++) bb = bb * bb % p;
        M_val = i;
        c = bb * bb % p;
        t = t * c % p;
        R = R * bb % p;
    }
}

/* ==================== Knuth-Schroeppel Multiplier ==================== */
static int choose_multiplier(mpz_t N) {
    static const int ks[] = {1,2,3,5,6,7,10,11,13,14,15,17,19,21,23,29,31,37,41,43,0};
    double best = -1e30;
    int best_k = 1;
    for (int ki = 0; ks[ki]; ki++) {
        int k = ks[ki];
        mpz_t kN; mpz_init(kN); mpz_mul_ui(kN, N, k);
        double s = -0.5 * log((double)k);
        unsigned long m8 = mpz_fdiv_ui(kN, 8);
        if (m8 == 1) s += 2*log(2.0);
        else if (m8 == 5) s += log(2.0);
        else if (m8 == 3 || m8 == 7) s += 0.5*log(2.0);
        int ps[] = {3,5,7,11,13,17,19,23,29,31,37,41,43,47};
        for (int i = 0; i < 14; i++) {
            if (k % ps[i] == 0) { s += log((double)ps[i]); continue; }
            if (sqrt_mod(mpz_fdiv_ui(kN, ps[i]), ps[i]))
                s += 2.0*log(ps[i])/(ps[i]-1);
        }
        if (s > best) { best = s; best_k = k; }
        mpz_clear(kN);
    }
    return best_k;
}

/* ==================== Factor Base ==================== */
typedef struct {
    unsigned int *prime;
    unsigned int *sqrtN;   /* sqrt(kN) mod p */
    unsigned char *logp;
    int size;
} fb_t;

static fb_t *fb_create(mpz_t kN, int target) {
    fb_t *fb = malloc(sizeof(fb_t));
    int alloc = target + 100;
    fb->prime = malloc(alloc * sizeof(unsigned int));
    fb->sqrtN = malloc(alloc * sizeof(unsigned int));
    fb->logp  = malloc(alloc * sizeof(unsigned char));

    fb->prime[0] = 2; fb->sqrtN[0] = 1; fb->logp[0] = 1; fb->size = 1;

    int bound = target * 30 + 100000;
    char *sv = calloc(bound + 1, 1);
    for (int i = 2; (long)i*i <= bound; i++)
        if (!sv[i]) for (int j = i*i; j <= bound; j += i) sv[j] = 1;

    for (int i = 3; i <= bound && fb->size < target; i += 2) {
        if (sv[i]) continue;
        unsigned long nm = mpz_fdiv_ui(kN, i);
        if (nm == 0) {
            fb->prime[fb->size] = i;
            fb->sqrtN[fb->size] = 0;
            fb->logp[fb->size] = (unsigned char)(log2(i) + 0.5);
            fb->size++;
            continue;
        }
        unsigned int r = sqrt_mod((unsigned int)nm, i);
        if (!r) continue;
        fb->prime[fb->size] = i;
        fb->sqrtN[fb->size] = r;
        fb->logp[fb->size] = (unsigned char)(log2(i) + 0.5);
        fb->size++;
    }
    free(sv);
    return fb;
}

/* ==================== Large Prime Hash Table ==================== */
#define LP_HASH_BITS 22
#define LP_HASH_SIZE (1 << LP_HASH_BITS)

typedef struct lp_entry {
    unsigned long lp;
    int rel_idx;
    struct lp_entry *next;
} lp_entry_t;

typedef struct {
    lp_entry_t **buckets;
    lp_entry_t *pool;
    int used, max;
} lp_hash_t;

static lp_hash_t *lp_create(int max) {
    lp_hash_t *h = calloc(1, sizeof(lp_hash_t));
    h->buckets = calloc(LP_HASH_SIZE, sizeof(lp_entry_t*));
    h->pool = calloc(max, sizeof(lp_entry_t));
    h->max = max;
    return h;
}

static inline unsigned int lp_hashfn(unsigned long lp) {
    return (unsigned int)((lp * 0x9E3779B97F4A7C15ULL) >> (64 - LP_HASH_BITS));
}

static int lp_find(lp_hash_t *h, unsigned long lp) {
    unsigned int idx = lp_hashfn(lp);
    for (lp_entry_t *e = h->buckets[idx]; e; e = e->next)
        if (e->lp == lp) return e->rel_idx;
    return -1;
}

static void lp_insert(lp_hash_t *h, unsigned long lp, int rel_idx) {
    if (h->used >= h->max) return;
    unsigned int idx = lp_hashfn(lp);
    lp_entry_t *e = &h->pool[h->used++];
    e->lp = lp; e->rel_idx = rel_idx;
    e->next = h->buckets[idx];
    h->buckets[idx] = e;
}

/* ==================== Relation Storage ==================== */
typedef struct {
    mpz_t *Y;       /* ax+b values */
    short **exps;   /* exponent vectors (fb_size + 1 entries, last = sign) */
    unsigned long *lp;
    int count, alloc, fb_size, neg_col;
} rels_t;

static rels_t *rels_create(int alloc, int fb_size) {
    rels_t *r = malloc(sizeof(rels_t));
    r->Y    = malloc(alloc * sizeof(mpz_t));
    r->exps = malloc(alloc * sizeof(short*));
    r->lp   = calloc(alloc, sizeof(unsigned long));
    for (int i = 0; i < alloc; i++) {
        mpz_init(r->Y[i]);
        r->exps[i] = calloc(fb_size + 2, sizeof(short));
    }
    r->count = 0; r->alloc = alloc; r->fb_size = fb_size;
    r->neg_col = fb_size;
    return r;
}

/* ==================== Bucket Sieve (per polynomial) ==================== */
/*
 * For BATCH_POLYS simultaneous polynomials, each large prime generates up to
 * BATCH_POLYS bucket entries per hit (one per poly). We store them in separate
 * per-poly bucket arrays to allow independent application during block processing.
 *
 * Each entry is packed: (fb_idx << 16) | block_offset
 * logp is looked up from fb->logp[fb_idx] during application.
 */

typedef struct {
    unsigned int *data;
    int count, alloc;
} bucket_t;

/* bucket_array[poly_idx][block_idx] */
typedef struct {
    bucket_t **poly_blocks;  /* [BATCH_POLYS][num_blocks] */
    int num_polys;
    int num_blocks;
} multi_bucket_t;

static multi_bucket_t *mbucket_create(int num_polys, int num_blocks) {
    multi_bucket_t *mb = malloc(sizeof(multi_bucket_t));
    mb->num_polys  = num_polys;
    mb->num_blocks = num_blocks;
    mb->poly_blocks = malloc(num_polys * sizeof(bucket_t*));
    for (int p = 0; p < num_polys; p++) {
        mb->poly_blocks[p] = calloc(num_blocks, sizeof(bucket_t));
        for (int b = 0; b < num_blocks; b++) {
            mb->poly_blocks[p][b].alloc = BUCKET_INIT_ALLOC;
            mb->poly_blocks[p][b].data  = malloc(BUCKET_INIT_ALLOC * sizeof(unsigned int));
            mb->poly_blocks[p][b].count = 0;
        }
    }
    return mb;
}

static void mbucket_reset(multi_bucket_t *mb, int num_polys) {
    for (int p = 0; p < num_polys; p++)
        for (int b = 0; b < mb->num_blocks; b++)
            mb->poly_blocks[p][b].count = 0;
}

static inline void mbucket_add(multi_bucket_t *mb, int poly_idx, int block_idx,
                                unsigned int fb_idx, unsigned int offset) {
    bucket_t *bkt = &mb->poly_blocks[poly_idx][block_idx];
    if (bkt->count >= bkt->alloc) {
        bkt->alloc *= 2;
        bkt->data = realloc(bkt->data, bkt->alloc * sizeof(unsigned int));
    }
    bkt->data[bkt->count++] = (fb_idx << 16) | offset;
}

/* ==================== GF(2) Linear Algebra ==================== */
typedef unsigned long long u64;
typedef struct { u64 **rows; int nr, nc, fbw, idw, wprow; } gf2_t;

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
    m->rows[r][c/64] ^= (1ULL << (c % 64));
}

static int gf2_solve(gf2_t *m, int ***deps_out, int **dlen_out, int max_deps) {
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
    *deps_out = malloc(max_deps * sizeof(int*));
    *dlen_out = malloc(max_deps * sizeof(int));
    for (int r = piv; r < m->nr && nd < max_deps; r++) {
        int zero = 1;
        for (int w = 0; w < m->fbw && zero; w++) {
            u64 mask = (w < m->fbw-1) ? ~0ULL : (m->nc%64==0 ? ~0ULL : (1ULL << (m->nc%64))-1);
            if (m->rows[r][w] & mask) zero = 0;
        }
        if (!zero) continue;
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
        if (dl > 0) { (*deps_out)[nd] = d; (*dlen_out)[nd] = dl; nd++; }
        else free(d);
    }
    return nd;
}

/* ==================== Trial Division (sieve-informed) ==================== */
/*
 * Given a candidate at x_global (global x in [-M, M)), trial divide Q(x).
 * Uses sieve roots to skip non-matching primes (fast path).
 * Returns 1 for full relation, 2 for partial (single large prime), 0 for neither.
 */
static int trial_divide(
    mpz_t Q, mpz_t Y,
    short *exponents,
    fb_t *fb, int fb_size,
    unsigned int *root1, unsigned int *root2,
    int x_global,
    unsigned long lp_bound, unsigned long *lp_out,
    mpz_t kN, mpz_t a, mpz_t b,
    int neg_col)
{
    /* Y = a*x + b */
    mpz_mul_si(Y, a, x_global);
    mpz_add(Y, Y, b);
    /* Q = Y^2 - kN */
    mpz_mul(Q, Y, Y);
    mpz_sub(Q, Q, kN);

    memset(exponents, 0, (fb_size + 2) * sizeof(short));

    if (mpz_sgn(Q) < 0) {
        exponents[neg_col] = 1;
        mpz_neg(Q, Q);
    }

    /* p=2 */
    while (mpz_even_p(Q)) {
        exponents[0]++;
        mpz_tdiv_q_2exp(Q, Q, 1);
    }

    /* Try fast __int128 path */
    int use_fast = (mpz_sizeinbase(Q, 2) <= 127);
    __int128 Q128 = 0;
    if (use_fast) {
        if (mpz_fits_ulong_p(Q)) {
            Q128 = mpz_get_ui(Q);
        } else {
            mpz_t hi; mpz_init(hi);
            mpz_tdiv_q_2exp(hi, Q, 64);
            Q128 = ((__int128)mpz_get_ui(hi) << 64) | mpz_getlimbn(Q, 0);
            mpz_clear(hi);
        }
    }

    for (int i = 1; i < fb_size; i++) {
        unsigned int p = fb->prime[i];
        if (p < 3) continue;

        if (root1[i] != 0xFFFFFFFF) {
            int xp = x_global % (int)p;
            unsigned int xmod = (unsigned int)(xp < 0 ? xp + (int)p : xp);
            if (xmod != root1[i] && xmod != root2[i]) {
                if (use_fast) { if (Q128 % p != 0) continue; }
                else { if (!mpz_divisible_ui_p(Q, p)) continue; }
            }
        } else {
            if (use_fast) { if (Q128 % p != 0) continue; }
            else { if (!mpz_divisible_ui_p(Q, p)) continue; }
        }

        if (use_fast && Q128 > 0) {
            while (Q128 % p == 0) { exponents[i]++; Q128 /= p; }
            if (Q128 <= (unsigned long long)(-1ULL)) {
                mpz_set_ui(Q, (unsigned long long)Q128);
            } else {
                mpz_set_ui(Q, (unsigned long long)(Q128 >> 64));
                mpz_mul_2exp(Q, Q, 64);
                mpz_add_ui(Q, Q, (unsigned long long)Q128);
            }
        } else {
            do { exponents[i]++; mpz_divexact_ui(Q, Q, p); }
            while (mpz_divisible_ui_p(Q, p));
            if (mpz_sizeinbase(Q, 2) <= 127) {
                use_fast = 1;
                if (mpz_fits_ulong_p(Q)) Q128 = mpz_get_ui(Q);
                else {
                    mpz_t hi; mpz_init(hi);
                    mpz_tdiv_q_2exp(hi, Q, 64);
                    Q128 = ((__int128)mpz_get_ui(hi) << 64) | mpz_getlimbn(Q, 0);
                    mpz_clear(hi);
                }
            }
        }
    }

    /* Sync Q128 back */
    if (use_fast) {
        if (Q128 <= (unsigned long long)(-1ULL))
            mpz_set_ui(Q, (unsigned long long)Q128);
        else {
            mpz_set_ui(Q, (unsigned long long)(Q128 >> 64));
            mpz_mul_2exp(Q, Q, 64);
            mpz_add_ui(Q, Q, (unsigned long long)Q128);
        }
    }

    if (mpz_cmp_ui(Q, 1) == 0) { *lp_out = 0; return 1; }

    if (mpz_fits_ulong_p(Q)) {
        unsigned long cof = mpz_get_ui(Q);
        if (cof > 1 && cof <= lp_bound) {
            mpz_t tmp; mpz_init_set_ui(tmp, cof);
            int is_prime = mpz_probab_prime_p(tmp, 5);
            mpz_clear(tmp);
            if (is_prime) { *lp_out = cof; return 2; }
        }
    }
    return 0;
}

/* ==================== MAIN ==================== */
int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }
    clock_gettime(CLOCK_MONOTONIC, &g_start);

    mpz_t N, kN;
    mpz_inits(N, kN, NULL);
    if (mpz_set_str(N, argv[1], 10) != 0) {
        fprintf(stderr, "Invalid number: %s\n", argv[1]); return 1;
    }

    int digits = (int)mpz_sizeinbase(N, 10);
    int bits   = (int)mpz_sizeinbase(N, 2);

    /* Trial division for small factors */
    for (unsigned long p = 2; p < 100000; p++) {
        if (p > 2 && p % 2 == 0) continue;
        if (mpz_divisible_ui_p(N, p)) {
            mpz_t q; mpz_init(q);
            mpz_divexact_ui(q, N, p);
            if (mpz_cmp_ui(q, 1) > 0) {
                if (mpz_cmp_ui(q, p) < 0) gmp_printf("%lu\n%Zd\n", p, q);
                else gmp_printf("%Zd\n%lu\n", q, p);
            } else {
                printf("%lu\n", p);
            }
            mpz_clear(q);
            return 0;
        }
    }

    /* Choose multiplier */
    int mult = choose_multiplier(N);
    mpz_mul_ui(kN, N, mult);
    int kN_bits = (int)mpz_sizeinbase(kN, 2);

    /* Parameters */
    params_t P = get_params(kN_bits);

    /* Factor base */
    fb_t *fb = fb_create(kN, P.fb_size);
    int M = BLOCKSIZE * P.num_blocks;
    unsigned long lp_bound = (unsigned long)fb->prime[fb->size-1] * P.lp_mult;
    int target = fb->size + P.extra_rels;
    int total_blocks = 2 * P.num_blocks;

    /* Sieve threshold */
    double log2_kN  = kN_bits;
    double log2_M   = log2((double)M);
    double log2_a   = log2_kN / 2.0;
    double log2_Qmax = log2_a + log2_M;
    int threshold = (int)(log2_Qmax * P.thresh_adj);

    /* Bucket sieve cutoff: primes >= BLOCKSIZE use buckets */
    int bucket_start = fb->size; /* default: no bucket sieve */
    for (int i = 0; i < fb->size; i++) {
        if (fb->prime[i] >= BLOCKSIZE) { bucket_start = i; break; }
    }

    /* Sieve start: skip p < 5 */
    int sieve_start = 0;
    for (int i = 0; i < fb->size; i++) {
        if (fb->prime[i] >= 5) { sieve_start = i; break; }
    }

    fprintf(stderr,
        "HybridSIQS: %dd (%db), k=%d, FB=%d, blocks=%d, M=%d, thresh=%d, "
        "LP=%lu, target=%d, bucket_start=%d, batch=%d\n",
        digits, bits, mult, fb->size, total_blocks, M, threshold,
        lp_bound, target, bucket_start, BATCH_POLYS);

    /* Sieve arrays: one per polynomial in batch */
    unsigned char *sieve[BATCH_POLYS];
    for (int b = 0; b < BATCH_POLYS; b++)
        sieve[b] = malloc(BLOCKSIZE);

    /* Bucket sieve: BATCH_POLYS sets of per-block buckets */
    multi_bucket_t *mbuckets = mbucket_create(BATCH_POLYS, total_blocks);

    /* Sieve roots per polynomial */
    unsigned int *root1[BATCH_POLYS], *root2[BATCH_POLYS];
    for (int b = 0; b < BATCH_POLYS; b++) {
        root1[b] = malloc(fb->size * sizeof(unsigned int));
        root2[b] = malloc(fb->size * sizeof(unsigned int));
    }

    /* Polynomial state */
    mpz_t a, bs[BATCH_POLYS], cs[BATCH_POLYS], B_vals[MAX_A_FACTORS];
    mpz_init(a);
    for (int b = 0; b < BATCH_POLYS; b++) { mpz_init(bs[b]); mpz_init(cs[b]); }
    for (int j = 0; j < MAX_A_FACTORS; j++) mpz_init(B_vals[j]);

    /* ainv_data[j * fb_size + i] = 2*B_j * (2*a)^{-1} mod p_i (for Gray code incremental updates) */
    unsigned int *ainv_data = malloc((size_t)MAX_A_FACTORS * fb->size * sizeof(unsigned int));

    /* Relation storage */
    rels_t *full_rels = rels_create(MAX_RELS,     fb->size);
    rels_t *part_rels = rels_create(MAX_PARTIALS, fb->size);
    lp_hash_t *lp_hash = lp_create(MAX_PARTIALS);

    /* RNG */
    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, SEED);

    /* Temporaries */
    mpz_t Q_val, Y_val, tmp, tmp2;
    mpz_inits(Q_val, Y_val, tmp, tmp2, NULL);
    short *tmp_exps = calloc(fb->size + 2, sizeof(short));

    int total_polys  = 0;
    int combined_rels = 0;
    int a_count      = 0;
    int a_idx[MAX_A_FACTORS];
    int num_a_factors = 0;

    /* ================================================================
     * MAIN SIEVING LOOP
     * ================================================================ */
    while (full_rels->count < target) {
        double t = elapsed();
        if (t > 300) {
            fprintf(stderr, "TIMEOUT at %.1fs with %d/%d relations\n", t, full_rels->count, target);
            break;
        }

        /* ---- Generate new 'a' coefficient ---- */
        {
            mpz_t tgt; mpz_init(tgt);
            mpz_mul_ui(tgt, kN, 2);
            mpz_sqrt(tgt, tgt);
            mpz_tdiv_q_ui(tgt, tgt, M);

            double log_tgt = mpz_sizeinbase(tgt, 2) * log(2.0);

            int lo = fb->size / 4, hi = 3 * fb->size / 4;
            if (lo < 2) lo = 2;
            if (hi <= lo + 3) hi = fb->size - 1;

            double avg_logp = 0; int cnt = 0;
            for (int i = lo; i < hi; i++) {
                if (fb->sqrtN[i] == 0) continue;
                avg_logp += log((double)fb->prime[i]);
                cnt++;
            }
            if (cnt == 0) { mpz_clear(tgt); break; }
            avg_logp /= cnt;

            int s = (int)(log_tgt / avg_logp + 0.5);
            if (s < 3) s = 3;
            if (s > MAX_A_FACTORS) s = MAX_A_FACTORS;
            if (s > hi - lo) s = hi - lo;
            num_a_factors = s;

            double best_ratio = 1e30;
            int best[MAX_A_FACTORS];

            for (int att = 0; att < 50; att++) {
                mpz_set_ui(a, 1);
                int idx[MAX_A_FACTORS]; int ok = 1;
                for (int i = 0; i < s && ok; i++) {
                    int tries = 0, good;
                    do {
                        idx[i] = lo + (int)gmp_urandomm_ui(rng, hi - lo);
                        good = 1;
                        for (int j = 0; j < i; j++) if (idx[j] == idx[i]) { good = 0; break; }
                        if (fb->sqrtN[idx[i]] == 0) good = 0;
                        tries++;
                    } while (!good && tries < 100);
                    if (!good) { ok = 0; break; }
                    mpz_mul_ui(a, a, fb->prime[idx[i]]);
                }
                if (!ok) continue;

                double ratio;
                if (mpz_cmp(a, tgt) > 0) { mpz_tdiv_q(tmp, a, tgt); ratio = mpz_get_d(tmp); }
                else                      { mpz_tdiv_q(tmp, tgt, a); ratio = mpz_get_d(tmp); }
                if (ratio < best_ratio) { best_ratio = ratio; memcpy(best, idx, s * sizeof(int)); }
                if (ratio < 1.5) break;
            }

            memcpy(a_idx, best, s * sizeof(int));
            mpz_set_ui(a, 1);
            for (int i = 0; i < s; i++) mpz_mul_ui(a, a, fb->prime[a_idx[i]]);
            mpz_clear(tgt);
            a_count++;

            /* Compute B values for SIQS self-initialization */
            for (int j = 0; j < s; j++) {
                unsigned int qj = fb->prime[a_idx[j]];
                unsigned int rj = fb->sqrtN[a_idx[j]];
                mpz_t a_q; mpz_init(a_q);
                mpz_divexact_ui(a_q, a, qj);
                unsigned long aqmod = mpz_fdiv_ui(a_q, qj);
                unsigned int iv = mod_inverse((unsigned int)aqmod, qj);
                mpz_mul_ui(B_vals[j], a_q, ((unsigned long)rj * iv) % qj);
                mpz_clear(a_q);
            }

            /* Precompute ainv[j][i] = 2*B_j * (2*a)^{-1} mod p_i (Gray code delta step) */
            for (int j = 0; j < s; j++) {
                for (int i = 0; i < fb->size; i++) {
                    unsigned int p = fb->prime[i];
                    if (p < 3) { ainv_data[j * fb->size + i] = 0; continue; }
                    unsigned long am = mpz_fdiv_ui(a, p);
                    if (am == 0) { ainv_data[j * fb->size + i] = 0; continue; }
                    unsigned int ai = mod_inverse((unsigned int)((2UL * am) % p), p);
                    if (ai == 0) { ainv_data[j * fb->size + i] = 0; continue; }
                    unsigned long Bm = mpz_fdiv_ui(B_vals[j], p);
                    ainv_data[j * fb->size + i] = (unsigned int)((2ULL * ai % p * Bm) % p);
                }
            }
        }

        /* ---- Enumerate b-values via Gray code, batching BATCH_POLYS at a time ---- */
        int num_b = 1 << (num_a_factors - 1);

        for (int b_start = 0; b_start < num_b && full_rels->count < target; b_start += BATCH_POLYS) {
            int batch = BATCH_POLYS;
            if (b_start + batch > num_b) batch = num_b - b_start;

            /* Progress reporting */
            if (total_polys > 0 && total_polys % 500 == 0) {
                fprintf(stderr,
                    "  poly=%d, rels=%d/%d (full=%d+%d), part=%d, t=%.1fs\n",
                    total_polys, full_rels->count, target,
                    full_rels->count - combined_rels, combined_rels,
                    part_rels->count, elapsed());
            }

            /* Compute b-values and sieve roots for each polynomial in the batch */
            int valid_batch = 0; /* actual number with valid b */
            for (int bi = 0; bi < batch; bi++) {
                int b_idx = b_start + bi;
                int gray  = b_idx ^ (b_idx >> 1);

                /* b = sum(±B_j) per Gray code bits */
                mpz_set_ui(bs[bi], 0);
                for (int j = 0; j < num_a_factors; j++) {
                    if (gray & (1 << j)) mpz_add(bs[bi], bs[bi], B_vals[j]);
                    else                  mpz_sub(bs[bi], bs[bi], B_vals[j]);
                }

                /* Verify b^2 ≡ kN (mod a) */
                mpz_mul(tmp, bs[bi], bs[bi]);
                mpz_sub(tmp, tmp, kN);
                mpz_mod(tmp, tmp, a);
                if (mpz_sgn(tmp) != 0) {
                    mpz_neg(bs[bi], bs[bi]);
                    mpz_mul(tmp, bs[bi], bs[bi]);
                    mpz_sub(tmp, tmp, kN);
                    mpz_mod(tmp, tmp, a);
                    if (mpz_sgn(tmp) != 0) { batch = bi; break; }
                }

                /* c = (b^2 - kN) / a */
                mpz_mul(cs[bi], bs[bi], bs[bi]);
                mpz_sub(cs[bi], cs[bi], kN);
                mpz_divexact(cs[bi], cs[bi], a);

                /* Compute sieve roots for this polynomial */
                for (int i = 0; i < fb->size; i++) {
                    unsigned int p = fb->prime[i];
                    if (p < 3 || fb->sqrtN[i] == 0) {
                        root1[bi][i] = 0xFFFFFFFF;
                        root2[bi][i] = 0xFFFFFFFF;
                        continue;
                    }
                    unsigned long am = mpz_fdiv_ui(a, p);
                    if (am == 0) {
                        /* p | a: Q(x) ≡ 2bx + c (mod p), single root */
                        unsigned long bm = mpz_fdiv_ui(bs[bi], p);
                        unsigned long cm = mpz_fdiv_ui(cs[bi], p);
                        if (bm == 0) {
                            root1[bi][i] = 0xFFFFFFFF;
                            root2[bi][i] = 0xFFFFFFFF;
                            continue;
                        }
                        unsigned int inv2b = mod_inverse((unsigned int)((2UL * bm) % p), p);
                        unsigned int root  = (unsigned int)((unsigned long)(p - cm) % p * inv2b % p);
                        root1[bi][i] = root;
                        root2[bi][i] = root;
                        continue;
                    }
                    unsigned int ai = mod_inverse((unsigned int)am, p);
                    if (ai == 0) {
                        root1[bi][i] = 0xFFFFFFFF;
                        root2[bi][i] = 0xFFFFFFFF;
                        continue;
                    }
                    unsigned long bm = mpz_fdiv_ui(bs[bi], p);
                    unsigned int r   = fb->sqrtN[i];
                    root1[bi][i] = (unsigned int)((unsigned long)ai * ((r + p - bm) % p) % p);
                    root2[bi][i] = (unsigned int)((unsigned long)ai * ((p - r + p - bm) % p) % p);
                }
                valid_batch = bi + 1;
            }
            if (valid_batch == 0) continue;
            batch = valid_batch;
            total_polys++;

            /* ---- Fill buckets for large primes (one entry per poly per hit) ---- */
            mbucket_reset(mbuckets, batch);

            for (int i = bucket_start; i < fb->size; i++) {
                unsigned int p = fb->prime[i];

                for (int bi = 0; bi < batch; bi++) {
                    if (root1[bi][i] == 0xFFFFFFFF) continue;

                    /* Two roots; scatter across sieve interval [-M, M) */
                    for (int ri = 0; ri < 2; ri++) {
                        unsigned int root = (ri == 0) ? root1[bi][i] : root2[bi][i];
                        if (ri == 1 && root1[bi][i] == root2[bi][i]) break; /* single root */

                        /* First hit >= -M with x ≡ root (mod p) */
                        long off = ((long)root - (-(long)M)) % (long)p;
                        if (off < 0) off += p;
                        long pos = -(long)M + off;

                        while (pos < (long)M) {
                            int blk_idx = (int)((pos + M) >> BLOCKBITS);
                            int blk_off = (int)((pos + M) & BLOCKMASK);
                            if (blk_idx >= 0 && blk_idx < total_blocks)
                                mbucket_add(mbuckets, bi, blk_idx, (unsigned int)i, (unsigned int)blk_off);
                            pos += p;
                        }
                    }
                }
            }

            /* ---- Process each sieve block ---- */
            for (int blk = 0; blk < total_blocks; blk++) {
                int block_start_x = -(int)M + blk * BLOCKSIZE;

                /* Initialize all sieve arrays */
                for (int bi = 0; bi < batch; bi++)
                    memset(sieve[bi], 0, BLOCKSIZE);

                /* Small/medium primes: sieve ALL batch polynomials */
                for (int i = sieve_start; i < bucket_start && i < fb->size; i++) {
                    unsigned int p  = fb->prime[i];
                    unsigned char lp = fb->logp[i];

                    /* Iterate over batch polynomials; inner loop is fast per-prime */
                    for (int bi = 0; bi < batch; bi++) {
                        if (root1[bi][i] == 0xFFFFFFFF) continue;

                        /* Root 1 */
                        long off1 = ((long)root1[bi][i] - (long)block_start_x) % (long)p;
                        if (off1 < 0) off1 += p;
                        for (int j = (int)off1; j < BLOCKSIZE; j += (int)p)
                            sieve[bi][j] += lp;

                        /* Root 2 (if distinct) */
                        if (root1[bi][i] != root2[bi][i]) {
                            long off2 = ((long)root2[bi][i] - (long)block_start_x) % (long)p;
                            if (off2 < 0) off2 += p;
                            for (int j = (int)off2; j < BLOCKSIZE; j += (int)p)
                                sieve[bi][j] += lp;
                        }
                    }
                }

                /* Large primes: apply per-polynomial bucket entries */
                for (int bi = 0; bi < batch; bi++) {
                    bucket_t *bkt = &mbuckets->poly_blocks[bi][blk];
                    for (int k = 0; k < bkt->count; k++) {
                        unsigned int packed = bkt->data[k];
                        unsigned int fb_idx = packed >> 16;
                        unsigned int offset = packed & 0xFFFF;
                        sieve[bi][offset] += fb->logp[fb_idx];
                    }
                }

                /* ---- Scan candidates across all batch polynomials ---- */
                for (int bi = 0; bi < batch; bi++) {
                    unsigned char *sv = sieve[bi];

                    for (int j = 0; j < BLOCKSIZE; j += 4) {
                        /* Fast 4-byte threshold check */
                        unsigned int v = *(unsigned int*)(sv + j);
                        if (((v & 0xFF)         < (unsigned)threshold) &&
                            (((v >> 8)  & 0xFF)  < (unsigned)threshold) &&
                            (((v >> 16) & 0xFF)  < (unsigned)threshold) &&
                            (((v >> 24) & 0xFF)  < (unsigned)threshold))
                            continue;

                        for (int jj = j; jj < j + 4 && jj < BLOCKSIZE; jj++) {
                            if (sv[jj] < (unsigned char)threshold) continue;

                            int x_global = block_start_x + jj;
                            unsigned long lp_val = 0;

                            int result = trial_divide(
                                Q_val, Y_val, tmp_exps,
                                fb, fb->size,
                                root1[bi], root2[bi],
                                x_global,
                                lp_bound, &lp_val,
                                kN, a, bs[bi],
                                full_rels->neg_col);

                            if (result == 1) {
                                int idx = full_rels->count;
                                if (idx < full_rels->alloc) {
                                    mpz_set(full_rels->Y[idx], Y_val);
                                    memcpy(full_rels->exps[idx], tmp_exps, (fb->size + 2) * sizeof(short));
                                    full_rels->lp[idx] = 0;
                                    full_rels->count++;
                                }
                            } else if (result == 2) {
                                int match = lp_find(lp_hash, lp_val);
                                if (match >= 0) {
                                    int fidx = full_rels->count;
                                    if (fidx < full_rels->alloc) {
                                        mpz_mul(full_rels->Y[fidx], Y_val, part_rels->Y[match]);
                                        mpz_mod(full_rels->Y[fidx], full_rels->Y[fidx], N);
                                        for (int e = 0; e <= fb->size; e++)
                                            full_rels->exps[fidx][e] =
                                                tmp_exps[e] + part_rels->exps[match][e];
                                        full_rels->lp[fidx] = lp_val;
                                        full_rels->count++;
                                        combined_rels++;
                                    }
                                } else {
                                    int pidx = part_rels->count;
                                    if (pidx < part_rels->alloc) {
                                        mpz_set(part_rels->Y[pidx], Y_val);
                                        memcpy(part_rels->exps[pidx], tmp_exps, (fb->size + 2) * sizeof(short));
                                        part_rels->lp[pidx] = lp_val;
                                        part_rels->count++;
                                        lp_insert(lp_hash, lp_val, pidx);
                                    }
                                }
                            }
                        }
                    }
                }
            } /* end block loop */
        } /* end batch loop */
    } /* end main sieving loop */

    fprintf(stderr,
        "Sieving done: %d relations (%d full + %d combined) in %.2fs\n",
        full_rels->count, full_rels->count - combined_rels, combined_rels, elapsed());

    /* ================================================================
     * LINEAR ALGEBRA
     * ================================================================ */
    if (full_rels->count < fb->size + 1) {
        fprintf(stderr, "Not enough relations: %d < %d\n", full_rels->count, fb->size + 1);
        return 1;
    }

    int nrels = full_rels->count;
    int ncols = fb->size + 1; /* +1 for sign column */

    fprintf(stderr, "Building %d x %d GF(2) matrix...\n", nrels, ncols);
    gf2_t *mat = gf2_create(nrels, ncols);
    for (int r = 0; r < nrels; r++)
        for (int c = 0; c < ncols; c++)
            if (full_rels->exps[r][c] & 1)
                gf2_set(mat, r, c);

    fprintf(stderr, "Solving GF(2) system...\n");
    int **deps; int *dlen;
    int ndeps = gf2_solve(mat, &deps, &dlen, MAX_DEPS);
    fprintf(stderr, "Found %d dependencies, trying square root...\n", ndeps);

    /* ================================================================
     * SQUARE ROOT EXTRACTION (exponent-based, from spqs.c / fast_siqs.c)
     * ================================================================ */
    mpz_t X, Y2, g;
    mpz_inits(X, Y2, g, NULL);
    int found = 0;

    for (int d = 0; d < ndeps && !found; d++) {
        mpz_set_ui(X, 1);
        int *exps = calloc(ncols + 2, sizeof(int));

        for (int i = 0; i < dlen[d]; i++) {
            int ri = deps[d][i];
            mpz_mul(X, X, full_rels->Y[ri]);
            mpz_mod(X, X, N);
            for (int c = 0; c <= fb->size; c++)
                exps[c] += full_rels->exps[ri][c];
        }

        /* All exponents must be even (includes sign column) */
        int all_even = 1;
        for (int c = 0; c <= fb->size; c++)
            if (exps[c] & 1) { all_even = 0; break; }
        if (!all_even) { free(exps); continue; }

        /* Y = product of p_i^(exp_i/2) mod N */
        mpz_set_ui(Y2, 1);
        for (int c = 0; c < fb->size; c++) {
            if (exps[c] <= 0) continue;
            int half = exps[c] / 2;
            mpz_set_ui(tmp, fb->prime[c]);
            mpz_powm_ui(tmp, tmp, half, N);
            mpz_mul(Y2, Y2, tmp);
            mpz_mod(Y2, Y2, N);
        }

        /* For combined LP relations, each contributes lp^2 to Q product => lp^1 to sqrt */
        for (int i = 0; i < dlen[d]; i++) {
            int ri = deps[d][i];
            if (full_rels->lp[ri] > 1) {
                mpz_set_ui(tmp, full_rels->lp[ri]);
                mpz_mul(Y2, Y2, tmp);
                mpz_mod(Y2, Y2, N);
            }
        }

        /* gcd(X ± Y2, N) */
        mpz_sub(tmp, X, Y2);
        mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t cofactor; mpz_init(cofactor);
            mpz_divexact(cofactor, N, g);
            if (mpz_cmp(g, cofactor) > 0) mpz_swap(g, cofactor);
            gmp_printf("%Zd\n%Zd\n", g, cofactor);
            mpz_clear(cofactor);
            found = 1;
        }
        if (!found) {
            mpz_add(tmp, X, Y2);
            mpz_gcd(g, tmp, N);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                mpz_t cofactor; mpz_init(cofactor);
                mpz_divexact(cofactor, N, g);
                if (mpz_cmp(g, cofactor) > 0) mpz_swap(g, cofactor);
                gmp_printf("%Zd\n%Zd\n", g, cofactor);
                mpz_clear(cofactor);
                found = 1;
            }
        }

        free(exps);
    }

    fprintf(stderr, "HybridSIQS: %s in %.3fs\n", found ? "factored" : "FAILED", elapsed());
    if (!found) { fprintf(stderr, "FAILED: no non-trivial factor found\n"); return 1; }
    return 0;
}
