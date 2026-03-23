/*
<<<<<<< HEAD
 * SIQS with Bucket Sieve - Optimized Implementation
 *
 * Key optimizations over spqs.c:
 * 1. BUCKET SIEVE for large primes (p >= BLOCK_SIZE): pre-sorts sieve hits
 *    by block, converting random memory access to sequential access
 * 2. Gray code polynomial switching with incremental root updates
 * 3. Root-aware trial division (only test primes whose roots match x)
 * 4. SLP (Single Large Prime) with hash table matching
 * 5. Cache-optimized 32KB block sieve
=======
 * siqs_bucket.c - Optimized SIQS with bucket sieve + DLP
 *
 * Based on the proven SPQS code with:
 * 1. Bucket sieve for large FB primes (p > SIEVE_BLOCK)
 * 2. DLP (Double Large Prime) with Pollard rho cofactorization
 * 3. Sieve-informed trial division (only check primes whose roots match)
 * 4. Improved parameter tuning
 * 5. Multi-polynomial batch sieving (from SPQS)
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
 *
 * Compile: gcc -O3 -march=native -o siqs_bucket library/siqs_bucket.c -lgmp -lm
 * Usage: ./siqs_bucket <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
<<<<<<< HEAD
#include <stdint.h>
#include <gmp.h>

#define SEED 42
#define BLOCK_SIZE 32768       /* 32KB = L1 cache line friendly */
#define BLOCK_BITS 15          /* log2(BLOCK_SIZE) */
#define MAX_FB 80000
#define MAX_A_FACTORS 20
#define MAX_RELS 300000
#define MAX_PARTIALS 1500000
#define BUCKET_ALLOC 2048      /* max entries per bucket per slice */
#define SMALL_PRIME_CUTOFF 256 /* primes below this: skip sieve (small logp) */
=======
#include <gmp.h>

#define SEED 42
#define SIEVE_BLOCK 32768
#define MAX_FB 100000
#define MAX_A_FACTORS 25
#define MAX_RELS 500000
#define MAX_PARTIALS 2000000
#define BATCH_POLYS 4
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)

static struct timespec g_start;
static double elapsed(void) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec - g_start.tv_sec) + (now.tv_nsec - g_start.tv_nsec) / 1e9;
}

/* ==================== Modular Arithmetic ==================== */
static unsigned int mod_inverse(unsigned int a, unsigned int m) {
    int old_r = (int)a, r = (int)m, old_s = 1, s = 0;
<<<<<<< HEAD
    while (r) { int q = old_r / r, t = r; r = old_r - q * r; old_r = t; t = s; s = old_s - q * s; old_s = t; }
=======
    while (r) { int q = old_r / r, t; t = r; r = old_r - q * r; old_r = t; t = s; s = old_s - q * s; old_s = t; }
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
    if (old_r != 1) return 0;
    return (unsigned int)(((long long)old_s % m + m) % m);
}

static unsigned int sqrt_mod(unsigned int n, unsigned int p) {
    if (n == 0) return 0;
    if (p == 2) return n & 1;
    unsigned long long b, e, m = p, r;
    b = n % p; e = (p - 1) / 2; r = 1;
    { unsigned long long bb = b, ee = e; while (ee) { if (ee & 1) r = (r * bb) % m; bb = (bb * bb) % m; ee >>= 1; } }
    if (r != 1) return 0;
    if (p % 4 == 3) { b = n % p; e = (p + 1) / 4; r = 1; while (e) { if (e & 1) r = (r * b) % m; b = (b * b) % m; e >>= 1; } return (unsigned int)r; }
    unsigned int Q = p - 1, S = 0;
    while (Q % 2 == 0) { Q /= 2; S++; }
    unsigned int z = 2;
    while (1) { b = z; e = (p - 1) / 2; r = 1; while (e) { if (e & 1) r = (r * b) % m; b = (b * b) % m; e >>= 1; } if (r == (unsigned long long)(p - 1)) break; z++; }
    unsigned long long M_val = S;
    b = z; e = Q; unsigned long long c = 1; while (e) { if (e & 1) c = (c * b) % m; b = (b * b) % m; e >>= 1; }
    b = n % p; e = Q; unsigned long long t = 1; while (e) { if (e & 1) t = (t * b) % m; b = (b * b) % m; e >>= 1; }
    b = n % p; e = (Q + 1) / 2; unsigned long long R = 1; while (e) { if (e & 1) R = (R * b) % m; b = (b * b) % m; e >>= 1; }
<<<<<<< HEAD
    while (1) { if (t == 1) return (unsigned int)R; int i = 0; unsigned long long tt = t; while (tt != 1) { tt = (tt * tt) % p; i++; } unsigned long long bb2 = c; for (int j = 0; j < (int)M_val - i - 1; j++) bb2 = (bb2 * bb2) % p; M_val = i; c = (bb2 * bb2) % p; t = (t * c) % p; R = (R * bb2) % p; }
}

/* ==================== Multiplier (Knuth-Schroeppel) ==================== */
=======
    while (1) { if (t == 1) return (unsigned int)R; int i = 0; unsigned long long tt = t; while (tt != 1) { tt = (tt * tt) % p; i++; } unsigned long long bb = c; for (int j = 0; j < (int)M_val - i - 1; j++) bb = (bb * bb) % p; M_val = i; c = (bb * bb) % p; t = (t * c) % p; R = (R * bb) % p; }
}

/* ==================== Multiplier ==================== */
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
static int choose_multiplier(mpz_t N) {
    static const int ks[] = {1,2,3,5,6,7,10,11,13,14,15,17,19,21,23,29,31,37,41,43,0};
    double best = -1e30; int best_k = 1;
    for (int ki = 0; ks[ki]; ki++) {
        int k = ks[ki]; mpz_t kN; mpz_init(kN); mpz_mul_ui(kN, N, k);
        double s = -0.5 * log((double)k);
        unsigned long m8 = mpz_fdiv_ui(kN, 8);
<<<<<<< HEAD
        if (m8 == 1) s += 2*log(2.0); else if (m8 == 5) s += log(2.0);
        else if (m8 == 3 || m8 == 7) s += 0.5*log(2.0);
        int ps[] = {3,5,7,11,13,17,19,23,29,31,37,41,43,47};
        for (int i = 0; i < 14; i++) {
            if (k % ps[i] == 0) { s += log(ps[i]); continue; }
            if (sqrt_mod(mpz_fdiv_ui(kN, ps[i]), ps[i])) s += 2.0*log(ps[i])/(ps[i]-1);
        }
=======
        if (m8 == 1) s += 2*log(2.0); else if (m8 == 5) s += log(2.0); else if (m8 == 3 || m8 == 7) s += 0.5*log(2.0);
        int ps[] = {3,5,7,11,13,17,19,23,29,31,37,41,43,47};
        for (int i = 0; i < 14; i++) { if (k % ps[i] == 0) continue; if (sqrt_mod(mpz_fdiv_ui(kN, ps[i]), ps[i])) s += 2.0*log(ps[i])/(ps[i]-1); }
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
        if (s > best) { best = s; best_k = k; }
        mpz_clear(kN);
    }
    return best_k;
}

/* ==================== Factor Base ==================== */
<<<<<<< HEAD
typedef struct {
    unsigned int *prime;
    unsigned int *root1;     /* Tonelli-Shanks sqrt of kN mod p */
    unsigned char *logp;
    int size;
    int med_start;           /* index where p >= SMALL_PRIME_CUTOFF */
    int large_start;         /* index where p >= BLOCK_SIZE */
} fb_t;

static fb_t *fb_create(mpz_t kN, int target) {
    fb_t *fb = calloc(1, sizeof(fb_t));
    int alloc = target + 100;
    fb->prime = malloc(alloc * sizeof(unsigned int));
    fb->root1 = malloc(alloc * sizeof(unsigned int));
    fb->logp = malloc(alloc * sizeof(unsigned char));

    fb->prime[0] = 2; fb->root1[0] = 1; fb->logp[0] = 1; fb->size = 1;
    int bound = target * 30 + 50000;
=======
typedef struct { unsigned int *prime; unsigned int *root; unsigned char *logp; int size; } fb_t;

static fb_t *fb_create(mpz_t kN, int target) {
    fb_t *fb = malloc(sizeof(fb_t));
    int alloc = target + 20;
    fb->prime = malloc(alloc * sizeof(unsigned int));
    fb->root = malloc(alloc * sizeof(unsigned int));
    fb->logp = malloc(alloc * sizeof(unsigned char));
    fb->prime[0] = 2; fb->root[0] = 1; fb->logp[0] = 1; fb->size = 1;
    int bound = target * 30 + 100000;
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
    char *sv = calloc(bound + 1, 1);
    for (int i = 2; (long)i*i <= bound; i++) if (!sv[i]) for (int j = i*i; j <= bound; j += i) sv[j] = 1;
    for (int i = 3; i <= bound && fb->size < target; i += 2) {
        if (sv[i]) continue;
        unsigned long nm = mpz_fdiv_ui(kN, i);
<<<<<<< HEAD
        if (nm == 0) {
            fb->prime[fb->size] = i; fb->root1[fb->size] = 0;
            fb->logp[fb->size] = (unsigned char)(log2(i)+0.5); fb->size++; continue;
        }
        unsigned int r = sqrt_mod((unsigned int)nm, i);
        if (!r) continue;
        fb->prime[fb->size] = i; fb->root1[fb->size] = r;
        fb->logp[fb->size] = (unsigned char)(log2(i)+0.5); fb->size++;
    }
    free(sv);

    /* Mark boundaries */
    fb->med_start = fb->size;
    fb->large_start = fb->size;
    for (int i = 0; i < fb->size; i++) {
        if (fb->med_start == fb->size && fb->prime[i] >= SMALL_PRIME_CUTOFF)
            fb->med_start = i;
        if (fb->large_start == fb->size && fb->prime[i] >= BLOCK_SIZE)
            fb->large_start = i;
    }
    return fb;
}

/* ==================== Bucket Sieve Structures ==================== */
/*
 * For large primes (p >= BLOCK_SIZE), each prime hits at most 2 positions
 * per block (one per root). We pre-sort hits by block number.
 * Each entry: upper 16 bits = FB index offset, lower 16 bits = sieve offset.
 */
typedef struct {
    uint32_t *data;       /* flat array: [block0 entries][block1 entries]... */
    uint32_t *count;      /* number of entries per block */
    uint8_t *slice_logp;  /* logp per slice */
    int num_blocks;
    int num_slices;
    int max_slices;
    uint32_t *slice_bound; /* FB index start of each slice */
} bucket_t;

static bucket_t *bucket_create(int nblocks, int max_slices) {
    bucket_t *b = calloc(1, sizeof(bucket_t));
    b->num_blocks = nblocks;
    b->max_slices = max_slices;
    /* Each slice has nblocks buckets, each with BUCKET_ALLOC entries */
    b->data = malloc((size_t)max_slices * nblocks * BUCKET_ALLOC * sizeof(uint32_t));
    b->count = calloc((size_t)max_slices * nblocks, sizeof(uint32_t));
    b->slice_logp = malloc(max_slices * sizeof(uint8_t));
    b->slice_bound = malloc(max_slices * sizeof(uint32_t));
    b->num_slices = 0;
    return b;
}

static void bucket_reset(bucket_t *b) {
    memset(b->count, 0, (size_t)b->max_slices * b->num_blocks * sizeof(uint32_t));
    b->num_slices = 0;
}

static void bucket_free(bucket_t *b) {
    free(b->data); free(b->count); free(b->slice_logp); free(b->slice_bound); free(b);
}

/* Fill buckets for large primes. soln1/soln2 are the two sieve roots. */
static void bucket_fill(bucket_t *bkt, fb_t *fb, unsigned int *soln1, unsigned int *soln2,
                         int nblocks, int sieve_offset) {
    bucket_reset(bkt);

    int fb_start = fb->large_start;
    int fb_end = fb->size;
    int total_blocks = nblocks;

    int slice = 0;
    int bound_val = fb_start;
    bkt->slice_bound[0] = fb_start;
    bkt->slice_logp[0] = fb->logp[fb_start < fb_end ? fb_start : 0];

    for (int i = fb_start; i < fb_end; i++) {
        unsigned int p = fb->prime[i];
        if (soln1[i] == 0xFFFFFFFF) continue;

        /* Check if we need a new slice (bucket getting full or 65K FB indices) */
        if (i - bound_val >= 65535 || (slice < bkt->max_slices - 1)) {
            /* Check max bucket fullness for current slice */
            int max_count = 0;
            for (int bl = 0; bl < total_blocks; bl++) {
                int idx = slice * total_blocks + bl;
                if ((int)bkt->count[idx] > max_count) max_count = bkt->count[idx];
            }
            if (max_count > (int)(BUCKET_ALLOC * 3 / 4) || i - bound_val >= 65535) {
                slice++;
                if (slice >= bkt->max_slices) { slice--; break; }
                bound_val = i;
                bkt->slice_bound[slice] = i;
                bkt->slice_logp[slice] = fb->logp[i];
            }
        }

        /* For each root, compute which block it falls in */
        for (int rt = 0; rt < 2; rt++) {
            unsigned int root = (rt == 0) ? soln1[i] : soln2[i];
            if (rt == 1 && soln1[i] == soln2[i]) continue;

            /* root is the position in [0, sieve_interval) where this prime divides Q(x).
             * We need to adjust for the sieve offset. */
            long pos = ((long)root - sieve_offset);
            /* Reduce to [0, p) range first, then find positions in [0, nblocks*BLOCK_SIZE) */
            pos = pos % (long)p;
            if (pos < 0) pos += p;

            /* Since p >= BLOCK_SIZE, at most one hit per block for this root */
            while (pos < (long)total_blocks * BLOCK_SIZE) {
                int bnum = (int)(pos >> BLOCK_BITS);
                int offset = (int)(pos & (BLOCK_SIZE - 1));
                int idx = slice * total_blocks + bnum;
                if (bkt->count[idx] < BUCKET_ALLOC) {
                    bkt->data[(size_t)idx * BUCKET_ALLOC + bkt->count[idx]] =
                        ((uint32_t)(i - bound_val) << 16) | (uint32_t)offset;
                    bkt->count[idx]++;
                }
                pos += p;
            }
        }
    }
    bkt->num_slices = slice + 1;
}

/* Apply bucket sieve hits to a single block */
static void bucket_sieve_block(unsigned char *sieve, bucket_t *bkt, int bnum) {
    int total_blocks = bkt->num_blocks;
    for (int s = 0; s < bkt->num_slices; s++) {
        int idx = s * total_blocks + bnum;
        uint32_t n = bkt->count[idx];
        uint32_t *entries = &bkt->data[(size_t)idx * BUCKET_ALLOC];
        uint8_t logp = bkt->slice_logp[s];

        /* Unrolled loop for performance */
        uint32_t i = 0;
        for (; i + 8 <= n; i += 8) {
            sieve[entries[i]   & 0xFFFF] += logp;
            sieve[entries[i+1] & 0xFFFF] += logp;
            sieve[entries[i+2] & 0xFFFF] += logp;
            sieve[entries[i+3] & 0xFFFF] += logp;
            sieve[entries[i+4] & 0xFFFF] += logp;
            sieve[entries[i+5] & 0xFFFF] += logp;
            sieve[entries[i+6] & 0xFFFF] += logp;
            sieve[entries[i+7] & 0xFFFF] += logp;
        }
        for (; i < n; i++) {
            sieve[entries[i] & 0xFFFF] += logp;
        }
    }
}

/* ==================== Large Prime Hash ==================== */
#define LP_HASH_BITS 21
=======
        if (nm == 0) { fb->prime[fb->size] = i; fb->root[fb->size] = 0; fb->logp[fb->size] = (unsigned char)(log2(i)+0.5); fb->size++; continue; }
        unsigned int r = sqrt_mod((unsigned int)nm, i);
        if (!r) continue;
        fb->prime[fb->size] = i; fb->root[fb->size] = r; fb->logp[fb->size] = (unsigned char)(log2(i)+0.5); fb->size++;
    }
    free(sv);
    return fb;
}

/* ==================== Large Prime Hash ==================== */
#define LP_HASH_BITS 20
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
#define LP_HASH_SIZE (1 << LP_HASH_BITS)
typedef struct lp_e { unsigned long lp; int idx; struct lp_e *next; } lp_e_t;
typedef struct { lp_e_t **b; lp_e_t *pool; int used, max; } lp_t;
static lp_t *lp_create(int m) { lp_t *t = calloc(1, sizeof(lp_t)); t->b = calloc(LP_HASH_SIZE, sizeof(lp_e_t*)); t->pool = calloc(m, sizeof(lp_e_t)); t->max = m; return t; }
static int lp_find(lp_t *t, unsigned long lp) { unsigned int h = (unsigned int)((lp * 0x9E3779B97F4A7C15ULL) >> (64-LP_HASH_BITS)); for (lp_e_t *e = t->b[h]; e; e = e->next) if (e->lp == lp) return e->idx; return -1; }
static void lp_insert(lp_t *t, unsigned long lp, int idx) { if (t->used >= t->max) return; unsigned int h = (unsigned int)((lp * 0x9E3779B97F4A7C15ULL) >> (64-LP_HASH_BITS)); lp_e_t *e = &t->pool[t->used++]; e->lp = lp; e->idx = idx; e->next = t->b[h]; t->b[h] = e; }

/* ==================== Relation Storage ==================== */
<<<<<<< HEAD
typedef struct {
    mpz_t *ax_b, *Qx;
    unsigned long *lp;
    int count, alloc;
} rels_t;

static rels_t *rels_create(int n) {
    rels_t *r = malloc(sizeof(rels_t));
    r->ax_b = malloc(n*sizeof(mpz_t)); r->Qx = malloc(n*sizeof(mpz_t));
    r->lp = calloc(n, sizeof(unsigned long));
    for (int i = 0; i < n; i++) { mpz_init(r->ax_b[i]); mpz_init(r->Qx[i]); }
    r->count = 0; r->alloc = n;
    return r;
}

/* ==================== GF(2) Matrix / Gauss ==================== */
=======
typedef struct { mpz_t *ax_b, *Qx; unsigned long *lp; int count, alloc; } rels_t;
static rels_t *rels_create(int n) { rels_t *r = malloc(sizeof(rels_t)); r->ax_b = malloc(n*sizeof(mpz_t)); r->Qx = malloc(n*sizeof(mpz_t)); r->lp = calloc(n, sizeof(unsigned long)); for (int i = 0; i < n; i++) { mpz_init(r->ax_b[i]); mpz_init(r->Qx[i]); } r->count = 0; r->alloc = n; return r; }

/* ==================== GF(2) Matrix ==================== */
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
typedef unsigned long long u64;
typedef struct { u64 **rows; int nr, nc, fbw, idw, wprow; } gf2_t;
static gf2_t *gf2_create(int nr, int nc) {
    gf2_t *m = malloc(sizeof(gf2_t)); m->nr = nr; m->nc = nc;
    m->fbw = (nc+63)/64; m->idw = (nr+63)/64; m->wprow = m->fbw + m->idw;
    m->rows = malloc(nr * sizeof(u64*));
    for (int i = 0; i < nr; i++) { m->rows[i] = calloc(m->wprow, sizeof(u64)); m->rows[i][m->fbw + i/64] |= (1ULL << (i%64)); }
    return m;
}
static void gf2_set(gf2_t *m, int r, int c) { m->rows[r][c/64] |= (1ULL << (c%64)); }
static int gf2_solve(gf2_t *m, int ***deps, int **dlen, int max) {
    int piv = 0;
    for (int c = 0; c < m->nc && piv < m->nr; c++) {
<<<<<<< HEAD
        int pr = -1;
        for (int r = piv; r < m->nr; r++) if ((m->rows[r][c/64] >> (c%64)) & 1) { pr = r; break; }
        if (pr < 0) continue;
        if (pr != piv) { u64 *t = m->rows[pr]; m->rows[pr] = m->rows[piv]; m->rows[piv] = t; }
        for (int r = 0; r < m->nr; r++) {
            if (r == piv) continue;
            if ((m->rows[r][c/64] >> (c%64)) & 1)
                for (int w = 0; w < m->wprow; w++) m->rows[r][w] ^= m->rows[piv][w];
        }
=======
        int pr = -1; for (int r = piv; r < m->nr; r++) if ((m->rows[r][c/64] >> (c%64)) & 1) { pr = r; break; }
        if (pr < 0) continue;
        if (pr != piv) { u64 *t = m->rows[pr]; m->rows[pr] = m->rows[piv]; m->rows[piv] = t; }
        for (int r = 0; r < m->nr; r++) { if (r == piv) continue; if ((m->rows[r][c/64] >> (c%64)) & 1) for (int w = 0; w < m->wprow; w++) m->rows[r][w] ^= m->rows[piv][w]; }
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
        piv++;
    }
    int nd = 0; *deps = malloc(max * sizeof(int*)); *dlen = malloc(max * sizeof(int));
    for (int r = piv; r < m->nr && nd < max; r++) {
<<<<<<< HEAD
        int z = 1;
        for (int w = 0; w < m->fbw && z; w++) {
            u64 mask = (w < m->fbw-1) ? ~0ULL : (m->nc%64==0 ? ~0ULL : (1ULL << (m->nc%64))-1);
            if (m->rows[r][w] & mask) z = 0;
        }
        if (!z) continue;
        int *d = malloc(m->nr * sizeof(int)); int dl = 0;
        for (int w = 0; w < m->idw; w++) {
            u64 bits = m->rows[r][m->fbw+w];
            while (bits) { int bit = __builtin_ctzll(bits); int idx2 = w*64+bit;
                if (idx2 < m->nr) d[dl++] = idx2; bits &= bits-1; }
        }
=======
        int z = 1; for (int w = 0; w < m->fbw && z; w++) { u64 mask = (w < m->fbw-1) ? ~0ULL : (m->nc%64==0 ? ~0ULL : (1ULL << (m->nc%64))-1); if (m->rows[r][w] & mask) z = 0; }
        if (!z) continue;
        int *d = malloc(m->nr * sizeof(int)); int dl = 0;
        for (int w = 0; w < m->idw; w++) { u64 bits = m->rows[r][m->fbw+w]; while (bits) { int bit = __builtin_ctzll(bits); int idx = w*64+bit; if (idx < m->nr) d[dl++] = idx; bits &= bits-1; } }
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
        if (dl > 0) { (*deps)[nd] = d; (*dlen)[nd] = dl; nd++; } else free(d);
    }
    return nd;
}

<<<<<<< HEAD
/* ==================== Parameters ==================== */
typedef struct { int fb_size, nblocks, lp_mult, extra; double thresh; } params_t;
static params_t get_params(int bits) {
    if (bits <= 100) return (params_t){120, 1, 40, 40, 0.72};
    if (bits <= 110) return (params_t){180, 2, 40, 45, 0.73};
    if (bits <= 120) return (params_t){250, 3, 45, 50, 0.75};
    if (bits <= 130) return (params_t){350, 4, 50, 55, 0.77};
    if (bits <= 140) return (params_t){500, 5, 55, 60, 0.78};
    if (bits <= 150) return (params_t){700, 7, 60, 70, 0.79};
    if (bits <= 160) return (params_t){1000, 10, 65, 80, 0.80};
    if (bits <= 170) return (params_t){1400, 14, 70, 90, 0.81};
    if (bits <= 180) return (params_t){2000, 18, 75, 100, 0.82};
    if (bits <= 190) return (params_t){2800, 24, 80, 110, 0.83};
    if (bits <= 200) return (params_t){4000, 30, 85, 120, 0.84};
    if (bits <= 210) return (params_t){5500, 38, 90, 140, 0.85};
    if (bits <= 220) return (params_t){7500, 46, 95, 160, 0.86};
    if (bits <= 230) return (params_t){10000, 54, 100, 180, 0.87};
    if (bits <= 240) return (params_t){14000, 64, 105, 200, 0.875};
    if (bits <= 250) return (params_t){18000, 76, 110, 220, 0.88};
    if (bits <= 260) return (params_t){24000, 90, 115, 250, 0.885};
    if (bits <= 270) return (params_t){32000, 104, 120, 280, 0.89};
    if (bits <= 280) return (params_t){42000, 120, 125, 320, 0.895};
    if (bits <= 290) return (params_t){56000, 140, 130, 360, 0.90};
    if (bits <= 300) return (params_t){72000, 160, 135, 400, 0.905};
    return (params_t){80000, 180, 140, 450, 0.91};
}

/* ==================== Main SIQS ==================== */
=======
/* ==================== Bucket Sieve Structures ==================== */
typedef struct {
    unsigned short offset;
    unsigned char logp;
    unsigned char pad;
} bucket_hit_t;

typedef struct {
    bucket_hit_t *entries;
    int count, alloc;
} bucket_t;

/* ==================== DLP: Pollard Rho for cofactor splitting ==================== */
static unsigned long pollard_rho_64(unsigned long n) {
    if (n % 2 == 0) return 2;
    if (n < 4) return n;
    mpz_t N, x, y, d, c, tmp;
    mpz_init_set_ui(N, n); mpz_init(x); mpz_init(y); mpz_init(d); mpz_init(c); mpz_init(tmp);
    unsigned long result = n;
    for (int att = 0; att < 30 && result == n; att++) {
        mpz_set_ui(c, att + 1);
        mpz_set_ui(x, 2 + att);
        mpz_set(y, x);
        for (int iter = 0; iter < 100000; iter++) {
            mpz_mul(x, x, x); mpz_add(x, x, c); mpz_mod(x, x, N);
            mpz_mul(y, y, y); mpz_add(y, y, c); mpz_mod(y, y, N);
            mpz_mul(y, y, y); mpz_add(y, y, c); mpz_mod(y, y, N);
            mpz_sub(tmp, x, y); mpz_abs(tmp, tmp);
            mpz_gcd(d, tmp, N);
            if (mpz_cmp_ui(d, 1) > 0) {
                if (mpz_cmp(d, N) < 0) { result = mpz_get_ui(d); break; }
                else break;  /* d == N, try different c */
            }
        }
    }
    mpz_clear(N); mpz_clear(x); mpz_clear(y); mpz_clear(d); mpz_clear(c); mpz_clear(tmp);
    return result;
}

/* ==================== Parameters ==================== */
/* Calibrated against YAFU's proven parameter table - digit based */
typedef struct { int fb_size, nblocks, lp_mult, extra; double thresh; int use_dlp; int dlp_mult; } params_t;
static params_t get_params_by_digits(int digits) {
    /*                     FB    nblk LP   extra thresh DLP  DLP_m */
    if (digits <= 30) return (params_t){100,  1,  30,  40, 0.73, 0, 0};
    if (digits <= 33) return (params_t){150,  1,  30,  40, 0.74, 0, 0};
    if (digits <= 36) return (params_t){250,  1,  35,  50, 0.76, 0, 0};
    if (digits <= 39) return (params_t){450,  1,  40,  50, 0.78, 0, 0};
    if (digits <= 42) return (params_t){550,  1,  40,  60, 0.79, 0, 0};
    if (digits <= 45) return (params_t){700,  2,  40,  60, 0.80, 0, 0};
    if (digits <= 48) return (params_t){1000, 2,  40,  60, 0.80, 1, 80};
    if (digits <= 51) return (params_t){1400, 2,  40,  80, 0.81, 1, 80};
    if (digits <= 55) return (params_t){2000, 2,  50,  80, 0.82, 1, 100};
    if (digits <= 60) return (params_t){3000, 2,  50, 100, 0.83, 1, 100};
    if (digits <= 64) return (params_t){5400, 3,  80, 120, 0.84, 1, 150};
    if (digits <= 70) return (params_t){10000,3, 100, 150, 0.86, 1, 200};
    if (digits <= 75) return (params_t){27000,4, 100, 200, 0.875, 1, 250};
    if (digits <= 80) return (params_t){50000,4, 100, 250, 0.885, 1, 250};
    if (digits <= 85) return (params_t){55000,3,  80, 300, 0.89, 1, 300};
    if (digits <= 90) return (params_t){60000,9,  80, 350, 0.90, 1, 300};
    return (params_t){70000, 12, 80, 400, 0.91, 1, 300};
}

/* ==================== Main ==================== */
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }
    clock_gettime(CLOCK_MONOTONIC, &g_start);

    mpz_t N, kN;
    mpz_inits(N, kN, NULL);
    mpz_set_str(N, argv[1], 10);

    int digits = (int)mpz_sizeinbase(N, 10);
    int bits = (int)mpz_sizeinbase(N, 2);

    /* Trial division */
    for (int p = 2; p < 10000; p++) {
        if (mpz_divisible_ui_p(N, p)) {
            mpz_t c; mpz_init(c); mpz_divexact_ui(c, N, p);
<<<<<<< HEAD
            gmp_printf("%Zd\n", c);
            fprintf(stderr, "trivial factor %d in %.3fs\n", p, elapsed());
            return 0;
=======
            gmp_printf("%d\n", p); return 0;
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
        }
    }

    int mult = choose_multiplier(N);
    mpz_mul_ui(kN, N, mult);
    int kN_bits = (int)mpz_sizeinbase(kN, 2);
<<<<<<< HEAD
    params_t P = get_params(kN_bits);

    fb_t *fb = fb_create(kN, P.fb_size);
    int M = BLOCK_SIZE * P.nblocks;      /* half sieve interval */
    int total_blocks = 2 * P.nblocks;     /* blocks on both sides */
    unsigned long lp_bound = (unsigned long)fb->prime[fb->size-1] * P.lp_mult;
=======
    params_t P = get_params_by_digits(digits);

    fb_t *fb = fb_create(kN, P.fb_size);
    int M = SIEVE_BLOCK * P.nblocks;
    unsigned long lp_bound = (unsigned long)fb->prime[fb->size-1] * P.lp_mult;
    unsigned long dlp_bound = P.use_dlp ? (unsigned long)fb->prime[fb->size-1] * P.dlp_mult : 0;
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
    int target = fb->size + P.extra;

    double log2_Qmax = kN_bits / 2.0 + 0.5 + log2(M);
    int threshold = (int)(log2_Qmax * P.thresh);
<<<<<<< HEAD
    threshold -= 2;  /* slightly more aggressive to catch more candidates */

    fprintf(stderr, "SIQS-Bucket: %dd (%db), k=%d, FB=%d (med=%d,large=%d), M=%d, blocks=%d, thresh=%d, LP=%lu, target=%d\n",
            digits, bits, mult, fb->size, fb->med_start, fb->large_start,
            M, total_blocks, threshold, lp_bound, target);

    /* Sieve array */
    unsigned char *sieve = aligned_alloc(64, BLOCK_SIZE);

    /* Bucket sieve for large primes */
    int max_slices = 64;
    bucket_t *buckets = bucket_create(total_blocks, max_slices);

    /* Solution arrays */
    unsigned int *soln1 = malloc(fb->size * sizeof(unsigned int));
    unsigned int *soln2 = malloc(fb->size * sizeof(unsigned int));

    /* Relations */
=======
    threshold -= 3;
    if (P.use_dlp) threshold -= 6;  /* Lower threshold to catch more DLP candidates */

    fprintf(stderr, "SIQS-Bucket: %dd (%db), k=%d, FB=%d, M=%d, thresh=%d, LP=%lu",
            digits, bits, mult, fb->size, M, threshold, lp_bound);
    if (P.use_dlp) fprintf(stderr, ", DLP=%lu", dlp_bound);
    fprintf(stderr, ", target=%d\n", target);

    /* Determine bucket sieve threshold */
    int small_fb_end = fb->size;
    for (int i = 0; i < fb->size; i++) {
        if (fb->prime[i] > (unsigned int)SIEVE_BLOCK) { small_fb_end = i; break; }
    }
    int large_fb_start = small_fb_end;
    fprintf(stderr, "FB split: small=%d (primes<%d), large=%d (primes>%d)\n",
            small_fb_end, SIEVE_BLOCK, fb->size - large_fb_start, SIEVE_BLOCK);

    /* Allocate sieve arrays */
    unsigned char *sieves[BATCH_POLYS];
    for (int b = 0; b < BATCH_POLYS; b++) sieves[b] = malloc(SIEVE_BLOCK);

    /* Allocate bucket arrays: per-block, per-polynomial */
    int n_blocks = 2 * P.nblocks;  /* -nblocks to +nblocks */
    int n_bucket_arrays = n_blocks * BATCH_POLYS;
    bucket_t *buckets = calloc(n_bucket_arrays, sizeof(bucket_t));
    for (int i = 0; i < n_bucket_arrays; i++) {
        buckets[i].alloc = 256;
        buckets[i].entries = malloc(256 * sizeof(bucket_hit_t));
    }

>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
    rels_t *full = rels_create(MAX_RELS);
    rels_t *part = rels_create(MAX_PARTIALS);
    lp_t *lpt = lp_create(MAX_PARTIALS);

<<<<<<< HEAD
    /* Polynomial state */
    mpz_t a, b_val, c_val, B_vals[MAX_A_FACTORS];
    mpz_inits(a, b_val, c_val, NULL);
    for (int j = 0; j < MAX_A_FACTORS; j++) mpz_init(B_vals[j]);

    /* ainv[j][i] = 2*a^{-1}*B_j mod p_i, for Gray code root updates */
    unsigned int **ainv = malloc(MAX_A_FACTORS * sizeof(unsigned int *));
    for (int j = 0; j < MAX_A_FACTORS; j++)
        ainv[j] = malloc(fb->size * sizeof(unsigned int));
=======
    /* DLP storage */
    rels_t *dlp_part = P.use_dlp ? rels_create(MAX_PARTIALS) : NULL;
    lp_t *dlp_lpt = P.use_dlp ? lp_create(MAX_PARTIALS) : NULL;

    mpz_t a, bs[BATCH_POLYS], cs[BATCH_POLYS], B_vals[MAX_A_FACTORS];
    mpz_init(a);
    for (int b = 0; b < BATCH_POLYS; b++) { mpz_init(bs[b]); mpz_init(cs[b]); }
    for (int j = 0; j < MAX_A_FACTORS; j++) mpz_init(B_vals[j]);

    unsigned int *soln1[BATCH_POLYS], *soln2[BATCH_POLYS];
    for (int b = 0; b < BATCH_POLYS; b++) {
        soln1[b] = malloc(fb->size * sizeof(unsigned int));
        soln2[b] = malloc(fb->size * sizeof(unsigned int));
    }
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)

    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, SEED);

    mpz_t ax_b, Qx, residue, tmp;
    mpz_inits(ax_b, Qx, residue, tmp, NULL);

<<<<<<< HEAD
    int total_polys = 0, a_count = 0, combined = 0;
    int a_idx[MAX_A_FACTORS];
    int num_a_factors = 0;

    /* ===== Main sieving loop ===== */
    while (full->count < target) {
        if (elapsed() > 280) {
            fprintf(stderr, "TIMEOUT approaching, stopping sieve at %.1fs\n", elapsed());
            break;
        }

        /* ===== Generate new 'a' value ===== */
=======
    int total_polys = 0, a_count = 0, combined = 0, dlp_combined = 0;
    int a_idx[MAX_A_FACTORS];
    int num_a_factors = 0;
    unsigned int *ainv_data = malloc(MAX_A_FACTORS * fb->size * sizeof(unsigned int));

    while (full->count < target) {
        if (total_polys > 0 && total_polys % (2000/BATCH_POLYS) == 0) {
            double t = elapsed();
            if (t > 280) { fprintf(stderr, "TIMEOUT at %.1fs\n", t); break; }
            if (total_polys % (4000/BATCH_POLYS) == 0)
                fprintf(stderr, "  p=%d r=%d/%d (f=%d+slp=%d+dlp=%d) part=%d dlp_part=%d t=%.1fs\n",
                        total_polys * BATCH_POLYS, full->count, target,
                        full->count - combined - dlp_combined, combined, dlp_combined,
                        part->count, dlp_part ? dlp_part->count : 0, t);
        }

        /* Generate new 'a' */
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
        {
            mpz_t tgt; mpz_init(tgt);
            mpz_mul_ui(tgt, kN, 2); mpz_sqrt(tgt, tgt); mpz_tdiv_q_ui(tgt, tgt, M);
            double log_tgt = mpz_sizeinbase(tgt, 2) * log(2.0);

            int lo = fb->size / 3, hi = 2 * fb->size / 3;
            if (lo < 2) lo = 2; if (hi <= lo + 3) hi = fb->size - 1;

            double avg = 0; int cnt = 0;
<<<<<<< HEAD
            for (int i = lo; i < hi; i++) { if (fb->root1[i] == 0) continue; avg += log(fb->prime[i]); cnt++; }
=======
            for (int i = lo; i < hi; i++) { if (fb->root[i] == 0) continue; avg += log(fb->prime[i]); cnt++; }
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
            if (cnt == 0) { mpz_clear(tgt); break; }
            avg /= cnt;

            int s = (int)(log_tgt / avg + 0.5);
            if (s < 3) s = 3; if (s > MAX_A_FACTORS) s = MAX_A_FACTORS; if (s > hi - lo) s = hi - lo;
            num_a_factors = s;

            double best_ratio = 1e30;
            int best[MAX_A_FACTORS];

<<<<<<< HEAD
            for (int att = 0; att < 50; att++) {
                mpz_set_ui(a, 1);
                int idx[MAX_A_FACTORS]; int ok = 1;
                for (int i2 = 0; i2 < s && ok; i2++) {
                    int tries = 0, good;
                    do { idx[i2] = lo + gmp_urandomm_ui(rng, hi-lo); good = 1;
                         for (int j = 0; j < i2; j++) if (idx[j]==idx[i2]) {good=0; break;}
                         if (fb->root1[idx[i2]]==0) good=0; tries++;
                    } while (!good && tries < 100);
                    if (!good) { ok=0; break; }
                    mpz_mul_ui(a, a, fb->prime[idx[i2]]);
                }
                if (!ok) continue;
                double ratio;
                if (mpz_cmp(a, tgt) > 0) { mpz_t q2; mpz_init(q2); mpz_tdiv_q(q2, a, tgt); ratio = mpz_get_d(q2); mpz_clear(q2); }
                else { mpz_t q2; mpz_init(q2); mpz_tdiv_q(q2, tgt, a); ratio = mpz_get_d(q2); mpz_clear(q2); }
                if (ratio < best_ratio) { best_ratio = ratio; memcpy(best, idx, s*sizeof(int)); }
                if (ratio < 1.5) break;
=======
            for (int att = 0; att < 40; att++) {
                mpz_set_ui(a, 1);
                int idx[MAX_A_FACTORS]; int ok = 1;
                for (int i = 0; i < s && ok; i++) {
                    int tries = 0, good;
                    do { idx[i] = lo + gmp_urandomm_ui(rng, hi-lo); good = 1;
                         for (int j = 0; j < i; j++) if (idx[j]==idx[i]) {good=0; break;}
                         if (fb->root[idx[i]]==0) good=0; tries++;
                    } while (!good && tries < 100);
                    if (!good) { ok=0; break; }
                    mpz_mul_ui(a, a, fb->prime[idx[i]]);
                }
                if (!ok) continue;
                double ratio;
                if (mpz_cmp(a, tgt) > 0) { mpz_t q; mpz_init(q); mpz_tdiv_q(q, a, tgt); ratio = mpz_get_d(q); mpz_clear(q); }
                else { mpz_t q; mpz_init(q); mpz_tdiv_q(q, tgt, a); ratio = mpz_get_d(q); mpz_clear(q); }
                if (ratio < best_ratio) { best_ratio = ratio; memcpy(best, idx, s*sizeof(int)); }
                if (ratio < 2.0) break;
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
            }

            memcpy(a_idx, best, s * sizeof(int));
            mpz_set_ui(a, 1);
<<<<<<< HEAD
            for (int i2 = 0; i2 < s; i2++) mpz_mul_ui(a, a, fb->prime[a_idx[i2]]);
            mpz_clear(tgt);
            a_count++;

            /* Compute B values: B_j = a/q_j * (sqrt(kN) mod q_j) * (a/q_j)^{-1} mod q_j */
            for (int j = 0; j < s; j++) {
                int idx2 = a_idx[j];
                unsigned int qj = fb->prime[idx2], rj = fb->root1[idx2];
                mpz_t a_q, mod_q, inv2;
                mpz_inits(a_q, mod_q, inv2, NULL);
                mpz_divexact_ui(a_q, a, qj); mpz_set_ui(mod_q, qj);
                mpz_invert(inv2, a_q, mod_q);
                unsigned long iv = mpz_get_ui(inv2);
                mpz_mul_ui(B_vals[j], a_q, (unsigned long)rj * iv % qj);
                mpz_clears(a_q, mod_q, inv2, NULL);
            }

            /* Precompute ainv[j][i] = 2 * a^{-1} * B_j mod p_i for each FB prime */
            for (int j = 0; j < s; j++) {
                for (int i2 = 0; i2 < fb->size; i2++) {
                    unsigned int p = fb->prime[i2];
                    unsigned long am = mpz_fdiv_ui(a, p);
                    if (am == 0 || fb->root1[i2] == 0) { ainv[j][i2] = 0; continue; }
                    unsigned int ai = mod_inverse((unsigned int)am, p);
                    unsigned long Bm = mpz_fdiv_ui(B_vals[j], p);
                    ainv[j][i2] = (unsigned int)((2ULL * ai % p * Bm) % p);
=======
            for (int i = 0; i < s; i++) mpz_mul_ui(a, a, fb->prime[a_idx[i]]);
            mpz_clear(tgt);
            a_count++;

            /* Compute B values */
            for (int j = 0; j < s; j++) {
                int idx = a_idx[j];
                unsigned int qj = fb->prime[idx], rj = fb->root[idx];
                mpz_t a_q, mod_q, inv; mpz_inits(a_q, mod_q, inv, NULL);
                mpz_divexact_ui(a_q, a, qj); mpz_set_ui(mod_q, qj);
                mpz_invert(inv, a_q, mod_q);
                unsigned long iv = mpz_get_ui(inv);
                mpz_mul_ui(B_vals[j], a_q, (rj * iv) % qj);
                mpz_clears(a_q, mod_q, inv, NULL);
            }

            /* Precompute ainv for Gray code */
            for (int j = 0; j < s; j++) {
                for (int i = 0; i < fb->size; i++) {
                    unsigned int p = fb->prime[i];
                    unsigned long am = mpz_fdiv_ui(a, p);
                    if (am == 0 || fb->root[i] == 0) { ainv_data[j*fb->size+i] = 0; continue; }
                    unsigned int ai = mod_inverse((unsigned int)am, p);
                    unsigned long Bm = mpz_fdiv_ui(B_vals[j], p);
                    ainv_data[j*fb->size+i] = (unsigned int)((2UL * ai % p * Bm) % p);
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
                }
            }
        }

<<<<<<< HEAD
        /* ===== Iterate through b-values using Gray code ===== */
        int num_b = 1 << (num_a_factors - 1);

        /* Compute initial b (Gray code index 0) */
        mpz_set_ui(b_val, 0);
        for (int j = 0; j < num_a_factors; j++)
            mpz_add(b_val, b_val, B_vals[j]);

        /* c = (b^2 - kN) / a */
        mpz_mul(c_val, b_val, b_val); mpz_sub(c_val, c_val, kN);
        if (!mpz_divisible_p(c_val, a)) {
            mpz_neg(b_val, b_val);
            mpz_mul(c_val, b_val, b_val); mpz_sub(c_val, c_val, kN);
        }
        mpz_divexact(c_val, c_val, a);

        /* Compute initial sieve solutions */
        for (int i = 0; i < fb->size; i++) {
            unsigned int p = fb->prime[i];
            unsigned long am = mpz_fdiv_ui(a, p);
            if (am == 0 || fb->root1[i] == 0) { soln1[i] = soln2[i] = 0xFFFFFFFF; continue; }
            unsigned int ai = mod_inverse((unsigned int)am, p);
            if (ai == 0) { soln1[i] = soln2[i] = 0xFFFFFFFF; continue; }
            unsigned long bm = mpz_fdiv_ui(b_val, p);
            unsigned int r = fb->root1[i];
            soln1[i] = (unsigned int)((unsigned long)ai * ((r + p - bm) % p) % p);
            soln2[i] = (unsigned int)((unsigned long)ai * ((p - r + p - bm) % p) % p);
        }

        for (int b_idx = 0; b_idx < num_b && full->count < target; b_idx++) {
            total_polys++;

            if (total_polys % 500 == 0) {
                double t = elapsed();
                if (t > 280) break;
                if (total_polys % 2000 == 0)
                    fprintf(stderr, "  polys=%d rels=%d/%d (full=%d+%d) part=%d t=%.1fs\n",
                            total_polys, full->count, target,
                            full->count - combined, combined, part->count, t);
            }

            /* Gray code update for b_idx > 0 */
            if (b_idx > 0) {
                /* Find which B_j to flip: position of lowest set bit in b_idx */
                int j = __builtin_ctz(b_idx);
                int sign = ((b_idx >> j) & 2) ? 1 : -1;
                /* Note: Gray code changes sign of B_j */
                /* Update b: b_new = b_old ± 2*B_j */
                if (sign > 0) {
                    mpz_addmul_ui(b_val, B_vals[j], 2);
                } else {
                    mpz_submul_ui(b_val, B_vals[j], 2);
                }
                /* Update c = (b^2 - kN) / a */
                mpz_mul(c_val, b_val, b_val); mpz_sub(c_val, c_val, kN);
                mpz_divexact(c_val, c_val, a);

                /* Update sieve solutions: soln_new = soln_old ± ainv[j][i] */
                for (int i = 0; i < fb->size; i++) {
                    if (soln1[i] == 0xFFFFFFFF) continue;
                    unsigned int p = fb->prime[i];
                    unsigned int delta = ainv[j][i];
                    if (sign > 0) {
                        soln1[i] = (soln1[i] >= delta) ? soln1[i] - delta : soln1[i] + p - delta;
                        soln2[i] = (soln2[i] >= delta) ? soln2[i] - delta : soln2[i] + p - delta;
                    } else {
                        soln1[i] = soln1[i] + delta;
                        if (soln1[i] >= p) soln1[i] -= p;
                        soln2[i] = soln2[i] + delta;
                        if (soln2[i] >= p) soln2[i] -= p;
                    }
                }
            }

            /* Fill buckets for large primes */
            int sieve_start = -M;
            bucket_fill(buckets, fb, soln1, soln2, total_blocks, sieve_start);

            /* ===== Sieve each block ===== */
            for (int bnum = 0; bnum < total_blocks; bnum++) {
                int block_start = sieve_start + bnum * BLOCK_SIZE;

                /* Initialize sieve */
                memset(sieve, 0, BLOCK_SIZE);

                /* Sieve with small/medium primes (direct sieve) */
                for (int i = 1; i < fb->large_start; i++) {
                    unsigned int p = fb->prime[i];
                    if (p < 5) continue;
                    if (soln1[i] == 0xFFFFFFFF) continue;
                    unsigned char lp = fb->logp[i];

                    /* Compute starting position in this block for root1 */
                    long off1 = ((long)soln1[i] - block_start) % (long)p;
                    if (off1 < 0) off1 += p;
                    /* Sieve root1 */
                    for (unsigned int j = (unsigned int)off1; j < BLOCK_SIZE; j += p)
                        sieve[j] += lp;

                    /* Sieve root2 (if different) */
                    if (soln1[i] != soln2[i]) {
                        long off2 = ((long)soln2[i] - block_start) % (long)p;
                        if (off2 < 0) off2 += p;
                        for (unsigned int j = (unsigned int)off2; j < BLOCK_SIZE; j += p)
                            sieve[j] += lp;
                    }
                }

                /* Apply bucket sieve for large primes */
                bucket_sieve_block(sieve, buckets, bnum);

                /* ===== Scan for smooth candidates ===== */
                for (int j = 0; j < BLOCK_SIZE; j++) {
                    if (sieve[j] < threshold) continue;
                    long x = (long)(block_start + j);
                    if (x == 0) continue;

                    /* Compute Q(x) = a*x^2 + 2*b*x + c */
                    mpz_set_si(tmp, x);
                    mpz_mul(Qx, a, tmp);
                    mpz_add(Qx, Qx, b_val); mpz_add(Qx, Qx, b_val);
                    mpz_mul(Qx, Qx, tmp);
                    mpz_add(Qx, Qx, c_val);

                    /* ax + b */
                    mpz_mul_si(ax_b, a, x);
                    mpz_add(ax_b, ax_b, b_val);

                    if (mpz_sgn(Qx) == 0) continue;
                    mpz_abs(residue, Qx);

                    /* Trial divide by factor base - root-aware */
                    while (mpz_even_p(residue)) mpz_tdiv_q_2exp(residue, residue, 1);

                    for (int i = 1; i < fb->size; i++) {
                        unsigned int p = fb->prime[i];
                        if (soln1[i] == 0xFFFFFFFF) continue;
                        /* Check if this prime divides Q(x) via sieve root match */
                        long xmod = ((x % (long)p) + p) % p;
                        if (xmod != (long)soln1[i] && xmod != (long)soln2[i]) continue;
                        if (mpz_divisible_ui_p(residue, p))
                            do { mpz_divexact_ui(residue, residue, p); } while (mpz_divisible_ui_p(residue, p));
                    }
                    /* Also check small primes that may have been skipped */
                    for (int i = 0; i < fb->size && fb->prime[i] < 5; i++) {
                        unsigned int p = fb->prime[i]; if (p <= 2) continue;
                        while (mpz_divisible_ui_p(residue, p)) mpz_divexact_ui(residue, residue, p);
                    }

                    /* Store relation */
                    mpz_t aQx; mpz_init(aQx); mpz_mul(aQx, Qx, a);

                    if (mpz_cmp_ui(residue, 1) == 0) {
                        int ri = full->count;
                        if (ri < full->alloc) {
                            mpz_set(full->ax_b[ri], ax_b);
                            mpz_set(full->Qx[ri], aQx);
                            full->lp[ri] = 0;
                            full->count++;
                        }
                    } else if (mpz_fits_ulong_p(residue) && mpz_get_ui(residue) <= lp_bound) {
                        unsigned long lp_val = mpz_get_ui(residue);
                        int match = lp_find(lpt, lp_val);
                        if (match >= 0) {
                            int ri = full->count;
                            if (ri < full->alloc) {
                                mpz_mul(full->ax_b[ri], ax_b, part->ax_b[match]);
                                mpz_mod(full->ax_b[ri], full->ax_b[ri], N);
                                mpz_mul(full->Qx[ri], aQx, part->Qx[match]);
                                full->lp[ri] = lp_val;
                                full->count++;
                                combined++;
                            }
                        } else {
                            int pi = part->count;
                            if (pi < part->alloc) {
                                mpz_set(part->ax_b[pi], ax_b);
                                mpz_set(part->Qx[pi], aQx);
                                part->lp[pi] = lp_val;
                                lp_insert(lpt, lp_val, pi);
                                part->count++;
                            }
                        }
                    }
                    mpz_clear(aQx);
                }
            }
        }
    }

    double sieve_time = elapsed();
    fprintf(stderr, "Sieving done: %d rels (%d full + %d combined) in %.2fs, %d polys, %d a-values\n",
            full->count, full->count - combined, combined, sieve_time, total_polys, a_count);

    if (full->count < fb->size + 1) {
        fprintf(stderr, "FAIL: not enough relations (%d < %d)\n", full->count, fb->size + 1);
        printf("FAIL\n"); return 1;
    }

    /* ===== Linear Algebra (Gaussian Elimination over GF(2)) ===== */
=======
        int num_b = 1 << (num_a_factors - 1);

        /* Persistent state for Gray code self-initialization */
        mpz_t cur_b; mpz_init(cur_b);
        unsigned int *cur_s1 = malloc(fb->size * sizeof(unsigned int));
        unsigned int *cur_s2 = malloc(fb->size * sizeof(unsigned int));

        /* Compute base solutions for b_idx=0: Gray(0)=0, all Bj negative */
        {
            mpz_set_ui(cur_b, 0);
            for (int j = 0; j < num_a_factors; j++)
                mpz_sub(cur_b, cur_b, B_vals[j]);
            mpz_mul(tmp, cur_b, cur_b); mpz_sub(tmp, tmp, kN); mpz_mod(tmp, tmp, a);
            if (mpz_sgn(tmp) != 0) mpz_neg(cur_b, cur_b);

            /* Compute base solutions from scratch */
            for (int i = 0; i < fb->size; i++) {
                unsigned int p = fb->prime[i];
                unsigned long am = mpz_fdiv_ui(a, p);
                if (am == 0 || fb->root[i] == 0) { cur_s1[i] = cur_s2[i] = 0xFFFFFFFF; continue; }
                unsigned int ai = mod_inverse((unsigned int)am, p);
                if (ai == 0) { cur_s1[i] = cur_s2[i] = 0xFFFFFFFF; continue; }
                unsigned long bm = mpz_fdiv_ui(cur_b, p);
                unsigned int r = fb->root[i];
                cur_s1[i] = (unsigned int)((unsigned long)ai * ((r + p - bm) % p) % p);
                cur_s2[i] = (unsigned int)((unsigned long)ai * ((p - r + p - bm) % p) % p);
            }
        }

        for (int b_start = 0; b_start < num_b && full->count < target; b_start += BATCH_POLYS) {
            int batch = BATCH_POLYS;
            if (b_start + batch > num_b) batch = num_b - b_start;

            /* Compute b-values and solutions for this batch using Gray code */
            for (int bi = 0; bi < batch; bi++) {
                int b_idx = b_start + bi;

                if (b_idx == 0) {
                    /* Copy base solutions */
                    mpz_set(bs[0], cur_b);
                    mpz_mul(cs[0], cur_b, cur_b); mpz_sub(cs[0], cs[0], kN);
                    mpz_divexact(cs[0], cs[0], a);
                    memcpy(soln1[0], cur_s1, fb->size * sizeof(unsigned int));
                    memcpy(soln2[0], cur_s2, fb->size * sizeof(unsigned int));
                } else {
                    /* Gray code update: find which bit changed */
                    int changed_bit = __builtin_ctz(b_idx);
                    int gray = b_idx ^ (b_idx >> 1);
                    int bit_is_set = (gray >> changed_bit) & 1;

                    /* Update b value */
                    if (bit_is_set)
                        mpz_addmul_ui(cur_b, B_vals[changed_bit], 2);
                    else
                        mpz_submul_ui(cur_b, B_vals[changed_bit], 2);

                    mpz_set(bs[bi], cur_b);
                    mpz_mul(cs[bi], cur_b, cur_b); mpz_sub(cs[bi], cs[bi], kN);
                    mpz_divexact(cs[bi], cs[bi], a);

                    /* Update sieve solutions using precomputed ainv_data */
                    unsigned int *delta = &ainv_data[changed_bit * fb->size];
                    for (int i = 0; i < fb->size; i++) {
                        if (cur_s1[i] == 0xFFFFFFFF) {
                            soln1[bi][i] = soln2[bi][i] = 0xFFFFFFFF;
                            continue;
                        }
                        unsigned int p = fb->prime[i];
                        unsigned int d = delta[i];
                        if (bit_is_set) {
                            cur_s1[i] = cur_s1[i] >= d ? cur_s1[i] - d : cur_s1[i] + p - d;
                            cur_s2[i] = cur_s2[i] >= d ? cur_s2[i] - d : cur_s2[i] + p - d;
                        } else {
                            cur_s1[i] = cur_s1[i] + d;
                            if (cur_s1[i] >= p) cur_s1[i] -= p;
                            cur_s2[i] = cur_s2[i] + d;
                            if (cur_s2[i] >= p) cur_s2[i] -= p;
                        }
                        soln1[bi][i] = cur_s1[i];
                        soln2[bi][i] = cur_s2[i];
                    }
                }
            }

            if (batch == 0) continue;
            total_polys++;

            /* Sieve all BATCH polynomials over each block */
            for (int block = -P.nblocks; block < P.nblocks; block++) {
                int block_start = block * SIEVE_BLOCK;
                int block_idx = block + P.nblocks;  /* 0 to n_blocks-1 */

                /* Initialize sieve arrays */
                for (int bi = 0; bi < batch; bi++)
                    memset(sieves[bi], 0, SIEVE_BLOCK);

                /* Clear buckets for this block */
                for (int bi = 0; bi < batch; bi++) {
                    int bkt_idx = block_idx * BATCH_POLYS + bi;
                    if (bkt_idx < n_bucket_arrays)
                        buckets[bkt_idx].count = 0;
                }

                /* Sieve small FB primes (direct sieve, cache-friendly) */
                for (int i = 1; i < small_fb_end; i++) {
                    unsigned int p = fb->prime[i];
                    if (p < 5) continue;
                    unsigned char lp = fb->logp[i];

                    for (int bi = 0; bi < batch; bi++) {
                        if (soln1[bi][i] == 0xFFFFFFFF) continue;

                        long off1 = ((long)soln1[bi][i] - block_start) % (long)p;
                        if (off1 < 0) off1 += p;
                        for (int j = (int)off1; j < SIEVE_BLOCK; j += p)
                            sieves[bi][j] += lp;

                        if (soln1[bi][i] != soln2[bi][i]) {
                            long off2 = ((long)soln2[bi][i] - block_start) % (long)p;
                            if (off2 < 0) off2 += p;
                            for (int j = (int)off2; j < SIEVE_BLOCK; j += p)
                                sieves[bi][j] += lp;
                        }
                    }
                }

                /* Large FB primes: add directly to sieve (they hit at most a few positions per block) */
                for (int i = large_fb_start; i < fb->size; i++) {
                    unsigned int p = fb->prime[i];
                    unsigned char lp = fb->logp[i];

                    for (int bi = 0; bi < batch; bi++) {
                        if (soln1[bi][i] == 0xFFFFFFFF) continue;

                        long off1 = ((long)soln1[bi][i] - block_start) % (long)p;
                        if (off1 < 0) off1 += p;
                        if (off1 < SIEVE_BLOCK)
                            sieves[bi][(int)off1] += lp;

                        if (soln1[bi][i] != soln2[bi][i]) {
                            long off2 = ((long)soln2[bi][i] - block_start) % (long)p;
                            if (off2 < 0) off2 += p;
                            if (off2 < SIEVE_BLOCK)
                                sieves[bi][(int)off2] += lp;
                        }
                    }
                }

                /* Scan for smooth candidates */
                for (int bi = 0; bi < batch; bi++) {
                    for (int j = 0; j < SIEVE_BLOCK; j++) {
                        if (sieves[bi][j] < threshold) continue;
                        long x = (long)(block_start + j);
                        if (x == 0) continue;

                        /* Compute Q(x) and ax+b */
                        mpz_set_si(tmp, x);
                        mpz_mul(Qx, a, tmp);
                        mpz_add(Qx, Qx, bs[bi]); mpz_add(Qx, Qx, bs[bi]);
                        mpz_mul(Qx, Qx, tmp);
                        mpz_add(Qx, Qx, cs[bi]);

                        mpz_mul_si(ax_b, a, x);
                        mpz_add(ax_b, ax_b, bs[bi]);

                        if (mpz_sgn(Qx) == 0) continue;
                        mpz_abs(residue, Qx);

                        /* Sieve-informed trial division: only check primes whose roots match */
                        while (mpz_even_p(residue)) mpz_tdiv_q_2exp(residue, residue, 1);
                        for (int i = 1; i < fb->size; i++) {
                            unsigned int p = fb->prime[i];
                            if (soln1[bi][i] == 0xFFFFFFFF) continue;
                            /* Check if this prime divides Q(x) by checking root match */
                            long xmod = ((x % (long)p) + p) % p;
                            if (xmod != (long)soln1[bi][i] && xmod != (long)soln2[bi][i]) continue;
                            if (mpz_divisible_ui_p(residue, p))
                                do { mpz_divexact_ui(residue, residue, p); } while (mpz_divisible_ui_p(residue, p));
                        }
                        /* Handle small primes that might have been skipped */
                        for (int i = 0; i < fb->size && fb->prime[i] < 5; i++) {
                            unsigned int p = fb->prime[i]; if (p <= 2) continue;
                            while (mpz_divisible_ui_p(residue, p)) mpz_divexact_ui(residue, residue, p);
                        }

                        /* Store a*Q(x) for proper congruence: (ax+b)^2 = a*Q(x) + kN */
                        mpz_t aQx; mpz_init(aQx); mpz_mul(aQx, Qx, a);

                        if (mpz_cmp_ui(residue, 1) == 0) {
                            /* Fully smooth */
                            int ri = full->count;
                            if (ri < full->alloc) {
                                mpz_set(full->ax_b[ri], ax_b);
                                mpz_set(full->Qx[ri], aQx);
                                full->lp[ri] = 0;
                                full->count++;
                            }
                        } else if (mpz_fits_ulong_p(residue)) {
                            unsigned long cof = mpz_get_ui(residue);
                            if (cof <= lp_bound) {
                                /* Single large prime partial */
                                int match = lp_find(lpt, cof);
                                if (match >= 0) {
                                    int ri = full->count;
                                    if (ri < full->alloc) {
                                        mpz_mul(full->ax_b[ri], ax_b, part->ax_b[match]);
                                        mpz_mod(full->ax_b[ri], full->ax_b[ri], N);
                                        mpz_mul(full->Qx[ri], aQx, part->Qx[match]);
                                        full->lp[ri] = cof;
                                        full->count++;
                                        combined++;
                                    }
                                } else {
                                    int pi = part->count;
                                    if (pi < part->alloc) {
                                        mpz_set(part->ax_b[pi], ax_b);
                                        mpz_set(part->Qx[pi], aQx);
                                        part->lp[pi] = cof;
                                        lp_insert(lpt, cof, pi);
                                        part->count++;
                                    }
                                }
                            } else if (P.use_dlp && cof <= dlp_bound * (unsigned long)fb->prime[fb->size-1]) {
                                /* Try DLP: split cofactor */
                                if (!mpz_probab_prime_p(residue, 1)) {
                                    unsigned long f1 = pollard_rho_64(cof);
                                    if (f1 > 1 && f1 < cof) {
                                        unsigned long f2 = cof / f1;
                                        if (f1 > f2) { unsigned long t = f1; f1 = f2; f2 = t; }
                                        if (f1 <= lp_bound && f2 <= lp_bound) {
                                            /* Valid DLP relation */
                                            unsigned long key = f1 ^ (f2 * 0x9E3779B9UL);
                                            int match = lp_find(dlp_lpt, key);
                                            if (match >= 0) {
                                                int ri = full->count;
                                                if (ri < full->alloc) {
                                                    mpz_mul(full->ax_b[ri], ax_b, dlp_part->ax_b[match]);
                                                    mpz_mod(full->ax_b[ri], full->ax_b[ri], N);
                                                    mpz_mul(full->Qx[ri], aQx, dlp_part->Qx[match]);
                                                    full->lp[ri] = cof;
                                                    full->count++;
                                                    dlp_combined++;
                                                }
                                            } else {
                                                int pi = dlp_part->count;
                                                if (pi < dlp_part->alloc) {
                                                    mpz_set(dlp_part->ax_b[pi], ax_b);
                                                    mpz_set(dlp_part->Qx[pi], aQx);
                                                    dlp_part->lp[pi] = cof;
                                                    lp_insert(dlp_lpt, key, pi);
                                                    dlp_part->count++;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        mpz_clear(aQx);
                    }
                }
            }
        }
        mpz_clear(cur_b);
        free(cur_s1);
        free(cur_s2);
    }

    double sieve_time = elapsed();
    fprintf(stderr, "Sieving: %d rels (%d full + %d SLP + %d DLP) in %.2fs, %d polys\n",
            full->count, full->count - combined - dlp_combined, combined, dlp_combined,
            sieve_time, total_polys * BATCH_POLYS);

    if (full->count < fb->size + 1) {
        fprintf(stderr, "FAIL: not enough relations\n"); printf("FAIL\n"); return 1;
    }

    /* Linear algebra */
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
    int nrels = full->count;
    if (nrels > target) nrels = target;
    int ncols = fb->size + 1;
    gf2_t *mat = gf2_create(nrels, ncols);

    for (int r = 0; r < nrels; r++) {
        mpz_t Qval; mpz_init(Qval); mpz_set(Qval, full->Qx[r]);
        if (mpz_sgn(Qval) < 0) { gf2_set(mat, r, 0); mpz_neg(Qval, Qval); }
        int e2 = 0; while (mpz_even_p(Qval)) { mpz_tdiv_q_2exp(Qval, Qval, 1); e2++; }
        if (e2 & 1) gf2_set(mat, r, 1);
        for (int i = 1; i < fb->size; i++) {
            unsigned int p = fb->prime[i]; int e = 0;
            while (mpz_divisible_ui_p(Qval, p)) { mpz_divexact_ui(Qval, Qval, p); e++; }
            if (e & 1) gf2_set(mat, r, i+1);
        }
        mpz_clear(Qval);
    }

    int **deps; int *dlen;
    int ndeps = gf2_solve(mat, &deps, &dlen, 64);
<<<<<<< HEAD
    fprintf(stderr, "LA: %d deps from %dx%d matrix in %.2fs\n", ndeps, nrels, ncols, elapsed());

    /* ===== Square Root ===== */
    for (int d = 0; d < ndeps; d++) {
        mpz_t X, Y, g, prod, rem;
        mpz_inits(X, Y, g, prod, rem, NULL);
        mpz_set_ui(X, 1);
        for (int k = 0; k < dlen[d]; k++) { mpz_mul(X, X, full->ax_b[deps[d][k]]); mpz_mod(X, X, N); }

        mpz_set_ui(prod, 1);
        for (int k = 0; k < dlen[d]; k++) {
            mpz_t aq2; mpz_init(aq2); mpz_abs(aq2, full->Qx[deps[d][k]]);
            mpz_mul(prod, prod, aq2); mpz_clear(aq2);
=======
    fprintf(stderr, "LA: %d deps from %dx%d (%.2fs)\n", ndeps, nrels, ncols, elapsed());

    /* Square root */
    for (int d = 0; d < ndeps; d++) {
        mpz_t X, Y, g, prod, rem; mpz_inits(X, Y, g, prod, rem, NULL);
        mpz_set_ui(X, 1);
        for (int k = 0; k < dlen[d]; k++) {
            mpz_mul(X, X, full->ax_b[deps[d][k]]);
            mpz_mod(X, X, N);
        }

        mpz_set_ui(Y, 1); mpz_set_ui(prod, 1);
        for (int k = 0; k < dlen[d]; k++) {
            mpz_t aq; mpz_init(aq);
            mpz_abs(aq, full->Qx[deps[d][k]]);
            mpz_mul(prod, prod, aq);
            mpz_clear(aq);
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
        }

        mpz_set(rem, prod);
        int e2 = 0; while (mpz_even_p(rem)) { mpz_tdiv_q_2exp(rem, rem, 1); e2++; }
        if (e2 & 1) goto next;
<<<<<<< HEAD
        mpz_set_ui(Y, 1);
=======
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
        if (e2/2 > 0) { mpz_set_ui(tmp, 2); mpz_powm_ui(tmp, tmp, e2/2, N); mpz_mul(Y, Y, tmp); mpz_mod(Y, Y, N); }

        { int valid = 1;
          for (int i = 1; i < fb->size; i++) {
              unsigned int p = fb->prime[i]; int e = 0;
              while (mpz_divisible_ui_p(rem, p)) { mpz_divexact_ui(rem, rem, p); e++; }
              if (e & 1) { valid = 0; break; }
              if (e/2 > 0) { mpz_set_ui(tmp, p); mpz_powm_ui(tmp, tmp, e/2, N); mpz_mul(Y, Y, tmp); mpz_mod(Y, Y, N); }
          }
          if (!valid) goto next;
        }
        if (mpz_cmp_ui(rem, 1) != 0) {
            if (mpz_perfect_square_p(rem)) {
                mpz_sqrt(tmp, rem); mpz_mod(tmp, tmp, N);
                mpz_mul(Y, Y, tmp); mpz_mod(Y, Y, N);
            } else goto next;
        }

        mpz_sub(tmp, X, Y); mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t o; mpz_init(o); mpz_divexact(o, N, g);
            if (mpz_cmp(g, o) > 0) mpz_swap(g, o);
            gmp_printf("%Zd\n", g);
<<<<<<< HEAD
            fprintf(stderr, "SIQS-Bucket: factored %dd in %.3fs\n", digits, elapsed());
=======
            fprintf(stderr, "SIQS-Bucket: factored in %.3fs\n", elapsed());
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
            mpz_clear(o); return 0;
        }
        mpz_add(tmp, X, Y); mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t o; mpz_init(o); mpz_divexact(o, N, g);
            if (mpz_cmp(g, o) > 0) mpz_swap(g, o);
            gmp_printf("%Zd\n", g);
<<<<<<< HEAD
            fprintf(stderr, "SIQS-Bucket: factored %dd in %.3fs\n", digits, elapsed());
=======
            fprintf(stderr, "SIQS-Bucket: factored in %.3fs\n", elapsed());
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
            mpz_clear(o); return 0;
        }
        next: mpz_clears(X, Y, g, prod, rem, NULL);
    }

<<<<<<< HEAD
    fprintf(stderr, "FAIL: no factor found from %d dependencies\n", ndeps);
=======
    fprintf(stderr, "SIQS-Bucket: FAILED after %.3fs\n", elapsed());
>>>>>>> c671567 (Add SIQS-Bucket with Gray code self-init and YAFU-calibrated params)
    printf("FAIL\n");
    return 1;
}
