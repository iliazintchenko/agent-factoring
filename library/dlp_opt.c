/*
 * DLP-Opt: SIQS with Double Large Prime and LP-column matrix
 *
 * Key features:
 * 1. Bucket sieve for large primes (p >= BLOCK_SIZE)
 * 2. SLP + DLP: accept cofactors with 1 or 2 large primes
 * 3. Pollard rho for fast cofactor splitting (up to 64 bits)
 * 4. LP columns in GF(2) matrix (no graph cycle detection needed)
 * 5. All relations (full, SLP, DLP) go directly into the matrix
 * 6. More aggressive sieve threshold to catch DLP candidates
 *
 * The LP-column approach: instead of combining partials by matching LPs,
 * we add a column for each distinct LP. A partial with LP l has an odd
 * exponent in column l. When Gaussian elimination finds a dependency,
 * it automatically ensures all LP columns have even exponents = LPs cancel.
 *
 * Compile: gcc -O3 -march=native -o dlp_opt library/dlp_opt.c -lgmp -lm
 * Usage: ./dlp_opt <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <gmp.h>

#define SEED 42
#define BLOCK_SIZE 32768
#define BLOCK_BITS 15
#define MAX_FB 80000
#define MAX_A_FACTORS 20
#define MAX_RELS 300000
#define MAX_PARTIALS 1500000
#define BUCKET_ALLOC 2048

static struct timespec g_start;
static double elapsed(void) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec - g_start.tv_sec) + (now.tv_nsec - g_start.tv_nsec) / 1e9;
}

/* ==================== Modular Arithmetic ==================== */
static unsigned int mod_inverse(unsigned int a, unsigned int m) {
    int old_r = (int)a, r = (int)m, old_s = 1, s = 0;
    while (r) { int q2 = old_r / r, t = r; r = old_r - q2 * r; old_r = t; t = s; s = old_s - q2 * s; old_s = t; }
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
    unsigned long long M2 = S;
    b = z; e = Q; unsigned long long c = 1; while (e) { if (e & 1) c = (c * b) % m; b = (b * b) % m; e >>= 1; }
    b = n % p; e = Q; unsigned long long t = 1; while (e) { if (e & 1) t = (t * b) % m; b = (b * b) % m; e >>= 1; }
    b = n % p; e = (Q + 1) / 2; unsigned long long R2 = 1; while (e) { if (e & 1) R2 = (R2 * b) % m; b = (b * b) % m; e >>= 1; }
    while (1) { if (t == 1) return (unsigned int)R2; int i = 0; unsigned long long tt = t; while (tt != 1) { tt = (tt * tt) % p; i++; } unsigned long long bb2 = c; for (int j = 0; j < (int)M2 - i - 1; j++) bb2 = (bb2 * bb2) % p; M2 = i; c = (bb2 * bb2) % p; t = (t * c) % p; R2 = (R2 * bb2) % p; }
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
        int ps[] = {3,5,7,11,13,17,19,23,29,31,37,41,43,47};
        for (int i = 0; i < 14; i++) {
            if (k % ps[i] == 0) { s += log(ps[i]); continue; }
            if (sqrt_mod(mpz_fdiv_ui(kN, ps[i]), ps[i])) s += 2.0*log(ps[i])/(ps[i]-1);
        }
        if (s > best) { best = s; best_k = k; }
        mpz_clear(kN);
    }
    return best_k;
}

/* ==================== Factor Base ==================== */
typedef struct {
    unsigned int *prime;
    unsigned int *root;
    unsigned char *logp;
    int size;
    int large_start;  /* index where p >= BLOCK_SIZE */
} fb_t;

static fb_t *fb_create(mpz_t kN, int target) {
    fb_t *fb = calloc(1, sizeof(fb_t));
    int alloc = target + 100;
    fb->prime = malloc(alloc * sizeof(unsigned int));
    fb->root = malloc(alloc * sizeof(unsigned int));
    fb->logp = malloc(alloc * sizeof(unsigned char));
    fb->prime[0] = 2; fb->root[0] = 1; fb->logp[0] = 1; fb->size = 1;
    int bound = target * 30 + 50000;
    char *sv = calloc(bound + 1, 1);
    for (int i = 2; (long)i*i <= bound; i++) if (!sv[i]) for (int j = i*i; j <= bound; j += i) sv[j] = 1;
    for (int i = 3; i <= bound && fb->size < target; i += 2) {
        if (sv[i]) continue;
        unsigned long nm = mpz_fdiv_ui(kN, i);
        if (nm == 0) { fb->prime[fb->size] = i; fb->root[fb->size] = 0; fb->logp[fb->size] = (unsigned char)(log2(i)+0.5); fb->size++; continue; }
        unsigned int r = sqrt_mod((unsigned int)nm, i);
        if (!r) continue;
        fb->prime[fb->size] = i; fb->root[fb->size] = r; fb->logp[fb->size] = (unsigned char)(log2(i)+0.5); fb->size++;
    }
    free(sv);
    fb->large_start = fb->size;
    for (int i = 0; i < fb->size; i++) {
        if (fb->prime[i] >= BLOCK_SIZE) { fb->large_start = i; break; }
    }
    return fb;
}

/* ==================== Bucket Sieve ==================== */
typedef struct {
    uint32_t *data;
    uint32_t *count;
    uint8_t *slice_logp;
    uint32_t *slice_bound;
    int num_blocks, num_slices, max_slices;
} bucket_t;

static bucket_t *bucket_create(int nblocks, int max_slices) {
    bucket_t *b = calloc(1, sizeof(bucket_t));
    b->num_blocks = nblocks;
    b->max_slices = max_slices;
    b->data = malloc((size_t)max_slices * nblocks * BUCKET_ALLOC * sizeof(uint32_t));
    b->count = calloc((size_t)max_slices * nblocks, sizeof(uint32_t));
    b->slice_logp = malloc(max_slices);
    b->slice_bound = malloc(max_slices * sizeof(uint32_t));
    return b;
}

static void bucket_fill(bucket_t *bkt, fb_t *fb, unsigned int *soln1, unsigned int *soln2,
                         int nblocks, int sieve_start) {
    memset(bkt->count, 0, (size_t)bkt->max_slices * nblocks * sizeof(uint32_t));
    bkt->num_slices = 0;

    int fb_start = fb->large_start;
    int fb_end = fb->size;
    if (fb_start >= fb_end) { bkt->num_slices = 0; return; }

    int slice = 0;
    int bound_val = fb_start;
    bkt->slice_bound[0] = fb_start;
    bkt->slice_logp[0] = fb->logp[fb_start];
    long interval_size = (long)nblocks * BLOCK_SIZE;
    int check_interval = 64; /* check bucket fullness every N primes */

    for (int i = fb_start; i < fb_end; i++) {
        unsigned int p = fb->prime[i];
        if (soln1[i] == 0xFFFFFFFF) continue;

        /* Check if we need a new slice: either FB index overflow or buckets getting full */
        if (i - bound_val >= 65535 || ((i - fb_start) % check_interval == 0 && i > fb_start)) {
            int need_new = (i - bound_val >= 65535);
            if (!need_new) {
                /* Check max bucket fullness */
                uint32_t max_count = 0;
                for (int bl = 0; bl < nblocks; bl++) {
                    int idx = slice * nblocks + bl;
                    if (bkt->count[idx] > max_count) max_count = bkt->count[idx];
                }
                if (max_count > BUCKET_ALLOC * 3 / 4) need_new = 1;
            }
            if (need_new) {
                slice++;
                if (slice >= bkt->max_slices) { slice--; break; }
                bound_val = i;
                bkt->slice_bound[slice] = i;
                bkt->slice_logp[slice] = fb->logp[i];
            }
        }

        for (int rt = 0; rt < 2; rt++) {
            unsigned int root_val = (rt == 0) ? soln1[i] : soln2[i];
            if (rt == 1 && soln1[i] == soln2[i]) continue;

            long pos = ((long)root_val - sieve_start) % (long)p;
            if (pos < 0) pos += p;

            while (pos < interval_size) {
                int bnum = (int)(pos >> BLOCK_BITS);
                int offset = (int)(pos & (BLOCK_SIZE - 1));
                int idx = slice * nblocks + bnum;
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

static void bucket_apply(unsigned char *sieve, bucket_t *bkt, int bnum) {
    for (int s = 0; s < bkt->num_slices; s++) {
        int idx = s * bkt->num_blocks + bnum;
        uint32_t n = bkt->count[idx];
        uint32_t *entries = &bkt->data[(size_t)idx * BUCKET_ALLOC];
        uint8_t lp = bkt->slice_logp[s];
        for (uint32_t i = 0; i < n; i++)
            sieve[entries[i] & 0xFFFF] += lp;
    }
}

/* ==================== Pollard Rho for cofactor splitting ==================== */
static int pollard_rho_64(uint64_t n, uint64_t *factor) {
    if (n <= 1) return 0;
    if (n % 2 == 0) { *factor = 2; return 1; }
    /* Simple Brent rho */
    for (uint64_t c = 1; c < 100; c++) {
        __uint128_t x = 2, y = 2, d = 1;
        uint64_t q = 1;
        int steps = 0;
        while (d == 1 && steps < 1000000) {
            for (int batch = 0; batch < 128 && d == 1; batch++) {
                x = ((__uint128_t)x * x + c) % n;
                x = ((__uint128_t)x * x + c) % n;
                y = ((__uint128_t)y * y + c) % n;
                uint64_t diff = ((uint64_t)x > (uint64_t)y) ? (uint64_t)x - (uint64_t)y : (uint64_t)y - (uint64_t)x;
                q = (uint64_t)(((__uint128_t)q * diff) % n);
                steps++;
            }
            /* GCD */
            uint64_t a2 = q, b2 = n;
            while (b2) { uint64_t t = b2; b2 = a2 % b2; a2 = t; }
            d = a2;
        }
        if (d > 1 && d < n) { *factor = d; return 1; }
    }
    return 0;
}

/* Simple primality test for 64-bit numbers */
static int is_prime_64(uint64_t n) {
    if (n < 2) return 0;
    if (n < 4) return 1;
    if (n % 2 == 0 || n % 3 == 0) return 0;
    /* Miller-Rabin with deterministic bases for n < 2^64 */
    uint64_t d = n - 1; int r = 0;
    while (d % 2 == 0) { d /= 2; r++; }
    uint64_t witnesses[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (int w = 0; w < 12; w++) {
        uint64_t a = witnesses[w];
        if (a >= n) continue;
        /* Compute a^d mod n */
        __uint128_t x = 1, base = a;
        uint64_t exp = d;
        while (exp > 0) {
            if (exp & 1) x = (x * base) % n;
            base = (base * base) % n;
            exp >>= 1;
        }
        uint64_t xv = (uint64_t)x;
        if (xv == 1 || xv == n - 1) continue;
        int composite = 1;
        for (int i = 0; i < r - 1; i++) {
            xv = (uint64_t)(((__uint128_t)xv * xv) % n);
            if (xv == n - 1) { composite = 0; break; }
        }
        if (composite) return 0;
    }
    return 1;
}

/* ==================== LP Index Map ==================== */
/* Maps large primes to column indices in the GF(2) matrix */
#define LP_MAP_BITS 22
#define LP_MAP_SIZE (1 << LP_MAP_BITS)
typedef struct lp_node { uint64_t lp; int col_idx; struct lp_node *next; } lp_node_t;
typedef struct {
    lp_node_t **buckets;
    lp_node_t *pool;
    int used, max;
    int next_col;  /* next available LP column index */
} lp_map_t;

static lp_map_t *lp_map_create(int max_lps, int first_col) {
    lp_map_t *m = calloc(1, sizeof(lp_map_t));
    m->buckets = calloc(LP_MAP_SIZE, sizeof(lp_node_t *));
    m->pool = calloc(max_lps, sizeof(lp_node_t));
    m->max = max_lps;
    m->next_col = first_col;
    return m;
}

/* Get or create a column index for this LP */
static int lp_map_get(lp_map_t *m, uint64_t lp) {
    uint32_t h = (uint32_t)((lp * 0x9E3779B97F4A7C15ULL) >> (64 - LP_MAP_BITS));
    for (lp_node_t *n = m->buckets[h]; n; n = n->next)
        if (n->lp == lp) return n->col_idx;
    /* Not found: create new */
    if (m->used >= m->max) return -1;
    lp_node_t *n = &m->pool[m->used++];
    n->lp = lp;
    n->col_idx = m->next_col++;
    n->next = m->buckets[h];
    m->buckets[h] = n;
    return n->col_idx;
}

/* ==================== Relation Storage ==================== */
#define MAX_LP_PER_REL 2
typedef struct {
    mpz_t *ax_b, *Qx;
    uint64_t lps[MAX_RELS][MAX_LP_PER_REL];  /* large primes */
    int nlps[MAX_RELS];                        /* 0, 1, or 2 LPs */
    int count, alloc;
} rels_t;

static rels_t *rels_create(int n) {
    rels_t *r = calloc(1, sizeof(rels_t));
    r->ax_b = malloc(n * sizeof(mpz_t));
    r->Qx = malloc(n * sizeof(mpz_t));
    for (int i = 0; i < n; i++) { mpz_init(r->ax_b[i]); mpz_init(r->Qx[i]); }
    r->alloc = n;
    return r;
}

/* ==================== GF(2) Matrix ==================== */
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
        int pr = -1;
        for (int r = piv; r < m->nr; r++) if ((m->rows[r][c/64] >> (c%64)) & 1) { pr = r; break; }
        if (pr < 0) continue;
        if (pr != piv) { u64 *t = m->rows[pr]; m->rows[pr] = m->rows[piv]; m->rows[piv] = t; }
        for (int r = 0; r < m->nr; r++) {
            if (r == piv) continue;
            if ((m->rows[r][c/64] >> (c%64)) & 1)
                for (int w = 0; w < m->wprow; w++) m->rows[r][w] ^= m->rows[piv][w];
        }
        piv++;
    }
    int nd = 0; *deps = malloc(max * sizeof(int*)); *dlen = malloc(max * sizeof(int));
    for (int r = piv; r < m->nr && nd < max; r++) {
        int z = 1;
        for (int w = 0; w < m->fbw && z; w++) {
            u64 mask = (w < m->fbw-1) ? ~0ULL : (m->nc%64==0 ? ~0ULL : (1ULL << (m->nc%64))-1);
            if (m->rows[r][w] & mask) z = 0;
        }
        if (!z) continue;
        int *d = malloc(m->nr * sizeof(int)); int dl = 0;
        for (int w = 0; w < m->idw; w++) {
            u64 bits = m->rows[r][m->fbw+w];
            while (bits) { int bit = __builtin_ctzll(bits); int ii = w*64+bit;
                if (ii < m->nr) d[dl++] = ii; bits &= bits-1; }
        }
        if (dl > 0) { (*deps)[nd] = d; (*dlen)[nd] = dl; nd++; } else free(d);
    }
    return nd;
}

/* ==================== Parameters ==================== */
/* DLP enabled at 70+ digits. dlp_exp: DLP bound = FB_max ^ dlp_exp */
typedef struct { int fb_size, nblocks, lp_mult, extra; double thresh, dlp_exp; } params_t;
static params_t get_params(int bits) {
    if (bits <= 100) return (params_t){120, 1, 40, 40, 0.72, 0};
    if (bits <= 110) return (params_t){180, 2, 40, 45, 0.73, 0};
    if (bits <= 120) return (params_t){250, 3, 45, 50, 0.75, 0};
    if (bits <= 130) return (params_t){350, 4, 50, 55, 0.77, 0};
    if (bits <= 140) return (params_t){500, 5, 55, 60, 0.78, 0};
    if (bits <= 150) return (params_t){700, 7, 60, 70, 0.79, 0};
    if (bits <= 160) return (params_t){1000, 10, 65, 80, 0.80, 0};
    if (bits <= 170) return (params_t){1400, 14, 70, 90, 0.81, 0};
    if (bits <= 180) return (params_t){2000, 18, 75, 100, 0.82, 0};
    if (bits <= 190) return (params_t){2800, 24, 80, 110, 0.83, 0};
    if (bits <= 200) return (params_t){4000, 30, 85, 120, 0.84, 0};
    if (bits <= 210) return (params_t){5500, 38, 90, 140, 0.85, 0};
    if (bits <= 220) return (params_t){7500, 46, 95, 160, 0.86, 1.75};
    if (bits <= 232) return (params_t){10000, 54, 100, 180, 0.85, 1.75};
    if (bits <= 248) return (params_t){16000, 64, 105, 200, 0.85, 1.75};
    if (bits <= 265) return (params_t){30000, 76, 110, 250, 0.86, 1.75};
    if (bits <= 281) return (params_t){45000, 90, 115, 300, 0.87, 1.8};
    if (bits <= 298) return (params_t){60000, 100, 120, 350, 0.88, 1.8};
    return (params_t){72000, 120, 135, 400, 0.89, 1.8};
}

/* ==================== Main ==================== */
int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }
    clock_gettime(CLOCK_MONOTONIC, &g_start);

    mpz_t N, kN;
    mpz_inits(N, kN, NULL);
    mpz_set_str(N, argv[1], 10);

    int digits = (int)mpz_sizeinbase(N, 10);
    int bits = (int)mpz_sizeinbase(N, 2);

    /* Trial division */
    for (int p2 = 2; p2 < 10000; p2++) {
        if (mpz_divisible_ui_p(N, p2)) {
            mpz_t c; mpz_init(c); mpz_divexact_ui(c, N, p2);
            gmp_printf("%Zd\n", c);
            return 0;
        }
    }

    int mult = choose_multiplier(N);
    mpz_mul_ui(kN, N, mult);
    int kN_bits = (int)mpz_sizeinbase(kN, 2);
    params_t P = get_params(kN_bits);

    fb_t *fb = fb_create(kN, P.fb_size);
    int M = BLOCK_SIZE * P.nblocks;
    int total_blocks = 2 * P.nblocks;
    unsigned long lp_bound = (unsigned long)fb->prime[fb->size-1] * P.lp_mult;
    /* DLP bound: each large prime must be <= dlp_bound */
    unsigned long dlp_bound = 0;
    int use_dlp = (P.dlp_exp > 0);
    if (use_dlp) {
        double fb_max_d = (double)fb->prime[fb->size-1];
        dlp_bound = (unsigned long)pow(fb_max_d, P.dlp_exp);
        if (dlp_bound < lp_bound) dlp_bound = lp_bound;
    }
    /* Target: for SLP-only (combined by matching), need FB+extra rels.
     * For DLP (LP columns), need more since LP columns add to matrix. */
    int target = fb->size + P.extra;
    if (use_dlp) target = (int)(fb->size * 2.0) + P.extra;

    double log2_Qmax = kN_bits / 2.0 + 0.5 + log2(M);
    int threshold = (int)(log2_Qmax * P.thresh);
    /* DLP: lower threshold by a few bits to get more candidates */
    if (use_dlp) threshold -= 4; else threshold -= 2;

    fprintf(stderr, "DLP-Opt: %dd (%db), k=%d, FB=%d (large@%d), M=%d, blocks=%d, thresh=%d, SLP=%lu%s, target=%d\n",
            digits, bits, mult, fb->size, fb->large_start, M, total_blocks, threshold, lp_bound,
            use_dlp ? "" : " (no DLP)", target);
    if (use_dlp) fprintf(stderr, "  DLP bound: %lu (FB_max^%.2f)\n", dlp_bound, P.dlp_exp);

    unsigned char *sieve = aligned_alloc(64, BLOCK_SIZE);
    bucket_t *buckets = bucket_create(total_blocks, 64);

    unsigned int *soln1 = malloc(fb->size * sizeof(unsigned int));
    unsigned int *soln2 = malloc(fb->size * sizeof(unsigned int));

    rels_t *rels = rels_create(MAX_RELS);  /* ALL relations: full, combined SLP, DLP */
    lp_map_t *lp_map = lp_map_create(500000, fb->size + 1); /* LP column mapping for DLP */

    /* SLP matching hash table (same as siqs_opt) */
    lp_node_t **slp_buckets = calloc(LP_MAP_SIZE, sizeof(lp_node_t *));
    lp_node_t *slp_pool = calloc(MAX_PARTIALS, sizeof(lp_node_t));
    int slp_pool_used = 0;
    /* Partial storage for SLP matching */
    mpz_t *slp_ax_b = malloc(MAX_PARTIALS * sizeof(mpz_t));
    mpz_t *slp_Qx = malloc(MAX_PARTIALS * sizeof(mpz_t));
    for (int i = 0; i < MAX_PARTIALS; i++) { mpz_init(slp_ax_b[i]); mpz_init(slp_Qx[i]); }
    int slp_count = 0;

    mpz_t a, b_val, c_val, B_vals[MAX_A_FACTORS];
    mpz_inits(a, b_val, c_val, NULL);
    for (int j = 0; j < MAX_A_FACTORS; j++) mpz_init(B_vals[j]);

    unsigned int **ainv_arr = malloc(MAX_A_FACTORS * sizeof(unsigned int *));
    for (int j = 0; j < MAX_A_FACTORS; j++)
        ainv_arr[j] = malloc(fb->size * sizeof(unsigned int));

    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, SEED);

    mpz_t ax_b, Qx, residue, tmp;
    mpz_inits(ax_b, Qx, residue, tmp, NULL);

    int total_polys = 0, a_count = 0;
    int n_full = 0, n_slp = 0, n_dlp = 0;
    int a_idx[MAX_A_FACTORS];
    int num_a_factors = 0;
    int sieve_start = -M;

    while (rels->count < target) {
        if (elapsed() > 280) {
            fprintf(stderr, "TIMEOUT at %.1fs\n", elapsed());
            break;
        }

        /* Generate new 'a' */
        {
            mpz_t tgt; mpz_init(tgt);
            mpz_mul_ui(tgt, kN, 2); mpz_sqrt(tgt, tgt); mpz_tdiv_q_ui(tgt, tgt, M);
            double log_tgt = mpz_sizeinbase(tgt, 2) * log(2.0);

            int lo = fb->size / 3, hi = 2 * fb->size / 3;
            if (lo < 2) lo = 2; if (hi <= lo + 3) hi = fb->size - 1;

            double avg = 0; int cnt = 0;
            for (int i = lo; i < hi; i++) { if (fb->root[i] == 0) continue; avg += log(fb->prime[i]); cnt++; }
            if (cnt == 0) { mpz_clear(tgt); break; }
            avg /= cnt;

            int s = (int)(log_tgt / avg + 0.5);
            if (s < 3) s = 3; if (s > MAX_A_FACTORS) s = MAX_A_FACTORS; if (s > hi - lo) s = hi - lo;
            num_a_factors = s;

            double best_ratio = 1e30;
            int best[MAX_A_FACTORS];
            for (int att = 0; att < 50; att++) {
                mpz_set_ui(a, 1);
                int idx[MAX_A_FACTORS]; int ok = 1;
                for (int i2 = 0; i2 < s && ok; i2++) {
                    int tries = 0, good;
                    do { idx[i2] = lo + gmp_urandomm_ui(rng, hi-lo); good = 1;
                         for (int j = 0; j < i2; j++) if (idx[j]==idx[i2]) {good=0; break;}
                         if (fb->root[idx[i2]]==0) good=0; tries++;
                    } while (!good && tries < 100);
                    if (!good) { ok=0; break; }
                    mpz_mul_ui(a, a, fb->prime[idx[i2]]);
                }
                if (!ok) continue;
                double ratio;
                if (mpz_cmp(a, tgt) > 0) { mpz_t q; mpz_init(q); mpz_tdiv_q(q, a, tgt); ratio = mpz_get_d(q); mpz_clear(q); }
                else { mpz_t q; mpz_init(q); mpz_tdiv_q(q, tgt, a); ratio = mpz_get_d(q); mpz_clear(q); }
                if (ratio < best_ratio) { best_ratio = ratio; memcpy(best, idx, s*sizeof(int)); }
                if (ratio < 1.5) break;
            }
            memcpy(a_idx, best, s * sizeof(int));
            mpz_set_ui(a, 1);
            for (int i2 = 0; i2 < s; i2++) mpz_mul_ui(a, a, fb->prime[a_idx[i2]]);
            mpz_clear(tgt);
            a_count++;

            for (int j = 0; j < s; j++) {
                int idx2 = a_idx[j];
                unsigned int qj = fb->prime[idx2], rj = fb->root[idx2];
                mpz_t a_q, mod_q, inv2; mpz_inits(a_q, mod_q, inv2, NULL);
                mpz_divexact_ui(a_q, a, qj); mpz_set_ui(mod_q, qj);
                mpz_invert(inv2, a_q, mod_q);
                unsigned long iv = mpz_get_ui(inv2);
                mpz_mul_ui(B_vals[j], a_q, (unsigned long)rj * iv % qj);
                mpz_clears(a_q, mod_q, inv2, NULL);
            }

            for (int j = 0; j < s; j++) {
                for (int i = 0; i < fb->size; i++) {
                    unsigned int p = fb->prime[i];
                    unsigned long am = mpz_fdiv_ui(a, p);
                    if (am == 0 || fb->root[i] == 0) { ainv_arr[j][i] = 0; continue; }
                    unsigned int ai = mod_inverse((unsigned int)am, p);
                    unsigned long Bm = mpz_fdiv_ui(B_vals[j], p);
                    ainv_arr[j][i] = (unsigned int)((2ULL * ai % p * Bm) % p);
                }
            }
        }

        int num_b = 1 << (num_a_factors - 1);

        /* Initial b */
        mpz_set_ui(b_val, 0);
        for (int j = 0; j < num_a_factors; j++) mpz_add(b_val, b_val, B_vals[j]);
        mpz_mul(c_val, b_val, b_val); mpz_sub(c_val, c_val, kN);
        if (!mpz_divisible_p(c_val, a)) {
            mpz_neg(b_val, b_val);
            mpz_mul(c_val, b_val, b_val); mpz_sub(c_val, c_val, kN);
        }
        mpz_divexact(c_val, c_val, a);

        /* Initial solutions */
        for (int i = 0; i < fb->size; i++) {
            unsigned int p = fb->prime[i];
            unsigned long am = mpz_fdiv_ui(a, p);
            if (am == 0 || fb->root[i] == 0) { soln1[i] = soln2[i] = 0xFFFFFFFF; continue; }
            unsigned int ai = mod_inverse((unsigned int)am, p);
            if (ai == 0) { soln1[i] = soln2[i] = 0xFFFFFFFF; continue; }
            unsigned long bm = mpz_fdiv_ui(b_val, p);
            unsigned int r = fb->root[i];
            soln1[i] = (unsigned int)((unsigned long)ai * ((r + p - bm) % p) % p);
            soln2[i] = (unsigned int)((unsigned long)ai * ((p - r + p - bm) % p) % p);
        }

        for (int b_idx = 0; b_idx < num_b && rels->count < target; b_idx++) {
            total_polys++;

            if (total_polys % 500 == 0) {
                double t = elapsed();
                if (t > 280) break;
                if (total_polys % 2000 == 0)
                    fprintf(stderr, "  polys=%d rels=%d/%d (full=%d slp=%d dlp=%d) t=%.1fs\n",
                            total_polys, rels->count, target,
                            n_full, n_slp, n_dlp, t);
            }

            /* Gray code update */
            if (b_idx > 0) {
                int j = __builtin_ctz(b_idx);
                int sign = ((b_idx >> j) & 2) ? 1 : -1;
                if (sign > 0) mpz_addmul_ui(b_val, B_vals[j], 2);
                else mpz_submul_ui(b_val, B_vals[j], 2);
                mpz_mul(c_val, b_val, b_val); mpz_sub(c_val, c_val, kN);
                mpz_divexact(c_val, c_val, a);

                for (int i = 0; i < fb->size; i++) {
                    if (soln1[i] == 0xFFFFFFFF) continue;
                    unsigned int p = fb->prime[i];
                    unsigned int delta = ainv_arr[j][i];
                    if (sign > 0) {
                        soln1[i] = (soln1[i] >= delta) ? soln1[i] - delta : soln1[i] + p - delta;
                        soln2[i] = (soln2[i] >= delta) ? soln2[i] - delta : soln2[i] + p - delta;
                    } else {
                        soln1[i] += delta; if (soln1[i] >= p) soln1[i] -= p;
                        soln2[i] += delta; if (soln2[i] >= p) soln2[i] -= p;
                    }
                }
            }

            /* Fill buckets for large primes */
            if (fb->large_start < fb->size)
                bucket_fill(buckets, fb, soln1, soln2, total_blocks, sieve_start);

            /* ===== Sieve each block ===== */
            for (int bnum = 0; bnum < total_blocks; bnum++) {
                int block_start = sieve_start + bnum * BLOCK_SIZE;

                memset(sieve, 0, BLOCK_SIZE);

                /* Small/medium primes: per-block offset computation */
                for (int i = 1; i < fb->large_start; i++) {
                    unsigned int p = fb->prime[i];
                    if (p < 5) continue;
                    if (soln1[i] == 0xFFFFFFFF) continue;
                    unsigned char lp = fb->logp[i];

                    long o1 = ((long)soln1[i] - block_start) % (long)p;
                    if (o1 < 0) o1 += p;
                    for (unsigned int j = (unsigned int)o1; j < BLOCK_SIZE; j += p)
                        sieve[j] += lp;

                    if (soln1[i] != soln2[i]) {
                        long o2 = ((long)soln2[i] - block_start) % (long)p;
                        if (o2 < 0) o2 += p;
                        for (unsigned int j = (unsigned int)o2; j < BLOCK_SIZE; j += p)
                            sieve[j] += lp;
                    }
                }

                /* Large primes: bucket sieve */
                if (fb->large_start < fb->size)
                    bucket_apply(sieve, buckets, bnum);

                /* ===== Scan for smooth candidates ===== */
                /* Use 64-bit word scan: check if any byte in word >= threshold */
                uint64_t thresh_mask = 0;
                for (int k = 0; k < 8; k++) thresh_mask |= ((uint64_t)(threshold - 1) << (k * 8));
                /* We look for bytes > threshold-1, i.e. >= threshold */

                for (int j = 0; j < BLOCK_SIZE; j += 8) {
                    /* Quick check: any byte in this 8-byte word above threshold? */
                    uint64_t word = *(uint64_t *)(sieve + j);
                    /* Check if any byte >= threshold using saturating subtraction trick */
                    /* A byte b >= threshold iff (b | 0x80) after subtracting threshold-1 has high bit set */
                    /* Simpler: just check each byte if the word has promising values */
                    uint64_t high_bits = word & 0x8080808080808080ULL;
                    if (!high_bits && threshold > 127) continue;
                    /* For threshold <= 127, need more careful check */
                    /* Just check byte by byte for now */
                    for (int k = 0; k < 8; k++) {
                        if (sieve[j + k] < threshold) continue;
                        long x = (long)(block_start + j + k);
                        if (x == 0) continue;

                        mpz_set_si(tmp, x);
                        mpz_mul(Qx, a, tmp); mpz_add(Qx, Qx, b_val); mpz_add(Qx, Qx, b_val);
                        mpz_mul(Qx, Qx, tmp); mpz_add(Qx, Qx, c_val);

                        mpz_mul_si(ax_b, a, x); mpz_add(ax_b, ax_b, b_val);

                        if (mpz_sgn(Qx) == 0) continue;
                        mpz_abs(residue, Qx);

                        /* Trial divide by FB primes */
                        while (mpz_even_p(residue)) mpz_tdiv_q_2exp(residue, residue, 1);
                        for (int i = 1; i < fb->size; i++) {
                            unsigned int p = fb->prime[i];
                            if (soln1[i] == 0xFFFFFFFF) continue;
                            long xmod = ((x % (long)p) + p) % p;
                            if (xmod != (long)soln1[i] && xmod != (long)soln2[i]) continue;
                            if (mpz_divisible_ui_p(residue, p))
                                do { mpz_divexact_ui(residue, residue, p); } while (mpz_divisible_ui_p(residue, p));
                        }
                        for (int i = 0; i < fb->size && fb->prime[i] < 5; i++) {
                            unsigned int p = fb->prime[i]; if (p <= 2) continue;
                            while (mpz_divisible_ui_p(residue, p)) mpz_divexact_ui(residue, residue, p);
                        }

                        mpz_t aQx; mpz_init(aQx); mpz_mul(aQx, Qx, a);

                        if (rels->count < rels->alloc) {
                            if (mpz_cmp_ui(residue, 1) == 0) {
                                /* Full relation */
                                int ri = rels->count;
                                mpz_set(rels->ax_b[ri], ax_b);
                                mpz_set(rels->Qx[ri], aQx);
                                rels->nlps[ri] = 0;
                                rels->count++;
                                n_full++;
                            } else if (mpz_fits_ulong_p(residue)) {
                                uint64_t cof = mpz_get_ui(residue);
                                /* Determine effective LP bound for matching */
                                unsigned long eff_lp_bound = use_dlp ? dlp_bound : lp_bound;
                                if (cof <= eff_lp_bound && is_prime_64(cof)) {
                                    /* SLP: try matching (extended to DLP bound when DLP active) */
                                    uint32_t sh = (uint32_t)((cof * 0x9E3779B97F4A7C15ULL) >> (64 - LP_MAP_BITS));
                                    int matched = 0;
                                    for (lp_node_t *n = slp_buckets[sh]; n; n = n->next) {
                                        if (n->lp == cof) {
                                            int ri = rels->count;
                                            int pi = n->col_idx;
                                            mpz_mul(rels->ax_b[ri], ax_b, slp_ax_b[pi]);
                                            mpz_mod(rels->ax_b[ri], rels->ax_b[ri], N);
                                            mpz_mul(rels->Qx[ri], aQx, slp_Qx[pi]);
                                            rels->nlps[ri] = 0;
                                            rels->count++;
                                            n_slp++;
                                            matched = 1;
                                            break;
                                        }
                                    }
                                    if (!matched && slp_count < MAX_PARTIALS) {
                                        int pi = slp_count;
                                        mpz_set(slp_ax_b[pi], ax_b);
                                        mpz_set(slp_Qx[pi], aQx);
                                        lp_node_t *n = &slp_pool[slp_pool_used++];
                                        n->lp = cof;
                                        n->col_idx = pi;
                                        n->next = slp_buckets[sh];
                                        slp_buckets[sh] = n;
                                        slp_count++;
                                    }
                                } else if (use_dlp && !is_prime_64(cof) && cof <= (uint64_t)dlp_bound * dlp_bound) {
                                    /* DLP: try to split composite cofactor into two primes */
                                    uint64_t f1 = 0;
                                    if (pollard_rho_64(cof, &f1)) {
                                        uint64_t f2 = cof / f1;
                                        if (f1 > f2) { uint64_t t2 = f1; f1 = f2; f2 = t2; }
                                        if (f2 <= dlp_bound && is_prime_64(f1) && is_prime_64(f2)) {
                                            /* DLP relation with LP columns */
                                            int ri = rels->count;
                                            mpz_set(rels->ax_b[ri], ax_b);
                                            mpz_set(rels->Qx[ri], aQx);
                                            rels->lps[ri][0] = f1;
                                            rels->lps[ri][1] = f2;
                                            rels->nlps[ri] = 2;
                                            rels->count++;
                                            n_dlp++;
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
    }

    double sieve_time = elapsed();
    int num_lp_cols = lp_map->next_col - (fb->size + 1);
    fprintf(stderr, "Sieving: %d rels (full=%d slp=%d dlp=%d) LP_cols=%d in %.2fs, %d polys\n",
            rels->count, n_full, n_slp, n_dlp, num_lp_cols, sieve_time, total_polys);

    /* First pass: register all LPs to determine matrix dimensions */
    int nrels = rels->count;
    for (int r = 0; r < nrels; r++) {
        for (int lpi = 0; lpi < rels->nlps[r]; lpi++)
            lp_map_get(lp_map, rels->lps[r][lpi]);
    }

    int ncols = lp_map->next_col; /* fb->size + 1 + num_lp_cols */
    int min_rels = ncols + 10;
    if (nrels < min_rels) {
        fprintf(stderr, "FAIL: not enough relations (%d < %d cols)\n", nrels, ncols);
        printf("FAIL\n"); return 1;
    }
    if (nrels > ncols + P.extra) nrels = ncols + P.extra;

    /* Linear Algebra with LP columns */
    gf2_t *mat = gf2_create(nrels, ncols);

    for (int r = 0; r < nrels; r++) {
        mpz_t Qval; mpz_init(Qval); mpz_set(Qval, rels->Qx[r]);
        if (mpz_sgn(Qval) < 0) { gf2_set(mat, r, 0); mpz_neg(Qval, Qval); }
        int e2 = 0; while (mpz_even_p(Qval)) { mpz_tdiv_q_2exp(Qval, Qval, 1); e2++; }
        if (e2 & 1) gf2_set(mat, r, 1);
        for (int i = 1; i < fb->size; i++) {
            unsigned int p = fb->prime[i]; int e = 0;
            while (mpz_divisible_ui_p(Qval, p)) { mpz_divexact_ui(Qval, Qval, p); e++; }
            if (e & 1) gf2_set(mat, r, i+1);
        }
        /* Add LP columns (already registered in first pass) */
        for (int lpi = 0; lpi < rels->nlps[r]; lpi++) {
            int col = lp_map_get(lp_map, rels->lps[r][lpi]);
            if (col >= 0 && col < ncols) gf2_set(mat, r, col);
        }
        mpz_clear(Qval);
    }

    int **deps; int *dlen;
    int ndeps = gf2_solve(mat, &deps, &dlen, 64);
    fprintf(stderr, "LA: %d deps from %dx%d in %.2fs\n", ndeps, nrels, ncols, elapsed());

    /* Square Root - need to handle LP factors in the product */
    for (int d = 0; d < ndeps; d++) {
        mpz_t X, Y, g, prod, rem;
        mpz_inits(X, Y, g, prod, rem, NULL);
        mpz_set_ui(X, 1);
        for (int k = 0; k < dlen[d]; k++) {
            mpz_mul(X, X, rels->ax_b[deps[d][k]]);
            mpz_mod(X, X, N);
        }
        mpz_set_ui(prod, 1);
        for (int k = 0; k < dlen[d]; k++) {
            mpz_t aq; mpz_init(aq); mpz_abs(aq, rels->Qx[deps[d][k]]);
            mpz_mul(prod, prod, aq); mpz_clear(aq);
        }
        mpz_set(rem, prod);
        int e2 = 0; while (mpz_even_p(rem)) { mpz_tdiv_q_2exp(rem, rem, 1); e2++; }
        if (e2 & 1) goto next;
        mpz_set_ui(Y, 1);
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
        /* Handle remaining factor (should be a perfect square of LP products) */
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
            fprintf(stderr, "DLP-Opt: factored %dd in %.3fs\n", digits, elapsed());
            mpz_clear(o); return 0;
        }
        mpz_add(tmp, X, Y); mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t o; mpz_init(o); mpz_divexact(o, N, g);
            if (mpz_cmp(g, o) > 0) mpz_swap(g, o);
            gmp_printf("%Zd\n", g);
            fprintf(stderr, "DLP-Opt: factored %dd in %.3fs\n", digits, elapsed());
            mpz_clear(o); return 0;
        }
        next: mpz_clears(X, Y, g, prod, rem, NULL);
    }

    fprintf(stderr, "FAIL: no factor found\n");
    printf("FAIL\n");
    return 1;
}
