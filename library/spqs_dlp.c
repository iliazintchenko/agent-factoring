/*
 * SPQS-DLP: Multi-Polynomial Batch Sieve SIQS with Double Large Primes
 *
 * Based on SPQS (best custom implementation). Key enhancement:
 * Double Large Prime (DLP) variation with SQUFOF cofactorization.
 *
 * DLP should improve scaling at larger sizes because:
 * - More partial relations are generated per sieve pass
 * - LP graph becomes denser (birthday paradox) as size increases
 * - Cycle-finding yields "free" full relations from partials
 * - The DLP advantage grows with digit size
 *
 * Compile: gcc -O3 -march=native -o spqs_dlp library/spqs_dlp.c -lgmp -lm
 * Usage: ./spqs_dlp <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <stdint.h>

#define SEED 42
#define SIEVE_BLOCK 49152  /* 48KB = L1d cache on EPYC 9R45 */
#define MAX_FB 80000
#define MAX_A_FACTORS 20
#define MAX_RELS 500000
#define MAX_PARTIALS 2000000
#define BATCH_POLYS 4  /* Number of polynomials to sieve simultaneously */

static struct timespec g_start;
static double elapsed(void) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec - g_start.tv_sec) + (now.tv_nsec - g_start.tv_nsec) / 1e9;
}

/* ==================== Modular Arithmetic ==================== */
static unsigned int mod_inverse(unsigned int a, unsigned int m) {
    int old_r = (int)a, r = (int)m, old_s = 1, s = 0;
    while (r) { int q = old_r / r, t = r; r = old_r - q * r; old_r = t; t = s; s = old_s - q * s; old_s = t; }
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
    while (1) { if (t == 1) return (unsigned int)R; int i = 0; unsigned long long tt = t; while (tt != 1) { tt = (tt * tt) % p; i++; } unsigned long long bb = c; for (int j = 0; j < (int)M_val - i - 1; j++) bb = (bb * bb) % p; M_val = i; c = (bb * bb) % p; t = (t * c) % p; R = (R * bb) % p; }
}

/* ==================== Multiplier ==================== */
static int choose_multiplier(mpz_t N) {
    static const int ks[] = {1,2,3,5,6,7,10,11,13,14,15,17,19,21,23,29,31,37,41,43,0};
    double best = -1e30; int best_k = 1;
    for (int ki = 0; ks[ki]; ki++) {
        int k = ks[ki]; mpz_t kN; mpz_init(kN); mpz_mul_ui(kN, N, k);
        double s = -0.5 * log((double)k);
        unsigned long m8 = mpz_fdiv_ui(kN, 8);
        if (m8 == 1) s += 2*log(2.0); else if (m8 == 5) s += log(2.0); else if (m8 == 3 || m8 == 7) s += 0.5*log(2.0);
        int ps[] = {3,5,7,11,13,17,19,23,29,31,37,41,43,47};
        for (int i = 0; i < 14; i++) { if (k % ps[i] == 0) continue; if (sqrt_mod(mpz_fdiv_ui(kN, ps[i]), ps[i])) s += 2.0*log(ps[i])/(ps[i]-1); }
        if (s > best) { best = s; best_k = k; }
        mpz_clear(kN);
    }
    return best_k;
}

/* ==================== Factor Base ==================== */
typedef struct { unsigned int *prime; unsigned int *root; unsigned char *logp; int size; } fb_t;

static fb_t *fb_create(mpz_t kN, int target) {
    fb_t *fb = malloc(sizeof(fb_t));
    int alloc = target + 20;
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
    return fb;
}

/* ==================== Primality / SQUFOF cofactorization ==================== */
static uint64_t gcd64(uint64_t a, uint64_t b) {
    while (b) { uint64_t t = b; b = a % b; a = t; }
    return a;
}
/* Miller-Rabin with bases sufficient for n < 2^64 */
static uint64_t mulmod64(uint64_t a, uint64_t b, uint64_t m) {
    return (unsigned __int128)a * b % m;
}
static uint64_t powmod64(uint64_t base, uint64_t exp, uint64_t mod) {
    uint64_t result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) result = mulmod64(result, base, mod);
        base = mulmod64(base, base, mod);
        exp >>= 1;
    }
    return result;
}
static int is_prime64(uint64_t n) {
    if (n < 2) return 0;
    if (n < 4) return 1;
    if (n % 2 == 0 || n % 3 == 0) return 0;
    /* Deterministic Miller-Rabin for n < 3,317,044,064,679,887,385,961,981 */
    static const uint64_t witnesses[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    uint64_t d = n - 1; int r = 0;
    while (d % 2 == 0) { d /= 2; r++; }
    for (int i = 0; i < 12; i++) {
        uint64_t a = witnesses[i];
        if (a >= n) continue;
        uint64_t x = powmod64(a, d, n);
        if (x == 1 || x == n - 1) continue;
        int found = 0;
        for (int j = 0; j < r - 1; j++) {
            x = mulmod64(x, x, n);
            if (x == n - 1) { found = 1; break; }
        }
        if (!found) return 0;
    }
    return 1;
}
static uint64_t isqrt64(uint64_t n) {
    if (n == 0) return 0;
    uint64_t x = (uint64_t)sqrt((double)n);
    while (x > 0 && x * x > n) x--;
    while ((x+1)*(x+1) <= n) x++;
    return x;
}
static uint64_t squfof_factor(uint64_t n) {
    if (n <= 1) return n;
    if (n % 2 == 0) return 2;
    if (n % 3 == 0) return 3;
    { uint64_t s = isqrt64(n); if (s*s == n) return s; }
    static const uint32_t mults[] = {1,3,5,7,11,15,21,33,35,55,77,105,165,231,385,1155};
    uint32_t sqN[16], cut[16], Q0[16], P1[16], Q1[16], ns[16];
    uint16_t sv[16][50];
    uint8_t fail[16];
    int nm = 0;
    for (int i = 0; i < 16; i++) {
        __uint128_t sn = (__uint128_t)n * mults[i];
        if (sn >> 62) break;
        uint64_t scaled = (uint64_t)sn;
        uint64_t sq = isqrt64(scaled);
        sqN[i] = (uint32_t)sq; cut[i] = (uint32_t)sqrt(2.0*(double)sq);
        Q0[i] = 1; P1[i] = (uint32_t)sq;
        Q1[i] = (uint32_t)(scaled - sq*sq); ns[i] = 0; fail[i] = 0;
        if (Q1[i] == 0) { uint64_t g = gcd64(n,sq); if (g>1 && g<n) return g; fail[i]=1; }
        nm = i+1;
    }
    if (nm == 0) return 0;
    uint32_t ti = 0;
    while (ti < 100000) {
        for (int mi = nm-1; mi >= 0; mi--) {
            if (fail[mi]) continue;
            uint32_t sn = sqN[mi], ct = cut[mi], mult = 2*mults[mi], ccut = ct*mult;
            uint32_t q0=Q0[mi], p1=P1[mi], q1=Q1[mi], nsv=ns[mi];
            uint32_t p0=0, sqq=0; int found=0;
            for (uint32_t it = 0; it < 300; it++) {
                uint32_t tmp=sn+p1-q1, q=1;
                if (tmp>=q1) q+=tmp/q1;
                p0=q*q1-p1; q0=q0+(p1-p0)*q;
                if (q1<ccut) { uint32_t t=q1/(uint32_t)gcd64(q1,mult); if (t<ct) { if (nsv>=50){fail[mi]=1;break;} sv[mi][nsv++]=(uint16_t)t; } }
                uint32_t bits=0; tmp=q0; while(tmp&&!(tmp&1)){bits++;tmp>>=1;}
                if (!(bits&1)&&((tmp&7)==1)) { sqq=(uint32_t)sqrt((double)q0); while(sqq>0&&sqq*sqq>q0)sqq--; while((sqq+1)*(sqq+1)<=q0)sqq++;
                    if (sqq*sqq==q0) { uint32_t j; for(j=0;j<nsv;j++) if(sv[mi][j]==sqq) break; if(j==nsv){found=1;break;} } }
                tmp=sn+p0-q0; q=1; if(tmp>=q0)q+=tmp/q0; p1=q*q0-p0; q1=q1+(p0-p1)*q;
                if(q0<ccut){uint32_t t=q0/(uint32_t)gcd64(q0,mult);if(t<ct){if(nsv>=50){fail[mi]=1;break;}sv[mi][nsv++]=(uint16_t)t;}}
                ti++;
            }
            if (fail[mi]) continue;
            if (!found) { Q0[mi]=q0;P1[mi]=p1;Q1[mi]=q1;ns[mi]=nsv; continue; }
            if (sqq<=1) { fail[mi]=1; continue; }
            uint64_t scaledn=(uint64_t)n*mults[mi];
            q0=sqq; p1=p0+sqq*((sn-p0)/sqq); q1=(uint32_t)((scaledn-(uint64_t)p1*p1)/q0);
            for(int inv=0;inv<10000;inv++) {
                uint32_t tmp=sn+p1-q1,q=1; if(tmp>=q1)q+=tmp/q1;
                p0=q*q1-p1; q0=q0+(p1-p0)*q;
                if(p0==p1){q0=q1;break;}
                tmp=sn+p0-q0;q=1;if(tmp>=q0)q+=tmp/q0;
                p1=q*q0-p0;q1=q1+(p0-p1)*q;
                if(p0==p1)break;
            }
            q0=q0/(uint32_t)gcd64(q0,mult);
            if (q0>1&&q0<n&&n%q0==0) return q0;
        }
    }
    return 0;
}

/* ==================== Bucket Sieve for Large Primes ==================== */
typedef struct { uint16_t pos; uint8_t logp; } bucket_hit_t;

static bucket_hit_t **g_buckets[BATCH_POLYS]; /* g_buckets[poly][block] */
static int **g_bucket_cnt;                     /* g_bucket_cnt[poly][block] -- shared alloc tracking */
static int g_bucket_alloc_per_block;
static int g_num_blocks_total;
static int g_fb_bucket_start; /* first FB index with prime > SIEVE_BLOCK */

static void bucket_init(int nblocks_total, int fb_size) {
    g_num_blocks_total = nblocks_total;
    /* Estimate max bucket entries per block per poly */
    g_bucket_alloc_per_block = 512;
    for (int bi = 0; bi < BATCH_POLYS; bi++) {
        g_buckets[bi] = malloc(nblocks_total * sizeof(bucket_hit_t*));
        for (int b = 0; b < nblocks_total; b++)
            g_buckets[bi][b] = malloc(g_bucket_alloc_per_block * sizeof(bucket_hit_t));
    }
    g_bucket_cnt = malloc(BATCH_POLYS * sizeof(int*));
    for (int bi = 0; bi < BATCH_POLYS; bi++)
        g_bucket_cnt[bi] = calloc(nblocks_total, sizeof(int));
}

static inline void bucket_add_hit(int poly, int block, uint16_t pos, uint8_t logp) {
    int c = g_bucket_cnt[poly][block];
    if (c >= g_bucket_alloc_per_block) {
        g_bucket_alloc_per_block *= 2;
        for (int bi = 0; bi < BATCH_POLYS; bi++)
            for (int b = 0; b < g_num_blocks_total; b++)
                g_buckets[bi][b] = realloc(g_buckets[bi][b],
                    g_bucket_alloc_per_block * sizeof(bucket_hit_t));
    }
    g_buckets[poly][block][c].pos = pos;
    g_buckets[poly][block][c].logp = logp;
    g_bucket_cnt[poly][block]++;
}

static void bucket_clear_all(void) {
    for (int bi = 0; bi < BATCH_POLYS; bi++)
        memset(g_bucket_cnt[bi], 0, g_num_blocks_total * sizeof(int));
}

/* ==================== Large Prime Hash ==================== */
#define LP_HASH_BITS 20
#define LP_HASH_SIZE (1 << LP_HASH_BITS)
typedef struct lp_e { unsigned long lp; int idx; struct lp_e *next; } lp_e_t;
typedef struct { lp_e_t **b; lp_e_t *pool; int used, max; } lp_t;
static lp_t *lp_create(int m) { lp_t *t = calloc(1, sizeof(lp_t)); t->b = calloc(LP_HASH_SIZE, sizeof(lp_e_t*)); t->pool = calloc(m, sizeof(lp_e_t)); t->max = m; return t; }
static int lp_find(lp_t *t, unsigned long lp) { unsigned int h = (unsigned int)((lp * 0x9E3779B97F4A7C15ULL) >> (64-LP_HASH_BITS)); for (lp_e_t *e = t->b[h]; e; e = e->next) if (e->lp == lp) return e->idx; return -1; }
static void lp_insert(lp_t *t, unsigned long lp, int idx) { if (t->used >= t->max) return; unsigned int h = (unsigned int)((lp * 0x9E3779B97F4A7C15ULL) >> (64-LP_HASH_BITS)); lp_e_t *e = &t->pool[t->used++]; e->lp = lp; e->idx = idx; e->next = t->b[h]; t->b[h] = e; }

/* ==================== Relation Storage ==================== */
typedef struct { mpz_t *ax_b, *Qx; unsigned long *lp; int count, alloc; } rels_t;
static rels_t *rels_create(int n) { rels_t *r = malloc(sizeof(rels_t)); r->ax_b = malloc(n*sizeof(mpz_t)); r->Qx = malloc(n*sizeof(mpz_t)); r->lp = calloc(n, sizeof(unsigned long)); for (int i = 0; i < n; i++) { mpz_init(r->ax_b[i]); mpz_init(r->Qx[i]); } r->count = 0; r->alloc = n; return r; }

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
        int pr = -1; for (int r = piv; r < m->nr; r++) if ((m->rows[r][c/64] >> (c%64)) & 1) { pr = r; break; }
        if (pr < 0) continue;
        if (pr != piv) { u64 *t = m->rows[pr]; m->rows[pr] = m->rows[piv]; m->rows[piv] = t; }
        for (int r = 0; r < m->nr; r++) { if (r == piv) continue; if ((m->rows[r][c/64] >> (c%64)) & 1) for (int w = 0; w < m->wprow; w++) m->rows[r][w] ^= m->rows[piv][w]; }
        piv++;
    }
    int nd = 0; *deps = malloc(max * sizeof(int*)); *dlen = malloc(max * sizeof(int));
    for (int r = piv; r < m->nr && nd < max; r++) {
        int z = 1; for (int w = 0; w < m->fbw && z; w++) { u64 mask = (w < m->fbw-1) ? ~0ULL : (m->nc%64==0 ? ~0ULL : (1ULL << (m->nc%64))-1); if (m->rows[r][w] & mask) z = 0; }
        if (!z) continue;
        int *d = malloc(m->nr * sizeof(int)); int dl = 0;
        for (int w = 0; w < m->idw; w++) { u64 bits = m->rows[r][m->fbw+w]; while (bits) { int bit = __builtin_ctzll(bits); int idx = w*64+bit; if (idx < m->nr) d[dl++] = idx; bits &= bits-1; } }
        if (dl > 0) { (*deps)[nd] = d; (*dlen)[nd] = dl; nd++; } else free(d);
    }
    return nd;
}

/* ==================== Parameters ==================== */
typedef struct { int fb_size, nblocks, lp_mult, extra; double thresh; int dlp_mult; } params_t;
static params_t get_params(int bits) {
    /* Balanced nblocks: keep enough blocks for consistency at small sizes,
     * reduce blocks at larger sizes for faster poly switching + smaller Q(x).
     * dlp_mult: 0 = no DLP. DLP starts at FB >= 400 (140+ bits). */
    if (bits <= 100) return (params_t){100, 1, 30, 40, 0.73, 0};
    if (bits <= 110) return (params_t){150, 1, 30, 40, 0.74, 0};
    if (bits <= 120) return (params_t){200, 2, 35, 50, 0.76, 0};
    if (bits <= 130) return (params_t){300, 3, 40, 50, 0.78, 0};
    if (bits <= 140) return (params_t){400, 4, 40, 60, 0.79, 50};
    if (bits <= 150) return (params_t){600, 6, 50, 60, 0.80, 80};
    if (bits <= 160) return (params_t){900, 6, 50, 80, 0.81, 100};
    if (bits <= 170) return (params_t){1200, 4, 60, 80, 0.82, 150};
    if (bits <= 180) return (params_t){1800, 3, 60, 80, 0.83, 200};
    if (bits <= 190) return (params_t){2500, 2, 60, 100, 0.84, 300};
    if (bits <= 200) return (params_t){3500, 2, 70, 100, 0.85, 400};
    if (bits <= 210) return (params_t){5000, 2, 70, 120, 0.86, 500};
    if (bits <= 220) return (params_t){7000, 3, 80, 120, 0.87, 600};
    if (bits <= 230) return (params_t){9000, 3, 80, 150, 0.875, 800};
    if (bits <= 240) return (params_t){12000, 4, 80, 150, 0.88, 1000};
    if (bits <= 250) return (params_t){16000, 5, 90, 200, 0.885, 1200};
    if (bits <= 260) return (params_t){22000, 6, 90, 200, 0.89, 1500};
    if (bits <= 270) return (params_t){30000, 8, 100, 250, 0.895, 2000};
    if (bits <= 280) return (params_t){40000, 10, 100, 300, 0.90, 2500};
    if (bits <= 290) return (params_t){55000, 12, 110, 350, 0.905, 3000};
    return (params_t){75000, 160, 120, 400, 0.91, 4000};
}

/* ==================== SIQS with Multi-Polynomial Sieving ==================== */

/*
 * Novel: For each 'a' value, instead of sieving one b-value at a time,
 * we sieve MULTIPLE b-values simultaneously over the same sieve block.
 *
 * This works because switching b doesn't change the FB prime divisibility
 * positions much (they shift by ±delta for each b-change). By keeping
 * multiple sieve arrays and doing the FB prime sieve loop once, we amortize
 * the prime iteration overhead.
 *
 * With k polynomials per block, the overhead is:
 * - k * SIEVE_BLOCK bytes of memory (k * 32KB)
 * - k sieve updates per FB prime position (instead of 1)
 * - BUT: the OUTER loop over FB primes is done once
 *
 * The benefit: k times more smooth candidates per sieve pass.
 * For small FB (where the outer loop cost dominates), this is a win.
 *
 * We batch BATCH_POLYS polynomials per sieve pass.
 */

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
            gmp_printf("%d\n", p); return 0;
        }
    }

    int mult = choose_multiplier(N);
    mpz_mul_ui(kN, N, mult);
    int kN_bits = (int)mpz_sizeinbase(kN, 2);
    params_t P = get_params(kN_bits);

    fb_t *fb = fb_create(kN, P.fb_size);
    int M = SIEVE_BLOCK * P.nblocks;
    unsigned long lp_bound = (unsigned long)fb->prime[fb->size-1] * P.lp_mult;
    /* DLP bound = LP_bound^2 (both factors must be < LP_bound) */
    uint64_t dlp_bound = (P.dlp_mult > 0) ? (uint64_t)lp_bound * (uint64_t)lp_bound : 0;
    /* Cap at 2^62 for SQUFOF */
    if (dlp_bound > (1ULL << 62)) dlp_bound = (1ULL << 62);
    int target = fb->size + P.extra;

    double log2_Qmax = kN_bits / 2.0 + 0.5 + log2(M);
    int threshold = (int)(log2_Qmax * P.thresh);
    threshold -= 3;
    /* Lower threshold when DLP is enabled to capture more candidates.
     * Only for larger sizes where DLP has enough candidates. */
    if (P.dlp_mult > 0 && fb->size >= 2000) {
        int dlp_fudge = (int)(log2((double)lp_bound) * 0.3);
        if (dlp_fudge > 7) dlp_fudge = 7;
        threshold -= dlp_fudge;
        if (threshold < 20) threshold = 20;
    }

    /* Determine bucket sieve threshold */
    g_fb_bucket_start = fb->size;
    for (int i = 0; i < fb->size; i++) {
        if (fb->prime[i] > SIEVE_BLOCK) { g_fb_bucket_start = i; break; }
    }

    /* Initialize bucket sieve */
    int nblocks_total = 2 * P.nblocks;
    bucket_init(nblocks_total, fb->size);

    fprintf(stderr, "SPQS-DLP: %dd (%db), k=%d, FB=%d (bucket@%d), M=%d, thresh=%d, LP=%lu, DLP=%lu, target=%d\n",
            digits, bits, mult, fb->size, g_fb_bucket_start, M, threshold, lp_bound, (unsigned long)dlp_bound, target);

    /* Allocate BATCH_POLYS sieve arrays */
    unsigned char *sieves[BATCH_POLYS];
    for (int b = 0; b < BATCH_POLYS; b++)
        sieves[b] = malloc(SIEVE_BLOCK);

    rels_t *full = rels_create(MAX_RELS);
    rels_t *part = rels_create(MAX_PARTIALS);
    lp_t *lpt = lp_create(MAX_PARTIALS);

    /* Polynomial state - we maintain BATCH_POLYS polynomials at once */
    mpz_t a, bs[BATCH_POLYS], cs[BATCH_POLYS], B_vals[MAX_A_FACTORS];
    mpz_init(a);
    for (int b = 0; b < BATCH_POLYS; b++) { mpz_init(bs[b]); mpz_init(cs[b]); }
    for (int j = 0; j < MAX_A_FACTORS; j++) mpz_init(B_vals[j]);

    unsigned int *soln1[BATCH_POLYS], *soln2[BATCH_POLYS];
    for (int b = 0; b < BATCH_POLYS; b++) {
        soln1[b] = malloc(fb->size * sizeof(unsigned int));
        soln2[b] = malloc(fb->size * sizeof(unsigned int));
    }

    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, SEED);

    mpz_t ax_b, Qx, residue, tmp;
    mpz_inits(ax_b, Qx, residue, tmp, NULL);

    int total_polys = 0, a_count = 0, combined = 0, combined_dlp = 0;
    int dlp_candidates = 0, dlp_composite = 0, dlp_split_ok = 0, dlp_both_ok = 0;
    int a_idx[MAX_A_FACTORS];
    int num_a_factors = 0;
    unsigned int *ainv_data = malloc(MAX_A_FACTORS * fb->size * sizeof(unsigned int));

    /* Main loop */
    while (full->count < target) {
        if (total_polys > 0 && total_polys % (2000/BATCH_POLYS) == 0) {
            double t = elapsed();
            if (t > 280) { fprintf(stderr, "TIMEOUT at %.1fs\n", t); break; }
            if (total_polys % (8000/BATCH_POLYS) == 0)
                fprintf(stderr, "  p=%d r=%d/%d (full=%d+slp%d+dlp%d) part=%d t=%.1fs\n",
                        total_polys * BATCH_POLYS, full->count, target,
                        full->count - combined - combined_dlp, combined, combined_dlp,
                        part->count, t);
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
            }

            memcpy(a_idx, best, s * sizeof(int));
            mpz_set_ui(a, 1);
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
                }
            }
        }

        /* Generate BATCH_POLYS b-values using Gray code positions */
        int num_b = 1 << (num_a_factors - 1);

        for (int b_start = 0; b_start < num_b && full->count < target; b_start += BATCH_POLYS) {
            int batch = BATCH_POLYS;
            if (b_start + batch > num_b) batch = num_b - b_start;

            /* Compute b-values for this batch */
            /* Start from the base b-value and apply Gray code changes */
            for (int bi = 0; bi < batch; bi++) {
                int b_idx = b_start + bi;
                int gray = b_idx ^ (b_idx >> 1);

                /* b = sum of (±1) * B_j based on gray code bits */
                mpz_set_ui(bs[bi], 0);
                for (int j = 0; j < num_a_factors; j++) {
                    if (gray & (1 << j))
                        mpz_add(bs[bi], bs[bi], B_vals[j]);
                    else
                        mpz_sub(bs[bi], bs[bi], B_vals[j]);
                }

                /* Verify b^2 ≡ kN (mod a) */
                mpz_mul(tmp, bs[bi], bs[bi]); mpz_sub(tmp, tmp, kN); mpz_mod(tmp, tmp, a);
                if (mpz_sgn(tmp) != 0) {
                    /* Try negating */
                    mpz_neg(bs[bi], bs[bi]);
                    mpz_mul(tmp, bs[bi], bs[bi]); mpz_sub(tmp, tmp, kN); mpz_mod(tmp, tmp, a);
                    if (mpz_sgn(tmp) != 0) { batch = bi; break; }
                }

                /* c = (b^2 - kN) / a */
                mpz_mul(cs[bi], bs[bi], bs[bi]); mpz_sub(cs[bi], cs[bi], kN);
                mpz_divexact(cs[bi], cs[bi], a);

                /* Compute sieve solutions */
                for (int i = 0; i < fb->size; i++) {
                    unsigned int p = fb->prime[i];
                    unsigned long am = mpz_fdiv_ui(a, p);
                    if (am == 0 || fb->root[i] == 0) { soln1[bi][i] = soln2[bi][i] = 0xFFFFFFFF; continue; }
                    unsigned int ai = mod_inverse((unsigned int)am, p);
                    if (ai == 0) { soln1[bi][i] = soln2[bi][i] = 0xFFFFFFFF; continue; }
                    unsigned long bm = mpz_fdiv_ui(bs[bi], p);
                    unsigned int r = fb->root[i];
                    soln1[bi][i] = (unsigned int)((unsigned long)ai * ((r + p - bm) % p) % p);
                    soln2[bi][i] = (unsigned int)((unsigned long)ai * ((p - r + p - bm) % p) % p);
                }
            }

            if (batch == 0) continue;
            total_polys++;

            /* Fill buckets for large primes across all blocks */
            if (g_fb_bucket_start < fb->size) {
                bucket_clear_all();
                int total_sieve = 2 * P.nblocks * SIEVE_BLOCK;
                for (int bi = 0; bi < batch; bi++) {
                    for (int i = g_fb_bucket_start; i < fb->size; i++) {
                        if (soln1[bi][i] == 0xFFFFFFFF) continue;
                        unsigned int p = fb->prime[i];
                        unsigned char lp = fb->logp[i];
                        /* Root 1 */
                        long s1 = ((long)soln1[bi][i] + P.nblocks * SIEVE_BLOCK);
                        for (long pos = s1 % (long)p; pos < total_sieve; pos += p) {
                            int blk = (int)(pos / SIEVE_BLOCK);
                            bucket_add_hit(bi, blk, (uint16_t)(pos % SIEVE_BLOCK), lp);
                        }
                        /* Root 2 */
                        if (soln1[bi][i] != soln2[bi][i]) {
                            long s2 = ((long)soln2[bi][i] + P.nblocks * SIEVE_BLOCK);
                            for (long pos = s2 % (long)p; pos < total_sieve; pos += p) {
                                int blk = (int)(pos / SIEVE_BLOCK);
                                bucket_add_hit(bi, blk, (uint16_t)(pos % SIEVE_BLOCK), lp);
                            }
                        }
                    }
                }
            }

            /* Sieve all BATCH polynomials over each block */
            for (int block = -P.nblocks; block < P.nblocks; block++) {
                int block_start = block * SIEVE_BLOCK;
                int block_idx = block + P.nblocks; /* 0-based block index */

                /* Initialize all sieve arrays */
                for (int bi = 0; bi < batch; bi++)
                    memset(sieves[bi], 0, SIEVE_BLOCK);

                /* Sieve with small/medium FB primes only (p <= SIEVE_BLOCK) */
                int fb_end = (g_fb_bucket_start < fb->size) ? g_fb_bucket_start : fb->size;
                for (int i = 1; i < fb_end; i++) {
                    unsigned int p = fb->prime[i];
                    if (p < 5) continue;
                    unsigned char lp = fb->logp[i];

                    for (int bi = 0; bi < batch; bi++) {
                        if (soln1[bi][i] == 0xFFFFFFFF) continue;

                        long off1 = ((long)soln1[bi][i] - block_start) % (long)p;
                        if (off1 < 0) off1 += p;
                        /* Unrolled inner loop for small primes */
                        if (p < 64) {
                            int j = (int)off1;
                            for (; j + 3*(int)p < SIEVE_BLOCK; j += 4*(int)p) {
                                sieves[bi][j] += lp;
                                sieves[bi][j + p] += lp;
                                sieves[bi][j + 2*p] += lp;
                                sieves[bi][j + 3*p] += lp;
                            }
                            for (; j < SIEVE_BLOCK; j += p) sieves[bi][j] += lp;
                        } else {
                            for (int j = (int)off1; j < SIEVE_BLOCK; j += p)
                                sieves[bi][j] += lp;
                        }

                        if (soln1[bi][i] != soln2[bi][i]) {
                            long off2 = ((long)soln2[bi][i] - block_start) % (long)p;
                            if (off2 < 0) off2 += p;
                            if (p < 64) {
                                int j = (int)off2;
                                for (; j + 3*(int)p < SIEVE_BLOCK; j += 4*(int)p) {
                                    sieves[bi][j] += lp;
                                    sieves[bi][j + p] += lp;
                                    sieves[bi][j + 2*p] += lp;
                                    sieves[bi][j + 3*p] += lp;
                                }
                                for (; j < SIEVE_BLOCK; j += p) sieves[bi][j] += lp;
                            } else {
                                for (int j = (int)off2; j < SIEVE_BLOCK; j += p)
                                    sieves[bi][j] += lp;
                            }
                        }
                    }
                }

                /* Apply bucket sieve hits for large primes */
                if (g_fb_bucket_start < fb->size) {
                    for (int bi = 0; bi < batch; bi++) {
                        int cnt = g_bucket_cnt[bi][block_idx];
                        bucket_hit_t *hits = g_buckets[bi][block_idx];
                        for (int j = 0; j < cnt; j++)
                            sieves[bi][hits[j].pos] += hits[j].logp;
                    }
                }

                /* Scan for smooth candidates in all polynomial sieves */
                for (int bi = 0; bi < batch; bi++) {
                    for (int j = 0; j < SIEVE_BLOCK; j++) {
                        if (sieves[bi][j] < threshold) continue;
                        long x = (long)(block_start + j);
                        if (x == 0) continue;

                        /* Compute Q(x) for polynomial bi */
                        mpz_set_si(tmp, x);
                        mpz_mul(Qx, a, tmp);
                        mpz_add(Qx, Qx, bs[bi]); mpz_add(Qx, Qx, bs[bi]);
                        mpz_mul(Qx, Qx, tmp);
                        mpz_add(Qx, Qx, cs[bi]);

                        mpz_mul_si(ax_b, a, x);
                        mpz_add(ax_b, ax_b, bs[bi]);

                        if (mpz_sgn(Qx) == 0) continue;
                        mpz_abs(residue, Qx);

                        /* Trial divide */
                        while (mpz_even_p(residue)) mpz_tdiv_q_2exp(residue, residue, 1);
                        for (int i = 1; i < fb->size; i++) {
                            unsigned int p = fb->prime[i];
                            if (soln1[bi][i] == 0xFFFFFFFF) continue;
                            long xmod = ((x % (long)p) + p) % p;
                            if (xmod != (long)soln1[bi][i] && xmod != (long)soln2[bi][i]) continue;
                            if (mpz_divisible_ui_p(residue, p))
                                do { mpz_divexact_ui(residue, residue, p); } while (mpz_divisible_ui_p(residue, p));
                        }
                        for (int i = 0; i < fb->size && fb->prime[i] < 5; i++) {
                            unsigned int p = fb->prime[i]; if (p <= 2) continue;
                            while (mpz_divisible_ui_p(residue, p)) mpz_divexact_ui(residue, residue, p);
                        }

                        mpz_t aQx; mpz_init(aQx); mpz_mul(aQx, Qx, a);

                        if (mpz_cmp_ui(residue, 1) == 0) {
                            int ri = full->count;
                            if (ri < full->alloc) {
                                mpz_set(full->ax_b[ri], ax_b);
                                mpz_set(full->Qx[ri], aQx);
                                full->lp[ri] = 0;
                                full->count++;
                            }
                        } else if (mpz_fits_ulong_p(residue)) {
                            uint64_t cofactor = mpz_get_ui(residue);
                            if (cofactor <= (uint64_t)lp_bound) {
                                /* Single Large Prime */
                                int match = lp_find(lpt, cofactor);
                                if (match >= 0) {
                                    int ri = full->count;
                                    if (ri < full->alloc) {
                                        mpz_mul(full->ax_b[ri], ax_b, part->ax_b[match]);
                                        mpz_mod(full->ax_b[ri], full->ax_b[ri], N);
                                        mpz_mul(full->Qx[ri], aQx, part->Qx[match]);
                                        full->lp[ri] = cofactor;
                                        full->count++;
                                        combined++;
                                    }
                                } else {
                                    int pi = part->count;
                                    if (pi < part->alloc) {
                                        mpz_set(part->ax_b[pi], ax_b);
                                        mpz_set(part->Qx[pi], aQx);
                                        part->lp[pi] = cofactor;
                                        lp_insert(lpt, cofactor, pi);
                                        part->count++;
                                    }
                                }
                            } else if (dlp_bound > 0 && cofactor <= dlp_bound) {
                                dlp_candidates++;
                                if (is_prime64(cofactor)) goto skip_dlp;
                                dlp_composite++;
                                /* DLP candidate: composite cofactor, try SQUFOF to split */
                                uint64_t f1 = squfof_factor(cofactor);
                                if (f1 > 1 && f1 < cofactor) {
                                    dlp_split_ok++;
                                    uint64_t f2 = cofactor / f1;
                                    if (f1 > f2) { uint64_t t2=f1; f1=f2; f2=t2; }
                                    if (f1 <= (uint64_t)lp_bound && f2 <= (uint64_t)lp_bound) {
                                        dlp_both_ok++;
                                        /* Both primes within LP bound */
                                        int m1 = lp_find(lpt, f1);
                                        int m2 = lp_find(lpt, f2);
                                        if (m1 >= 0 && m2 >= 0) {
                                            /* Triple combination -> full relation */
                                            int ri = full->count;
                                            if (ri < full->alloc) {
                                                mpz_mul(full->ax_b[ri], ax_b, part->ax_b[m1]);
                                                mpz_mul(full->ax_b[ri], full->ax_b[ri], part->ax_b[m2]);
                                                mpz_mod(full->ax_b[ri], full->ax_b[ri], N);
                                                mpz_mul(full->Qx[ri], aQx, part->Qx[m1]);
                                                mpz_mul(full->Qx[ri], full->Qx[ri], part->Qx[m2]);
                                                full->lp[ri] = cofactor;
                                                full->count++;
                                                combined_dlp++;
                                            }
                                        } else if (m1 >= 0) {
                                            /* f1 matched, store partial with f2 */
                                            int pi = part->count;
                                            if (pi < part->alloc) {
                                                mpz_mul(part->ax_b[pi], ax_b, part->ax_b[m1]);
                                                mpz_mod(part->ax_b[pi], part->ax_b[pi], N);
                                                mpz_mul(part->Qx[pi], aQx, part->Qx[m1]);
                                                part->lp[pi] = f2;
                                                lp_insert(lpt, f2, pi);
                                                part->count++;
                                            }
                                        } else if (m2 >= 0) {
                                            int pi = part->count;
                                            if (pi < part->alloc) {
                                                mpz_mul(part->ax_b[pi], ax_b, part->ax_b[m2]);
                                                mpz_mod(part->ax_b[pi], part->ax_b[pi], N);
                                                mpz_mul(part->Qx[pi], aQx, part->Qx[m2]);
                                                part->lp[pi] = f1;
                                                lp_insert(lpt, f1, pi);
                                                part->count++;
                                            }
                                        } else {
                                            /* Store both LPs as partials */
                                            int pi = part->count;
                                            if (pi + 1 < part->alloc) {
                                                mpz_set(part->ax_b[pi], ax_b);
                                                mpz_set(part->Qx[pi], aQx);
                                                part->lp[pi] = f1;
                                                lp_insert(lpt, f1, pi);
                                                part->count++;
                                                pi = part->count;
                                                mpz_set(part->ax_b[pi], ax_b);
                                                mpz_set(part->Qx[pi], aQx);
                                                part->lp[pi] = f2;
                                                lp_insert(lpt, f2, pi);
                                                part->count++;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        skip_dlp:
                        mpz_clear(aQx);
                    }
                }
            }
        }
    }

    double sieve_time = elapsed();
    fprintf(stderr, "Sieving: %d rels (%d full + %d SLP + %d DLP) in %.2fs\n",
            full->count, full->count - combined - combined_dlp, combined, combined_dlp, sieve_time);
    fprintf(stderr, "DLP stats: candidates=%d composite=%d split=%d both_ok=%d\n",
            dlp_candidates, dlp_composite, dlp_split_ok, dlp_both_ok);

    if (full->count < fb->size + 1) {
        fprintf(stderr, "FAIL: not enough relations\n"); printf("FAIL\n"); return 1;
    }

    /* Linear algebra */
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
    int ndeps = gf2_solve(mat, &deps, &dlen, 256);
    fprintf(stderr, "LA: %d deps from %dx%d\n", ndeps, nrels, ncols);

    /* Square root */
    for (int d = 0; d < ndeps; d++) {
        mpz_t X, Y, g, prod, rem; mpz_inits(X, Y, g, prod, rem, NULL);
        mpz_set_ui(X, 1);
        for (int k = 0; k < dlen[d]; k++) { mpz_mul(X, X, full->ax_b[deps[d][k]]); mpz_mod(X, X, N); }

        mpz_set_ui(Y, 1); mpz_set_ui(prod, 1);
        for (int k = 0; k < dlen[d]; k++) { mpz_t aq; mpz_init(aq); mpz_abs(aq, full->Qx[deps[d][k]]); mpz_mul(prod, prod, aq); mpz_clear(aq); }

        mpz_set(rem, prod);
        int e2 = 0; while (mpz_even_p(rem)) { mpz_tdiv_q_2exp(rem, rem, 1); e2++; }
        if (e2 & 1) goto next;
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
            if (mpz_perfect_square_p(rem)) { mpz_sqrt(tmp, rem); mpz_mod(tmp, tmp, N); mpz_mul(Y, Y, tmp); mpz_mod(Y, Y, N); }
            else goto next;
        }

        mpz_sub(tmp, X, Y); mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t o; mpz_init(o); mpz_divexact(o, N, g);
            if (mpz_cmp(g, o) > 0) mpz_swap(g, o);
            gmp_printf("%Zd\n", g);
            fprintf(stderr, "SPQS-DLP: factored %dd in %.3fs\n", digits, elapsed());
            mpz_clear(o); return 0;
        }
        mpz_add(tmp, X, Y); mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t o; mpz_init(o); mpz_divexact(o, N, g);
            if (mpz_cmp(g, o) > 0) mpz_swap(g, o);
            gmp_printf("%Zd\n", g);
            fprintf(stderr, "SPQS-DLP: factored %dd in %.3fs\n", digits, elapsed());
            mpz_clear(o); return 0;
        }
        next: mpz_clears(X, Y, g, prod, rem, NULL);
    }

    fprintf(stderr, "SPQS-DLP: FAILED\n"); printf("FAIL\n"); return 1;
}
