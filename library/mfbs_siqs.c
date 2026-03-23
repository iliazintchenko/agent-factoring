/*
 * MFBS-SIQS: Multi-Factor-Base Sieve SIQS
 *
 * Novel approach: Use TWO factor bases:
 * - Small FB (primes <= B1): used for sieving (determines sieve cost)
 * - Extended FB (primes B1 < p <= B2): used ONLY for trial division on candidates
 *
 * The sieve only iterates over B1 primes (cheap). Candidates that pass the
 * sieve threshold are trial-divided by the extended FB too. This gives
 * "free" additional prime factors without sieve overhead.
 *
 * Each relation records which extended primes divide Q(x). These get columns
 * in the GF(2) matrix. The matrix is larger (B2 columns total) but relations
 * are richer (more prime factors identified per Q(x)).
 *
 * Key insight: at 70d+, the sieve cost grows as O(B1 * M * npolys).
 * By using B1 = B_optimal / 2, sieve cost halves. The extended primes
 * [B1, B2] are tested cheaply via trial division only on candidates
 * (maybe 1 per 1000 sieve positions), so their cost is negligible.
 *
 * The tradeoff: smaller B1 → fewer sieve hits → higher threshold needed →
 * fewer candidates. But extended TD recovers primes that the sieve missed,
 * making more candidates viable.
 *
 * Compile: gcc -O3 -march=native -o mfbs_siqs library/mfbs_siqs.c -lgmp -lm
 * Usage: ./mfbs_siqs <N>
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
#define MAX_RELS 400000
#define MAX_PARTIALS 1500000
#define BATCH_POLYS 4
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

/* ==================== Factor Base (includes extended) ==================== */
typedef struct {
    unsigned int *prime;
    unsigned int *root;
    unsigned char *logp;
    int size;          /* total FB size (sieve + extended) */
    int sieve_size;    /* primes used for sieving */
    int large_start;   /* index where p >= BLOCK_SIZE (within sieve primes) */
} fb_t;

static fb_t *fb_create(mpz_t kN, int sieve_target, int extended_target) {
    fb_t *fb = calloc(1, sizeof(fb_t));
    int total = sieve_target + extended_target;
    int alloc = total + 100;
    fb->prime = malloc(alloc * sizeof(unsigned int));
    fb->root = malloc(alloc * sizeof(unsigned int));
    fb->logp = malloc(alloc * sizeof(unsigned char));
    fb->prime[0] = 2; fb->root[0] = 1; fb->logp[0] = 1; fb->size = 1;
    int bound = total * 30 + 100000;
    char *sv = calloc(bound + 1, 1);
    for (int i = 2; (long)i*i <= bound; i++) if (!sv[i]) for (int j = i*i; j <= bound; j += i) sv[j] = 1;
    for (int i = 3; i <= bound && fb->size < total; i += 2) {
        if (sv[i]) continue;
        unsigned long nm = mpz_fdiv_ui(kN, i);
        if (nm == 0) { fb->prime[fb->size] = i; fb->root[fb->size] = 0; fb->logp[fb->size] = (unsigned char)(log2(i)+0.5); fb->size++; continue; }
        unsigned int r = sqrt_mod((unsigned int)nm, i);
        if (!r) continue;
        fb->prime[fb->size] = i; fb->root[fb->size] = r; fb->logp[fb->size] = (unsigned char)(log2(i)+0.5); fb->size++;
    }
    free(sv);
    fb->sieve_size = (fb->size < sieve_target) ? fb->size : sieve_target;
    fb->large_start = fb->sieve_size;
    for (int i = 0; i < fb->sieve_size; i++) {
        if (fb->prime[i] >= BLOCK_SIZE) { fb->large_start = i; break; }
    }
    return fb;
}

/* ==================== LP Hash ==================== */
#define LP_HASH_BITS 21
#define LP_HASH_SIZE (1 << LP_HASH_BITS)
typedef struct lp_e { unsigned long lp; int idx; struct lp_e *next; } lp_e_t;
typedef struct { lp_e_t **b; lp_e_t *pool; int used, max; } lp_t;
static lp_t *lp_create(int m) { lp_t *t = calloc(1, sizeof(lp_t)); t->b = calloc(LP_HASH_SIZE, sizeof(lp_e_t*)); t->pool = calloc(m, sizeof(lp_e_t)); t->max = m; return t; }
static int lp_find(lp_t *t, unsigned long lp) { uint32_t h = (uint32_t)((lp * 0x9E3779B97F4A7C15ULL) >> (64-LP_HASH_BITS)); for (lp_e_t *e = t->b[h]; e; e = e->next) if (e->lp == lp) return e->idx; return -1; }
static void lp_insert(lp_t *t, unsigned long lp, int idx) { if (t->used >= t->max) return; uint32_t h = (uint32_t)((lp * 0x9E3779B97F4A7C15ULL) >> (64-LP_HASH_BITS)); lp_e_t *e = &t->pool[t->used++]; e->lp = lp; e->idx = idx; e->next = t->b[h]; t->b[h] = e; }

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
typedef struct {
    int sieve_fb;     /* primes used for sieving */
    int extended_fb;  /* additional primes for trial division only */
    int nblocks;
    int lp_mult;
    int extra;
    double thresh;
} params_t;

static params_t get_params(int bits) {
    /* sieve_fb is SMALLER than usual, extended_fb adds cheap extra primes */
    if (bits <= 100) return (params_t){100, 50, 1, 40, 40, 0.73};
    if (bits <= 120) return (params_t){200, 100, 2, 45, 50, 0.76};
    if (bits <= 140) return (params_t){400, 200, 4, 50, 60, 0.78};
    if (bits <= 160) return (params_t){800, 400, 8, 55, 70, 0.80};
    if (bits <= 180) return (params_t){1500, 750, 14, 65, 80, 0.82};
    if (bits <= 200) return (params_t){2500, 1500, 22, 70, 100, 0.84};
    if (bits <= 210) return (params_t){3500, 2000, 28, 75, 120, 0.85};
    if (bits <= 220) return (params_t){5000, 3000, 36, 80, 140, 0.86};
    if (bits <= 232) return (params_t){7000, 4000, 44, 80, 160, 0.87};
    if (bits <= 248) return (params_t){12000, 6000, 56, 90, 200, 0.88};
    if (bits <= 265) return (params_t){20000, 10000, 70, 90, 250, 0.89};
    if (bits <= 281) return (params_t){35000, 15000, 90, 100, 300, 0.895};
    return (params_t){50000, 25000, 110, 120, 350, 0.90};
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
    for (int p = 2; p < 10000; p++) {
        if (mpz_divisible_ui_p(N, p)) {
            mpz_t c; mpz_init(c); mpz_divexact_ui(c, N, p);
            gmp_printf("%Zd\n", c); return 0;
        }
    }

    int mult = choose_multiplier(N);
    mpz_mul_ui(kN, N, mult);
    int kN_bits = (int)mpz_sizeinbase(kN, 2);
    params_t P = get_params(kN_bits);

    fb_t *fb = fb_create(kN, P.sieve_fb, P.extended_fb);
    int M = BLOCK_SIZE * P.nblocks;
    int total_blocks = 2 * P.nblocks;
    unsigned long lp_bound = (unsigned long)fb->prime[fb->size-1] * P.lp_mult;
    /* Target: need rels > total FB size (sieve + extended) + extra */
    int target = fb->size + P.extra;

    double log2_Qmax = kN_bits / 2.0 + 0.5 + log2(M);
    /* Lower threshold: sieve only has sieve_fb primes, not the full FB */
    double sieve_logp_total = 0;
    for (int i = 0; i < fb->sieve_size; i++) sieve_logp_total += fb->logp[i] * 2.0 / fb->prime[i];
    int threshold = (int)(log2_Qmax * P.thresh);
    /* Adjust: sieve accumulates less because extended primes aren't sieved */
    threshold -= 5;  /* accept more candidates for extended TD */

    fprintf(stderr, "MFBS-SIQS: %dd (%db), k=%d, sieve_FB=%d, ext_FB=%d (total=%d), M=%d, thresh=%d, LP=%lu, target=%d\n",
            digits, bits, mult, fb->sieve_size, fb->size - fb->sieve_size,
            fb->size, M, threshold, lp_bound, target);

    /* Sieve arrays for batch polys */
    unsigned char *sieves[BATCH_POLYS];
    for (int b = 0; b < BATCH_POLYS; b++) sieves[b] = aligned_alloc(64, BLOCK_SIZE);

    unsigned int *soln1[BATCH_POLYS], *soln2[BATCH_POLYS];
    for (int b = 0; b < BATCH_POLYS; b++) {
        soln1[b] = malloc(fb->size * sizeof(unsigned int));
        soln2[b] = malloc(fb->size * sizeof(unsigned int));
    }

    rels_t *full = rels_create(MAX_RELS);
    rels_t *part = rels_create(MAX_PARTIALS);
    lp_t *lpt = lp_create(MAX_PARTIALS);

    mpz_t a, B_vals[MAX_A_FACTORS], bs[BATCH_POLYS], cs[BATCH_POLYS];
    mpz_init(a);
    for (int b = 0; b < BATCH_POLYS; b++) { mpz_init(bs[b]); mpz_init(cs[b]); }
    for (int j = 0; j < MAX_A_FACTORS; j++) mpz_init(B_vals[j]);
    unsigned int *ainv_data = malloc(MAX_A_FACTORS * fb->sieve_size * sizeof(unsigned int));

    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, SEED);

    mpz_t ax_b, Qx, residue, tmp;
    mpz_inits(ax_b, Qx, residue, tmp, NULL);

    int total_polys = 0, a_count = 0, combined = 0;
    int a_idx[MAX_A_FACTORS];
    int num_a_factors = 0;

    while (full->count < target) {
        if (elapsed() > 280) { fprintf(stderr, "TIMEOUT\n"); break; }

        /* Generate new 'a' from SIEVE FB primes */
        {
            mpz_t tgt; mpz_init(tgt);
            mpz_mul_ui(tgt, kN, 2); mpz_sqrt(tgt, tgt); mpz_tdiv_q_ui(tgt, tgt, M);
            double log_tgt = mpz_sizeinbase(tgt, 2) * log(2.0);

            int lo = fb->sieve_size / 3, hi = 2 * fb->sieve_size / 3;
            if (lo < 2) lo = 2; if (hi <= lo + 3) hi = fb->sieve_size - 1;

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
                if (ratio < 1.5) break;
            }
            memcpy(a_idx, best, s * sizeof(int));
            mpz_set_ui(a, 1);
            for (int i = 0; i < s; i++) mpz_mul_ui(a, a, fb->prime[a_idx[i]]);
            mpz_clear(tgt);
            a_count++;

            for (int j = 0; j < s; j++) {
                int idx2 = a_idx[j];
                unsigned int qj = fb->prime[idx2], rj = fb->root[idx2];
                mpz_t a_q, mod_q, inv; mpz_inits(a_q, mod_q, inv, NULL);
                mpz_divexact_ui(a_q, a, qj); mpz_set_ui(mod_q, qj);
                mpz_invert(inv, a_q, mod_q);
                unsigned long iv = mpz_get_ui(inv);
                mpz_mul_ui(B_vals[j], a_q, (unsigned long)rj * iv % qj);
                mpz_clears(a_q, mod_q, inv, NULL);
            }

            /* ainv only for sieve FB */
            for (int j = 0; j < s; j++) {
                for (int i = 0; i < fb->sieve_size; i++) {
                    unsigned int p = fb->prime[i];
                    unsigned long am = mpz_fdiv_ui(a, p);
                    if (am == 0 || fb->root[i] == 0) { ainv_data[j*fb->sieve_size+i] = 0; continue; }
                    unsigned int ai = mod_inverse((unsigned int)am, p);
                    unsigned long Bm = mpz_fdiv_ui(B_vals[j], p);
                    ainv_data[j*fb->sieve_size+i] = (unsigned int)((2ULL * ai % p * Bm) % p);
                }
            }
        }

        int num_b = 1 << (num_a_factors - 1);

        for (int b_start = 0; b_start < num_b && full->count < target; b_start += BATCH_POLYS) {
            int batch = BATCH_POLYS;
            if (b_start + batch > num_b) batch = num_b - b_start;

            /* Compute b-values for this batch (Gray code) */
            for (int bi = 0; bi < batch; bi++) {
                int b_idx = b_start + bi;
                int gray = b_idx ^ (b_idx >> 1);
                mpz_set_ui(bs[bi], 0);
                for (int j = 0; j < num_a_factors; j++) {
                    if (gray & (1 << j)) mpz_add(bs[bi], bs[bi], B_vals[j]);
                    else mpz_sub(bs[bi], bs[bi], B_vals[j]);
                }
                mpz_mul(tmp, bs[bi], bs[bi]); mpz_sub(tmp, tmp, kN); mpz_mod(tmp, tmp, a);
                if (mpz_sgn(tmp) != 0) {
                    mpz_neg(bs[bi], bs[bi]);
                    mpz_mul(tmp, bs[bi], bs[bi]); mpz_sub(tmp, tmp, kN); mpz_mod(tmp, tmp, a);
                    if (mpz_sgn(tmp) != 0) { batch = bi; break; }
                }
                mpz_mul(cs[bi], bs[bi], bs[bi]); mpz_sub(cs[bi], cs[bi], kN);
                mpz_divexact(cs[bi], cs[bi], a);

                /* Solutions for SIEVE primes */
                for (int i = 0; i < fb->sieve_size; i++) {
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
                /* Solutions for EXTENDED primes (only used for root-aware TD) */
                for (int i = fb->sieve_size; i < fb->size; i++) {
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

            if (total_polys % (500/BATCH_POLYS) == 0) {
                double t = elapsed();
                if (t > 280) break;
                if (total_polys % (2000/BATCH_POLYS) == 0)
                    fprintf(stderr, "  p=%d r=%d/%d (full=%d+%d) t=%.1fs\n",
                            total_polys * BATCH_POLYS, full->count, target,
                            full->count - combined, combined, t);
            }

            /* Sieve blocks - only use sieve_fb primes */
            for (int block = -P.nblocks; block < P.nblocks; block++) {
                int block_start = block * BLOCK_SIZE;

                for (int bi = 0; bi < batch; bi++) memset(sieves[bi], 0, BLOCK_SIZE);

                /* Sieve with SIEVE FB primes (amortized across batch polys) */
                for (int i = 1; i < fb->sieve_size; i++) {
                    unsigned int p = fb->prime[i];
                    if (p < 5) continue;
                    unsigned char lp = fb->logp[i];

                    for (int bi = 0; bi < batch; bi++) {
                        if (soln1[bi][i] == 0xFFFFFFFF) continue;
                        long off1 = ((long)soln1[bi][i] - block_start) % (long)p;
                        if (off1 < 0) off1 += p;
                        for (int j = (int)off1; j < BLOCK_SIZE; j += p) sieves[bi][j] += lp;
                        if (soln1[bi][i] != soln2[bi][i]) {
                            long off2 = ((long)soln2[bi][i] - block_start) % (long)p;
                            if (off2 < 0) off2 += p;
                            for (int j = (int)off2; j < BLOCK_SIZE; j += p) sieves[bi][j] += lp;
                        }
                    }
                }

                /* Scan and trial divide (including EXTENDED primes) */
                for (int bi = 0; bi < batch; bi++) {
                    for (int j = 0; j < BLOCK_SIZE; j++) {
                        if (sieves[bi][j] < threshold) continue;
                        long x = (long)(block_start + j);
                        if (x == 0) continue;

                        mpz_set_si(tmp, x);
                        mpz_mul(Qx, a, tmp); mpz_add(Qx, Qx, bs[bi]); mpz_add(Qx, Qx, bs[bi]);
                        mpz_mul(Qx, Qx, tmp); mpz_add(Qx, Qx, cs[bi]);
                        mpz_mul_si(ax_b, a, x); mpz_add(ax_b, ax_b, bs[bi]);

                        if (mpz_sgn(Qx) == 0) continue;
                        mpz_abs(residue, Qx);

                        /* Trial divide by ALL FB primes (sieve + extended) */
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
    }

    double sieve_time = elapsed();
    fprintf(stderr, "Sieving: %d rels (%d full + %d combined) in %.2fs, %d polys\n",
            full->count, full->count - combined, combined, sieve_time, total_polys * BATCH_POLYS);

    if (full->count < fb->size + 1) {
        fprintf(stderr, "FAIL: not enough relations (%d/%d)\n", full->count, fb->size + 1);
        printf("FAIL\n"); return 1;
    }

    /* Linear Algebra - matrix has columns for ALL FB primes (sieve + extended) */
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
    fprintf(stderr, "LA: %d deps from %dx%d in %.2fs\n", ndeps, nrels, ncols, elapsed());

    /* Square Root */
    for (int d = 0; d < ndeps; d++) {
        mpz_t X, Y, g, prod, rem;
        mpz_inits(X, Y, g, prod, rem, NULL);
        mpz_set_ui(X, 1);
        for (int k = 0; k < dlen[d]; k++) { mpz_mul(X, X, full->ax_b[deps[d][k]]); mpz_mod(X, X, N); }
        mpz_set_ui(prod, 1);
        for (int k = 0; k < dlen[d]; k++) {
            mpz_t aq; mpz_init(aq); mpz_abs(aq, full->Qx[deps[d][k]]);
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
            fprintf(stderr, "MFBS-SIQS: factored %dd in %.3fs\n", digits, elapsed());
            mpz_clear(o); return 0;
        }
        mpz_add(tmp, X, Y); mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t o; mpz_init(o); mpz_divexact(o, N, g);
            if (mpz_cmp(g, o) > 0) mpz_swap(g, o);
            gmp_printf("%Zd\n", g);
            fprintf(stderr, "MFBS-SIQS: factored %dd in %.3fs\n", digits, elapsed());
            mpz_clear(o); return 0;
        }
        next: mpz_clears(X, Y, g, prod, rem, NULL);
    }

    fprintf(stderr, "FAIL: no factor found\n");
    printf("FAIL\n");
    return 1;
}
