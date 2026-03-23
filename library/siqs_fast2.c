/*
 * siqs_fast2.c - Highly Optimized SIQS with Small Prime Variation
 *
 * Key optimizations over siqs_ub.c:
 * 1. Small Prime Variation (SPV): skip sieving with primes < 256, adjust threshold
 * 2. Multi-polynomial batching (BATCH_POLYS=8) with Gray code updates
 * 3. Faster trial division: sieve-informed for large primes, unconditional for small
 * 4. Better threshold computation with SPV adjustment
 * 5. Larger LP bound (lp_mult 80-120) for more combined relations
 *
 * Compile: gcc -O3 -march=native -o siqs_fast2 library/siqs_fast2.c -lgmp -lm
 * Usage: ./siqs_fast2 <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

#define SEED 42
#define BLOCK_SIZE 32768
#define MAX_FB 80000
#define MAX_A_FACTORS 20
#define MAX_RELS 400000
#define MAX_PARTIALS 2000000
#define MAX_DEPS 128
#define BATCH_POLYS 8
#define SPV_CUTOFF 256

static struct timespec g_start;
static double elapsed(void) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec - g_start.tv_sec) + (now.tv_nsec - g_start.tv_nsec) / 1e9;
}

/* ==================== Modular Arithmetic ==================== */
static unsigned int mod_inverse(unsigned int a, unsigned int m) {
    int old_r = (int)a, r = (int)m, old_s = 1, s = 0;
    while (r) { int q = old_r / r; int t = r; r = old_r - q * r; old_r = t; t = s; s = old_s - q * s; old_s = t; }
    if (old_r != 1) return 0;
    return (unsigned int)(((long long)old_s % m + m) % m);
}

static unsigned int sqrt_mod(unsigned int n, unsigned int p) {
    if (n == 0) return 0;
    if (p == 2) return n & 1;
    unsigned long long nn = n % p, r = 1, b = nn, e = (p - 1) / 2, m = p;
    while (e) { if (e & 1) r = r * b % m; b = b * b % m; e >>= 1; }
    if (r != 1) return 0;
    if (p % 4 == 3) { r = 1; b = nn; e = (p + 1) / 4; while (e) { if (e & 1) r = r * b % m; b = b * b % m; e >>= 1; } return (unsigned int)r; }
    unsigned int Q = p - 1, S = 0;
    while (Q % 2 == 0) { Q /= 2; S++; }
    unsigned int z = 2;
    for (;;) { r = 1; b = z; e = (p-1)/2; while (e) { if (e&1) r=r*b%m; b=b*b%m; e>>=1; } if (r==(unsigned long long)(p-1)) break; z++; }
    unsigned long long M_val = S;
    b = z; e = Q; unsigned long long c = 1; while (e) { if (e&1) c=c*b%m; b=b*b%m; e>>=1; }
    b = nn; e = Q; unsigned long long t = 1; while (e) { if (e&1) t=t*b%m; b=b*b%m; e>>=1; }
    b = nn; e = (Q+1)/2; unsigned long long R = 1; while (e) { if (e&1) R=R*b%m; b=b*b%m; e>>=1; }
    for (;;) { if (t==1) return (unsigned int)R; int i=0; unsigned long long tt=t; while (tt!=1){tt=tt*tt%p;i++;} unsigned long long bb=c; for (int j=0;j<(int)M_val-i-1;j++) bb=bb*bb%p; M_val=i; c=bb*bb%p; t=t*c%p; R=R*bb%p; }
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
        for (int i = 0; i < 14; i++) { if (k%ps[i]==0) continue; if (sqrt_mod(mpz_fdiv_ui(kN,ps[i]),ps[i])) s += 2.0*log(ps[i])/(ps[i]-1); }
        if (s > best) { best = s; best_k = k; }
        mpz_clear(kN);
    }
    return best_k;
}

/* ==================== Factor Base ==================== */
typedef struct { unsigned int *prime, *root; unsigned char *logp; int size; } fb_t;

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
        if (nm == 0) { fb->prime[fb->size]=i; fb->root[fb->size]=0; fb->logp[fb->size]=(unsigned char)(log2(i)+0.5); fb->size++; continue; }
        unsigned int r = sqrt_mod((unsigned int)nm, i);
        if (!r) continue;
        fb->prime[fb->size]=i; fb->root[fb->size]=r; fb->logp[fb->size]=(unsigned char)(log2(i)+0.5); fb->size++;
    }
    free(sv);
    return fb;
}

/* ==================== Large Prime Hash ==================== */
#define LP_HASH_BITS 22
#define LP_HASH_SIZE (1 << LP_HASH_BITS)
typedef struct lp_e { unsigned long lp; int idx; struct lp_e *next; } lp_e_t;
typedef struct { lp_e_t **b; lp_e_t *pool; int used, max; } lp_t;
static lp_t *lp_create(int m) { lp_t *t=calloc(1,sizeof(lp_t)); t->b=calloc(LP_HASH_SIZE,sizeof(lp_e_t*)); t->pool=calloc(m,sizeof(lp_e_t)); t->max=m; return t; }
static int lp_find(lp_t *t, unsigned long lp) { unsigned int h=(unsigned int)((lp*0x9E3779B97F4A7C15ULL)>>(64-LP_HASH_BITS)); for (lp_e_t *e=t->b[h]; e; e=e->next) if (e->lp==lp) return e->idx; return -1; }
static void lp_insert(lp_t *t, unsigned long lp, int idx) { if (t->used>=t->max) return; unsigned int h=(unsigned int)((lp*0x9E3779B97F4A7C15ULL)>>(64-LP_HASH_BITS)); lp_e_t *e=&t->pool[t->used++]; e->lp=lp; e->idx=idx; e->next=t->b[h]; t->b[h]=e; }

/* ==================== Relation Storage ==================== */
typedef struct { mpz_t *ax_b, *Qx; unsigned long *lp; int count, alloc; } rels_t;
static rels_t *rels_create(int n) { rels_t *r=malloc(sizeof(rels_t)); r->ax_b=malloc(n*sizeof(mpz_t)); r->Qx=malloc(n*sizeof(mpz_t)); r->lp=calloc(n,sizeof(unsigned long)); for (int i=0;i<n;i++){mpz_init(r->ax_b[i]);mpz_init(r->Qx[i]);} r->count=0; r->alloc=n; return r; }

/* ==================== GF(2) Matrix ==================== */
typedef unsigned long long u64;
typedef struct { u64 **rows; int nr,nc,fbw,idw,wprow; } gf2_t;
static gf2_t *gf2_create(int nr, int nc) {
    gf2_t *m=malloc(sizeof(gf2_t)); m->nr=nr; m->nc=nc;
    m->fbw=(nc+63)/64; m->idw=(nr+63)/64; m->wprow=m->fbw+m->idw;
    m->rows=malloc(nr*sizeof(u64*));
    for (int i=0;i<nr;i++){m->rows[i]=calloc(m->wprow,sizeof(u64)); m->rows[i][m->fbw+i/64]|=(1ULL<<(i%64));}
    return m;
}
static void gf2_set(gf2_t *m, int r, int c) { m->rows[r][c/64]|=(1ULL<<(c%64)); }
static int gf2_solve(gf2_t *m, int ***deps, int **dlen, int max) {
    int piv=0;
    for (int c=0;c<m->nc&&piv<m->nr;c++){
        int pr=-1; for(int r=piv;r<m->nr;r++) if((m->rows[r][c/64]>>(c%64))&1){pr=r;break;}
        if(pr<0) continue;
        if(pr!=piv){u64*t=m->rows[pr];m->rows[pr]=m->rows[piv];m->rows[piv]=t;}
        for(int r=0;r<m->nr;r++){if(r==piv||!((m->rows[r][c/64]>>(c%64))&1))continue; for(int w=0;w<m->wprow;w++) m->rows[r][w]^=m->rows[piv][w];}
        piv++;
    }
    int nd=0; *deps=malloc(max*sizeof(int*)); *dlen=malloc(max*sizeof(int));
    for(int r=piv;r<m->nr&&nd<max;r++){
        int z=1; for(int w=0;w<m->fbw&&z;w++){u64 mask=(w<m->fbw-1)?~0ULL:(m->nc%64==0?~0ULL:(1ULL<<(m->nc%64))-1); if(m->rows[r][w]&mask) z=0;}
        if(!z) continue;
        int *d=malloc(m->nr*sizeof(int)); int dl=0;
        for(int w=0;w<m->idw;w++){u64 bits=m->rows[r][m->fbw+w]; while(bits){int bit=__builtin_ctzll(bits); int idx=w*64+bit; if(idx<m->nr) d[dl++]=idx; bits&=bits-1;}}
        if(dl>0){(*deps)[nd]=d;(*dlen)[nd]=dl;nd++;} else free(d);
    }
    return nd;
}

/* ==================== Parameters ==================== */
typedef struct { int fb_size, nblocks, lp_mult, extra; double thresh; } params_t;
static params_t get_params(int bits) {
    /* Tuned for balanced semiprimes with SPV optimization */
    if (bits <= 100) return (params_t){100, 1, 80, 30, 0.73};
    if (bits <= 110) return (params_t){140, 1, 85, 35, 0.73};
    if (bits <= 120) return (params_t){200, 1, 90, 40, 0.74};
    if (bits <= 130) return (params_t){250, 2, 90, 50, 0.75};
    if (bits <= 140) return (params_t){350, 3, 95, 60, 0.76};
    if (bits <= 150) return (params_t){500, 4, 100, 70, 0.77};
    if (bits <= 160) return (params_t){700, 5, 100, 80, 0.78};
    if (bits <= 170) return (params_t){950, 7, 105, 90, 0.79};
    if (bits <= 180) return (params_t){1300, 10, 110, 100, 0.80};
    if (bits <= 190) return (params_t){1800, 14, 110, 120, 0.81};
    if (bits <= 200) return (params_t){2500, 16, 115, 150, 0.82};
    if (bits <= 210) return (params_t){3500, 22, 115, 180, 0.83};
    if (bits <= 220) return (params_t){5000, 30, 120, 220, 0.84};
    if (bits <= 230) return (params_t){6500, 36, 120, 260, 0.85};
    if (bits <= 240) return (params_t){9000, 46, 120, 320, 0.86};
    if (bits <= 250) return (params_t){12000, 58, 120, 400, 0.87};
    if (bits <= 260) return (params_t){16000, 72, 120, 500, 0.875};
    if (bits <= 270) return (params_t){22000, 90, 120, 600, 0.88};
    if (bits <= 280) return (params_t){30000, 110, 120, 700, 0.885};
    if (bits <= 290) return (params_t){40000, 130, 130, 800, 0.89};
    if (bits <= 300) return (params_t){55000, 150, 130, 900, 0.895};
    return (params_t){75000, 170, 140, 1000, 0.90};
}

/* ==================== Main ==================== */
int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }
    clock_gettime(CLOCK_MONOTONIC, &g_start);

    mpz_t N, kN; mpz_inits(N, kN, NULL);
    mpz_set_str(N, argv[1], 10);
    int digits = (int)mpz_sizeinbase(N, 10);
    int bits = (int)mpz_sizeinbase(N, 2);

    /* Trial division */
    for (unsigned long p = 2; p < 100000; p++) {
        if (mpz_divisible_ui_p(N, p)) { printf("%lu\n", p); return 0; }
    }

    int mult = choose_multiplier(N);
    mpz_mul_ui(kN, N, mult);
    int kN_bits = (int)mpz_sizeinbase(kN, 2);
    params_t P = get_params(kN_bits);

    fb_t *fb = fb_create(kN, P.fb_size);
    int M = BLOCK_SIZE * P.nblocks;
    unsigned long lp_bound = (unsigned long)fb->prime[fb->size-1] * P.lp_mult;
    int target = fb->size + P.extra;

    /* ===== Small Prime Variation (SPV) ===== */
    /* Find cutoff index: first FB prime >= SPV_CUTOFF */
    int spv_cutoff_idx = fb->size; /* default: no SPV if all primes < cutoff */
    for (int i = 1; i < fb->size; i++) {
        if (fb->prime[i] >= SPV_CUTOFF) {
            spv_cutoff_idx = i;
            break;
        }
    }

    /* Compute expected sieve contribution of small primes */
    double spv_adj = 0.0;
    for (int i = 1; i < spv_cutoff_idx; i++) {
        unsigned int p = fb->prime[i];
        if (fb->root[i] == 0) {
            /* Divides kN, hits every sieve position */
            spv_adj += fb->logp[i];
            continue;
        }
        /* Two roots, each hits 1/p of positions */
        spv_adj += fb->logp[i] * 2.0 / p;
    }

    double log2_Qmax = kN_bits / 2.0 + 0.5 + log2((double)M);
    int threshold = (int)(log2_Qmax * P.thresh);
    threshold -= (int)(spv_adj + 0.5);

    /* Safety: threshold must stay positive and reasonable */
    if (threshold < 10) threshold = 10;

    fprintf(stderr, "SIQS-FAST2: %dd (%db), k=%d, FB=%d, M=%d, thresh=%d (spv_adj=%.1f, cutoff_idx=%d), LP=%lu, target=%d, batch=%d\n",
            digits, bits, mult, fb->size, M, threshold, spv_adj, spv_cutoff_idx, lp_bound, target, BATCH_POLYS);

    /* Allocate batch sieve arrays */
    unsigned char *sieves[BATCH_POLYS];
    for (int b = 0; b < BATCH_POLYS; b++) sieves[b] = malloc(BLOCK_SIZE);

    rels_t *full = rels_create(MAX_RELS);
    rels_t *part = rels_create(MAX_PARTIALS);
    lp_t *lpt = lp_create(MAX_PARTIALS);

    /* Polynomial state */
    mpz_t a, b_cur, c_cur, B_vals[MAX_A_FACTORS];
    mpz_init(a); mpz_init(b_cur); mpz_init(c_cur);
    for (int j = 0; j < MAX_A_FACTORS; j++) mpz_init(B_vals[j]);

    /* Solutions for current base polynomial */
    unsigned int *base_soln1 = malloc(fb->size * sizeof(unsigned int));
    unsigned int *base_soln2 = malloc(fb->size * sizeof(unsigned int));

    /* Batch solutions */
    unsigned int *soln1[BATCH_POLYS], *soln2[BATCH_POLYS];
    for (int b = 0; b < BATCH_POLYS; b++) {
        soln1[b] = malloc(fb->size * sizeof(unsigned int));
        soln2[b] = malloc(fb->size * sizeof(unsigned int));
    }

    /* Gray code deltas: delta[j][i] = 2 * a^{-1} * B_j mod prime[i] */
    unsigned int **gray_delta = malloc(MAX_A_FACTORS * sizeof(unsigned int*));
    for (int j = 0; j < MAX_A_FACTORS; j++)
        gray_delta[j] = malloc(fb->size * sizeof(unsigned int));

    /* Batch b and c values */
    mpz_t batch_b[BATCH_POLYS], batch_c[BATCH_POLYS];
    for (int b = 0; b < BATCH_POLYS; b++) { mpz_init(batch_b[b]); mpz_init(batch_c[b]); }

    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, SEED);

    mpz_t ax_b, Qx, residue, tmp;
    mpz_inits(ax_b, Qx, residue, tmp, NULL);

    int total_polys = 0, a_count = 0, combined = 0;
    int a_idx[MAX_A_FACTORS];
    int num_a_factors = 0;

    /* ==================== Main Loop ==================== */
    while (full->count < target) {
        if (elapsed() > 280) {
            fprintf(stderr, "TIMEOUT at %.1fs, rels=%d/%d\n", elapsed(), full->count, target);
            break;
        }

        /* ===== Generate new 'a' ===== */
        {
            mpz_t tgt; mpz_init(tgt);
            mpz_mul_ui(tgt, kN, 2); mpz_sqrt(tgt, tgt); mpz_tdiv_q_ui(tgt, tgt, M);
            double log_tgt = mpz_sizeinbase(tgt, 2) * log(2.0);

            int lo = fb->size / 4, hi = 3 * fb->size / 4;
            if (lo < 2) lo = 2; if (hi <= lo + 3) hi = fb->size - 1;

            double avg = 0; int cnt = 0;
            for (int i = lo; i < hi; i++) { if (fb->root[i]==0) continue; avg += log(fb->prime[i]); cnt++; }
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
        }

        /* ===== Compute B values ===== */
        for (int j = 0; j < num_a_factors; j++) {
            int idx = a_idx[j];
            unsigned int qj = fb->prime[idx], rj = fb->root[idx];
            mpz_t a_q, mod_q, inv; mpz_inits(a_q, mod_q, inv, NULL);
            mpz_divexact_ui(a_q, a, qj); mpz_set_ui(mod_q, qj);
            mpz_invert(inv, a_q, mod_q);
            unsigned long iv = mpz_get_ui(inv);
            mpz_mul_ui(B_vals[j], a_q, (unsigned long)rj * iv % qj);
            mpz_clears(a_q, mod_q, inv, NULL);
        }

        /* ===== Compute Gray code deltas ===== */
        for (int j = 0; j < num_a_factors; j++) {
            for (int i = 0; i < fb->size; i++) {
                unsigned int p = fb->prime[i];
                if (p <= 2 || fb->root[i] == 0) { gray_delta[j][i] = 0; continue; }
                unsigned long am = mpz_fdiv_ui(a, p);
                if (am == 0) { gray_delta[j][i] = 0; continue; }
                unsigned int ai = mod_inverse((unsigned int)am, p);
                unsigned long Bm = mpz_fdiv_ui(B_vals[j], p);
                gray_delta[j][i] = (unsigned int)((2ULL * ai % p * Bm) % p);
            }
        }

        /* ===== Initial b (Gray index 0): b = -sum(B_j) ===== */
        mpz_set_ui(b_cur, 0);
        for (int j = 0; j < num_a_factors; j++) mpz_sub(b_cur, b_cur, B_vals[j]);

        /* Compute initial solutions using GMP */
        for (int i = 0; i < fb->size; i++) {
            unsigned int p = fb->prime[i];
            if (p <= 2 || fb->root[i] == 0) { base_soln1[i] = base_soln2[i] = 0xFFFFFFFF; continue; }
            unsigned long am = mpz_fdiv_ui(a, p);
            if (am == 0) { base_soln1[i] = base_soln2[i] = 0xFFFFFFFF; continue; }
            unsigned int ai = mod_inverse((unsigned int)am, p);
            if (ai == 0) { base_soln1[i] = base_soln2[i] = 0xFFFFFFFF; continue; }
            mpz_t bmod_gmp; mpz_init(bmod_gmp);
            mpz_fdiv_r_ui(bmod_gmp, b_cur, p);
            unsigned long bmod = mpz_get_ui(bmod_gmp);
            mpz_clear(bmod_gmp);
            unsigned int r = fb->root[i];
            base_soln1[i] = (unsigned int)((unsigned long long)ai * ((r + p - bmod) % p) % p);
            base_soln2[i] = (unsigned int)((unsigned long long)ai * ((p - r + p - bmod) % p) % p);
        }

        /* ===== Iterate through polynomials in batches ===== */
        int num_polys = 1 << (num_a_factors - 1);
        int cur_gray = 0;

        for (int batch_start = 0; batch_start < num_polys && full->count < target; batch_start += BATCH_POLYS) {
            int batch = BATCH_POLYS;
            if (batch_start + batch > num_polys) batch = num_polys - batch_start;

            /* Advance base solutions to batch_start using Gray code */
            while (cur_gray < batch_start) {
                int next = cur_gray + 1;
                int g_prev = cur_gray ^ (cur_gray >> 1);
                int g_next = next ^ (next >> 1);
                int changed = __builtin_ctz(g_prev ^ g_next);
                int sign = (g_next >> changed) & 1;

                for (int i = 0; i < fb->size; i++) {
                    if (base_soln1[i] == 0xFFFFFFFF) continue;
                    unsigned int p = fb->prime[i];
                    unsigned int d = gray_delta[changed][i];
                    if (d == 0) continue;
                    if (sign) {
                        base_soln1[i] = (base_soln1[i] >= d) ? base_soln1[i] - d : base_soln1[i] + p - d;
                        base_soln2[i] = (base_soln2[i] >= d) ? base_soln2[i] - d : base_soln2[i] + p - d;
                    } else {
                        base_soln1[i] = base_soln1[i] + d; if (base_soln1[i] >= p) base_soln1[i] -= p;
                        base_soln2[i] = base_soln2[i] + d; if (base_soln2[i] >= p) base_soln2[i] -= p;
                    }
                }
                if (sign) mpz_addmul_ui(b_cur, B_vals[changed], 2);
                else mpz_submul_ui(b_cur, B_vals[changed], 2);
                cur_gray = next;
            }

            /* Poly 0 in batch = base */
            memcpy(soln1[0], base_soln1, fb->size * sizeof(unsigned int));
            memcpy(soln2[0], base_soln2, fb->size * sizeof(unsigned int));
            mpz_set(batch_b[0], b_cur);

            /* Polys 1..batch-1: apply Gray code from base */
            unsigned int *tmp_s1 = malloc(fb->size * sizeof(unsigned int));
            unsigned int *tmp_s2 = malloc(fb->size * sizeof(unsigned int));
            memcpy(tmp_s1, base_soln1, fb->size * sizeof(unsigned int));
            memcpy(tmp_s2, base_soln2, fb->size * sizeof(unsigned int));
            mpz_t tmp_b; mpz_init_set(tmp_b, b_cur);

            for (int bi = 1; bi < batch; bi++) {
                int prev_idx = batch_start + bi - 1;
                int curr_idx = batch_start + bi;
                int g_prev = prev_idx ^ (prev_idx >> 1);
                int g_curr = curr_idx ^ (curr_idx >> 1);
                int changed = __builtin_ctz(g_prev ^ g_curr);
                int sign = (g_curr >> changed) & 1;

                for (int i = 0; i < fb->size; i++) {
                    if (tmp_s1[i] == 0xFFFFFFFF) continue;
                    unsigned int p = fb->prime[i];
                    unsigned int d = gray_delta[changed][i];
                    if (d == 0) continue;
                    if (sign) {
                        tmp_s1[i] = (tmp_s1[i] >= d) ? tmp_s1[i] - d : tmp_s1[i] + p - d;
                        tmp_s2[i] = (tmp_s2[i] >= d) ? tmp_s2[i] - d : tmp_s2[i] + p - d;
                    } else {
                        tmp_s1[i] = tmp_s1[i] + d; if (tmp_s1[i] >= p) tmp_s1[i] -= p;
                        tmp_s2[i] = tmp_s2[i] + d; if (tmp_s2[i] >= p) tmp_s2[i] -= p;
                    }
                }
                if (sign) mpz_addmul_ui(tmp_b, B_vals[changed], 2);
                else mpz_submul_ui(tmp_b, B_vals[changed], 2);

                memcpy(soln1[bi], tmp_s1, fb->size * sizeof(unsigned int));
                memcpy(soln2[bi], tmp_s2, fb->size * sizeof(unsigned int));
                mpz_set(batch_b[bi], tmp_b);
            }
            free(tmp_s1); free(tmp_s2);

            /* Compute c for each batch polynomial */
            for (int bi = 0; bi < batch; bi++) {
                mpz_mul(batch_c[bi], batch_b[bi], batch_b[bi]);
                mpz_sub(batch_c[bi], batch_c[bi], kN);
                mpz_divexact(batch_c[bi], batch_c[bi], a);
            }

            total_polys += batch;

            if (total_polys % 5000 < BATCH_POLYS) {
                double t = elapsed();
                if (t > 280) break;
                if (total_polys % 20000 < BATCH_POLYS)
                    fprintf(stderr, "  polys=%d rels=%d/%d (full=%d+%d) part=%d t=%.1fs\n",
                            total_polys, full->count, target,
                            full->count - combined, combined, part->count, t);
            }

            /* ===== Sieve all blocks ===== */
            for (int block = -P.nblocks; block < P.nblocks; block++) {
                int block_start = block * BLOCK_SIZE;

                /* Initialize all batch sieve arrays */
                for (int bi = 0; bi < batch; bi++)
                    memset(sieves[bi], 0, BLOCK_SIZE);

                /* SPV: Start sieving from spv_cutoff_idx, skip small primes */
                for (int i = spv_cutoff_idx; i < fb->size; i++) {
                    unsigned int p = fb->prime[i];
                    unsigned char lp_val = fb->logp[i];

                    for (int bi = 0; bi < batch; bi++) {
                        if (soln1[bi][i] == 0xFFFFFFFF) continue;
                        unsigned char *sv = sieves[bi];

                        long off1 = ((long)soln1[bi][i] - block_start) % (long)p;
                        if (off1 < 0) off1 += p;
                        for (int j = (int)off1; j < BLOCK_SIZE; j += p)
                            sv[j] += lp_val;

                        if (soln1[bi][i] != soln2[bi][i]) {
                            long off2 = ((long)soln2[bi][i] - block_start) % (long)p;
                            if (off2 < 0) off2 += p;
                            for (int j = (int)off2; j < BLOCK_SIZE; j += p)
                                sv[j] += lp_val;
                        }
                    }
                }

                /* ===== Scan for smooth candidates ===== */
                for (int bi = 0; bi < batch; bi++) {
                    for (int j = 0; j < BLOCK_SIZE; j++) {
                        if (sieves[bi][j] < threshold) continue;
                        long x = (long)(block_start + j);
                        if (x == 0) continue;

                        /* Q(x) = a*x^2 + 2*b*x + c */
                        mpz_set_si(tmp, x);
                        mpz_mul(Qx, a, tmp);
                        mpz_add(Qx, Qx, batch_b[bi]);
                        mpz_add(Qx, Qx, batch_b[bi]);
                        mpz_mul(Qx, Qx, tmp);
                        mpz_add(Qx, Qx, batch_c[bi]);

                        mpz_mul_si(ax_b, a, x);
                        mpz_add(ax_b, ax_b, batch_b[bi]);

                        if (mpz_sgn(Qx) == 0) continue;
                        mpz_abs(residue, Qx);

                        /* Trial divide by 2 first */
                        while (mpz_even_p(residue)) mpz_tdiv_q_2exp(residue, residue, 1);

                        /* Trial divide by small primes unconditionally (cheap, and they
                         * were not sieved so we must check them all) */
                        for (int i = 1; i < spv_cutoff_idx; i++) {
                            unsigned int p = fb->prime[i];
                            while (mpz_divisible_ui_p(residue, p))
                                mpz_divexact_ui(residue, residue, p);
                        }

                        /* Trial divide by sieved primes using root information */
                        for (int i = spv_cutoff_idx; i < fb->size; i++) {
                            unsigned int p = fb->prime[i];
                            if (soln1[bi][i] == 0xFFFFFFFF) continue;
                            long xmod = ((x % (long)p) + p) % p;
                            if (xmod != (long)soln1[bi][i] && xmod != (long)soln2[bi][i]) continue;
                            while (mpz_divisible_ui_p(residue, p))
                                mpz_divexact_ui(residue, residue, p);
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
                            unsigned long lp_val2 = mpz_get_ui(residue);
                            int match = lp_find(lpt, lp_val2);
                            if (match >= 0) {
                                int ri = full->count;
                                if (ri < full->alloc) {
                                    mpz_mul(full->ax_b[ri], ax_b, part->ax_b[match]);
                                    mpz_mod(full->ax_b[ri], full->ax_b[ri], N);
                                    mpz_mul(full->Qx[ri], aQx, part->Qx[match]);
                                    full->lp[ri] = lp_val2;
                                    full->count++;
                                    combined++;
                                }
                            } else {
                                int pi = part->count;
                                if (pi < part->alloc) {
                                    mpz_set(part->ax_b[pi], ax_b);
                                    mpz_set(part->Qx[pi], aQx);
                                    part->lp[pi] = lp_val2;
                                    lp_insert(lpt, lp_val2, pi);
                                    part->count++;
                                }
                            }
                        }
                        mpz_clear(aQx);
                    }
                }
            }

            /* Advance base to end of batch for next iteration */
            if (batch > 1) {
                memcpy(base_soln1, soln1[batch - 1], fb->size * sizeof(unsigned int));
                memcpy(base_soln2, soln2[batch - 1], fb->size * sizeof(unsigned int));
                mpz_set(b_cur, batch_b[batch - 1]);
                cur_gray = batch_start + batch - 1;
            }
            mpz_clear(tmp_b);
        }
    }

    double sieve_time = elapsed();
    fprintf(stderr, "Sieving: %d rels (%d full + %d combined) in %.2fs, %d polys\n",
            full->count, full->count - combined, combined, sieve_time, total_polys);

    if (full->count < fb->size + 1) {
        fprintf(stderr, "FAIL: not enough relations (%d < %d)\n", full->count, fb->size + 1);
        printf("FAIL\n"); return 1;
    }

    /* ===== Linear Algebra ===== */
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
    int ndeps = gf2_solve(mat, &deps, &dlen, MAX_DEPS);
    fprintf(stderr, "LA: %d deps from %dx%d in %.2fs\n", ndeps, nrels, ncols, elapsed());

    /* ===== Square Root ===== */
    for (int d = 0; d < ndeps; d++) {
        mpz_t X, Y, g, prod, rem; mpz_inits(X, Y, g, prod, rem, NULL);
        mpz_set_ui(X, 1);
        for (int k = 0; k < dlen[d]; k++) { mpz_mul(X, X, full->ax_b[deps[d][k]]); mpz_mod(X, X, N); }
        mpz_set_ui(prod, 1);
        for (int k = 0; k < dlen[d]; k++) { mpz_t aq; mpz_init(aq); mpz_abs(aq, full->Qx[deps[d][k]]); mpz_mul(prod, prod, aq); mpz_clear(aq); }

        mpz_set(rem, prod);
        int e2 = 0; while (mpz_even_p(rem)) { mpz_tdiv_q_2exp(rem, rem, 1); e2++; }
        if (e2 & 1) goto next;
        mpz_set_ui(Y, 1);
        if (e2/2 > 0) { mpz_set_ui(tmp, 2); mpz_powm_ui(tmp, tmp, e2/2, N); mpz_mul(Y, Y, tmp); mpz_mod(Y, Y, N); }

        { int valid = 1;
          for (int i = 1; i < fb->size && valid; i++) {
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
            gmp_printf("%Zd\n", g); fprintf(stderr, "Factored in %.3fs\n", elapsed());
            mpz_clear(o); return 0;
        }
        mpz_add(tmp, X, Y); mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
            mpz_t o; mpz_init(o); mpz_divexact(o, N, g);
            if (mpz_cmp(g, o) > 0) mpz_swap(g, o);
            gmp_printf("%Zd\n", g); fprintf(stderr, "Factored in %.3fs\n", elapsed());
            mpz_clear(o); return 0;
        }
        next: mpz_clears(X, Y, g, prod, rem, NULL);
    }

    fprintf(stderr, "FAIL: no factor found\n"); printf("FAIL\n"); return 1;
}
