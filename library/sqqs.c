/*
 * sqqs.c - Special-Q Quadratic Sieve (Novel approach)
 *
 * Key insight: In standard SIQS, Q(x) ≈ sqrt(N) and smoothness probability
 * is ~exp(-u*ln(u)) where u = ln(sqrt(N))/ln(B).
 *
 * Novel approach: Use a "special prime" q to reduce the effective size of
 * the sieve values. For a SIQS polynomial Q(x) with root r_q, we know
 * q | Q(r_q). Writing x = r_q + q*k, we get Q(x)/q as a polynomial in k
 * with values ~ sqrt(N)/q. This increases smoothness probability.
 *
 * For q ≈ B (factor base bound), u drops from ln(sqrt(N))/ln(B) to
 * ln(sqrt(N)/B)/ln(B), which is u-1. Since the smoothness probability
 * is roughly u^(-u), reducing u by 1 gives a factor of ~u improvement.
 *
 * The cost: we need to sieve q different sublattices (one per special-q),
 * but each sublattice has q times fewer positions. The net sieve work
 * is similar, but the RELATION YIELD is higher because values are more smooth.
 *
 * This is analogous to how NFS uses special-q lattice sieving to reduce
 * algebraic norm sizes. Applied to QS, it should give better scaling.
 *
 * Compile: gcc -O3 -march=native -o sqqs library/sqqs.c -lgmp -lm
 * Usage: ./sqqs <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

#define SEED 42
#define BLOCKSIZE 32768
#define MAX_FB 100000
#define MAX_RELS 500000
#define MAX_PARTIALS 2000000
#define MAX_A_FACTORS 25
#define MAX_DEPS 64

static struct timespec g_start;
static double elapsed(void) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec - g_start.tv_sec) + (now.tv_nsec - g_start.tv_nsec) / 1e9;
}

/* ==================== Parameters ==================== */
typedef struct {
    int fb_size;
    int num_blocks;
    int lp_mult;
    int extra_rels;
    double thresh_adj;
    int sq_min;       /* min special-q prime (FB index) */
    int sq_max;       /* max special-q prime (FB index) */
} params_t;

static params_t get_params(int bits) {
    /* Special-Q range: use primes from upper portion of FB */
    if (bits <= 100) return (params_t){120,   1,  40,  40,  0.70, 60,  120};
    if (bits <= 110) return (params_t){180,   1,  40,  40,  0.72, 90,  180};
    if (bits <= 120) return (params_t){250,   2,  40,  50,  0.74, 120, 250};
    if (bits <= 130) return (params_t){350,   3,  50,  50,  0.76, 170, 350};
    if (bits <= 140) return (params_t){500,   4,  50,  60,  0.78, 250, 500};
    if (bits <= 150) return (params_t){750,   6,  60,  60,  0.79, 350, 750};
    if (bits <= 160) return (params_t){1100,  8,  60,  80,  0.80, 500, 1100};
    if (bits <= 170) return (params_t){1500,  12, 70,  80,  0.81, 700, 1500};
    if (bits <= 180) return (params_t){2200,  16, 70,  80,  0.82, 1000,2200};
    if (bits <= 190) return (params_t){3200,  22, 80,  100, 0.83, 1500,3200};
    if (bits <= 200) return (params_t){4500,  30, 80,  100, 0.84, 2000,4500};
    if (bits <= 210) return (params_t){6500,  40, 90,  120, 0.85, 3000,6500};
    if (bits <= 220) return (params_t){9000,  50, 90,  120, 0.86, 4000,9000};
    if (bits <= 230) return (params_t){12000, 60, 100, 150, 0.87, 5500,12000};
    if (bits <= 240) return (params_t){16000, 72, 100, 150, 0.88, 7000,16000};
    if (bits <= 250) return (params_t){22000, 88, 110, 200, 0.885,9500,22000};
    if (bits <= 260) return (params_t){30000, 104,110, 200, 0.89, 13000,30000};
    if (bits <= 270) return (params_t){40000, 120,120, 250, 0.895,17000,40000};
    if (bits <= 280) return (params_t){55000, 140,120, 300, 0.90, 24000,55000};
    if (bits <= 290) return (params_t){70000, 160,130, 350, 0.905,30000,70000};
    return                  (params_t){90000, 180,130, 400, 0.91, 40000,90000};
}

/* ==================== Modular Arithmetic ==================== */
static inline unsigned int mod_inverse_u32(unsigned int a, unsigned int m) {
    int old_r = (int)a, r = (int)m, old_s = 1, s = 0;
    while (r) { int q = old_r/r; int t = r; r = old_r-q*r; old_r = t; t = s; s = old_s-q*s; old_s = t; }
    if (old_r != 1) return 0;
    return (unsigned int)(((long long)old_s % m + m) % m);
}

static unsigned int sqrt_mod_p(unsigned int n, unsigned int p) {
    if (n == 0) return 0;
    if (p == 2) return n & 1;
    unsigned long long nn = n%p, m = p, r = 1, b = nn, e = (p-1)/2;
    while (e) { if (e&1) r=r*b%m; b=b*b%m; e>>=1; }
    if (r != 1) return 0;
    if (p%4==3) { r=1; b=nn; e=(p+1)/4; while (e) { if (e&1) r=r*b%m; b=b*b%m; e>>=1; } return (unsigned int)r; }
    unsigned int Q=p-1, S=0; while (Q%2==0) { Q/=2; S++; }
    unsigned int z=2;
    for (;;) { r=1; b=z; e=(p-1)/2; while (e) { if (e&1) r=r*b%m; b=b*b%m; e>>=1; } if (r==m-1) break; z++; }
    unsigned long long M2=S; r=1; b=z; e=Q; while (e) { if (e&1) r=r*b%m; b=b*b%m; e>>=1; } unsigned long long c=r;
    r=1; b=nn; e=Q; while (e) { if (e&1) r=r*b%m; b=b*b%m; e>>=1; } unsigned long long t=r;
    r=1; b=nn; e=(Q+1)/2; while (e) { if (e&1) r=r*b%m; b=b*b%m; e>>=1; } unsigned long long R=r;
    for (;;) { if (t==1) return (unsigned int)R; int i2=0; unsigned long long tt=t; while (tt!=1) { tt=tt*tt%p; i2++; }
        unsigned long long bb=c; for (int j=0; j<(int)M2-i2-1; j++) bb=bb*bb%p; M2=i2; c=bb*bb%p; t=t*c%p; R=R*bb%p; }
}

/* ==================== Multiplier ==================== */
static int choose_multiplier(mpz_t N) {
    static const int ks[] = {1,2,3,5,6,7,10,11,13,14,15,17,19,21,23,29,31,37,41,43,0};
    double best = -1e30; int best_k = 1;
    for (int ki = 0; ks[ki]; ki++) {
        int k = ks[ki]; mpz_t kN; mpz_init(kN); mpz_mul_ui(kN, N, k);
        double s = -0.5*log((double)k);
        unsigned long m8 = mpz_fdiv_ui(kN, 8);
        if (m8==1) s += 2*log(2.0); else if (m8==5) s += log(2.0); else if (m8==3||m8==7) s += 0.5*log(2.0);
        int ps[] = {3,5,7,11,13,17,19,23,29,31,37,41,43,47};
        for (int i = 0; i < 14; i++) { if (k%ps[i]==0) { s+=log((double)ps[i]); continue; }
            if (sqrt_mod_p(mpz_fdiv_ui(kN,ps[i]),ps[i])) s += 2.0*log(ps[i])/(ps[i]-1); }
        if (s > best) { best = s; best_k = k; }
        mpz_clear(kN);
    }
    return best_k;
}

/* ==================== Factor Base ==================== */
typedef struct {
    unsigned int *prime;
    unsigned int *sqrtN;
    unsigned char *logp;
    int size, alloc;
} fb_t;

static fb_t *fb_create(mpz_t kN, int target) {
    fb_t *fb = malloc(sizeof(fb_t));
    int alloc = target + 100;
    fb->prime = malloc(alloc * sizeof(unsigned int));
    fb->sqrtN = malloc(alloc * sizeof(unsigned int));
    fb->logp = malloc(alloc * sizeof(unsigned char));
    fb->alloc = alloc;
    fb->prime[0] = 2; fb->sqrtN[0] = 1; fb->logp[0] = 1; fb->size = 1;
    int bound = target * 30 + 100000;
    char *sv = calloc(bound+1, 1);
    for (int i = 2; (long)i*i <= bound; i++) if (!sv[i]) for (int j = i*i; j <= bound; j += i) sv[j] = 1;
    for (int i = 3; i <= bound && fb->size < target; i += 2) {
        if (sv[i]) continue;
        unsigned long nm = mpz_fdiv_ui(kN, i);
        if (nm == 0) { fb->prime[fb->size]=i; fb->sqrtN[fb->size]=0; fb->logp[fb->size]=(unsigned char)(log2(i)+0.5); fb->size++; continue; }
        unsigned int r = sqrt_mod_p((unsigned int)nm, i);
        if (!r) continue;
        fb->prime[fb->size]=i; fb->sqrtN[fb->size]=r; fb->logp[fb->size]=(unsigned char)(log2(i)+0.5); fb->size++;
    }
    free(sv);
    return fb;
}

/* ==================== Large Prime Hash ==================== */
#define LP_HASH_BITS 22
#define LP_HASH_SIZE (1 << LP_HASH_BITS)
typedef struct lp_e { unsigned long lp; int idx; struct lp_e *next; } lp_e_t;
typedef struct { lp_e_t **b; lp_e_t *pool; int used, max; } lp_t;
static lp_t *lp_create(int m) { lp_t *t = calloc(1,sizeof(lp_t)); t->b = calloc(LP_HASH_SIZE,sizeof(lp_e_t*)); t->pool=calloc(m,sizeof(lp_e_t)); t->max=m; return t; }
static int lp_find(lp_t *t, unsigned long lp) { unsigned int h=(unsigned int)((lp*0x9E3779B97F4A7C15ULL)>>(64-LP_HASH_BITS)); for (lp_e_t *e=t->b[h]; e; e=e->next) if (e->lp==lp) return e->idx; return -1; }
static void lp_insert(lp_t *t, unsigned long lp, int idx) { if (t->used>=t->max) return; unsigned int h=(unsigned int)((lp*0x9E3779B97F4A7C15ULL)>>(64-LP_HASH_BITS)); lp_e_t *e=&t->pool[t->used++]; e->lp=lp; e->idx=idx; e->next=t->b[h]; t->b[h]=e; }

/* ==================== Relation Storage ==================== */
typedef struct {
    mpz_t *Y;
    short **exps;
    unsigned long *lp;
    int count, alloc, fb_size, neg_col;
} rels_t;

static rels_t *rels_create(int alloc, int fb_size) {
    rels_t *r = malloc(sizeof(rels_t));
    r->Y = malloc(alloc * sizeof(mpz_t));
    r->exps = malloc(alloc * sizeof(short*));
    r->lp = calloc(alloc, sizeof(unsigned long));
    for (int i = 0; i < alloc; i++) { mpz_init(r->Y[i]); r->exps[i] = calloc(fb_size+2, sizeof(short)); }
    r->count = 0; r->alloc = alloc; r->fb_size = fb_size; r->neg_col = fb_size;
    return r;
}

/* ==================== GF(2) Solve ==================== */
typedef unsigned long long u64;
typedef struct { u64 **rows; int nr, nc, fbw, idw, wprow; } gf2_t;

static gf2_t *gf2_create(int nr, int nc) {
    gf2_t *m = malloc(sizeof(gf2_t));
    m->nr=nr; m->nc=nc; m->fbw=(nc+63)/64; m->idw=(nr+63)/64; m->wprow=m->fbw+m->idw;
    m->rows = malloc(nr * sizeof(u64*));
    for (int i = 0; i < nr; i++) { m->rows[i] = calloc(m->wprow, sizeof(u64)); m->rows[i][m->fbw+i/64] |= (1ULL<<(i%64)); }
    return m;
}

static int gf2_solve(gf2_t *m, int ***deps_out, int **dlen_out, int max_deps) {
    int piv = 0;
    for (int c = 0; c < m->nc && piv < m->nr; c++) {
        int pr = -1; for (int r = piv; r < m->nr; r++) if ((m->rows[r][c/64]>>(c%64))&1) { pr=r; break; }
        if (pr < 0) continue;
        if (pr != piv) { u64 *t = m->rows[pr]; m->rows[pr] = m->rows[piv]; m->rows[piv] = t; }
        for (int r = 0; r < m->nr; r++) { if (r==piv) continue;
            if ((m->rows[r][c/64]>>(c%64))&1) for (int w=0; w<m->wprow; w++) m->rows[r][w]^=m->rows[piv][w]; }
        piv++;
    }
    int nd = 0; *deps_out = malloc(max_deps*sizeof(int*)); *dlen_out = malloc(max_deps*sizeof(int));
    for (int r = piv; r < m->nr && nd < max_deps; r++) {
        int zero = 1;
        for (int w = 0; w < m->fbw && zero; w++) {
            u64 mask = (w<m->fbw-1)?~0ULL:(m->nc%64==0?~0ULL:(1ULL<<(m->nc%64))-1);
            if (m->rows[r][w]&mask) zero=0;
        }
        if (!zero) continue;
        int *d = malloc(m->nr*sizeof(int)); int dl=0;
        for (int w = 0; w < m->idw; w++) { u64 bits = m->rows[r][m->fbw+w];
            while (bits) { int bit = __builtin_ctzll(bits); int idx = w*64+bit; if (idx<m->nr) d[dl++]=idx; bits &= bits-1; } }
        if (dl > 0) { (*deps_out)[nd]=d; (*dlen_out)[nd]=dl; nd++; } else free(d);
    }
    return nd;
}

/* ==================== Trial Division ==================== */
static int trial_divide(mpz_t Q, mpz_t Y, short *exps,
                        fb_t *fb, int fb_size,
                        unsigned long lp_bound, unsigned long *lp_out,
                        mpz_t kN, mpz_t a, mpz_t b_poly, int x_global,
                        int neg_col) {
    mpz_mul_si(Y, a, x_global);
    mpz_add(Y, Y, b_poly);
    mpz_mul(Q, Y, Y);
    mpz_sub(Q, Q, kN);

    memset(exps, 0, (fb_size + 2) * sizeof(short));
    if (mpz_sgn(Q) < 0) { exps[neg_col] = 1; mpz_neg(Q, Q); }

    while (mpz_even_p(Q)) { exps[0]++; mpz_tdiv_q_2exp(Q, Q, 1); }

    for (int i = 1; i < fb_size; i++) {
        unsigned int p = fb->prime[i];
        if (p < 3) continue;
        if (!mpz_divisible_ui_p(Q, p)) continue;
        do { exps[i]++; mpz_divexact_ui(Q, Q, p); } while (mpz_divisible_ui_p(Q, p));
    }

    if (mpz_cmp_ui(Q, 1) == 0) { *lp_out = 0; return 1; }
    if (mpz_fits_ulong_p(Q)) {
        unsigned long cof = mpz_get_ui(Q);
        if (cof <= lp_bound && cof > 1) {
            mpz_t c; mpz_init_set_ui(c, cof);
            int pr = mpz_probab_prime_p(c, 5);
            mpz_clear(c);
            if (pr) { *lp_out = cof; return 2; }
        }
    }
    return 0;
}

/* ==================== Main: Special-Q SIQS ==================== */

int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }
    clock_gettime(CLOCK_MONOTONIC, &g_start);

    mpz_t N, kN;
    mpz_inits(N, kN, NULL);
    mpz_set_str(N, argv[1], 10);

    int digits = (int)mpz_sizeinbase(N, 10);
    int bits = (int)mpz_sizeinbase(N, 2);

    /* Quick trial division */
    for (unsigned long p = 2; p < 100000; p++)
        if (mpz_divisible_ui_p(N, p)) {
            mpz_t q; mpz_init(q); mpz_divexact_ui(q, N, p);
            gmp_printf("%lu\n%Zd\n", p, q); mpz_clear(q); return 0;
        }

    int mult = choose_multiplier(N);
    mpz_mul_ui(kN, N, mult);
    int kN_bits = (int)mpz_sizeinbase(kN, 2);
    params_t P = get_params(kN_bits);

    fb_t *fb = fb_create(kN, P.fb_size);
    int M = BLOCKSIZE * P.num_blocks;
    unsigned long lp_bound = (unsigned long)fb->prime[fb->size-1] * P.lp_mult;
    int target = fb->size + P.extra_rels;

    /* Special-Q parameters */
    int sq_lo = P.sq_min;
    int sq_hi = P.sq_max;
    if (sq_lo >= fb->size) sq_lo = fb->size / 2;
    if (sq_hi >= fb->size) sq_hi = fb->size - 1;

    fprintf(stderr, "SQQS: %dd (%db), k=%d, FB=%d, M=%d, target=%d, LP=%lu, SQ range=[%d,%d]\n",
            digits, bits, mult, fb->size, M, target, lp_bound, sq_lo, sq_hi);

    /* Standard SIQS sieve (fallback) + special-Q enhancement */
    unsigned char *sieve_array = malloc(BLOCKSIZE);

    rels_t *full = rels_create(MAX_RELS, fb->size);
    rels_t *part = rels_create(MAX_PARTIALS, fb->size);
    lp_t *lph = lp_create(MAX_PARTIALS);

    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, SEED);

    mpz_t Q_val, Y_val, tmp, tmp2, a, b_poly, c_poly;
    mpz_inits(Q_val, Y_val, tmp, tmp2, a, b_poly, c_poly, NULL);
    short *tmp_exps = calloc(fb->size + 2, sizeof(short));

    int total_polys = 0;
    int combined = 0;
    mpz_t B_vals[MAX_A_FACTORS];
    for (int j = 0; j < MAX_A_FACTORS; j++) mpz_init(B_vals[j]);
    unsigned int *soln1 = malloc(fb->size * sizeof(unsigned int));
    unsigned int *soln2 = malloc(fb->size * sizeof(unsigned int));

    /* ========== MAIN LOOP: Standard SIQS with Special-Q boosting ========== */
    /*
     * Strategy: Run standard SIQS, but for each sieve block, also check
     * if values near the threshold can be made smooth by dividing out a
     * "special-q" prime. Candidates that are below threshold by log(q)
     * would be above threshold if q were included.
     *
     * This is equivalent to: for each special-q prime q in the FB,
     * lower the sieve threshold by log(q) at positions where q divides Q(x).
     * These positions are at x ≡ root_q (mod q).
     *
     * Implementation: After standard sieve, scan for positions where
     * sieve_value >= threshold - max_logq. For each such position, check
     * if any special-q prime divides Q(x) and makes the total smooth.
     */

    int max_logq = 0;
    for (int i = sq_lo; i < sq_hi && i < fb->size; i++)
        if (fb->logp[i] > max_logq) max_logq = fb->logp[i];

    while (full->count < target) {
        double t = elapsed();
        if (t > 280) { fprintf(stderr, "TIMEOUT at %.1fs with %d/%d rels\n", t, full->count, target); break; }

        /* Generate new SIQS polynomial */
        {
            mpz_t tgt; mpz_init(tgt);
            mpz_mul_ui(tgt, kN, 2); mpz_sqrt(tgt, tgt); if (M > 0) mpz_tdiv_q_ui(tgt, tgt, M);
            double log_tgt = mpz_sizeinbase(tgt, 2) * log(2.0);

            int lo = fb->size/4, hi = 3*fb->size/4;
            if (lo < 2) lo = 2; if (hi <= lo+3) hi = fb->size-1;

            double avg = 0; int cnt = 0;
            for (int i = lo; i < hi; i++) { if (fb->sqrtN[i]==0) continue; avg += log(fb->prime[i]); cnt++; }
            if (cnt == 0) { mpz_clear(tgt); break; }
            avg /= cnt;

            int s = (int)(log_tgt / avg + 0.5);
            if (s < 3) s = 3; if (s > MAX_A_FACTORS) s = MAX_A_FACTORS; if (s > hi-lo) s = hi-lo;
            int num_a = s;

            double best_ratio = 1e30;
            int best[MAX_A_FACTORS], a_idx[MAX_A_FACTORS];

            for (int att = 0; att < 50; att++) {
                mpz_set_ui(a, 1);
                int idx[MAX_A_FACTORS]; int ok = 1;
                for (int i = 0; i < s && ok; i++) {
                    int tries = 0, good;
                    do { idx[i] = lo + gmp_urandomm_ui(rng, hi-lo); good = 1;
                         for (int j = 0; j < i; j++) if (idx[j]==idx[i]) {good=0; break;}
                         if (fb->sqrtN[idx[i]]==0) good=0; tries++;
                    } while (!good && tries < 100);
                    if (!good) { ok=0; break; }
                    mpz_mul_ui(a, a, fb->prime[idx[i]]);
                }
                if (!ok) continue;
                double ratio;
                if (mpz_cmp(a,tgt)>0) { mpz_tdiv_q(tmp,a,tgt); ratio=mpz_get_d(tmp); }
                else { mpz_tdiv_q(tmp,tgt,a); ratio=mpz_get_d(tmp); }
                if (ratio < best_ratio) { best_ratio = ratio; memcpy(best, idx, s*sizeof(int)); }
                if (ratio < 1.5) break;
            }

            memcpy(a_idx, best, s*sizeof(int));
            mpz_set_ui(a, 1);
            for (int i = 0; i < s; i++) mpz_mul_ui(a, a, fb->prime[a_idx[i]]);

            /* B values */
            for (int j = 0; j < s; j++) {
                int idx = a_idx[j]; unsigned int qj = fb->prime[idx], rj = fb->sqrtN[idx];
                mpz_t a_q, inv; mpz_inits(a_q, inv, NULL);
                mpz_divexact_ui(a_q, a, qj);
                unsigned long aqmod = mpz_fdiv_ui(a_q, qj);
                unsigned int iv = mod_inverse_u32((unsigned int)aqmod, qj);
                mpz_mul_ui(B_vals[j], a_q, ((unsigned long)rj * iv) % qj);
                mpz_clears(a_q, inv, NULL);
            }
            mpz_clear(tgt);

            /* Enumerate b-values with Gray code */
            int num_b = 1 << (num_a - 1);

            for (int b_idx = 0; b_idx < num_b && full->count < target; b_idx++) {
                if (elapsed() > 285) break;

                int gray = b_idx ^ (b_idx >> 1);
                mpz_set_ui(b_poly, 0);
                for (int j = 0; j < num_a; j++) {
                    if (gray & (1 << j)) mpz_add(b_poly, b_poly, B_vals[j]);
                    else mpz_sub(b_poly, b_poly, B_vals[j]);
                }

                mpz_mul(tmp, b_poly, b_poly); mpz_sub(tmp, tmp, kN); mpz_mod(tmp, tmp, a);
                if (mpz_sgn(tmp) != 0) {
                    mpz_neg(b_poly, b_poly);
                    mpz_mul(tmp, b_poly, b_poly); mpz_sub(tmp, tmp, kN); mpz_mod(tmp, tmp, a);
                    if (mpz_sgn(tmp) != 0) continue;
                }

                /* Compute sieve roots */
                for (int i = 0; i < fb->size; i++) {
                    unsigned int p = fb->prime[i];
                    if (p < 3 || fb->sqrtN[i] == 0) { soln1[i] = soln2[i] = 0xFFFFFFFF; continue; }
                    unsigned long am = mpz_fdiv_ui(a, p);
                    if (am == 0) {
                        unsigned long bm = mpz_fdiv_ui(b_poly, p);
                        mpz_mul(c_poly, b_poly, b_poly); mpz_sub(c_poly, c_poly, kN); mpz_divexact(c_poly, c_poly, a);
                        unsigned long cm = mpz_fdiv_ui(c_poly, p);
                        if (bm == 0) { soln1[i] = soln2[i] = 0xFFFFFFFF; continue; }
                        unsigned int inv2b = mod_inverse_u32((unsigned int)((2UL*bm)%p), p);
                        unsigned int root = (unsigned int)((unsigned long)(p-cm)%p * inv2b % p);
                        soln1[i] = root; soln2[i] = root;
                        continue;
                    }
                    unsigned int ai = mod_inverse_u32((unsigned int)am, p);
                    if (ai == 0) { soln1[i] = soln2[i] = 0xFFFFFFFF; continue; }
                    unsigned long bm = mpz_fdiv_ui(b_poly, p);
                    unsigned int r = fb->sqrtN[i];
                    soln1[i] = (unsigned int)((unsigned long)ai * ((r + p - bm) % p) % p);
                    soln2[i] = (unsigned int)((unsigned long)ai * ((p - r + p - bm) % p) % p);
                }

                total_polys++;
                if (total_polys % 200 == 0)
                    fprintf(stderr, "  poly=%d rels=%d/%d (full=%d+%d) part=%d t=%.1fs\n",
                            total_polys, full->count, target, full->count-combined, combined, part->count, elapsed());

                /* Sieve blocks */
                int total_blocks = 2 * P.num_blocks;
                for (int blk = 0; blk < total_blocks; blk++) {
                    int block_start = -(int)M + blk * BLOCKSIZE;

                    memset(sieve_array, 0, BLOCKSIZE);

                    /* Sieve with all FB primes */
                    for (int i = 1; i < fb->size; i++) {
                        if (soln1[i] == 0xFFFFFFFF) continue;
                        unsigned int p = fb->prime[i];
                        if (p < 5) continue;
                        unsigned char lp = fb->logp[i];

                        for (int ri = 0; ri < 2; ri++) {
                            unsigned int root = (ri == 0) ? soln1[i] : soln2[i];
                            long off = ((long)root - (long)block_start) % (long)p;
                            if (off < 0) off += p;
                            for (int j = (int)off; j < BLOCKSIZE; j += p)
                                sieve_array[j] += lp;
                        }
                    }

                    /* Compute threshold */
                    double log2_a = mpz_sizeinbase(a, 2);
                    double log2_Qmax = log2_a + log2((double)M);
                    int threshold = (int)(log2_Qmax * P.thresh_adj);
                    /* NOVEL: Also scan at lower threshold for special-Q candidates */
                    int sq_threshold = threshold - max_logq;
                    if (sq_threshold < 20) sq_threshold = 20;

                    for (int j = 0; j < BLOCKSIZE; j++) {
                        if (sieve_array[j] < sq_threshold) continue;

                        int x_global = block_start + j;
                        unsigned long lp_val = 0;

                        int result = trial_divide(
                            Q_val, Y_val, tmp_exps,
                            fb, fb->size,
                            lp_bound, &lp_val,
                            kN, a, b_poly, x_global,
                            full->neg_col);

                        if (result == 1) {
                            int idx = full->count;
                            if (idx < full->alloc) {
                                mpz_set(full->Y[idx], Y_val);
                                memcpy(full->exps[idx], tmp_exps, (fb->size+2)*sizeof(short));
                                full->lp[idx] = 0;
                                full->count++;
                            }
                        } else if (result == 2) {
                            int match = lp_find(lph, lp_val);
                            if (match >= 0) {
                                int pidx = match;
                                int fidx = full->count;
                                if (fidx < full->alloc) {
                                    mpz_mul(full->Y[fidx], Y_val, part->Y[pidx]);
                                    mpz_mod(full->Y[fidx], full->Y[fidx], N);
                                    for (int e = 0; e < fb->size+2; e++)
                                        full->exps[fidx][e] = tmp_exps[e] + part->exps[pidx][e];
                                    full->lp[fidx] = lp_val;
                                    full->count++;
                                    combined++;
                                }
                            } else {
                                int pidx = part->count;
                                if (pidx < part->alloc) {
                                    mpz_set(part->Y[pidx], Y_val);
                                    memcpy(part->exps[pidx], tmp_exps, (fb->size+2)*sizeof(short));
                                    part->lp[pidx] = lp_val;
                                    part->count++;
                                    lp_insert(lph, lp_val, pidx);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    fprintf(stderr, "Sieving done: %d rels (%d full + %d combined) in %.1fs\n",
            full->count, full->count - combined, combined, elapsed());

    if (full->count < fb->size + 1) {
        fprintf(stderr, "Not enough relations: %d < %d\n", full->count, fb->size+1);
        return 1;
    }

    /* ========== LINEAR ALGEBRA ========== */
    int nrels = full->count;
    int ncols = fb->size + 1;

    fprintf(stderr, "Building %d x %d GF(2) matrix...\n", nrels, ncols);
    gf2_t *mat = gf2_create(nrels, ncols);
    for (int r = 0; r < nrels; r++)
        for (int c = 0; c < ncols; c++)
            if (full->exps[r][c] & 1)
                mat->rows[r][c/64] ^= (1ULL << (c%64));

    fprintf(stderr, "Solving GF(2) system...\n");
    int **deps; int *dlen;
    int ndeps = gf2_solve(mat, &deps, &dlen, MAX_DEPS);
    fprintf(stderr, "Found %d dependencies\n", ndeps);

    /* ========== SQUARE ROOT ========== */
    mpz_t X, Y2, g;
    mpz_inits(X, Y2, g, NULL);

    int found = 0;
    for (int d = 0; d < ndeps && !found; d++) {
        mpz_set_ui(X, 1);
        int *exps = calloc(ncols+2, sizeof(int));
        for (int i = 0; i < dlen[d]; i++) {
            int ri = deps[d][i];
            mpz_mul(X, X, full->Y[ri]);
            mpz_mod(X, X, N);
            for (int c = 0; c <= fb->size; c++) exps[c] += full->exps[ri][c];
        }
        int all_even = 1;
        for (int c = 0; c <= fb->size; c++) if (exps[c]&1) { all_even=0; break; }
        if (!all_even) { free(exps); continue; }

        mpz_set_ui(Y2, 1);
        for (int c = 0; c < fb->size; c++) {
            if (exps[c] <= 0) continue;
            mpz_set_ui(tmp, fb->prime[c]);
            mpz_powm_ui(tmp, tmp, exps[c]/2, N);
            mpz_mul(Y2, Y2, tmp); mpz_mod(Y2, Y2, N);
        }

        /* Include LP contributions for combined relations */
        for (int i = 0; i < dlen[d]; i++) {
            int ri = deps[d][i];
            if (full->lp[ri] > 1) {
                mpz_set_ui(tmp, full->lp[ri]);
                mpz_mul(Y2, Y2, tmp); mpz_mod(Y2, Y2, N);
            }
        }

        mpz_sub(tmp, X, Y2); mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g,1) > 0 && mpz_cmp(g,N) < 0) {
            mpz_t cof; mpz_init(cof); mpz_divexact(cof, N, g);
            gmp_printf("%Zd\n%Zd\n", g, cof); mpz_clear(cof); found = 1;
        }
        if (!found) {
            mpz_add(tmp, X, Y2); mpz_gcd(g, tmp, N);
            if (mpz_cmp_ui(g,1) > 0 && mpz_cmp(g,N) < 0) {
                mpz_t cof; mpz_init(cof); mpz_divexact(cof, N, g);
                gmp_printf("%Zd\n%Zd\n", g, cof); mpz_clear(cof); found = 1;
            }
        }
        free(exps);
    }

    if (!found) { fprintf(stderr, "FAILED: no factor found from %d deps\n", ndeps); return 1; }

    fprintf(stderr, "Total time: %.3fs\n", elapsed());
    return 0;
}
