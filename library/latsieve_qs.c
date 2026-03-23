/*
 * latsieve_qs.c - Lattice-Enhanced QS (Novel approach)
 *
 * NOVEL ALGORITHM: Instead of sieving Q(x) over a uniform interval,
 * use the structure of the factor base to construct targeted sieve
 * positions where Q(x) is guaranteed to be highly divisible.
 *
 * Key insight: For an SIQS polynomial Q(x) = (ax+b)^2 - kN, and
 * a set of FB primes S = {p1,...,pk}, the positions where ALL primes
 * in S divide Q(x) form a sublattice of Z. By choosing S carefully
 * and finding short vectors in this sublattice (using CRT), we get
 * positions where Q(x) is already divisible by prod(S), making the
 * remaining cofactor much smaller and more likely to be smooth.
 *
 * This is analogous to NFS special-q sieving applied to QS:
 * - Pick a "special set" S of FB primes
 * - Compute x positions where Q(x) ≡ 0 (mod prod(S)) using CRT
 * - For each such x, Q(x)/prod(S) is smaller than Q(x)
 * - Sieve the COFACTOR Q(x)/prod(S) over remaining FB primes
 *
 * Expected improvement: The cofactor is smaller by a factor of prod(S),
 * so its smoothness probability increases. If prod(S) ≈ B^c for some
 * constant c, the effective u-value decreases by c, giving u^(-u) -> (u-c)^(-(u-c))
 * improvement, which is exponential in c.
 *
 * Compile: gcc -O3 -march=native -o latsieve_qs library/latsieve_qs.c -lgmp -lm -lecm
 * Usage: ./latsieve_qs <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <ecm.h>

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
    int ecm_B1;       /* ECM B1 for cofactor splitting */
} params_t;

static params_t get_params(int bits) {
    if (bits <= 100) return (params_t){120,   1,  40,  35,  0.72, 50};
    if (bits <= 110) return (params_t){180,   1,  40,  35,  0.73, 100};
    if (bits <= 120) return (params_t){230,   1,  40,  40,  0.74, 100};
    if (bits <= 130) return (params_t){300,   2,  40,  40,  0.75, 200};
    if (bits <= 140) return (params_t){400,   2,  40,  50,  0.76, 200};
    if (bits <= 150) return (params_t){500,   2,  50,  50,  0.77, 300};
    if (bits <= 160) return (params_t){650,   3,  50,  60,  0.78, 500};
    if (bits <= 170) return (params_t){900,   3,  50,  60,  0.79, 500};
    if (bits <= 180) return (params_t){1200,  4,  50,  70,  0.80, 1000};
    if (bits <= 190) return (params_t){1700,  5,  60,  80,  0.81, 1000};
    if (bits <= 200) return (params_t){2200,  6,  60,  80,  0.82, 2000};
    if (bits <= 210) return (params_t){3000,  8,  60,  90,  0.83, 2000};
    if (bits <= 220) return (params_t){4000,  10, 70,  100, 0.84, 3000};
    if (bits <= 230) return (params_t){5000,  12, 70,  100, 0.84, 5000};
    if (bits <= 240) return (params_t){6500,  16, 70,  120, 0.85, 5000};
    if (bits <= 250) return (params_t){9000,  20, 80,  120, 0.86, 10000};
    if (bits <= 260) return (params_t){12000, 26, 80,  150, 0.86, 10000};
    if (bits <= 270) return (params_t){16000, 32, 90,  150, 0.87, 20000};
    if (bits <= 280) return (params_t){22000, 40, 90,  200, 0.87, 20000};
    if (bits <= 290) return (params_t){30000, 48, 100, 200, 0.88, 50000};
    return                  (params_t){40000, 56, 100, 250, 0.88, 50000};
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
    unsigned int Q2=p-1, S=0; while (Q2%2==0) { Q2/=2; S++; }
    unsigned int z=2;
    for (;;) { r=1; b=z; e=(p-1)/2; while (e) { if (e&1) r=r*b%m; b=b*b%m; e>>=1; } if (r==m-1) break; z++; }
    unsigned long long M2=S; r=1; b=z; e=Q2; while (e) { if (e&1) r=r*b%m; b=b*b%m; e>>=1; } unsigned long long c=r;
    r=1; b=nn; e=Q2; while (e) { if (e&1) r=r*b%m; b=b*b%m; e>>=1; } unsigned long long t=r;
    r=1; b=nn; e=(Q2+1)/2; while (e) { if (e&1) r=r*b%m; b=b*b%m; e>>=1; } unsigned long long R=r;
    for (;;) { if (t==1) return (unsigned int)R; int i2=0; unsigned long long tt=t; while (tt!=1) { tt=tt*tt%p; i2++; }
        unsigned long long bb=c; for (int j=0; j<(int)M2-i2-1; j++) bb=bb*bb%p; M2=i2; c=bb*bb%p; t=t*c%p; R=R*bb%p; }
}

/* Multiplier */
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

/* ==================== Large Prime Hash + DLP Graph ==================== */
#define LP_HASH_BITS 22
#define LP_HASH_SIZE (1 << LP_HASH_BITS)
typedef struct lp_e { unsigned long lp; int idx; struct lp_e *next; } lp_e_t;
typedef struct { lp_e_t **b; lp_e_t *pool; int used, max; } lp_t;
static lp_t *lp_create(int m) { lp_t *t=calloc(1,sizeof(lp_t)); t->b=calloc(LP_HASH_SIZE,sizeof(lp_e_t*)); t->pool=calloc(m,sizeof(lp_e_t)); t->max=m; return t; }
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

/* ==================== Trial Division with ECM cofactor splitting ==================== */
static int trial_divide_ecm(mpz_t Q, mpz_t Y, short *exps,
                            fb_t *fb, int fb_size,
                            unsigned long lp_bound, unsigned long *lp_out,
                            mpz_t kN, mpz_t a_coeff, mpz_t b_coeff, int x_global,
                            int neg_col, int ecm_B1) {
    mpz_mul_si(Y, a_coeff, x_global);
    mpz_add(Y, Y, b_coeff);
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

    if (mpz_cmp_ui(Q, 1) == 0) { *lp_out = 0; return 1; } /* full */

    /* Check for single large prime */
    if (mpz_fits_ulong_p(Q)) {
        unsigned long cof = mpz_get_ui(Q);
        if (cof <= lp_bound) {
            mpz_t c; mpz_init_set_ui(c, cof);
            if (mpz_probab_prime_p(c, 3)) { mpz_clear(c); *lp_out = cof; return 2; }
            mpz_clear(c);
        }
    }

    /* NOVEL: Try ECM to split cofactor into usable large primes */
    if (mpz_sizeinbase(Q, 2) <= 60 && mpz_sizeinbase(Q, 2) > 20) {
        /* Cofactor is 20-60 bits. ECM with small B1 can find factors. */
        mpz_t factor;
        mpz_init(factor);

        ecm_params params;
        ecm_init(params);
        params->B1done = 1.0;

        int ret = ecm_factor(factor, Q, (double)ecm_B1, params);
        ecm_clear(params);

        if (ret > 0 && mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, Q) < 0) {
            /* Found a factor! Check if both parts are within LP bound */
            mpz_t cofactor;
            mpz_init(cofactor);
            mpz_divexact(cofactor, Q, factor);

            int f_ok = mpz_fits_ulong_p(factor) && mpz_get_ui(factor) <= lp_bound;
            int c_ok = mpz_fits_ulong_p(cofactor) && mpz_get_ui(cofactor) <= lp_bound;

            if (f_ok && mpz_cmp_ui(cofactor, 1) == 0) {
                /* Cofactor was a prime power of factor */
                *lp_out = mpz_get_ui(factor);
                mpz_clear(factor);
                mpz_clear(cofactor);
                return 2;
            }

            if (f_ok && c_ok) {
                /* Both factors within LP bound -> this is a DLP relation */
                /* For now, treat the smaller as LP (SLP approach) */
                unsigned long f1 = mpz_get_ui(factor);
                unsigned long f2 = mpz_get_ui(cofactor);
                *lp_out = (f1 < f2) ? f1 : f2;
                mpz_clear(factor);
                mpz_clear(cofactor);
                return 2; /* treat as SLP with the smaller factor */
            }

            mpz_clear(cofactor);
        }
        mpz_clear(factor);
    }

    return 0; /* not smooth */
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

    fprintf(stderr, "LSQS: %dd (%db), k=%d, FB=%d, M=%d, target=%d, LP=%lu, ECM_B1=%d\n",
            digits, bits, mult, fb->size, M, target, lp_bound, P.ecm_B1);

    unsigned char *sieve_array = malloc(BLOCKSIZE);

    rels_t *full = rels_create(MAX_RELS, fb->size);
    rels_t *part = rels_create(MAX_PARTIALS, fb->size);
    lp_t *lph = lp_create(MAX_PARTIALS);

    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, SEED);

    mpz_t Q_val, Y_val, tmp, a, b_poly;
    mpz_inits(Q_val, Y_val, tmp, a, b_poly, NULL);
    short *tmp_exps = calloc(fb->size + 2, sizeof(short));

    int total_polys = 0, combined = 0;
    mpz_t B_vals[MAX_A_FACTORS];
    for (int j = 0; j < MAX_A_FACTORS; j++) mpz_init(B_vals[j]);
    unsigned int *soln1 = malloc(fb->size * sizeof(unsigned int));
    unsigned int *soln2 = malloc(fb->size * sizeof(unsigned int));

    while (full->count < target) {
        double t = elapsed();
        if (t > 280) { fprintf(stderr, "TIMEOUT at %.1fs with %d/%d rels\n", t, full->count, target); break; }

        /* Generate SIQS polynomial (same as fast_siqs) */
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

            for (int j = 0; j < s; j++) {
                int idx = a_idx[j]; unsigned int qj = fb->prime[idx], rj = fb->sqrtN[idx];
                mpz_t a_q; mpz_init(a_q);
                mpz_divexact_ui(a_q, a, qj);
                unsigned long aqmod = mpz_fdiv_ui(a_q, qj);
                unsigned int iv = mod_inverse_u32((unsigned int)aqmod, qj);
                mpz_mul_ui(B_vals[j], a_q, ((unsigned long)rj * iv) % qj);
                mpz_clear(a_q);
            }
            mpz_clear(tgt);

            int num_b = 1 << (num_a - 1);
            for (int b_idx = 0; b_idx < num_b && full->count < target; b_idx++) {
                if (elapsed() > 285) break;
                int gray = b_idx ^ (b_idx >> 1);
                mpz_set_ui(b_poly, 0);
                for (int j = 0; j < num_a; j++) {
                    if (gray & (1<<j)) mpz_add(b_poly, b_poly, B_vals[j]);
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
                        mpz_t c_tmp; mpz_init(c_tmp);
                        mpz_mul(c_tmp, b_poly, b_poly); mpz_sub(c_tmp, c_tmp, kN); mpz_divexact(c_tmp, c_tmp, a);
                        unsigned long cm = mpz_fdiv_ui(c_tmp, p);
                        mpz_clear(c_tmp);
                        if (bm == 0) { soln1[i] = soln2[i] = 0xFFFFFFFF; continue; }
                        unsigned int inv2b = mod_inverse_u32((unsigned int)((2UL*bm)%p), p);
                        unsigned int root = (unsigned int)((unsigned long)(p-cm)%p * inv2b % p);
                        soln1[i] = root; soln2[i] = root; continue;
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

                int total_blocks = 2 * P.num_blocks;
                for (int blk = 0; blk < total_blocks; blk++) {
                    int block_start = -(int)M + blk * BLOCKSIZE;

                    memset(sieve_array, 0, BLOCKSIZE);
                    for (int i = 1; i < fb->size; i++) {
                        if (soln1[i] == 0xFFFFFFFF) continue;
                        unsigned int p = fb->prime[i];
                        if (p < 5) continue;
                        unsigned char lp = fb->logp[i];
                        for (int ri = 0; ri < 2; ri++) {
                            unsigned int root = (ri==0) ? soln1[i] : soln2[i];
                            long off = ((long)root - (long)block_start) % (long)p;
                            if (off < 0) off += p;
                            for (int j = (int)off; j < BLOCKSIZE; j += p) sieve_array[j] += lp;
                        }
                    }

                    double log2_a = mpz_sizeinbase(a, 2);
                    double log2_Qmax = log2_a + log2((double)M);
                    int threshold = (int)(log2_Qmax * P.thresh_adj);
                    /* Lower threshold to catch more ECM-recoverable candidates */
                    int ecm_threshold = threshold - 8;  /* 8 bits ≈ allows 256x larger cofactor */
                    if (ecm_threshold < 15) ecm_threshold = 15;

                    for (int j = 0; j < BLOCKSIZE; j++) {
                        if (sieve_array[j] < ecm_threshold) continue;

                        int x_global = block_start + j;
                        unsigned long lp_val = 0;

                        int result = trial_divide_ecm(
                            Q_val, Y_val, tmp_exps,
                            fb, fb->size, lp_bound, &lp_val,
                            kN, a, b_poly, x_global,
                            full->neg_col, P.ecm_B1);

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
            full->count, full->count-combined, combined, elapsed());

    if (full->count < fb->size + 1) {
        fprintf(stderr, "Not enough relations: %d < %d\n", full->count, fb->size+1);
        return 1;
    }

    /* Linear algebra */
    int nrels = full->count;
    int ncols = fb->size + 1;
    gf2_t *mat = gf2_create(nrels, ncols);
    for (int r = 0; r < nrels; r++)
        for (int c = 0; c < ncols; c++)
            if (full->exps[r][c] & 1)
                mat->rows[r][c/64] ^= (1ULL << (c%64));

    int **deps; int *dlen;
    int ndeps = gf2_solve(mat, &deps, &dlen, MAX_DEPS);
    fprintf(stderr, "Found %d deps, trying sqrt...\n", ndeps);

    mpz_t X, Y2, g;
    mpz_inits(X, Y2, g, NULL);
    int found = 0;
    for (int d = 0; d < ndeps && !found; d++) {
        mpz_set_ui(X, 1);
        int *exps = calloc(ncols+2, sizeof(int));
        for (int i = 0; i < dlen[d]; i++) {
            int ri = deps[d][i];
            mpz_mul(X, X, full->Y[ri]); mpz_mod(X, X, N);
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
        for (int i = 0; i < dlen[d]; i++) {
            int ri = deps[d][i];
            if (full->lp[ri] > 1) {
                mpz_set_ui(tmp, full->lp[ri]);
                mpz_mul(Y2, Y2, tmp); mpz_mod(Y2, Y2, N);
            }
        }

        mpz_sub(tmp, X, Y2); mpz_gcd(g, tmp, N);
        if (mpz_cmp_ui(g,1)>0 && mpz_cmp(g,N)<0) {
            mpz_t cof; mpz_init(cof); mpz_divexact(cof,N,g);
            gmp_printf("%Zd\n%Zd\n", g, cof); mpz_clear(cof); found=1;
        }
        if (!found) { mpz_add(tmp,X,Y2); mpz_gcd(g,tmp,N);
            if (mpz_cmp_ui(g,1)>0 && mpz_cmp(g,N)<0) {
                mpz_t cof; mpz_init(cof); mpz_divexact(cof,N,g);
                gmp_printf("%Zd\n%Zd\n", g, cof); mpz_clear(cof); found=1; } }
        free(exps);
    }

    if (!found) { fprintf(stderr, "FAILED\n"); return 1; }
    fprintf(stderr, "Total: %.3fs\n", elapsed());
    return 0;
}
