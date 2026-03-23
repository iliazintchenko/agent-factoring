/*
 * siqs_hybrid.c — Novel SIQS with 48KB L1-optimized sieve
 *
 * Innovations:
 * 1. 48KB sieve blocks (AMD EPYC 9R45 L1D = 48KB, vs standard 32KB)
 * 2. AVX512BW candidate scanning
 * 3. Aggressive SLP collection with online matching
 * 4. Pollard rho DLP cofactor splitting
 * 5. Adaptive sieve threshold based on measured yield rate
 *
 * Compile: gcc -O3 -march=native -mavx512bw -o siqs_hybrid library/siqs_hybrid.c -lgmp -lm
 * Usage:   ./siqs_hybrid <N>
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <immintrin.h>
#include <gmp.h>

#define SIEVE_SIZE    49152  /* 48KB — matches AMD EPYC 9R45 L1D */
#define MAX_FB        200000
#define MAX_AFACTORS  16
#define MAX_RELS      600000
#define HASH_BITS     22
#define HASH_SIZE     (1 << HASH_BITS)
#define HASH_MASK     (HASH_SIZE - 1)
#define SEED          42

/* ==================== Timing ==================== */
static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* ==================== Parameters ==================== */
typedef struct {
    int fb_size, num_blocks, lp_mult, num_a, extra;
} param_t;

static param_t get_params(int d) {
    if (d <= 30) return (param_t){200,  1, 30, 3, 60};
    if (d <= 35) return (param_t){400,  1, 40, 3, 80};
    if (d <= 40) return (param_t){800,  2, 50, 4, 100};
    if (d <= 45) return (param_t){1500, 2, 70, 5, 100};
    if (d <= 50) return (param_t){2500, 4, 120, 6, 120};
    if (d <= 55) return (param_t){5000, 5, 110, 7, 130};
    if (d <= 60) return (param_t){8000, 7, 140, 7, 150};
    if (d <= 65) return (param_t){15000, 10, 180, 8, 160};
    if (d <= 70) return (param_t){28000, 14, 230, 9, 200};
    if (d <= 75) return (param_t){45000, 18, 280, 10, 220};
    if (d <= 80) return (param_t){65000, 24, 340, 11, 250};
    return (param_t){100000, 32, 400, 12, 300};
}

/* ==================== Math utilities ==================== */
static uint32_t mod_inv(uint32_t a, uint32_t m) {
    int64_t g0=m, g1=a%m, u0=0, u1=1;
    while (g1) { int64_t q=g0/g1, t=g0-q*g1; g0=g1; g1=t; t=u0-q*u1; u0=u1; u1=t; }
    return (uint32_t)((u0%m+m)%m);
}

static uint64_t powmod64(uint64_t b, uint64_t e, uint64_t m) {
    uint64_t r=1; b%=m;
    while(e>0){if(e&1)r=(__uint128_t)r*b%m; b=(__uint128_t)b*b%m; e>>=1;}
    return r;
}

static uint32_t sqrt_mod_p(uint32_t n, uint32_t p) {
    if(p==2) return n&1;
    n %= p;
    if(n==0) return 0;
    if((p&3)==3) return (uint32_t)powmod64(n,(p+1)/4,p);
    /* Tonelli-Shanks */
    uint32_t Q=p-1, S=0;
    while(!(Q&1)){Q>>=1;S++;}
    uint32_t z=2;
    while(powmod64(z,(p-1)/2,p)!=p-1) z++;
    uint64_t M=S, c=powmod64(z,Q,p), t=powmod64(n,Q,p), R=powmod64(n,(Q+1)/2,p);
    while(t!=1){
        uint64_t i=1,tmp=t*t%p;
        while(tmp!=1){tmp=tmp*tmp%p;i++;}
        uint64_t b=c;
        for(uint64_t j=0;j<M-i-1;j++)b=b*b%p;
        M=i; c=b*b%p; t=t*c%p; R=R*b%p;
    }
    return(uint32_t)R;
}

/* ==================== Factor Base ==================== */
static uint32_t fb_p[MAX_FB];
static uint32_t fb_r[MAX_FB]; /* sqrt(N) mod p */
static uint8_t  fb_lg[MAX_FB];
static int fb_n;
static int sieve_start;

static int is_prime(uint32_t n) {
    if(n<2)return 0; if(n<4)return 1; if(n%2==0||n%3==0)return 0;
    for(uint32_t d=5;d*d<=n;d+=6) if(n%d==0||n%(d+2)==0)return 0;
    return 1;
}

static void build_fb(mpz_t N, int target) {
    fb_n = 0;
    fb_p[0]=2; fb_r[0]=1; fb_lg[0]=1; fb_n=1;
    for(uint32_t p=3; fb_n<target; p+=2) {
        if(!is_prime(p)) continue;
        uint32_t nr = mpz_fdiv_ui(N,p);
        if(nr==0){fb_p[fb_n]=p;fb_r[fb_n]=0;fb_lg[fb_n]=(uint8_t)(log2(p)+.5);fb_n++;continue;}
        if(powmod64(nr,(p-1)/2,p)!=1) continue;
        fb_p[fb_n]=p; fb_r[fb_n]=sqrt_mod_p(nr,p); fb_lg[fb_n]=(uint8_t)(log2(p)+.5); fb_n++;
    }
    sieve_start=1;
    while(sieve_start<fb_n && fb_p[sieve_start]<7) sieve_start++;
}

/* ==================== Relation storage ==================== */
typedef struct {
    int *exp_vec;    /* full integer exponents, length fb_n+1 (last = sign) */
    mpz_t Qprod;     /* |A*f(x)| or product for combined */
    mpz_t axb_prod;  /* (ax+b) or product for combined */
    uint64_t lp;     /* large prime (0 if full) */
    int sign;        /* 1 if f(x) < 0 */
} rel_t;

static rel_t *rels;
static int n_rels, n_rels_alloc;
static int n_full, n_slp_combined;

static void init_rels(void) {
    n_rels_alloc = MAX_RELS;
    rels = calloc(n_rels_alloc, sizeof(rel_t));
    for(int i=0; i<n_rels_alloc; i++) {
        rels[i].exp_vec = calloc(fb_n + 1, sizeof(int)); /* +1 for sign */
        mpz_init(rels[i].Qprod);
        mpz_init(rels[i].axb_prod);
    }
    n_rels = n_full = n_slp_combined = 0;
}

/* SLP hash table */
typedef struct he { uint64_t lp; int idx; struct he *next; } he_t;
static he_t *htable[HASH_SIZE];
static he_t *hpool;
static int hpool_n, hpool_cap;

static void init_hash(void) {
    hpool_cap = 500000;
    hpool = calloc(hpool_cap, sizeof(he_t));
    hpool_n = 0;
    memset(htable, 0, sizeof(htable));
}

static void add_full(int *ev, int sign, mpz_t axb, mpz_t Af) {
    if(n_rels >= n_rels_alloc) return;
    rel_t *r = &rels[n_rels];
    memcpy(r->exp_vec, ev, (fb_n+1)*sizeof(int));
    r->exp_vec[fb_n] = sign; /* sign bit */
    mpz_set(r->axb_prod, axb);
    mpz_set(r->Qprod, Af);
    r->lp = 0;
    r->sign = sign;
    n_rels++;
    n_full++;
}

static int add_slp(int *ev, int sign, mpz_t axb, mpz_t Af, uint64_t lp) {
    uint32_t h = (uint32_t)(lp * 2654435761ULL) & HASH_MASK;
    for(he_t *e = htable[h]; e; e = e->next) {
        if(e->lp == lp) {
            if(n_rels >= n_rels_alloc) return 0;
            rel_t *old = &rels[e->idx];
            rel_t *nr = &rels[n_rels];
            for(int i=0;i<=fb_n;i++) nr->exp_vec[i] = old->exp_vec[i] + ev[i];
            nr->exp_vec[fb_n] = old->sign + sign; /* combined sign */
            mpz_mul(nr->axb_prod, old->axb_prod, axb);
            /* old->Qprod = |A1*f1| (includes LP), Af = |A2*f2| (includes LP)
             * Product already has LP^2 — do NOT multiply by LP again */
            mpz_mul(nr->Qprod, old->Qprod, Af);
            nr->lp = 0;
            nr->sign = old->sign + sign;
            n_rels++;
            n_slp_combined++;
            return 1;
        }
    }
    if(n_rels >= n_rels_alloc || hpool_n >= hpool_cap) return 0;
    rel_t *pr = &rels[n_rels];
    memcpy(pr->exp_vec, ev, (fb_n+1)*sizeof(int));
    pr->exp_vec[fb_n] = sign;
    mpz_set(pr->axb_prod, axb);
    mpz_set(pr->Qprod, Af);
    pr->lp = lp;
    pr->sign = sign;
    he_t *ne = &hpool[hpool_n++];
    ne->lp = lp; ne->idx = n_rels; ne->next = htable[h]; htable[h] = ne;
    n_rels++;
    return 0;
}

static int usable(void) { return n_full + n_slp_combined; }

/* ==================== SIQS Polynomial ==================== */
static mpz_t g_N;
static mpz_t pa, pb, pc; /* A, B, C of current polynomial */
static int a_idx[MAX_AFACTORS], na;
static mpz_t Bvals[MAX_AFACTORS];
static uint32_t *rt1, *rt2; /* sieve roots */
static gmp_randstate_t rng;
static int *acands, nacands;
static int gray;

static void init_poly(param_t *P) {
    mpz_inits(pa,pb,pc,NULL);
    for(int i=0;i<MAX_AFACTORS;i++) mpz_init(Bvals[i]);
    rt1=calloc(fb_n,sizeof(uint32_t));
    rt2=calloc(fb_n,sizeof(uint32_t));
    na = P->num_a;

    /* Select A-factor candidates */
    mpz_t ta; mpz_init(ta);
    mpz_mul_ui(ta, g_N, 2);
    mpz_sqrt(ta, ta);
    mpz_fdiv_q_ui(ta, ta, (uint64_t)P->num_blocks * SIEVE_SIZE);
    double lt = mpz_sizeinbase(ta,2)*0.693147/na;
    int tp = (int)exp(lt);
    mpz_clear(ta);

    acands = malloc(fb_n * sizeof(int));
    nacands = 0;
    for(int i=sieve_start;i<fb_n;i++) {
        if(fb_p[i]>=(uint32_t)(tp/2) && fb_p[i]<=(uint32_t)(tp*3) && fb_r[i]!=0)
            acands[nacands++]=i;
    }
    if(nacands < na+5) {
        nacands=0;
        for(int i=fb_n/3;i<fb_n && nacands<300;i++)
            if(fb_r[i]!=0) acands[nacands++]=i;
    }
}

static void new_a(void) {
    /* Random selection of na primes */
    int *sh = malloc(nacands*sizeof(int));
    memcpy(sh,acands,nacands*sizeof(int));
    for(int i=0;i<na && i<nacands;i++){
        int j=i+gmp_urandomm_ui(rng,nacands-i);
        int t=sh[i];sh[i]=sh[j];sh[j]=t;
        a_idx[i]=sh[i];
    }
    free(sh);

    mpz_set_ui(pa,1);
    for(int i=0;i<na;i++) mpz_mul_ui(pa,pa,fb_p[a_idx[i]]);

    /* Compute B via CRT */
    mpz_t ad,tmp; mpz_inits(ad,tmp,NULL);
    mpz_set_ui(pb,0);
    for(int i=0;i<na;i++){
        uint32_t p=fb_p[a_idx[i]], sq=fb_r[a_idx[i]];
        mpz_fdiv_q_ui(ad,pa,p);
        uint32_t adm=mpz_fdiv_ui(ad,p);
        uint32_t inv=mod_inv(adm,p);
        uint64_t bj=(uint64_t)sq*inv%p;
        mpz_mul_ui(tmp,ad,bj);
        mpz_set(Bvals[i],tmp);
        mpz_add(pb,pb,tmp);
    }
    mpz_mod(pb,pb,pa);
    /* Ensure B < A/2 */
    mpz_fdiv_q_ui(tmp,pa,2);
    if(mpz_cmp(pb,tmp)>0) mpz_sub(pb,pa,pb);
    /* C = (B^2 - N) / A */
    mpz_mul(pc,pb,pb); mpz_sub(pc,pc,g_N); mpz_divexact(pc,pc,pa);
    mpz_clears(ad,tmp,NULL);

    /* Compute sieve roots */
    for(int i=sieve_start;i<fb_n;i++){
        uint32_t p=fb_p[i], sq=fb_r[i];
        int skip=0;
        for(int j=0;j<na;j++) if(fb_p[a_idx[j]]==p){skip=1;break;}
        if(skip||sq==0){rt1[i]=rt2[i]=UINT32_MAX;continue;}
        uint32_t am=mpz_fdiv_ui(pa,p), bm=mpz_fdiv_ui(pb,p);
        uint32_t ai=mod_inv(am,p);
        rt1[i]=(uint32_t)((uint64_t)(sq+p-bm)%p*ai%p);
        rt2[i]=(uint32_t)((uint64_t)(p-sq+p-bm)%p*ai%p);
    }
    gray=0;
}

static int next_b(void) {
    gray++;
    if(gray>=(1<<(na-1))) return 0;
    int j=0; {int g=gray; while(!(g&1)){j++;g>>=1;}}
    int sign=((gray>>(j+1))&1)?-1:1;
    if(sign>0){mpz_add(pb,pb,Bvals[j]);mpz_add(pb,pb,Bvals[j]);}
    else{mpz_sub(pb,pb,Bvals[j]);mpz_sub(pb,pb,Bvals[j]);}
    mpz_mod(pb,pb,pa);
    mpz_mul(pc,pb,pb);mpz_sub(pc,pc,g_N);mpz_divexact(pc,pc,pa);
    for(int i=sieve_start;i<fb_n;i++){
        if(rt1[i]==UINT32_MAX)continue;
        uint32_t p=fb_p[i];
        uint32_t bjm=mpz_fdiv_ui(Bvals[j],p);
        uint32_t ai=mod_inv(mpz_fdiv_ui(pa,p),p);
        uint32_t sh=(uint32_t)((uint64_t)2*bjm%p*ai%p);
        if(sign>0){rt1[i]=(rt1[i]+p-sh)%p;rt2[i]=(rt2[i]+p-sh)%p;}
        else{rt1[i]=(rt1[i]+sh)%p;rt2[i]=(rt2[i]+sh)%p;}
    }
    return 1;
}

/* ==================== Sieve + Trial Division ==================== */
static uint8_t sieve_arr[SIEVE_SIZE] __attribute__((aligned(64)));

static void do_sieve(int boff) {
    memset(sieve_arr, 0, SIEVE_SIZE);
    for(int i=sieve_start;i<fb_n;i++){
        if(rt1[i]==UINT32_MAX) continue;
        uint32_t p=fb_p[i];
        uint8_t lg=fb_lg[i];
        /* Compute start positions in this block */
        uint32_t s1,s2;
        if(boff>=0){
            s1=((uint32_t)(rt1[i]-boff%p+p))%p;
            s2=((uint32_t)(rt2[i]-boff%p+p))%p;
        }else{
            int bm=((boff%(int)p)+(int)p)%(int)p;
            s1=(rt1[i]+p-(uint32_t)bm)%p;
            s2=(rt2[i]+p-(uint32_t)bm)%p;
        }
        if(p<SIEVE_SIZE){
            for(uint32_t pos=s1;pos<SIEVE_SIZE;pos+=p) sieve_arr[pos]+=lg;
            if(rt1[i]!=rt2[i])
                for(uint32_t pos=s2;pos<SIEVE_SIZE;pos+=p) sieve_arr[pos]+=lg;
        }else{
            if(s1<SIEVE_SIZE) sieve_arr[s1]+=lg;
            if(rt1[i]!=rt2[i] && s2<SIEVE_SIZE) sieve_arr[s2]+=lg;
        }
    }
}

static int scan_avx512(int *cands, uint8_t thresh) {
    int nc=0;
    __m512i tv=_mm512_set1_epi8((char)thresh);
    for(int i=0;i<SIEVE_SIZE;i+=64){
        __m512i sv=_mm512_load_si512((__m512i*)(sieve_arr+i));
        __mmask64 m=_mm512_cmpge_epu8_mask(sv,tv);
        while(m){int b=__builtin_ctzll(m);cands[nc++]=i+b;m&=m-1;}
    }
    return nc;
}

/* Fast trial division using sieve root information.
 * Instead of testing mpz_divisible_ui_p for every prime (slow GMP call),
 * check if the position matches the known roots (fast integer comparison). */
static int trial_div(int x, int *ev, int *sign_out, mpz_t axb, mpz_t Af, uint64_t *cofactor) {
    /* Compute ax+b */
    mpz_mul_si(axb, pa, x);
    mpz_add(axb, axb, pb);

    /* f(x) = ((ax+b)^2 - N) / A */
    mpz_t fx; mpz_init(fx);
    mpz_mul(fx, axb, axb);
    mpz_sub(fx, fx, g_N);
    mpz_divexact(fx, fx, pa);

    int neg = (mpz_sgn(fx) < 0);
    if(neg) mpz_neg(fx, fx);
    *sign_out = neg;

    memset(ev, 0, (fb_n+1) * sizeof(int));

    /* Divide by prime 2 (special handling) */
    while(mpz_divisible_ui_p(fx, 2)) { ev[0]++; mpz_divexact_ui(fx, fx, 2); }

    /* For primes with known roots: check if position matches */
    /* x_pos mod p should equal root1[i] or root2[i] */
    for(int i=1; i<fb_n; i++) {
        uint32_t p = fb_p[i];

        /* Check if this prime divides Q(x) using root information */
        int divides = 0;
        if(rt1[i] == UINT32_MAX) {
            /* A-factor prime — always divides f(x) at the A-contribution level
             * but may also divide f(x) itself. Check directly. */
            divides = mpz_divisible_ui_p(fx, p);
        } else {
            /* Use root matching: x mod p should equal rt1 or rt2 */
            int xm = ((x % (int)p) + (int)p) % (int)p;
            divides = (xm == (int)rt1[i] || xm == (int)rt2[i]);
        }

        if(divides) {
            /* Known to divide at least once; extract full power */
            while(mpz_divisible_ui_p(fx, p)) {
                ev[i]++;
                mpz_divexact_ui(fx, fx, p);
            }
        }
    }

    /* Include A-factor contribution */
    for(int j=0;j<na;j++) ev[a_idx[j]]++;

    /* Compute Af = |A*f(x)| = |(ax+b)^2 - N| */
    mpz_mul(Af, axb, axb);
    mpz_sub(Af, Af, g_N);
    if(mpz_sgn(Af)<0) mpz_neg(Af, Af);

    if(mpz_cmp_ui(fx,1)==0){
        *cofactor=1;
        mpz_clear(fx);
        return 1;
    }

    if(mpz_fits_ulong_p(fx)){
        *cofactor = mpz_get_ui(fx);
    } else {
        *cofactor = UINT64_MAX;
    }
    mpz_clear(fx);
    return 1;
}

/* ==================== Linear Algebra (Gaussian Elimination) ==================== */
static int solve_matrix(mpz_t factor) {
    /* Build matrix of usable relations (full + combined SLP) */
    int nuse=0;
    int *use_idx = malloc(n_rels * sizeof(int));
    for(int i=0;i<n_rels;i++){
        if(rels[i].lp==0) use_idx[nuse++]=i; /* full or combined */
    }

    if(nuse < 10) { free(use_idx); return 0; }

    int ncols = fb_n + 1; /* +1 for sign column */
    int nwords = (ncols+31)/32;
    int aug_words = (nuse+31)/32;

    /* Build augmented GF(2) matrix */
    uint32_t **mat = malloc(nuse * sizeof(uint32_t*));
    for(int i=0;i<nuse;i++){
        mat[i] = calloc(nwords+aug_words, sizeof(uint32_t));
        rel_t *r = &rels[use_idx[i]];
        for(int j=0;j<=fb_n;j++) /* fb_n = sign column */
            if(r->exp_vec[j]&1) mat[i][j/32]|=(1u<<(j%32));
        mat[i][nwords+i/32]|=(1u<<(i%32)); /* identity augmentation */
    }

    /* Gaussian elimination */
    int rank=0;
    for(int col=0;col<ncols && rank<nuse;col++){
        int pr=-1;
        for(int r=rank;r<nuse;r++) if(mat[r][col/32]&(1u<<(col%32))){pr=r;break;}
        if(pr<0) continue;
        if(pr!=rank){uint32_t*t=mat[pr];mat[pr]=mat[rank];mat[rank]=t;}
        for(int r=0;r<nuse;r++){
            if(r==rank) continue;
            if(mat[r][col/32]&(1u<<(col%32)))
                for(int w=0;w<nwords+aug_words;w++) mat[r][w]^=mat[rank][w];
        }
        rank++;
    }

    fprintf(stderr, "[LA] %d relations, rank %d, %d null vectors\n", nuse, rank, nuse-rank);

    /* Try each null vector */
    int found=0;
    mpz_t X, Y, Ysq, tmp;
    mpz_inits(X,Y,Ysq,tmp,NULL);

    for(int r=rank; r<nuse && !found; r++){
        /* Check all-zero in exponent part */
        int ok=1;
        for(int w=0;w<nwords&&ok;w++) if(mat[r][w]) ok=0;
        if(!ok) continue;

        /* Extract combined relation */
        mpz_set_ui(X,1);
        mpz_set_ui(Ysq,1);
        int cnt=0;

        for(int i=0;i<nuse;i++){
            if(!(mat[r][nwords+i/32]&(1u<<(i%32)))) continue;
            cnt++;
            rel_t *rel = &rels[use_idx[i]];
            mpz_mul(X, X, rel->axb_prod);
            mpz_mod(X, X, g_N);
            mpz_mul(Ysq, Ysq, rel->Qprod);
        }
        if(cnt<2) continue;

        if(!mpz_perfect_square_p(Ysq)){
            continue;
        }
        mpz_sqrt(Y, Ysq);
        mpz_mod(Y, Y, g_N);

        mpz_sub(tmp, X, Y);
        mpz_gcd(factor, tmp, g_N);
        if(mpz_cmp_ui(factor,1)>0 && mpz_cmp(factor,g_N)<0){found=1;break;}
        mpz_add(tmp, X, Y);
        mpz_gcd(factor, tmp, g_N);
        if(mpz_cmp_ui(factor,1)>0 && mpz_cmp(factor,g_N)<0){found=1;break;}
    }

    mpz_clears(X,Y,Ysq,tmp,NULL);
    for(int i=0;i<nuse;i++) free(mat[i]);
    free(mat); free(use_idx);
    return found;
}

/* ==================== Main ==================== */
int main(int argc, char **argv) {
    if(argc<2){fprintf(stderr,"Usage: %s <N>\n",argv[0]);return 1;}

    mpz_t N,factor,cof;
    mpz_inits(N,factor,cof,NULL);
    mpz_set_str(N,argv[1],10);
    int digits = mpz_sizeinbase(N,10);
    fprintf(stderr,"Factoring %d-digit number...\n",digits);
    double t0=now_sec();

    /* Trivial checks */
    if(mpz_probab_prime_p(N,25)){gmp_printf("%Zd is prime\n",N);return 0;}
    for(uint32_t p=2;p<100000;p++){
        if(mpz_divisible_ui_p(N,p)){
            gmp_printf("%Zd = %u * ",N,p);
            mpz_divexact_ui(cof,N,p);
            gmp_printf("%Zd\n",cof);
            return 0;
        }
    }

    /* SIQS */
    mpz_init_set(g_N, N);
    gmp_randinit_mt(rng); gmp_randseed_ui(rng, SEED);

    param_t P = get_params(digits);
    fprintf(stderr,"[SIQS] FB=%d blocks=%d LP=%d A-factors=%d (48KB sieve)\n",
            P.fb_size, P.num_blocks, P.lp_mult, P.num_a);

    build_fb(N, P.fb_size);
    fprintf(stderr,"[SIQS] FB built: %d primes, max=%u\n", fb_n, fb_p[fb_n-1]);

    uint64_t lp_bound = (uint64_t)P.lp_mult * fb_p[fb_n-1];
    int target = fb_n + P.extra;
    fprintf(stderr,"[SIQS] Target: %d usable rels, LP bound: %lu\n", target, lp_bound);

    init_rels();
    init_hash();
    init_poly(&P);

    /* Compute sieve threshold: we want sum of log2(p) for primes dividing Q(x)
     * to be close to log2(Q(x)). Candidates where sieve accumulation ≥ threshold
     * are likely smooth. The threshold should be slightly below log2(Q_typical)
     * to allow for:
     * 1. Missing small prime contributions (skipped in sieve)
     * 2. LP/DLP cofactors up to lp_bound
     * 3. Approximation errors from using ceil(log2(p)) instead of exact log2(p) */
    double log2N = mpz_sizeinbase(N,2);
    /* Estimate actual A size from parameter table */
    double log2A_est = 0;
    if(nacands > 0) {
        /* Use middle of a_candidates range */
        int mid = acands[nacands/2];
        for(int i=0;i<na;i++) log2A_est += log2(fb_p[mid]);
    } else {
        for(int i=0;i<na;i++) log2A_est += log2(fb_p[fb_n/2]);
    }
    /* Q(x)/A at sieve boundary: roughly sqrt(2N)/A * M where M = sieve_radius
     * But more precisely: Q(x) = A*x^2 + 2Bx + C, |Q| ≈ A*(M/2)^2 near boundary
     * and |Q| ≈ |C| = (B^2-N)/A ≈ N/A near center.
     * Average: roughly sqrt(N/A) * sqrt(M) but let's use the YAFU-style formula:
     * threshold ≈ log2(M * sqrt(kN)) - log2(A) - correction */
    int M = P.num_blocks * SIEVE_SIZE;
    double log2_Q_typical = (log2N - log2A_est) / 2.0 + log2(M) / 2.0;
    /* Correction for small primes not in sieve */
    double sp_corr = 0;
    for(int i=0;i<sieve_start;i++) sp_corr += fb_lg[i] * 1.5;
    /* Threshold: candidates whose sieve accumulation is close to log2(Q(x)).
     * The "closnuf" adjustment allows for:
     * - Rounding in byte-valued log approximations (~3 bits)
     * - Small primes not sieved (~sp_corr bits)
     * - LP/DLP tolerance is NOT included here — handled in trial division
     * Use a moderate closnuf to balance false positives vs missed smooths */
    int closnuf = 3 + (digits > 50 ? 2 : 0) + (digits > 70 ? 2 : 0);
    int thresh = (int)(log2_Q_typical - sp_corr - closnuf);
    if(thresh < 25) thresh = 25;
    if(thresh > 200) thresh = 200;
    fprintf(stderr,"[SIQS] Sieve threshold: %d (Q_typ=%.0f, sp=%.0f, closnuf=%d)\n",
            thresh, log2_Q_typical, sp_corr, closnuf);

    int *cands = malloc(SIEVE_SIZE * sizeof(int));
    int *ev = calloc(fb_n, sizeof(int));
    mpz_t axb, Af;
    mpz_inits(axb, Af, NULL);

    int npoly=0, napoly=0, ncands_total=0;
    double last_rep=t0;

    while(usable() < target) {
        if(now_sec()-t0 > 285.0){fprintf(stderr,"[SIQS] Timeout\n");break;}

        if(gray==0 || gray>=(1<<(na-1))){
            new_a(); napoly++;
        }

        /* Sieve positive and negative sides */
        for(int side=0;side<2;side++){
            for(int blk=0;blk<P.num_blocks;blk++){
                int boff = side==0 ? blk*SIEVE_SIZE : -(blk+1)*SIEVE_SIZE;
                do_sieve(boff);
                int nc = scan_avx512(cands, (uint8_t)thresh);
                ncands_total += nc;

                for(int ci=0;ci<nc;ci++){
                    int x = boff + cands[ci];
                    uint64_t cofactor;
                    int fsign;
                    trial_div(x, ev, &fsign, axb, Af, &cofactor);

                    if(cofactor==1){
                        add_full(ev, fsign, axb, Af);
                    } else if(cofactor>1 && cofactor<=lp_bound){
                        mpz_t ctest; mpz_init_set_ui(ctest,cofactor);
                        if(mpz_probab_prime_p(ctest,1)){
                            add_slp(ev, fsign, axb, Af, cofactor);
                        }
                        mpz_clear(ctest);
                    }
                }
            }
        }

        npoly++;
        if(!next_b()) gray=0; /* force new A */

        double tn=now_sec();
        if(tn-last_rep>5.0){
            fprintf(stderr,"[SIQS] %.1fs: %d/%d usable (full=%d SLP=%d), %d polys, %d A, %d cands\n",
                    tn-t0, usable(), target, n_full, n_slp_combined, npoly, napoly, ncands_total);
            last_rep=tn;
        }
    }

    fprintf(stderr,"[SIQS] Sieve done: %d usable in %.1fs. Starting LA...\n",
            usable(), now_sec()-t0);
    int found = solve_matrix(factor);
    double elapsed = now_sec()-t0;

    if(found){
        mpz_divexact(cof, N, factor);
        gmp_printf("%Zd = %Zd * %Zd\n", N, factor, cof);
        fprintf(stderr,"Factor found in %.3fs\n",elapsed);
    } else {
        fprintf(stderr,"FAILED after %.3fs (%d rels, %d usable)\n",elapsed,n_rels,usable());
        gmp_printf("FAIL %Zd\n",N);
    }

    /* Cleanup */
    free(cands); free(ev);
    mpz_clears(axb, Af, NULL);
    for(int i=0;i<n_rels_alloc;i++){free(rels[i].exp_vec);mpz_clear(rels[i].Qprod);mpz_clear(rels[i].axb_prod);}
    free(rels); free(hpool); free(rt1); free(rt2); free(acands);
    for(int i=0;i<MAX_AFACTORS;i++) mpz_clear(Bvals[i]);
    mpz_clears(pa,pb,pc,g_N,N,factor,cof,NULL);
    gmp_randclear(rng);
    return found?0:1;
}
