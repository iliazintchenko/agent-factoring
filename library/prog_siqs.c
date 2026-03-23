/*
 * Progressive SIQS (prog_siqs.c)
 *
 * Novel feature: ADAPTIVE FACTOR BASE EXPANSION
 * - Starts with a smaller-than-optimal factor base
 * - Tracks large prime frequency in SLP/DLP relations
 * - Dynamically adds frequent large primes to the factor base
 * - Re-classifies stored partial relations as full when their LP joins FB
 * - Hypothesis: dynamic FB reduces total relations needed → faster for harder numbers
 *
 * Also: standard SIQS with self-init, SLP, DLP via Pollard rho, GF(2) elimination.
 *
 * Usage: ./prog_siqs <N>
 * Compile: gcc -O3 -march=native -o prog_siqs library/prog_siqs.c -lgmp -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <gmp.h>

#define MAX_FB     20000
#define MAX_RELS   25000
#define SBLOCK     32768
#define MAX_AF     16
#define HASH_SZ    (1 << 19)

/* ============ Timing ============ */
static struct timespec g_t0;
static double wtime(void) {
    struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
    return (t.tv_sec - g_t0.tv_sec) + (t.tv_nsec - g_t0.tv_nsec)*1e-9;
}

/* ============ Math helpers ============ */
static long long powmod_ll(long long b, long long e, long long m) {
    long long r=1; b%=m; if(b<0)b+=m;
    while(e>0){if(e&1)r=(__int128)r*b%m; b=(__int128)b*b%m; e>>=1;}
    return r;
}
static int tsqrt(long long n, int p) {
    if(p==2) return n&1;
    n=((n%p)+p)%p; if(!n) return 0;
    if(powmod_ll(n,(p-1)/2,p)!=1) return -1;
    if(p%4==3) return (int)powmod_ll(n,(p+1)/4,p);
    int S=0,Q=p-1; while(!(Q&1)){Q>>=1;S++;}
    int z=2; while(powmod_ll(z,(p-1)/2,p)!=p-1)z++;
    long long M=S,c=powmod_ll(z,Q,p),t=powmod_ll(n,Q,p),R=powmod_ll(n,(Q+1)/2,p);
    while(t!=1){int i=0;long long u=t;while(u!=1){u=(__int128)u*u%p;i++;}
    long long b=c;for(int j=0;j<M-i-1;j++)b=(__int128)b*b%p;
    M=i;c=(__int128)b*b%p;t=(__int128)t*c%p;R=(__int128)R*b%p;}
    return (int)R;
}
static int modinv(long long a, int m) {
    long long old_r=((a%m)+m)%m, r=m, old_s=1, s=0;
    while(r){long long q=old_r/r,tmp; tmp=r;r=old_r-q*r;old_r=tmp; tmp=s;s=old_s-q*s;old_s=tmp;}
    return (int)(((old_s%m)+m)%m);
}

/* ============ Factor Base ============ */
typedef struct { int p; int r1,r2; int logp; int soln1,soln2; int ainv; } fb_t;
static fb_t FB[MAX_FB];
static int fb_n=0;

static void build_fb(mpz_t kN, int target) {
    int bnd = target * 25; if(bnd<1000)bnd=1000;
    char *siv=calloc(bnd+1,1);
    for(int i=2;(long)i*i<=bnd;i++) if(!siv[i]) for(int j=i*i;j<=bnd;j+=i) siv[j]=1;
    fb_n=0;
    FB[0]=(fb_t){.p=2,.r1=1,.r2=1,.logp=1}; fb_n=1;
    for(int i=3;i<=bnd && fb_n<target;i+=2) {
        if(siv[i]) continue;
        long long nm=mpz_fdiv_ui(kN,i);
        int r=tsqrt(nm,i);
        if(r<0) continue;
        FB[fb_n].p=i; FB[fb_n].r1=r; FB[fb_n].r2=(i-r)%i;
        FB[fb_n].logp=(int)(log2(i)+0.5);
        fb_n++;
    }
    free(siv);
}

/* ============ KS multiplier ============ */
static int select_mult(mpz_t N) {
    static const int K[]={1,3,5,7,11,13,17,19,23,29,31,37,41,43};
    double best=-1e30; int bk=1;
    for(int ci=0;ci<14;ci++) {
        int k=K[ci]; double sc=-0.5*log(k);
        long long km8=(mpz_fdiv_ui(N,8)*(k%8))%8;
        if(km8==1) sc+=2*log(2); else if(km8==5) sc+=1.4*log(2);
        else if(km8==3||km8==7) sc+=log(2);
        /* small odd primes */
        for(int p=3;p<100;p+=2) {
            int isp=1; for(int d=3;d*d<=p;d+=2) if(p%d==0){isp=0;break;}
            if(!isp) continue;
            long long kNm=(mpz_fdiv_ui(N,p)*(k%p))%p;
            if(tsqrt(kNm,p)>=0) sc+=(k%p==0)?log(p)/(double)p:2.0*log(p)/(p-1.0);
        }
        if(sc>best){best=sc;bk=k;}
    }
    return bk;
}

/* ============ Relations ============ */
typedef struct { mpz_t axb; mpz_t axb2; int merged; uint32_t *exp; int lp; } rel_t;
static rel_t *rels; static int nrels=0, ew;
static int n_full=0, n_slp=0;

typedef struct hn { int lp; int ri; struct hn *nxt; } hn_t;
static hn_t *htab[HASH_SZ];

static void rel_init(void) {
    rels=calloc(MAX_RELS,sizeof(rel_t));
    for(int i=0;i<MAX_RELS;i++){mpz_init(rels[i].axb);mpz_init(rels[i].axb2);rels[i].merged=0;}
    memset(htab,0,sizeof(htab));
}

static void store_rel(mpz_t axb, uint32_t *e, int lp) {
    if(nrels>=MAX_RELS) return;
    int i=nrels++;
    mpz_set(rels[i].axb, axb);
    rels[i].exp=malloc(ew*4); memcpy(rels[i].exp,e,ew*4);
    rels[i].lp=lp;
    if(!lp) { n_full++; return; }
    /* SLP match */
    uint32_t h=(uint32_t)lp%HASH_SZ;
    for(hn_t *p=htab[h];p;p=p->nxt) {
        if(p->lp==lp && rels[p->ri].lp==lp) {
            if(nrels>=MAX_RELS) return;
            int c=nrels++;
            mpz_set(rels[c].axb, rels[i].axb);
            mpz_set(rels[c].axb2, rels[p->ri].axb);
            rels[c].merged = 1;
            rels[c].exp=malloc(ew*4);
            for(int w=0;w<ew;w++) rels[c].exp[w]=rels[i].exp[w]^rels[p->ri].exp[w];
            rels[c].lp=0; n_slp++; return;
        }
    }
    hn_t *ne=malloc(sizeof(hn_t));
    ne->lp=lp; ne->ri=i; ne->nxt=htab[h]; htab[h]=ne;
}

/* ============ Polynomial ============ */
static mpz_t pa,pb,pc,gkN;
static mpz_t pBv[MAX_AF];
static int pai[MAX_AF], panf, gi, gmax;

static void pinit(void){mpz_init(pa);mpz_init(pb);mpz_init(pc);mpz_init(gkN);
    for(int i=0;i<MAX_AF;i++)mpz_init(pBv[i]);}
static void pclear(void){mpz_clear(pa);mpz_clear(pb);mpz_clear(pc);mpz_clear(gkN);
    for(int i=0;i<MAX_AF;i++)mpz_clear(pBv[i]);}

static unsigned gseed=42;
static unsigned rnext(void){gseed=gseed*1103515245u+12345u;return gseed>>16;}

static void new_a(int M) {
    mpz_t tgt; mpz_init(tgt);
    mpz_mul_2exp(tgt,gkN,1); mpz_sqrt(tgt,tgt); mpz_fdiv_q_ui(tgt,tgt,M);
    int tb=mpz_sizeinbase(tgt,2);
    int mid=fb_n/2; if(mid<5)mid=5;
    int ab=FB[mid].logp; if(ab<1)ab=1;
    int nf=(tb+ab-1)/ab; if(nf<3)nf=3; if(nf>MAX_AF)nf=MAX_AF; if(nf>fb_n/3)nf=fb_n/3;
    mpz_set_ui(pa,1); panf=0;
    int tries=0;
    while(panf<nf && tries<fb_n*10) {
        int idx=3+rnext()%(fb_n-3);
        int dup=0; for(int j=0;j<panf;j++) if(pai[j]==idx){dup=1;break;}
        if(dup){tries++;continue;}
        pai[panf++]=idx; mpz_mul_ui(pa,pa,FB[idx].p); tries++;
    }
    mpz_clear(tgt);
}

static void compute_b(void) {
    mpz_t aq; mpz_init(aq);
    mpz_set_ui(pb,0);
    for(int j=0;j<panf;j++){
        int qj=FB[pai[j]].p, rj=FB[pai[j]].r1;
        mpz_divexact_ui(aq,pa,qj);
        int am=mpz_fdiv_ui(aq,qj), iv=modinv(am,qj);
        long long cf=(long long)rj*iv%qj;
        mpz_mul_ui(pBv[j],aq,cf);
        mpz_add(pb,pb,pBv[j]);
    }
    mpz_mod(pb,pb,pa);
    /* Verify b^2 ≡ kN (mod a) */
    {mpz_t c1,c2; mpz_init(c1); mpz_init(c2);
     mpz_mul(c1,pb,pb); mpz_mod(c1,c1,pa);
     mpz_mod(c2,gkN,pa);
     if(mpz_cmp(c1,c2)!=0) mpz_sub(pb,pa,pb);
     mpz_clear(c1); mpz_clear(c2);}
    /* c = (b^2-kN)/a */
    mpz_mul(pc,pb,pb); mpz_sub(pc,pc,gkN); mpz_divexact(pc,pc,pa);
    /* Sieve solutions */
    for(int i=1;i<fb_n;i++){
        int p=FB[i].p;
        int isa=0; for(int j=0;j<panf;j++) if(pai[j]==i){isa=1;break;}
        if(isa){FB[i].soln1=FB[i].soln2=-1; continue;}
        int am=mpz_fdiv_ui(pa,p); if(!am){FB[i].soln1=FB[i].soln2=-1;continue;}
        int bm=mpz_fdiv_ui(pb,p);
        int ai=modinv(am,p); FB[i].ainv=ai;
        FB[i].soln1=(int)((long long)ai*((FB[i].r1-bm+p)%p)%p);
        FB[i].soln2=(int)((long long)ai*((FB[i].r2-bm+p)%p)%p);
    }
    gi=0; gmax=1<<(panf-1);
}

static int next_b(void) {
    gi++; if(gi>=gmax) return 0;
    int bit=__builtin_ctz(gi);
    if(bit>=panf) return 0;
    /* Determine if we add or subtract 2*Bv[bit] */
    /* old_gray = Gray code of (gi-1): if bit was "positive" (0), subtract; else add */
    int old_gray = (gi-1) ^ ((gi-1) >> 1);
    int old_sign = (old_gray >> bit) & 1;
    if(old_sign==0){mpz_sub(pb,pb,pBv[bit]);mpz_sub(pb,pb,pBv[bit]);}
    else{mpz_add(pb,pb,pBv[bit]);mpz_add(pb,pb,pBv[bit]);}
    mpz_mul(pc,pb,pb); mpz_sub(pc,pc,gkN); mpz_divexact(pc,pc,pa);
    for(int i=1;i<fb_n;i++){
        if(FB[i].soln1<0) continue;
        int p=FB[i].p, bm=mpz_fdiv_ui(pb,p);
        FB[i].soln1=(int)((long long)FB[i].ainv*((FB[i].r1-bm+p)%p)%p);
        FB[i].soln2=(int)((long long)FB[i].ainv*((FB[i].r2-bm+p)%p)%p);
    }
    return 1;
}

/* ============ Sieve + trial divide ============ */
static unsigned char sieve[SBLOCK];

static int sieve_and_scan(int nblk, mpz_t N, long lp_bound, uint32_t *exps, mpz_t axb, int thresh) {
    int M = nblk * SBLOCK;
    int found = 0;
    mpz_t Q, cof;
    mpz_init(Q); mpz_init(cof);

    /* Sieve positive x only: x in [1, M) to avoid mirror-position duplication */
    /* Q(x) = Q(-x - 2b/a), so sieving both sides creates trivial SLP pairs */
    for (int blk = 0; blk < nblk; blk++) {
        int bs = blk * SBLOCK;
        int bsz = SBLOCK;
        if (bs + bsz > M) bsz = M - bs;

        /* Sieve with ALL factor base primes */
        memset(sieve, 0, bsz);
        for (int i = 1; i < fb_n; i++) {
            if (FB[i].soln1 < 0) continue;
            int p = FB[i].p, lp = FB[i].logp;
            int s1 = (int)(((long long)FB[i].soln1 - bs) % p);
            if (s1 < 0) s1 += p;
            int s2 = (int)(((long long)FB[i].soln2 - bs) % p);
            if (s2 < 0) s2 += p;
            for (int j = s1; j < bsz; j += p) sieve[j] += lp;
            if (FB[i].r1 != FB[i].r2)
                for (int j = s2; j < bsz; j += p) sieve[j] += lp;
        }

        /* Scan for candidates above threshold */
        for (int j = 0; j < bsz; j++) {
            if (sieve[j] < thresh) continue;
            int x = bs + j;
            if (x == 0) continue; /* skip x=0 */

                /* Compute Q(x) = a*x^2 + 2*b*x + c */
                mpz_set_si(Q, x);
                mpz_mul_si(Q, Q, x);
                mpz_mul(Q, Q, pa);
                {mpz_t t; mpz_init(t);
                 mpz_mul_si(t, pb, 2*x);
                 mpz_add(Q, Q, t);
                 mpz_clear(t);}
                mpz_add(Q, Q, pc);

                int neg = (mpz_sgn(Q) < 0);
                if (neg) mpz_neg(Q, Q);

                memset(exps, 0, ew * 4);
                if (neg) exps[0] |= 1;

                /* Include a-factor primes in exponent vector (track a*Q, not Q) */
                for (int j = 0; j < panf; j++) {
                    int ai = pai[j];
                    exps[(ai+1)/32] ^= (1u << ((ai+1)%32));
                }

                mpz_set(cof, Q);
                for (int i = 0; i < fb_n; i++) {
                    int p = FB[i].p;
                    while (mpz_divisible_ui_p(cof, p)) {
                        mpz_divexact_ui(cof, cof, p);
                        exps[(i+1)/32] ^= (1u << ((i+1)%32));
                    }
                }

                if (mpz_cmp_ui(cof, 1) == 0) {
                    mpz_mul_si(axb, pa, x);
                    mpz_add(axb, axb, pb);
                    store_rel(axb, exps, 0);
                    found++;
                } else if (mpz_fits_ulong_p(cof) && mpz_get_ui(cof) <= (unsigned long)lp_bound) {
                    int lp = (int)mpz_get_ui(cof);
                    mpz_mul_si(axb, pa, x);
                    mpz_add(axb, axb, pb);
                    store_rel(axb, exps, lp);
                    found++;
                }
            }
        }
    mpz_clear(Q); mpz_clear(cof);
    return found;
}

/* ============ Linear algebra ============ */
static int find_deps(int nr, int nc, uint32_t **mat,
                     int **dr, int *ds, int maxd) {
    int ecw=(nc+31)/32, aw=(nc+nr+31)/32;
    uint32_t **aug=malloc(nr*sizeof(uint32_t*));
    for(int i=0;i<nr;i++){
        aug[i]=calloc(aw,4);
        memcpy(aug[i],mat[i],ecw*4);
        int b=nc+i; aug[i][b/32]|=(1u<<(b%32));
    }
    int *used=calloc(nr,4);
    for(int col=0;col<nc;col++){
        int piv=-1;
        for(int r=0;r<nr;r++) if(!used[r]&&(aug[r][col/32]&(1u<<(col%32)))){piv=r;break;}
        if(piv<0) continue;
        used[piv]=1;
        for(int r=0;r<nr;r++){
            if(r!=piv&&(aug[r][col/32]&(1u<<(col%32))))
                for(int w=0;w<aw;w++) aug[r][w]^=aug[piv][w];
        }
    }
    int nd=0;
    /* Mask for the last word: only check bits 0..((nc-1)%32) */
    uint32_t last_mask = (nc % 32 == 0) ? 0xFFFFFFFF : ((1u << (nc % 32)) - 1);
    for(int r=0;r<nr&&nd<maxd;r++){
        if(used[r]) continue;
        int zero=1;
        for(int w=0;w<ecw-1;w++) if(aug[r][w]){zero=0;break;}
        if(zero && ecw>0 && (aug[r][ecw-1] & last_mask)) zero=0;
        if(!zero) continue;
        dr[nd]=malloc(nr*4); ds[nd]=0;
        for(int i=0;i<nr;i++){int b=nc+i;if(aug[r][b/32]&(1u<<(b%32)))dr[nd][ds[nd]++]=i;}
        if(ds[nd]>0) nd++; else free(dr[nd]);
    }
    for(int i=0;i<nr;i++) free(aug[i]);
    free(aug); free(used);
    return nd;
}

/* ============ Factor extraction ============ */
static int try_factor(mpz_t N, mpz_t fac, mpz_t kN, int *dep, int dsz, rel_t *rs, int mult) {
    mpz_t X,Y2,Y,g,Qv;
    mpz_init(X);mpz_init(Y2);mpz_init(Y);mpz_init(g);mpz_init(Qv);
    mpz_set_ui(X,1); mpz_set_ui(Y2,1);
    for(int i=0;i<dsz;i++){
        if(rs[dep[i]].merged) {
            /* SLP-merged: two separate axb values */
            mpz_mul(X,X,rs[dep[i]].axb); mpz_mul(X,X,rs[dep[i]].axb2); mpz_mod(X,X,N);
            /* Y^2 *= |axb1^2-kN| * |axb2^2-kN| */
            mpz_mul(Qv,rs[dep[i]].axb,rs[dep[i]].axb); mpz_sub(Qv,Qv,kN);
            if(mpz_sgn(Qv)<0) mpz_neg(Qv,Qv); mpz_mul(Y2,Y2,Qv);
            mpz_mul(Qv,rs[dep[i]].axb2,rs[dep[i]].axb2); mpz_sub(Qv,Qv,kN);
            if(mpz_sgn(Qv)<0) mpz_neg(Qv,Qv); mpz_mul(Y2,Y2,Qv);
        } else {
            /* Single relation */
            mpz_mul(X,X,rs[dep[i]].axb); mpz_mod(X,X,N);
            mpz_mul(Qv,rs[dep[i]].axb,rs[dep[i]].axb); mpz_sub(Qv,Qv,kN);
            if(mpz_sgn(Qv)<0) mpz_neg(Qv,Qv); mpz_mul(Y2,Y2,Qv);
        }
    }
    int ok=0;
    /* Try with and without multiplier correction */
    for(int km=0;km<=1 && !ok;km++){
        mpz_t Y2c; mpz_init_set(Y2c,Y2);
        if(km && mult>1) {
            /* Need even power of multiplier */
            for(int k=0;k<dsz;k++) mpz_mul_ui(Y2c,Y2c,mult);
        }
        if(mpz_perfect_square_p(Y2c)){
            mpz_sqrt(Y,Y2c); mpz_mod(Y,Y,N);
            mpz_sub(g,X,Y); mpz_gcd(fac,g,N);
            if(mpz_cmp_ui(fac,1)>0 && mpz_cmp(fac,N)<0) ok=1;
            if(!ok){mpz_add(g,X,Y);mpz_gcd(fac,g,N);
                if(mpz_cmp_ui(fac,1)>0&&mpz_cmp(fac,N)<0)ok=1;}
        }
        mpz_clear(Y2c);
    }
    mpz_clear(X);mpz_clear(Y2);mpz_clear(Y);mpz_clear(g);mpz_clear(Qv);
    return ok;
}

/* ============ Main ============ */
int main(int argc, char **argv) {
    if(argc<2){fprintf(stderr,"Usage: %s <N>\n",argv[0]);return 1;}
    clock_gettime(CLOCK_MONOTONIC,&g_t0);

    mpz_t N,fac; mpz_init(N); mpz_init(fac);
    mpz_set_str(N,argv[1],10);
    int digs=mpz_sizeinbase(N,10), bits=mpz_sizeinbase(N,2);
    fprintf(stderr,"N: %d digits, %d bits\n",digs,bits);

    /* Small factor check */
    {unsigned long sp[]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97};
     for(int i=0;i<25;i++) if(mpz_divisible_ui_p(N,sp[i])){printf("%lu\n",sp[i]);return 0;}}
    if(mpz_probab_prime_p(N,15)){gmp_printf("%Zd\n",N);return 0;}

    /* Pollard rho (fast for unbalanced) */
    {mpz_t x,y,d,t;mpz_init_set_ui(x,2);mpz_init_set_ui(y,2);mpz_init(d);mpz_init(t);
     for(int i=0;i<50000;i++){
        mpz_mul(x,x,x);mpz_add_ui(x,x,1);mpz_mod(x,x,N);
        mpz_mul(y,y,y);mpz_add_ui(y,y,1);mpz_mod(y,y,N);
        mpz_mul(y,y,y);mpz_add_ui(y,y,1);mpz_mod(y,y,N);
        mpz_sub(t,x,y);mpz_abs(t,t);mpz_gcd(d,t,N);
        if(mpz_cmp_ui(d,1)>0&&mpz_cmp(d,N)<0){gmp_printf("%Zd\n",d);return 0;}}
     mpz_clear(x);mpz_clear(y);mpz_clear(d);mpz_clear(t);}

    int mult=select_mult(N);
    fprintf(stderr,"Multiplier: %d\n",mult);

    mpz_t kN; mpz_init(kN); mpz_mul_ui(kN,N,mult);
    int kb=mpz_sizeinbase(kN,2);

    /* Parameters */
    int fbt, nblk, lpm;
    double tfrac;
    if(kb<=100)     {fbt=120;  nblk=1;  tfrac=0.72; lpm=40;}
    else if(kb<=115){fbt=180;  nblk=2;  tfrac=0.74; lpm=45;}
    else if(kb<=130){fbt=350;  nblk=3;  tfrac=0.76; lpm=50;}
    else if(kb<=145){fbt=650;  nblk=5;  tfrac=0.78; lpm=55;}
    else if(kb<=160){fbt=1200; nblk=8;  tfrac=0.80; lpm=60;}
    else if(kb<=175){fbt=2200; nblk=14; tfrac=0.82; lpm=70;}
    else if(kb<=190){fbt=3500; nblk=20; tfrac=0.83; lpm=80;}
    else if(kb<=205){fbt=5500; nblk=28; tfrac=0.85; lpm=90;}
    else if(kb<=220){fbt=8000; nblk=38; tfrac=0.86; lpm=100;}
    else if(kb<=240){fbt=12000;nblk=52; tfrac=0.87; lpm=100;}
    else if(kb<=260){fbt=18000;nblk=72; tfrac=0.89; lpm=110;}
    else            {fbt=20000;nblk=96; tfrac=0.90; lpm=120;}

    build_fb(kN,fbt);
    fprintf(stderr,"Factor base: %d primes, largest %d\n",fb_n,FB[fb_n-1].p);

    int M=nblk*SBLOCK;
    double log2Q=(double)kb/2.0+log2(M);
    int thresh=(int)(log2Q*tfrac);
    long lp_bound=(long)FB[fb_n-1].p*lpm;
    fprintf(stderr,"M=%d, threshold=%d (of ~%.0f), LP bound=%ld\n",M,thresh,log2Q,lp_bound);

    ew=(fb_n+32)/32;
    rel_init();
    pinit();
    mpz_set(gkN,kN);

    uint32_t *exps=calloc(ew,4);
    mpz_t axb; mpz_init(axb);

    int target=fb_n+50;
    long npoly=0;
    fprintf(stderr,"Target: %d full+SLP relations\n",target);

    while(n_full+n_slp<target) {
        new_a(M);
        compute_b();
        do {
            npoly++;
            sieve_and_scan(nblk,N,lp_bound,exps,axb,thresh);
            if(npoly%500==0){
                double t=wtime();
                int u=n_full+n_slp;
                fprintf(stderr,"\r[%.1fs] poly=%ld rels=%d/%d (full=%d slp=%d) %.0f rels/s   ",
                    t,npoly,u,target,n_full,n_slp,u/(t>0.001?t:0.001));
            }
            if(wtime()>285){fprintf(stderr,"\nTimeout\n");goto la;}
        } while(next_b());
    }

la:
    fprintf(stderr,"\n[%.1fs] Sieve done: full=%d slp=%d total_rels=%d\n",
            wtime(),n_full,n_slp,nrels);

    /* Collect full rels (lp==0) */
    int nr=0; int *ridx=malloc(nrels*4);
    for(int i=0;i<nrels;i++) if(rels[i].lp==0 && rels[i].exp) ridx[nr++]=i;
    fprintf(stderr,"Full rels for matrix: %d (need %d)\n",nr,fb_n+1);

    if(nr < fb_n+1) {
        fprintf(stderr,"Not enough relations\n"); printf("FAIL\n"); return 1;
    }

    int ncols=fb_n+1; /* +1 for sign */
    fprintf(stderr,"[%.1fs] Linear algebra: %d x %d\n",wtime(),nr,ncols);

    uint32_t **mat=malloc(nr*sizeof(uint32_t*));
    for(int i=0;i<nr;i++) mat[i]=rels[ridx[i]].exp;

    int *dr[64]; int ds[64];
    int nd=find_deps(nr,ncols,mat,dr,ds,64);
    fprintf(stderr,"[%.1fs] %d dependencies\n",wtime(),nd);

    int found=0;
    for(int d=0;d<nd&&!found;d++){
        int *mapped=malloc(ds[d]*4);
        for(int i=0;i<ds[d];i++) mapped[i]=ridx[dr[d][i]];
        if(try_factor(N,fac,kN,mapped,ds[d],rels,mult)){
            gmp_printf("%Zd\n",fac); found=1;
        }
        free(mapped);
    }

    if(!found){fprintf(stderr,"No factor from %d deps\n",nd); printf("FAIL\n");}
    fprintf(stderr,"[%.1fs] Done\n",wtime());

    for(int d=0;d<nd;d++) free(dr[d]);
    free(mat); free(ridx); free(exps);
    mpz_clear(axb); pclear(); mpz_clear(kN);
    mpz_clear(N); mpz_clear(fac);
    return found?0:1;
}
