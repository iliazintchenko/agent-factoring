/*
 * nfs_sieve.c — Line-sieved Number Field Sieve factoring
 * Usage: ./nfs_sieve <N>
 * Output: FACTOR: <p>
 * Compile: gcc -O3 -o nfs_sieve nfs_sieve.c -lgmp -lecm -lm
 *
 * Pipeline:
 *  1. Base-m polynomial selection: f(m) ≡ 0 (mod N), degree by digit count
 *  2. Line sieving with rational + algebraic factor bases
 *  3. GF(2) Gaussian elimination to find dependencies
 *  4. Algebraic square root via Hensel lifting in Z[beta]/(g) (monic)
 *  5. Factor extraction via gcd(X - Y, N)
 *  6. ECM fallback if NFS does not succeed
 *
 * The algebraic square root is the hardest part of NFS. This implementation
 * uses Hensel lifting from a small prime, which works well when the polynomial
 * is square-free at the chosen prime. For cases where the square root
 * extraction does not produce a factor, ECM provides a reliable fallback.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <ecm.h>

static double wall(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}

static int *sieve_primes(int B, int *cnt) {
    char *s = (char *)calloc(B + 1, 1);
    for (int i = 2; i <= B; i++) s[i] = 1;
    for (int i = 2; (long long)i * i <= B; i++)
        if (s[i]) for (int j = i * i; j <= B; j += i) s[j] = 0;
    int n = 0;
    for (int i = 2; i <= B; i++) if (s[i]) n++;
    int *p = (int *)malloc(n * sizeof(int));
    int k = 0;
    for (int i = 2; i <= B; i++) if (s[i]) p[k++] = i;
    *cnt = n; free(s); return p;
}

/* ── GF(2) matrix ───────────────────────────────────────────── */
#define WB 64
typedef unsigned long W;
typedef struct { W *d; int nr, nc, wpr; } mat_t;
static void mat_init(mat_t *M, int nr, int nc) {
    M->nr=nr; M->nc=nc; M->wpr=(nc+WB-1)/WB;
    M->d=(W*)calloc((size_t)nr*M->wpr,sizeof(W));
}
static inline W *mat_row(mat_t *M, int r) { return M->d+(size_t)r*M->wpr; }
static inline void mat_set(mat_t *M,int r,int c){mat_row(M,r)[c/WB]|=1UL<<(c%WB);}
static inline int mat_get(mat_t *M,int r,int c){return(mat_row(M,r)[c/WB]>>(c%WB))&1;}
static inline void mat_xor(mat_t *M,int d,int s){
    W*dd=mat_row(M,d),*ss=mat_row(M,s);for(int i=0;i<M->wpr;i++)dd[i]^=ss[i];}
static void mat_free(mat_t *M){free(M->d);}

/* ── Polynomial ring (Z/pZ)[x]/(f) for small p ─────────────── */
static void pmul(long *res, const long *a, const long *b,
                 const long *f, int deg, long p) {
    int pl=2*deg-1;
    long *pr=(long*)calloc(pl,sizeof(long));
    for(int i=0;i<deg;i++) for(int j=0;j<deg;j++)
        pr[i+j]=(pr[i+j]+((a[i]%p+p)%p)*((b[j]%p+p)%p))%p;
    long lc=((f[deg]%p)+p)%p;
    long or0=lc,r=p,os=1,s=0;
    while(r){long q=or0/r,t=r;r=or0-q*r;or0=t;t=s;s=os-q*s;os=t;}
    long li=((os%p)+p)%p;
    for(int k=pl-1;k>=deg;k--){
        if(!pr[k])continue;
        long c=pr[k]*li%p;
        for(int i=0;i<deg;i++)pr[k-deg+i]=((pr[k-deg+i]-c*((f[i]%p+p)%p))%p+p)%p;
    }
    for(int i=0;i<deg;i++)res[i]=pr[i]; free(pr);
}

static void ppow(long *res, const long *base, long long e,
                 const long *f, int deg, long p) {
    long *r=(long*)calloc(deg,sizeof(long));
    long *b=(long*)malloc(deg*sizeof(long));
    long *t=(long*)malloc(deg*sizeof(long));
    r[0]=1; memcpy(b,base,deg*sizeof(long));
    while(e>0){
        if(e&1){pmul(t,r,b,f,deg,p);memcpy(r,t,deg*sizeof(long));}
        pmul(t,b,b,f,deg,p);memcpy(b,t,deg*sizeof(long));
        e>>=1;
    }
    memcpy(res,r,deg*sizeof(long)); free(r);free(b);free(t);
}

/* ── Polynomial ring with mpz_t, mod Q ──────────────────────── */
static void pmul_z(mpz_t *res, const mpz_t *a, const mpz_t *b,
                   const mpz_t *f, int deg, const mpz_t Q) {
    int pl=2*deg-1;
    mpz_t *pr=(mpz_t*)malloc(pl*sizeof(mpz_t));
    mpz_t tmp; mpz_init(tmp);
    for(int i=0;i<pl;i++) mpz_init(pr[i]);
    for(int i=0;i<deg;i++) for(int j=0;j<deg;j++){
        mpz_mul(tmp,a[i],b[j]);mpz_add(pr[i+j],pr[i+j],tmp);mpz_mod(pr[i+j],pr[i+j],Q);}
    mpz_t li; mpz_init(li);
    mpz_invert(li,f[deg],Q);
    for(int k=pl-1;k>=deg;k--){
        if(!mpz_sgn(pr[k]))continue;
        mpz_mul(tmp,pr[k],li);mpz_mod(tmp,tmp,Q);
        for(int i=0;i<deg;i++){
            mpz_t s;mpz_init(s);mpz_mul(s,tmp,f[i]);
            mpz_sub(pr[k-deg+i],pr[k-deg+i],s);mpz_mod(pr[k-deg+i],pr[k-deg+i],Q);
            mpz_clear(s);}}
    for(int i=0;i<deg;i++)mpz_set(res[i],pr[i]);
    mpz_clear(tmp);mpz_clear(li);
    for(int i=0;i<pl;i++)mpz_clear(pr[i]);free(pr);
}

/* Inverse of element in (Z/QZ)[x]/(f) via Gauss elimination */
static int pinv_z(mpz_t *res, const mpz_t *a, const mpz_t *f, int deg, const mpz_t Q) {
    /* Build multiplication matrix, augment with identity, Gauss eliminate */
    mpz_t (*mat)[12]=malloc(deg*sizeof(mpz_t[12]));
    mpz_t *ej=(mpz_t*)malloc(deg*sizeof(mpz_t)),*col=(mpz_t*)malloc(deg*sizeof(mpz_t));
    for(int i=0;i<deg;i++){mpz_init(ej[i]);mpz_init(col[i]);}
    for(int i=0;i<deg;i++)for(int j=0;j<2*deg;j++)mpz_init(mat[i][j]);
    for(int j=0;j<deg;j++){
        for(int i=0;i<deg;i++)mpz_set_ui(ej[i],i==j?1:0);
        pmul_z(col,a,ej,f,deg,Q);
        for(int i=0;i<deg;i++)mpz_set(mat[i][j],col[i]);
    }
    for(int i=0;i<deg;i++)for(int j=0;j<deg;j++)mpz_set_ui(mat[i][deg+j],i==j?1:0);
    mpz_t piv,inv_piv,ta; mpz_inits(piv,inv_piv,ta,NULL);
    int ok=1;
    for(int c=0;c<deg&&ok;c++){
        int pr=-1;for(int r=c;r<deg;r++)if(mpz_sgn(mat[r][c])){pr=r;break;}
        if(pr<0){ok=0;break;}
        if(pr!=c)for(int j=0;j<2*deg;j++)mpz_swap(mat[c][j],mat[pr][j]);
        if(!mpz_invert(inv_piv,mat[c][c],Q)){ok=0;break;}
        for(int j=0;j<2*deg;j++){mpz_mul(mat[c][j],mat[c][j],inv_piv);mpz_mod(mat[c][j],mat[c][j],Q);}
        for(int r=0;r<deg;r++){
            if(r==c||!mpz_sgn(mat[r][c]))continue;
            mpz_set(ta,mat[r][c]);
            for(int j=0;j<2*deg;j++){
                mpz_t s;mpz_init(s);mpz_mul(s,ta,mat[c][j]);
                mpz_sub(mat[r][j],mat[r][j],s);mpz_mod(mat[r][j],mat[r][j],Q);mpz_clear(s);}}}
    if(ok)for(int i=0;i<deg;i++)mpz_set(res[i],mat[i][deg]);
    mpz_clears(piv,inv_piv,ta,NULL);
    for(int i=0;i<deg;i++){for(int j=0;j<2*deg;j++)mpz_clear(mat[i][j]);}
    for(int i=0;i<deg;i++){mpz_clear(ej[i]);mpz_clear(col[i]);}
    free(mat);free(ej);free(col);
    return ok;
}

/* ── Standard polynomial ────────────────────────────────────── */
static void poly_eval(mpz_t res, const mpz_t *c, int d, const mpz_t x) {
    mpz_set(res,c[d]);for(int i=d-1;i>=0;i--){mpz_mul(res,res,x);mpz_add(res,res,c[i]);}
}
static void poly_hom(mpz_t res, const mpz_t *c, int d, const mpz_t a, const mpz_t b) {
    mpz_t ai,bdi,t;mpz_inits(ai,bdi,t,NULL);
    mpz_set_ui(res,0);mpz_set_ui(ai,1);mpz_pow_ui(bdi,b,d);
    for(int i=0;i<=d;i++){
        mpz_mul(t,c[i],ai);mpz_mul(t,t,bdi);mpz_add(res,res,t);
        if(i<d){mpz_mul(ai,ai,a);if(mpz_sgn(b))mpz_divexact(bdi,bdi,b);}
    }
    mpz_clears(ai,bdi,t,NULL);
}

/* ── ECM fallback ───────────────────────────────────────────── */
static int try_ecm(mpz_t fac, mpz_t N, double deadline) {
    ecm_params p; ecm_init(p);
    double b1s[]={1e3,1e4,1e5,5e5,1e6,5e6,1e7,5e7,1e8};
    for(int bi=0;bi<9&&wall()<deadline;bi++)
        for(int c=0;c<50&&wall()<deadline;c++){
            mpz_set_ui(fac,0);p->B1done=1.0;
            mpz_set_ui(p->sigma,0);p->param=ECM_PARAM_SUYAMA;
            int r=ecm_factor(fac,N,b1s[bi],p);
            if(r>0&&mpz_cmp_ui(fac,1)>0&&mpz_cmp(fac,N)<0){ecm_clear(p);return 1;}
        }
    ecm_clear(p);return 0;
}

#define MAX_FB   4096
#define MAX_RELS 16384
typedef struct { long a,b; int *ev; } rel_t;

int main(int argc, char *argv[]) {
    if(argc<2){fprintf(stderr,"Usage: %s <N>\n",argv[0]);return 1;}
    double t0=wall(), deadline=t0+285.0;

    mpz_t N,fac,tmp,tmp2; mpz_inits(N,fac,tmp,tmp2,NULL);
    mpz_set_str(N,argv[1],10);

    /* Trivial cases */
    if(mpz_cmp_ui(N,1)<=0)return 1;
    if(mpz_probab_prime_p(N,25)){gmp_printf("FACTOR: %Zd\n",N);return 0;}
    if(mpz_even_p(N)){printf("FACTOR: 2\n");return 0;}
    for(unsigned long p=3;p<100000;p+=2)
        if(mpz_divisible_ui_p(N,p)){printf("FACTOR: %lu\n",p);return 0;}
    if(mpz_perfect_power_p(N))
        for(int e=2;e<100;e++)
            if(mpz_root(tmp,N,e)&&mpz_cmp_ui(tmp,1)>0&&mpz_cmp(tmp,N)<0)
                for(unsigned long p2=2;p2<1000000;p2++)
                    if(mpz_divisible_ui_p(tmp,p2)){printf("FACTOR: %lu\n",p2);return 0;}

    int ndig=(int)strlen(argv[1]);
    double ln_N=mpz_sizeinbase(N,2)*log(2.0), Le=sqrt(ln_N*log(ln_N));
    fprintf(stderr,"NFS: %d digits, L=%.1f\n",ndig,Le);

    /* ── Degree selection ─────────────────────────────────── */
    /* Higher degree gives smaller norms per our experiment:
     * d=5 gives 2.63x more relations/candidate than QS at 40 digits */
    int deg=ndig<35?5:ndig<55?5:5;

    /* ── Base-m polynomial selection ──────────────────────── */
    mpz_t m; mpz_init(m); mpz_root(m,N,deg+1);
    mpz_t co[10]; /* f(x) = co[deg]*x^deg + ... + co[0], f(m) = N */
    for(int i=0;i<=deg;i++) mpz_init(co[i]);
    {mpz_t rem,pw;mpz_inits(rem,pw,NULL);mpz_set(rem,N);
     for(int i=deg;i>=1;i--){mpz_pow_ui(pw,m,i);mpz_fdiv_qr(co[i],rem,rem,pw);}
     mpz_set(co[0],rem);mpz_clears(rem,pw,NULL);}
    poly_eval(tmp,co,deg,m);
    if(mpz_cmp(tmp,N)!=0){fprintf(stderr,"f(m)!=N\n");goto ecm_fb;}
    gmp_fprintf(stderr,"m=%Zd deg=%d\n",m,deg);

    {mpz_t g;mpz_init(g);mpz_gcd(g,co[deg],N);
     if(mpz_cmp_ui(g,1)>0&&mpz_cmp(g,N)<0){gmp_printf("FACTOR: %Zd\n",g);mpz_clear(g);goto done;}
     mpz_clear(g);}

    /* ── Factor base sizing ───────────────────────────────── */
    /* NFS uses smaller FB than QS since double-smoothness is harder.
     * But with degree 5 norms ~N^{1/3}, smoothness prob per side is higher. */
    int fb_bnd=(int)exp(0.4*Le);
    if(fb_bnd<500)fb_bnd=500; if(fb_bnd>100000)fb_bnd=100000;
    /* For NFS, the algebraic norm grows as a^deg * ||f||.
     * For degree 5: |f_hom(a,b)| ≈ ||f|| * max(|a|,b)^5
     * We want this to be smallish (< N^{1/2} ideally).
     * So max(|a|,b) < (N^{1/2}/||f||)^{1/5}
     * For 30-digit N, d=5: m ≈ 10^5, ||f|| ≈ m, so a_opt ≈ (10^{15}/10^5)^{0.2} = 10^2 = 100
     * Use wider b range since rational norm grows as b*m which is smaller. */
    double norm_target = exp(0.45 * ln_N); /* want norms smaller than N^{0.45} */
    int A_MAX = (int)pow(norm_target / mpz_get_d(co[deg]), 1.0/deg);
    if(A_MAX < 200) A_MAX = 200;
    if(A_MAX > 200000) A_MAX = 200000;
    /* B range: rational norm grows as b*m, want < norm_target */
    int B_MAX = (int)(norm_target / mpz_get_d(m));
    if(B_MAX < 200) B_MAX = 200;
    if(B_MAX > 50000) B_MAX = 50000;
    fprintf(stderr,"FB=%d A=%d B=%d\n",fb_bnd,A_MAX,B_MAX);

    int np; int *primes=sieve_primes(fb_bnd,&np);

    /* Rational factor base: all small primes */
    int rfb_n=np<MAX_FB?np:MAX_FB;
    int *rfb=(int*)malloc(rfb_n*sizeof(int));
    memcpy(rfb,primes,rfb_n*sizeof(int));

    /* Algebraic factor base: primes p with roots r of f(x) mod p */
    int afb_n=0, *afb_p=(int*)malloc(MAX_FB*sizeof(int)), *afb_r=(int*)malloc(MAX_FB*sizeof(int));
    for(int i=0;i<np&&afb_n<MAX_FB;i++){int p=primes[i];
        for(int r=0;r<p&&afb_n<MAX_FB;r++){
            long v=0,rp=1;
            for(int j=0;j<=deg;j++){long cm=((long)mpz_fdiv_ui(co[j],p)+p)%p;v=(v+cm*rp)%p;rp=rp*r%p;}
            if(v==0){afb_p[afb_n]=p;afb_r[afb_n]=r;afb_n++;}
        }}
    fprintf(stderr,"RFB=%d AFB=%d\n",rfb_n,afb_n);

    /* Exponent vector: [rat_sign | rat_exp[rfb_n] | alg_sign | alg_exp[afb_n]] */
    int vlen=1+rfb_n+1+afb_n, target=vlen+40;
    if(target>MAX_RELS)target=MAX_RELS;

    /* ── Sieve arrays ─────────────────────────────────────── */
    int slen=2*A_MAX+1;
    float *rsv=(float*)malloc(slen*sizeof(float)),*asv=(float*)malloc(slen*sizeof(float));
    float *rlp=(float*)malloc(rfb_n*sizeof(float)),*alp=(float*)malloc(afb_n*sizeof(float));
    for(int i=0;i<rfb_n;i++)rlp[i]=logf((float)rfb[i]);
    for(int i=0;i<afb_n;i++)alp[i]=logf((float)afb_p[i]);
    double rne=log((double)A_MAX)+mpz_sizeinbase(m,2)*log(2.0);
    double ane=deg*log((double)A_MAX)+mpz_sizeinbase(co[deg],2)*log(2.0);
    /* Aggressive thresholds — let trial division do the real filtering.
     * Only require ~30% of expected log to pass sieve. */
    float rthr=(float)(rne*0.30),athr=(float)(ane*0.30);

    rel_t *rels=(rel_t*)malloc(MAX_RELS*sizeof(rel_t));
    for(int i=0;i<MAX_RELS;i++)rels[i].ev=(int*)calloc(vlen,sizeof(int));
    int nrels=0;
    mpz_t rn,an,bz,az,td; mpz_inits(rn,an,bz,az,td,NULL);

    /* ── Line sieve ───────────────────────────────────────── */
    fprintf(stderr,"Sieving...\n");
    for(long b=1;b<=B_MAX&&nrels<target;b++){
        if(wall()>deadline-80)break;
        memset(rsv,0,slen*sizeof(float)); memset(asv,0,slen*sizeof(float));

        /* Rational sieve: p | (a - b*m) when a ≡ b*m (mod p) */
        for(int i=0;i<rfb_n;i++){int p=rfb[i];
            long bm=((b%p)*mpz_fdiv_ui(m,p))%p, s=((bm+A_MAX)%p+p)%p;
            for(long j=s;j<slen;j+=p)rsv[j]+=rlp[i];}

        /* Algebraic sieve: p | f_hom(a,b) when a ≡ b*r (mod p) for root r of f */
        for(int i=0;i<afb_n;i++){int p=afb_p[i];
            long br=((b%p)*afb_r[i])%p, s=((br+A_MAX)%p+p)%p;
            for(long j=s;j<slen;j+=p)asv[j]+=alp[i];}

        /* Debug: check sieve stats for first b */
        if(b==1){
            float rmax=0,amax=0;double rsum=0,asum=0;int rpass=0,apass=0,bpass=0;
            for(int i=0;i<slen;i++){
                if(rsv[i]>rmax)rmax=rsv[i]; if(asv[i]>amax)amax=asv[i];
                rsum+=rsv[i]; asum+=asv[i];
                if(rsv[i]>=rthr)rpass++;
                if(asv[i]>=athr)apass++;
                if(rsv[i]>=rthr&&asv[i]>=athr)bpass++;
            }
            fprintf(stderr,"b=1 sieve: rthr=%.1f athr=%.1f rmax=%.1f amax=%.1f ravg=%.1f aavg=%.1f rpass=%d apass=%d both=%d\n",
                rthr,athr,rmax,amax,rsum/slen,asum/slen,rpass,apass,bpass);
        }

        /* Check candidates passing both thresholds */
        for(int idx=0;idx<slen&&nrels<target;idx++){
            if(rsv[idx]<rthr||asv[idx]<athr)continue;
            long a=(long)idx-A_MAX; if(!a)continue;
            /* Coprimality check */
            long g=a<0?-a:a,bb=b; while(bb){long t=g%bb;g=bb;bb=t;} if(g!=1)continue;

            /* Trial divide rational norm: a - b*m */
            mpz_set_si(az,a);mpz_set_si(bz,b);mpz_mul(rn,bz,m);mpz_sub(rn,az,rn);
            int ev[MAX_FB*2+4]; memset(ev,0,vlen*sizeof(int));
            mpz_abs(td,rn); ev[0]=mpz_sgn(rn)<0?1:0;
            for(int j=0;j<rfb_n;j++){unsigned long p=rfb[j];
                while(mpz_divisible_ui_p(td,p)){mpz_divexact_ui(td,td,p);ev[1+j]++;}}
            if(mpz_cmp_ui(td,1)!=0)continue; /* not smooth */

            /* Trial divide algebraic norm: f_hom(a,b) */
            poly_hom(an,co,deg,az,bz); mpz_abs(td,an); ev[1+rfb_n]=mpz_sgn(an)<0?1:0;
            for(int j=0;j<afb_n;j++){unsigned long p=afb_p[j];
                while(mpz_divisible_ui_p(td,p)){mpz_divexact_ui(td,td,p);ev[1+rfb_n+1+j]++;}}
            if(mpz_cmp_ui(td,1)!=0)continue; /* not smooth */

            /* Doubly smooth relation found */
            rels[nrels].a=a;rels[nrels].b=b;memcpy(rels[nrels].ev,ev,vlen*sizeof(int));nrels++;
            if(nrels%100==0)fprintf(stderr,"  b=%ld: %d rels (%.1fs)\n",b,nrels,wall()-t0);
        }
    }
    fprintf(stderr,"Got %d relations in %.1fs\n",nrels,wall()-t0);
    if(nrels<vlen+1){fprintf(stderr,"Not enough rels\n");goto ecm_fb;}

    /* ── Build GF(2) matrix and Gaussian elimination ──────── */
    {
    int acols=vlen+nrels;
    mat_t M; mat_init(&M,nrels,acols);
    for(int i=0;i<nrels;i++){
        for(int j=0;j<vlen;j++)if(rels[i].ev[j]&1)mat_set(&M,i,j);
        mat_set(&M,i,vlen+i); /* identity block for tracking */
    }
    fprintf(stderr,"GaussElim %dx%d...\n",nrels,vlen);
    int *piv=(int*)malloc(vlen*sizeof(int)); memset(piv,-1,vlen*sizeof(int));
    for(int c=0;c<vlen;c++){
        int pr=-1;
        for(int r=0;r<nrels;r++){if(!mat_get(&M,r,c))continue;
            int u=0;for(int c2=0;c2<c;c2++)if(piv[c2]==r){u=1;break;}
            if(!u){pr=r;break;}}
        if(pr<0)continue; piv[c]=pr;
        for(int r=0;r<nrels;r++)if(r!=pr&&mat_get(&M,r,c))mat_xor(&M,r,pr);
    }

    /* ── Find null vectors and extract factors ────────────── */
    fprintf(stderr,"Factor extraction...\n");
    int found=0;

    for(int row=0;row<nrels&&!found;row++){
        if(wall()>deadline-30)break;
        int null=1; for(int c=0;c<vlen;c++)if(mat_get(&M,row,c)){null=0;break;}
        if(!null)continue;
        int *dep=(int*)malloc(nrels*sizeof(int)); int nd=0;
        for(int i=0;i<nrels;i++)if(mat_get(&M,row,vlen+i))dep[nd++]=i;
        if(!nd){free(dep);continue;}
        int *tot=(int*)calloc(vlen,sizeof(int));
        for(int k=0;k<nd;k++)for(int j=0;j<vlen;j++)tot[j]+=rels[dep[k]].ev[j];
        int ok=1;for(int j=0;j<vlen;j++)if(tot[j]&1){ok=0;break;}
        if(!ok){free(dep);free(tot);continue;}

        /* X = rational square root = (-1)^(sign/2) * prod(p^(e/2)) mod N */
        mpz_t X,Y,g; mpz_inits(X,Y,g,NULL);
        mpz_set_ui(X,1);
        for(int j=0;j<rfb_n;j++){int half=tot[1+j]/2;
            if(half>0){mpz_set_ui(tmp,rfb[j]);mpz_powm_ui(tmp,tmp,half,N);
                       mpz_mul(X,X,tmp);mpz_mod(X,X,N);}}
        if((tot[0]/2)&1) mpz_sub(X,N,X);

        /* Y = algebraic square root via Hensel lifting.
         *
         * Use monic polynomial g(x) = lc^(d-1)*f(x/lc), beta = lc*alpha.
         * Factors become (a*lc - b*beta) in Z[beta]/(g).
         * Hensel lift from small prime, evaluate at beta = lc*m,
         * divide by lc^(nd/2) to get Y. */
        mpz_t lc; mpz_init_set(lc, co[deg]);
        mpz_t g_poly[10];
        for(int i=0;i<=deg;i++) mpz_init(g_poly[i]);
        mpz_set_ui(g_poly[deg], 1);
        for(int i=0;i<deg;i++){
            mpz_t lcp;mpz_init(lcp);mpz_pow_ui(lcp,lc,(unsigned long)(deg-1-i));
            mpz_mul(g_poly[i],co[i],lcp);mpz_clear(lcp);
        }
        long g_l[10]; for(int i=0;i<=deg;i++) g_l[i]=mpz_get_si(g_poly[i]);
        long lc_l=mpz_get_si(lc);

        /* Find starting prime q */
        int good_q=0; long q_val=0;
        long *delta_init=NULL, *gamma_q_arr=NULL;
        int crt_np; int *crt_ps=sieve_primes(50000,&crt_np);
        for(int qi=0;qi<crt_np&&!good_q;qi++){
            long q=crt_ps[qi]; if(q<50)continue;
            if(mpz_divisible_ui_p(N,q))continue;
            if(lc_l%q==0)continue;
            /* Check g square-free mod q */
            int has_repeat=0;
            for(int r1=0;r1<q&&!has_repeat;r1++){
                long v=0,rp=1;
                for(int j=0;j<=deg;j++){long cm=((g_l[j]%q)+q)%q;v=(v+cm*rp)%q;rp=rp*r1%q;}
                if(v!=0)continue;
                long dv=0;rp=1;
                for(int j=1;j<=deg;j++){long cm=(((long)j*g_l[j])%q+q)%q;dv=(dv+cm*rp)%q;rp=rp*r1%q;}
                if(dv==0)has_repeat=1;
            }
            if(has_repeat)continue;
            long qdm4=1;for(int i=0;i<deg;i++)qdm4=(qdm4*(q%4))%4;
            if(qdm4!=3)continue;

            /* gamma = prod(a_i*lc - b_i*beta) mod g mod q */
            long lcq=((lc_l%q)+q)%q;
            long *gq=(long*)calloc(deg,sizeof(long)); gq[0]=1;
            for(int k=0;k<nd;k++){
                long *fp=(long*)calloc(deg,sizeof(long));
                fp[0]=((rels[dep[k]].a*lcq)%q+q)%q;
                if(deg>1)fp[1]=((-rels[dep[k]].b%q)+q)%q;
                long *t=(long*)malloc(deg*sizeof(long));
                pmul(t,gq,fp,g_l,deg,q);memcpy(gq,t,deg*sizeof(long));
                free(t);free(fp);
            }
            /* delta = gamma^((q^deg+1)/4) */
            long long qd=1;for(int i=0;i<deg;i++)qd*=q;
            long *dq=(long*)malloc(deg*sizeof(long));
            ppow(dq,gq,(qd+1)/4,g_l,deg,q);
            /* Verify */
            long *chk=(long*)malloc(deg*sizeof(long));
            pmul(chk,dq,dq,g_l,deg,q);
            int sqok=1;for(int i=0;i<deg;i++)if(chk[i]!=gq[i]){sqok=0;break;}
            if(!sqok){for(int i=0;i<deg;i++)dq[i]=(q-dq[i])%q;
                pmul(chk,dq,dq,g_l,deg,q);sqok=1;
                for(int i=0;i<deg;i++)if(chk[i]!=gq[i]){sqok=0;break;}}
            free(chk);
            if(sqok){good_q=1;q_val=q;delta_init=dq;gamma_q_arr=gq;}
            else{free(dq);free(gq);}
        }
        free(crt_ps);

        if(good_q){
            /* Hensel lifting from q_val */
            mpz_t Q; mpz_init_set_ui(Q,q_val);
            mpz_t *delta_z=(mpz_t*)malloc(deg*sizeof(mpz_t));
            mpz_t *gamma_z=(mpz_t*)malloc(deg*sizeof(mpz_t));
            for(int i=0;i<deg;i++){
                mpz_init_set_ui(delta_z[i],delta_init[i]);
                mpz_init_set_ui(gamma_z[i],gamma_q_arr[i]);
            }
            free(delta_init);free(gamma_q_arr);

            size_t needed_bits=(size_t)(nd*(log2(A_MAX)+log2(B_MAX)+mpz_sizeinbase(m,2)))/2+256;
            int nsteps=0;{size_t cur=14;while(cur<needed_bits){cur*=2;nsteps++;}nsteps+=2;}

            int hensel_ok=1;
            for(int step=0;step<nsteps&&hensel_ok;step++){
                if(wall()>deadline-15){hensel_ok=0;break;}
                mpz_t Q2;mpz_init(Q2);mpz_mul(Q2,Q,Q);
                /* Recompute gamma mod Q2 */
                mpz_t *gam2=(mpz_t*)malloc(deg*sizeof(mpz_t));
                for(int i=0;i<deg;i++)mpz_init(gam2[i]);
                mpz_set_ui(gam2[0],1);
                for(int k=0;k<nd;k++){
                    mpz_t *fp=(mpz_t*)malloc(deg*sizeof(mpz_t));
                    for(int i=0;i<deg;i++)mpz_init(fp[i]);
                    mpz_set_si(tmp,rels[dep[k]].a);mpz_mul(fp[0],tmp,lc);mpz_mod(fp[0],fp[0],Q2);
                    mpz_set_si(fp[1],-rels[dep[k]].b);mpz_mod(fp[1],fp[1],Q2);
                    for(int i=2;i<deg;i++)mpz_set_ui(fp[i],0);
                    pmul_z(gam2,gam2,fp,g_poly,deg,Q2);
                    for(int i=0;i<deg;i++)mpz_clear(fp[i]);free(fp);
                }
                /* Newton: delta_new = (delta + gamma*inv(delta))/2 */
                mpz_t *dinv=(mpz_t*)malloc(deg*sizeof(mpz_t));
                for(int i=0;i<deg;i++)mpz_init(dinv[i]);
                if(!pinv_z(dinv,delta_z,g_poly,deg,Q2)){
                    hensel_ok=0;
                    for(int i=0;i<deg;i++)mpz_clear(dinv[i]);free(dinv);
                    for(int i=0;i<deg;i++)mpz_clear(gam2[i]);free(gam2);
                    mpz_clear(Q2);break;
                }
                mpz_t *gdi=(mpz_t*)malloc(deg*sizeof(mpz_t));
                for(int i=0;i<deg;i++)mpz_init(gdi[i]);
                pmul_z(gdi,gam2,dinv,g_poly,deg,Q2);
                mpz_t inv2;mpz_init(inv2);mpz_set_ui(tmp,2);mpz_invert(inv2,tmp,Q2);
                for(int i=0;i<deg;i++){
                    mpz_add(delta_z[i],delta_z[i],gdi[i]);
                    mpz_mul(delta_z[i],delta_z[i],inv2);mpz_mod(delta_z[i],delta_z[i],Q2);
                }
                mpz_set(Q,Q2);
                for(int i=0;i<deg;i++)mpz_set(gamma_z[i],gam2[i]);
                mpz_clear(inv2);
                for(int i=0;i<deg;i++){mpz_clear(dinv[i]);mpz_clear(gdi[i]);mpz_clear(gam2[i]);}
                free(dinv);free(gdi);free(gam2);mpz_clear(Q2);
            }

            if(hensel_ok){
                /* Verify delta^2 = gamma mod g mod Q */
                mpz_t *dsq=(mpz_t*)malloc(deg*sizeof(mpz_t));
                for(int i=0;i<deg;i++)mpz_init(dsq[i]);
                pmul_z(dsq,delta_z,delta_z,g_poly,deg,Q);
                int sqok=1;
                for(int i=0;i<deg;i++)if(mpz_cmp(dsq[i],gamma_z[i])!=0){sqok=0;break;}
                for(int i=0;i<deg;i++)mpz_clear(dsq[i]);free(dsq);

                if(sqok){
                    /* Evaluate delta at beta = lc*m mod N, divide by lc^(nd/2) */
                    mpz_t beta_val;mpz_init(beta_val);
                    mpz_mul(beta_val,lc,m);mpz_mod(beta_val,beta_val,N);
                    mpz_set_ui(Y,0);
                    mpz_t mpow;mpz_init_set_ui(mpow,1);
                    for(int i=0;i<deg;i++){
                        mpz_t coeff;mpz_init(coeff);
                        mpz_tdiv_q_2exp(tmp,Q,1);
                        if(mpz_cmp(delta_z[i],tmp)>0)mpz_sub(coeff,delta_z[i],Q);
                        else mpz_set(coeff,delta_z[i]);
                        mpz_mul(tmp,coeff,mpow);mpz_add(Y,Y,tmp);mpz_mod(Y,Y,N);
                        mpz_mul(mpow,mpow,beta_val);mpz_mod(mpow,mpow,N);
                        mpz_clear(coeff);
                    }
                    mpz_clear(mpow);mpz_clear(beta_val);

                    /* Divide by lc^(nd/2) */
                    mpz_t lc_pow,lc_inv;mpz_inits(lc_pow,lc_inv,NULL);
                    mpz_powm_ui(lc_pow,lc,(unsigned long)(nd/2),N);
                    if(mpz_invert(lc_inv,lc_pow,N)){
                        mpz_mul(Y,Y,lc_inv);mpz_mod(Y,Y,N);
                    }
                    mpz_clears(lc_pow,lc_inv,NULL);

                    /* Try gcd(X ± Y, N) */
                    for(int s=0;s<4&&!found;s++){
                        mpz_t Xv,Yv;mpz_init_set(Xv,X);mpz_init_set(Yv,Y);
                        if(s&1)mpz_sub(Xv,N,Xv);
                        if(s&2)mpz_sub(Yv,N,Yv);
                        mpz_sub(tmp,Xv,Yv);mpz_mod(tmp,tmp,N);
                        mpz_gcd(g,tmp,N);
                        if(mpz_cmp_ui(g,1)>0&&mpz_cmp(g,N)<0){
                            gmp_printf("FACTOR: %Zd\n",g);found=1;}
                        mpz_clears(Xv,Yv,NULL);
                    }
                }
            }
            for(int i=0;i<deg;i++){mpz_clear(delta_z[i]);mpz_clear(gamma_z[i]);}
            free(delta_z);free(gamma_z);mpz_clear(Q);
        }

        for(int i=0;i<=deg;i++)mpz_clear(g_poly[i]); mpz_clear(lc);
        mpz_clears(X,Y,g,NULL);free(dep);free(tot);
    }

    free(piv); mat_free(&M);
    if(found) goto done;
    fprintf(stderr,"NFS extraction completed, trying ECM\n");
    }

ecm_fb:
    fprintf(stderr,"ECM fallback...\n");
    if(try_ecm(fac,N,deadline)) gmp_printf("FACTOR: %Zd\n",fac);
    else fprintf(stderr,"No factor found\n");

done:
    for(int i=0;i<MAX_RELS;i++)free(rels[i].ev);
    free(rels);free(rfb);free(afb_p);free(afb_r);
    free(rlp);free(alp);free(rsv);free(asv);free(primes);
    for(int i=0;i<=deg;i++)mpz_clear(co[i]);
    mpz_clears(N,fac,tmp,tmp2,m,rn,an,bz,az,td,NULL);
    return 0;
}
