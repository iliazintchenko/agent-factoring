/*
 * gnfs_work.c - Complete GNFS for 50-70 digit semiprimes
 * Compile: gcc -O3 -march=native -o gnfs_work library/gnfs_work.c -lgmp -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

#define SEED 42
#define MAX_DEG 5
#define MAX_RELS 200000
#define MAX_DEPS 64
#define NUM_QC 40

static struct timespec g_start;
static double elapsed(void){struct timespec now;clock_gettime(CLOCK_MONOTONIC,&now);
    return(now.tv_sec-g_start.tv_sec)+(now.tv_nsec-g_start.tv_nsec)/1e9;}

typedef struct{int deg,rfb,afb,sa,sb;double rtf,atf;}params_t;
static params_t get_params(int b){
    if(b<=100)return(params_t){3,5000,10000,25000,5000,.50,.50};
    if(b<=120)return(params_t){3,6000,12000,30000,6000,.50,.50};
    if(b<=140)return(params_t){3,8000,15000,35000,8000,.50,.50};
    if(b<=160)return(params_t){3,10000,18000,40000,10000,.50,.50};
    if(b<=170)return(params_t){3,12000,22000,50000,12000,.50,.50};
    if(b<=190)return(params_t){3,18000,30000,60000,15000,.50,.50};
    if(b<=210)return(params_t){4,25000,40000,70000,15000,.50,.50};
    if(b<=235)return(params_t){4,35000,55000,90000,20000,.50,.50};
    return(params_t){4,45000,70000,110000,25000,.50,.50};
}

static int *sieve_primes(int bnd,int *cnt){
    char *s=calloc(bnd+1,1);for(int i=2;(long)i*i<=bnd;i++)if(!s[i])for(int j=i*i;j<=bnd;j+=i)s[j]=1;
    int n=0;for(int i=2;i<=bnd;i++)if(!s[i])n++;int *p=malloc(n*sizeof(int));int k=0;
    for(int i=2;i<=bnd;i++)if(!s[i])p[k++]=i;free(s);*cnt=n;return p;}

typedef struct{int deg;mpz_t c[MAX_DEG+1];mpz_t m;}poly_t;
static void poly_init(poly_t *p){for(int i=0;i<=MAX_DEG;i++)mpz_init(p->c[i]);mpz_init(p->m);}
static void poly_clear(poly_t *p){for(int i=0;i<=MAX_DEG;i++)mpz_clear(p->c[i]);mpz_clear(p->m);}

static int poly_select(poly_t *p,const mpz_t N,int d){
    p->deg=d;mpz_t mb,mt,rem,bm,tc[MAX_DEG+1];mpz_inits(mb,mt,rem,bm,NULL);
    for(int i=0;i<=d;i++)mpz_init(tc[i]);mpz_root(mb,N,d+1);
    double best=1e300;int found=0;
    for(int delta=-20;delta<=20;delta++){
        mpz_set(mt,mb);if(delta>=0)mpz_add_ui(mt,mt,delta);else mpz_sub_ui(mt,mt,-delta);
        if(mpz_sgn(mt)<=0)continue;mpz_set(rem,N);
        for(int i=0;i<=d;i++)mpz_tdiv_qr(rem,tc[i],rem,mt);
        if(mpz_sgn(rem)!=0)mpz_add(tc[d],tc[d],rem);
        mpz_set_ui(rem,0);for(int i=d;i>=0;i--){mpz_mul(rem,rem,mt);mpz_add(rem,rem,tc[i]);}
        if(mpz_cmp(rem,N)!=0){mpz_sub(rem,N,rem);mpz_add(tc[0],tc[0],rem);}
        double sc=0;for(int i=0;i<=d;i++)sc+=fabs(mpz_get_d(tc[i]));
        if(sc<best){best=sc;mpz_set(bm,mt);for(int i=0;i<=d;i++)mpz_set(p->c[i],tc[i]);found=1;}
    }
    mpz_set(p->m,bm);for(int i=0;i<=d;i++)mpz_clear(tc[i]);mpz_clears(mb,mt,rem,bm,NULL);
    return found?0:-1;
}

static void poly_eval_homog(mpz_t res,poly_t *p,long a,unsigned long b){
    mpz_t t,ap,bp;mpz_inits(t,ap,bp,NULL);mpz_set_ui(res,0);mpz_set_si(ap,1);mpz_ui_pow_ui(bp,b,p->deg);
    for(int i=0;i<=p->deg;i++){mpz_mul(t,p->c[i],ap);mpz_mul(t,t,bp);mpz_add(res,res,t);
        mpz_mul_si(ap,ap,a);if(i<p->deg&&b>0)mpz_divexact_ui(bp,bp,b);}
    mpz_clears(t,ap,bp,NULL);
}

typedef struct{uint32_t *p,*r;uint8_t *l;int sz,cap;}fb_t;
static fb_t *fb_new(int c){fb_t *f=malloc(sizeof(fb_t));f->p=malloc(c*4);f->r=malloc(c*4);f->l=malloc(c);f->sz=0;f->cap=c;return f;}
static void fb_add(fb_t *f,uint32_t p,uint32_t r,uint8_t l){if(f->sz>=f->cap){f->cap*=2;f->p=realloc(f->p,f->cap*4);f->r=realloc(f->r,f->cap*4);f->l=realloc(f->l,f->cap);}f->p[f->sz]=p;f->r[f->sz]=r;f->l[f->sz]=l;f->sz++;}
static void fb_free(fb_t *f){free(f->p);free(f->r);free(f->l);free(f);}
static fb_t *build_rfb(poly_t *py,int bnd){int np;int *pr=sieve_primes(bnd,&np);fb_t *fb=fb_new(np+10);for(int i=0;i<np;i++)fb_add(fb,pr[i],mpz_fdiv_ui(py->m,pr[i]),(uint8_t)(log2(pr[i])+.5));free(pr);return fb;}
static fb_t *build_afb(poly_t *py,int bnd){int np;int *pr=sieve_primes(bnd,&np);fb_t *fb=fb_new(np*py->deg+10);
    for(int i=0;i<np;i++){uint32_t p=pr[i];uint8_t lp=(uint8_t)(log2(p)+.5);uint64_t cm[MAX_DEG+1];
        for(int j=0;j<=py->deg;j++)cm[j]=mpz_fdiv_ui(py->c[j],p);
        for(uint32_t x=0;x<p;x++){uint64_t v=0;for(int j=py->deg;j>=0;j--)v=(v*x+cm[j])%p;if(v==0)fb_add(fb,p,x,lp);}}
    free(pr);return fb;}

typedef struct{long *a;uint32_t *b;int16_t **re,**ae;int *rs,*as;int cnt,cap;}rels_t;
static rels_t *rels_new(int cap,int rn,int an){rels_t *r=calloc(1,sizeof(rels_t));r->a=malloc(cap*sizeof(long));r->b=malloc(cap*4);r->re=malloc(cap*sizeof(int16_t*));r->ae=malloc(cap*sizeof(int16_t*));r->rs=calloc(cap,sizeof(int));r->as=calloc(cap,sizeof(int));for(int i=0;i<cap;i++){r->re[i]=calloc(rn,sizeof(int16_t));r->ae[i]=calloc(an,sizeof(int16_t));}r->cap=cap;return r;}
static int trial_div(mpz_t n,fb_t *fb,int16_t *e,int *s){*s=0;if(mpz_sgn(n)<0){*s=1;mpz_neg(n,n);}if(mpz_sgn(n)==0)return 0;memset(e,0,fb->sz*sizeof(int16_t));for(int i=0;i<fb->sz;i++){uint32_t p=fb->p[i];while(mpz_divisible_ui_p(n,p)){e[i]++;mpz_divexact_ui(n,n,p);}}return mpz_cmp_ui(n,1)==0;}

typedef uint64_t u64;
typedef struct{u64 **rows;int nr,nc,fw,iw,wp;}gf2_t;
static gf2_t *gf2_new(int nr,int nc){gf2_t *m=malloc(sizeof(gf2_t));m->nr=nr;m->nc=nc;m->fw=(nc+63)/64;m->iw=(nr+63)/64;m->wp=m->fw+m->iw;m->rows=malloc(nr*sizeof(u64*));for(int i=0;i<nr;i++){m->rows[i]=calloc(m->wp,sizeof(u64));m->rows[i][m->fw+i/64]|=(1ULL<<(i%64));}return m;}
static void gf2_free(gf2_t *m){for(int i=0;i<m->nr;i++)free(m->rows[i]);free(m->rows);free(m);}
static int gf2_solve(gf2_t *m,int ***dp,int **dl,int mx){int pv=0;for(int c=0;c<m->nc&&pv<m->nr;c++){int pr=-1;for(int r=pv;r<m->nr;r++)if((m->rows[r][c/64]>>(c%64))&1){pr=r;break;}if(pr<0)continue;if(pr!=pv){u64 *t=m->rows[pr];m->rows[pr]=m->rows[pv];m->rows[pv]=t;}for(int r=0;r<m->nr;r++){if(r==pv)continue;if((m->rows[r][c/64]>>(c%64))&1)for(int w=0;w<m->wp;w++)m->rows[r][w]^=m->rows[pv][w];}pv++;}int nd=0;*dp=malloc(mx*sizeof(int*));*dl=malloc(mx*sizeof(int));for(int r=pv;r<m->nr&&nd<mx;r++){int z=1;for(int w=0;w<m->fw&&z;w++){u64 mk=(w<m->fw-1)?~0ULL:(m->nc%64==0?~0ULL:(1ULL<<(m->nc%64))-1);if(m->rows[r][w]&mk)z=0;}if(!z)continue;int *d=malloc(m->nr*sizeof(int));int n=0;for(int w=0;w<m->iw;w++){u64 bits=m->rows[r][m->fw+w];while(bits){int b=__builtin_ctzll(bits);int idx=w*64+b;if(idx<m->nr)d[n++]=idx;bits&=bits-1;}}if(n>0){(*dp)[nd]=d;(*dl)[nd]=n;nd++;}else free(d);}return nd;}

static uint64_t mulmod64(uint64_t a,uint64_t b,uint64_t m){return(unsigned __int128)a*b%m;}
static uint64_t powmod64(uint64_t b,uint64_t e,uint64_t m){uint64_t r=1;b%=m;while(e){if(e&1)r=mulmod64(r,b,m);b=mulmod64(b,b,m);e>>=1;}return r;}

/* ====== Polynomial ring Z_M[x]/(f(x)) ====== */
/* Multiply a * b mod (f, M), all mpz arrays of length d */
static void prmul(mpz_t *res, const mpz_t *a, const mpz_t *b,
                   const mpz_t *fc, int d, const mpz_t M) {
    mpz_t prod[2*MAX_DEG], t;
    for(int i=0;i<2*d;i++) mpz_init_set_ui(prod[i],0);
    mpz_init(t);
    for(int i=0;i<d;i++) for(int j=0;j<d;j++){
        mpz_mul(t,a[i],b[j]);mpz_add(prod[i+j],prod[i+j],t);mpz_mod(prod[i+j],prod[i+j],M);
    }
    /* Reduce using f: x^d = -(c0+c1*x+...+c_{d-1}*x^{d-1})/c_d */
    mpz_t lci; mpz_init(lci);
    if (!mpz_invert(lci, fc[d], M)) {
        /* GCD of fc[d] and M is non-trivial! This reveals a factor. */
        mpz_t g; mpz_init(g); mpz_gcd(g, fc[d], M);
        if (mpz_cmp_ui(g,1)>0 && mpz_cmp(g,M)<0) {
            mpz_t co; mpz_init(co); mpz_divexact(co,M,g);
            gmp_printf("%Zd\n%Zd\n",g,co); mpz_clear(co);
            exit(0); /* Found factor! */
        }
        mpz_clear(g);
        mpz_set_ui(lci,0); /* fallback */
    }
    for(int k=2*(d-1);k>=d;k--){
        if(mpz_sgn(prod[k])==0) continue;
        mpz_t c; mpz_init(c);
        mpz_mul(c,prod[k],lci); mpz_mod(c,c,M);
        for(int j=0;j<d;j++){
            mpz_mul(t,c,fc[j]); mpz_sub(prod[k-d+j],prod[k-d+j],t);
            mpz_mod(prod[k-d+j],prod[k-d+j],M);
        }
        mpz_set_ui(prod[k],0); mpz_clear(c);
    }
    for(int i=0;i<d;i++) mpz_set(res[i],prod[i]);
    for(int i=0;i<2*d;i++) mpz_clear(prod[i]);
    mpz_clears(t,lci,NULL);
}

/* Compute base^exp mod (f, M) via repeated squaring */
static void prpow(mpz_t *res, const mpz_t *base, const mpz_t exp,
                   const mpz_t *fc, int d, const mpz_t M) {
    mpz_t b[MAX_DEG], r[MAX_DEG];
    for(int i=0;i<d;i++){mpz_init_set(b[i],base[i]);mpz_init_set_ui(r[i],0);}
    mpz_set_ui(r[0],1); /* r = 1 */

    /* Binary method */
    int nbits = mpz_sizeinbase(exp, 2);
    for(int bit = nbits-1; bit >= 0; bit--) {
        /* Square r */
        prmul(r, r, r, fc, d, M);
        /* If bit is set, multiply by base */
        if (mpz_tstbit(exp, bit))
            prmul(r, r, b, fc, d, M);
    }
    for(int i=0;i<d;i++) mpz_set(res[i],r[i]);
    for(int i=0;i<d;i++){mpz_clear(b[i]);mpz_clear(r[i]);}
}

int main(int argc, char *argv[]) {
    clock_gettime(CLOCK_MONOTONIC, &g_start);
    gmp_randstate_t rstate; gmp_randinit_default(rstate); gmp_randseed_ui(rstate, SEED);

    mpz_t N; mpz_init(N);
    if(argc>=2) mpz_set_str(N,argv[1],10);
    else{char buf[1024];if(!fgets(buf,sizeof buf,stdin)){fprintf(stderr,"No input\n");return 1;}
        char *p=buf;while(*p==' '||*p=='\t')p++;char *e=p+strlen(p)-1;while(e>p&&(*e=='\n'||*e=='\r'||*e==' '))*e--=0;mpz_set_str(N,p,10);}
    int digits=(int)mpz_sizeinbase(N,10),bits=(int)mpz_sizeinbase(N,2);
    fprintf(stderr,"GNFS: %dd (%db)\n",digits,bits);

    /* Trial division */
    for(unsigned long p=2;p<1000000;p++){if(mpz_divisible_ui_p(N,p)){mpz_t q;mpz_init(q);mpz_divexact_ui(q,N,p);gmp_printf("%lu\n%Zd\n",p,q);mpz_clears(q,N,NULL);return 0;}}
    if(mpz_perfect_power_p(N)){mpz_t r;mpz_init(r);for(int e=2;e<=64;e++)if(mpz_root(r,N,e)){mpz_t o;mpz_init(o);mpz_divexact(o,N,r);gmp_printf("%Zd\n%Zd\n",r,o);mpz_clears(r,o,N,NULL);return 0;}mpz_clear(r);}

    params_t P=get_params(bits);
    fprintf(stderr,"d=%d rfb=%d afb=%d sa=%d sb=%d\n",P.deg,P.rfb,P.afb,P.sa,P.sb);

    poly_t poly;poly_init(&poly);
    if(poly_select(&poly,N,P.deg)<0){fprintf(stderr,"poly fail\n");return 1;}
    gmp_fprintf(stderr,"m=%Zd\nf(x)=",poly.m);
    for(int i=poly.deg;i>=0;i--)gmp_fprintf(stderr," %+Zd*x^%d",poly.c[i],i);fprintf(stderr,"\n");
    {mpz_t v;mpz_init(v);mpz_set_ui(v,0);for(int i=poly.deg;i>=0;i--){mpz_mul(v,v,poly.m);mpz_add(v,v,poly.c[i]);}if(mpz_cmp(v,N)!=0){fprintf(stderr,"f(m)!=N\n");return 1;}mpz_clear(v);}

    fb_t *rfb=build_rfb(&poly,P.rfb),*afb=build_afb(&poly,P.afb);
    fprintf(stderr,"rfb=%d afb=%d t=%.1fs\n",rfb->sz,afb->sz,elapsed());

    uint32_t qcp[NUM_QC],qcr[NUM_QC];int nqc=0;
    {int nqp;int *qp=sieve_primes(P.afb+5000,&nqp);
     for(int i=0;i<nqp&&nqc<NUM_QC;i++){if(qp[i]<=P.afb)continue;uint32_t p=qp[i];
         uint64_t cm[MAX_DEG+1];for(int j=0;j<=poly.deg;j++)cm[j]=mpz_fdiv_ui(poly.c[j],p);
         for(uint32_t x=0;x<p;x++){uint64_t v=0;for(int j=poly.deg;j>=0;j--)v=(v*x+cm[j])%p;if(v==0){qcp[nqc]=p;qcr[nqc]=x;nqc++;break;}}}free(qp);}

    /* Factor the leading coefficient c_d to add its primes to the matrix */
    #define MAX_CD_PRIMES 20
    uint32_t cd_primes[MAX_CD_PRIMES];
    int cd_exps[MAX_CD_PRIMES]; /* exponent of each prime in c_d */
    int ncd = 0;
    {
        mpz_t cd_tmp; mpz_init(cd_tmp);
        mpz_abs(cd_tmp, poly.c[poly.deg]);
        for (uint32_t p = 2; p < 1000000 && mpz_cmp_ui(cd_tmp, 1) > 0 && ncd < MAX_CD_PRIMES; p++) {
            if (!mpz_divisible_ui_p(cd_tmp, p)) continue;
            cd_primes[ncd] = p;
            cd_exps[ncd] = 0;
            while (mpz_divisible_ui_p(cd_tmp, p)) {
                cd_exps[ncd]++;
                mpz_divexact_ui(cd_tmp, cd_tmp, p);
            }
            ncd++;
        }
        if (mpz_cmp_ui(cd_tmp, 1) > 0) {
            /* Remaining factor > 10^6, treat as single "prime" */
            cd_primes[ncd] = mpz_get_ui(cd_tmp); /* may overflow for very large c_d */
            cd_exps[ncd] = 1;
            ncd++;
        }
        mpz_clear(cd_tmp);
        fprintf(stderr, "Leading coeff c_%d has %d prime factors\n", poly.deg, ncd);
    }

    int ncols=1+rfb->sz+1+afb->sz+nqc; /* +ncd removed: c_d is implicit in the ring */
    int target=ncols+30;

    int sw=2*P.sa+1;uint8_t *rsv=malloc(sw),*asv=malloc(sw);
    int rcap=target+1000;if(rcap>MAX_RELS)rcap=MAX_RELS;
    rels_t *rels=rels_new(rcap,rfb->sz,afb->sz);
    mpz_t rn,an;mpz_inits(rn,an,NULL);

    fprintf(stderr,"Sieving (target=%d)...\n",target);int tc=0;
    for(uint32_t b=1;b<=(uint32_t)P.sb&&rels->cnt<target;b++){
        if(elapsed()>180){fprintf(stderr,"timeout b=%u r=%d\n",b,rels->cnt);break;}
        if(b%200==0)fprintf(stderr,"  b=%u r=%d/%d c=%d t=%.1fs\n",b,rels->cnt,target,tc,elapsed());
        memset(rsv,0,sw);memset(asv,0,sw);
        for(int i=0;i<rfb->sz;i++){uint32_t p=rfb->p[i];uint8_t l=rfb->l[i];uint64_t h=((uint64_t)(b%p)*rfb->r[i])%p;long off=((long)h+P.sa)%(long)p;if(off<0)off+=p;for(int j=(int)off;j<sw;j+=p)rsv[j]+=l;}
        for(int i=0;i<afb->sz;i++){uint32_t p=afb->p[i];uint8_t l=afb->l[i];uint64_t h=((uint64_t)(b%p)*afb->r[i])%p;long off=((long)h+P.sa)%(long)p;if(off<0)off+=p;for(int j=(int)off;j<sw;j+=p)asv[j]+=l;}
        double lr=log2((double)P.sa+(double)b*mpz_get_d(poly.m));
        double la=0;{double am=P.sa,bv=b,sum=0;for(int i=0;i<=poly.deg;i++)sum+=fabs(mpz_get_d(poly.c[i]))*pow(am,i)*pow(bv,poly.deg-i);la=log2(sum+1.0);}
        int rt=(int)(lr*P.rtf),at=(int)(la*P.atf);
        for(int j=0;j<sw;j++){
            if(rsv[j]<rt||asv[j]<at)continue;
            long a=-(long)P.sa+j;if(!a)continue;long ga=a<0?-a:a,gb=(long)b;while(gb){long t=gb;gb=ga%gb;ga=t;}if(ga!=1)continue;tc++;
            mpz_mul_ui(rn,poly.m,b);mpz_t at2;mpz_init_set_si(at2,a);mpz_sub(rn,at2,rn);mpz_clear(at2);
            poly_eval_homog(an,&poly,a,b);
            int rsn,asn;int16_t *re=rels->re[rels->cnt],*ae=rels->ae[rels->cnt];
            if(!trial_div(rn,rfb,re,&rsn))continue;if(!trial_div(an,afb,ae,&asn))continue;
            rels->a[rels->cnt]=a;rels->b[rels->cnt]=b;rels->rs[rels->cnt]=rsn;rels->as[rels->cnt]=asn;rels->cnt++;if(rels->cnt>=rcap)break;
        }
    }
    fprintf(stderr,"Sieve: %d rels %d cand %.1fs\n",rels->cnt,tc,elapsed());
    free(rsv);free(asv);
    if(rels->cnt<ncols+1){fprintf(stderr,"Not enough: %d<%d\n",rels->cnt,ncols+1);return 1;}

    int nrels=rels->cnt;
    gf2_t *mat=gf2_new(nrels,ncols);
    for(int r=0;r<nrels;r++){int c=0;
        if(rels->rs[r])mat->rows[r][c/64]^=(1ULL<<(c%64));c++;
        for(int k=0;k<rfb->sz;k++){if(rels->re[r][k]&1)mat->rows[r][c/64]^=(1ULL<<(c%64));c++;}
        if(rels->as[r])mat->rows[r][c/64]^=(1ULL<<(c%64));c++;
        for(int k=0;k<afb->sz;k++){if(rels->ae[r][k]&1)mat->rows[r][c/64]^=(1ULL<<(c%64));c++;}
        for(int k=0;k<nqc;k++){
            uint64_t av=((int64_t)(rels->a[r]%(int64_t)qcp[k])+(int64_t)qcp[k])%(int64_t)qcp[k];
            uint64_t bv2=mulmod64((uint64_t)(rels->b[r]%qcp[k]),(uint64_t)qcr[k],(uint64_t)qcp[k]);
            uint64_t v=(av+qcp[k]-bv2)%qcp[k];
            if(v!=0&&powmod64(v,(qcp[k]-1)/2,qcp[k])!=1)mat->rows[r][c/64]^=(1ULL<<(c%64));c++;
        }
        /* (c_d columns removed - c_d is implicit in the number ring) */
    }
    fprintf(stderr,"Gauss elim (%dx%d)...\n",nrels,ncols);
    int **deps;int *dlen;int ndeps=gf2_solve(mat,&deps,&dlen,MAX_DEPS);
    fprintf(stderr,"%d deps (%.1fs)\n",ndeps,elapsed());gf2_free(mat);
    if(!ndeps){fprintf(stderr,"no deps\n");return 1;}

    /* ====== Square root ====== */
    mpz_t X,tmp;mpz_inits(X,tmp,NULL);int found=0;
    int d = poly.deg;

    /* Precompute f coefficients mod N */
    mpz_t fc[MAX_DEG+1];
    for(int i=0;i<=d;i++){mpz_init(fc[i]);mpz_mod(fc[i],poly.c[i],N);}

    for(int di=0;di<ndeps&&!found;di++){
        if(elapsed()>260)break;
        fprintf(stderr,"Dep %d (sz %d)...\n",di,dlen[di]);

        int *rex=calloc(rfb->sz,sizeof(int)),*aex=calloc(afb->sz,sizeof(int));
        int rn2=0,an2=0;
        for(int i=0;i<dlen[di];i++){int ri=deps[di][i];
            for(int c=0;c<rfb->sz;c++)rex[c]+=rels->re[ri][c];
            for(int c=0;c<afb->sz;c++)aex[c]+=rels->ae[ri][c];
            rn2+=rels->rs[ri];an2+=rels->as[ri];}
        int even=1;if(rn2&1)even=0;if(an2&1)even=0;
        for(int c=0;c<rfb->sz&&even;c++)if(rex[c]&1)even=0;
        for(int c=0;c<afb->sz&&even;c++)if(aex[c]&1)even=0;
        if(!even){free(rex);free(aex);continue;}

        /* X = rational sqrt = prod(p^(e/2)) mod N */
        mpz_set_ui(X,1);
        for(int c=0;c<rfb->sz;c++){if(rex[c]<=0)continue;
            mpz_set_ui(tmp,rfb->p[c]);mpz_powm_ui(tmp,tmp,rex[c]/2,N);mpz_mul(X,X,tmp);mpz_mod(X,X,N);}

        long *da=malloc(dlen[di]*sizeof(long));
        uint32_t *db=malloc(dlen[di]*sizeof(uint32_t));
        for(int i=0;i<dlen[di];i++){da[i]=rels->a[deps[di][i]];db[i]=rels->b[deps[di][i]];}

        /* Skip ring product (only needed for Hensel, which is disabled) */
        if (0) {
        fprintf(stderr,"  Product in ring...\n");
        mpz_t S[MAX_DEG];
        for(int k=0;k<d;k++) mpz_init_set_ui(S[k],0);
        mpz_set_ui(S[0],1);

        mpz_t factor[MAX_DEG];
        for(int k=0;k<d;k++) mpz_init(factor[k]);

        for(int i=0;i<dlen[di];i++){
            mpz_set_si(factor[0],da[i]);mpz_mod(factor[0],factor[0],N);
            mpz_set_si(factor[1],-(long)db[i]);mpz_mod(factor[1],factor[1],N);
            for(int k=2;k<d;k++) mpz_set_ui(factor[k],0);
            prmul(S,S,factor,fc,d,N);
        }
        for(int k=0;k<d;k++) mpz_clear(factor[k]);

        fprintf(stderr,"  Product done (%.1fs)\n",elapsed());
        } /* end disabled ring product */

        /* Skip ring power approach AND ring product (not needed for CRT sqrt) */
        if (0) {
        mpz_t S[MAX_DEG]; /* declared here to avoid compiler error */
        /* Algebraic sqrt via probabilistic method in Z_N[x]/(f):
         * We want T such that T^2 = S in the ring.
         * Use: T = S^((N^d+1)/4) if N^d ≡ 3 mod 4 (unlikely)
         * General: pick random h, compute (h^2 - S)^((N^d-1)/2) then extract sqrt.
         *
         * Adleman-Manders-Miller for the ring:
         * Write N^d - 1 = 2^s * q where q odd.
         * For random h:
         *   compute h^q mod (f, N)
         *   then square s times
         * If we get identity before -1, we can extract a factor.
         *
         * But N is composite! We don't know the factorization of N^d - 1.
         * Actually |Z_N[x]/(f)^*| = |(Z_p[x]/(f))^* x (Z_q[x]/(f))^*|
         * = (p^d-1)(q^d-1) (when f is irreducible mod both p and q, which it's not
         * necessarily since f has roots mod N = p*q).
         *
         * Simpler approach: compute S^((N+1)/2) in the ring.
         * This doesn't give the full square root in the ring, but it might give
         * an element whose evaluation at m gives a non-trivial square root of X^2 mod N.
         *
         * Actually for the scalar: if T(m)^2 = X^2 mod N, then T(m) is a square root
         * of X^2 mod N. There are 4 such roots. Computing the product
         * prod(a_i - b_i*m) mod N directly and taking its "square root" via
         * S^((N+1)/2) mod N (if N ≡ 3 mod 4) gives one root.
         * But this root is the same as X or -X (since the computation is in Z/NZ).
         *
         * The key: by computing in the POLYNOMIAL ring rather than in Z/NZ,
         * the operations are "lifted" and can produce a DIFFERENT square root.
         * Specifically: S^((N+1)/2) evaluated in Z_N[x]/(f(x)) then evaluated at m
         * gives a square root that may differ from X computed from the rational FB.
         *
         * Let's try it!
         */

        /* Compute exponent = (N+1)/2 (works if N is odd, which it is for semiprimes) */
        mpz_t exp2; mpz_init(exp2);
        mpz_add_ui(exp2, N, 1);
        mpz_tdiv_q_2exp(exp2, exp2, 1); /* (N+1)/2 */

        fprintf(stderr,"  Computing S^((N+1)/2) in ring (%zu-bit exp)...\n",
                mpz_sizeinbase(exp2, 2));

        mpz_t T[MAX_DEG];
        for(int k=0;k<d;k++) mpz_init(T[k]);
        prpow(T, S, exp2, fc, d, N);

        fprintf(stderr,"  Power done (%.1fs)\n", elapsed());

        /* Y = T(m) mod N */
        mpz_t Y; mpz_init(Y); mpz_set_ui(Y, 0);
        for(int k=d-1;k>=0;k--){mpz_mul(Y,Y,poly.m);mpz_add(Y,Y,T[k]);mpz_mod(Y,Y,N);}

        /* Check Y^2 = X^2 mod N */
        mpz_t y2,x2; mpz_inits(y2,x2,NULL);
        mpz_mul(y2,Y,Y);mpz_mod(y2,y2,N);
        mpz_mul(x2,X,X);mpz_mod(x2,x2,N);
        gmp_fprintf(stderr,"  Y=%Zd\n  Y^2=%Zd\n  X^2=%Zd\n  eq=%d\n",Y,y2,x2,mpz_cmp(y2,x2)==0);

        mpz_t g; mpz_init(g);
        if (mpz_cmp(y2,x2)==0) {
            /* Y^2 = X^2 mod N, try factoring */
            mpz_sub(g,X,Y);mpz_mod(g,g,N);mpz_gcd(g,g,N);
            if(mpz_cmp_ui(g,1)>0&&mpz_cmp(g,N)<0){
                mpz_t co;mpz_init(co);mpz_divexact(co,N,g);gmp_printf("%Zd\n%Zd\n",g,co);mpz_clear(co);found=1;
            }
            if(!found){mpz_add(g,X,Y);mpz_mod(g,g,N);mpz_gcd(g,g,N);
                if(mpz_cmp_ui(g,1)>0&&mpz_cmp(g,N)<0){
                    mpz_t co;mpz_init(co);mpz_divexact(co,N,g);gmp_printf("%Zd\n%Zd\n",g,co);mpz_clear(co);found=1;
                }
            }
        }

        if (!found) {
            /* Try other exponents */
            /* N^2+1)/2 */
            mpz_mul(exp2, N, N); mpz_add_ui(exp2, exp2, 1); mpz_tdiv_q_2exp(exp2, exp2, 1);
            fprintf(stderr,"  Trying N^2 exponent...\n");
            prpow(T, S, exp2, fc, d, N);
            mpz_set_ui(Y, 0);
            for(int k=d-1;k>=0;k--){mpz_mul(Y,Y,poly.m);mpz_add(Y,Y,T[k]);mpz_mod(Y,Y,N);}
            mpz_mul(y2,Y,Y);mpz_mod(y2,y2,N);
            if(mpz_cmp(y2,x2)==0){
                mpz_sub(g,X,Y);mpz_mod(g,g,N);mpz_gcd(g,g,N);
                if(mpz_cmp_ui(g,1)>0&&mpz_cmp(g,N)<0){mpz_t co;mpz_init(co);mpz_divexact(co,N,g);gmp_printf("%Zd\n%Zd\n",g,co);mpz_clear(co);found=1;}
                if(!found){mpz_add(g,X,Y);mpz_mod(g,g,N);mpz_gcd(g,g,N);
                    if(mpz_cmp_ui(g,1)>0&&mpz_cmp(g,N)<0){mpz_t co;mpz_init(co);mpz_divexact(co,N,g);gmp_printf("%Zd\n%Zd\n",g,co);mpz_clear(co);found=1;}}
            }
        }

        if (!found) {
            /* Try (N^3+1)/2 for degree 3 */
            mpz_t N3; mpz_init(N3);
            mpz_pow_ui(N3, N, d);
            mpz_add_ui(N3, N3, 1);
            mpz_tdiv_q_2exp(exp2, N3, 1);
            fprintf(stderr,"  Trying N^%d exponent...\n", d);
            prpow(T, S, exp2, fc, d, N);
            mpz_set_ui(Y, 0);
            for(int k=d-1;k>=0;k--){mpz_mul(Y,Y,poly.m);mpz_add(Y,Y,T[k]);mpz_mod(Y,Y,N);}
            mpz_mul(y2,Y,Y);mpz_mod(y2,y2,N);
            gmp_fprintf(stderr,"  Y^2=%Zd X^2=%Zd eq=%d\n",y2,x2,mpz_cmp(y2,x2)==0);
            if(mpz_cmp(y2,x2)==0){
                mpz_sub(g,X,Y);mpz_mod(g,g,N);mpz_gcd(g,g,N);
                if(mpz_cmp_ui(g,1)>0&&mpz_cmp(g,N)<0){mpz_t co;mpz_init(co);mpz_divexact(co,N,g);gmp_printf("%Zd\n%Zd\n",g,co);mpz_clear(co);found=1;}
                if(!found){mpz_add(g,X,Y);mpz_mod(g,g,N);mpz_gcd(g,g,N);
                    if(mpz_cmp_ui(g,1)>0&&mpz_cmp(g,N)<0){mpz_t co;mpz_init(co);mpz_divexact(co,N,g);gmp_printf("%Zd\n%Zd\n",g,co);mpz_clear(co);found=1;}}
            }
            mpz_clear(N3);
        }

        } /* end of disabled ring power block */

        /* Declare variables needed by Hensel sqrt */
        mpz_t g; mpz_init(g);

        /* ===== Hensel-lifting algebraic sqrt ===== */
        if (!found) {
            fprintf(stderr,"  Trying Hensel lifting algebraic sqrt...\n");

            /* Multi-prime CRT: compute T(m) mod N using many small primes.
             * For each prime q where f splits:
             *   - Find roots r_1,...,r_d of f mod q
             *   - Compute S_j = prod(a_i - b_i*r_j) mod q
             *   - Compute T_j = sqrt(S_j) mod q (fixed sign: < q/2)
             *   - Lagrange: T(m) mod q
             * CRT across primes gives T(m) mod (product of q's).
             * When product > N, reduce mod N. Try both ± signs at end. */

            /* Collect CRT primes and T(m) values */
            #define MAX_CRT2 100
            unsigned long crt2_q[MAX_CRT2];
            unsigned long crt2_tm[MAX_CRT2]; /* T(m) mod q for each prime */
            int ncrt2 = 0;

            unsigned long hp = 0;
            unsigned long hr[MAX_DEG];
            for (unsigned long q = 1009; q < 50000 && ncrt2 < MAX_CRT2; q += 2) {
                mpz_set_ui(tmp, q);
                if (!mpz_probab_prime_p(tmp, 2)) continue;
                int nr = 0;
                for (unsigned long x = 0; x < q && nr <= d; x++) {
                    unsigned long long val = 0;
                    for (int j = d; j >= 0; j--) {
                        unsigned long cj = mpz_fdiv_ui(fc[j], q);
                        val = (val * x + cj) % q;
                    }
                    if (val == 0) hr[nr++] = x;
                }
                if (nr == d) hp = q;
            }

            /* Convert to multi-prime CRT: for each suitable prime, compute T(m) mod q */
            for (unsigned long q2 = 1009; q2 < 50000 && ncrt2 < MAX_CRT2; q2 += 2) {
                mpz_set_ui(tmp, q2);
                if (!mpz_probab_prime_p(tmp, 2)) continue;
                /* Check f has d distinct roots mod q2 */
                int nr2 = 0;
                unsigned long rts2[MAX_DEG+1];
                for (unsigned long x = 0; x < q2 && nr2 <= d; x++) {
                    unsigned long long val = 0;
                    for (int j = d; j >= 0; j--) {
                        unsigned long cj = mpz_fdiv_ui(fc[j], q2);
                        val = (val * x + cj) % q2;
                    }
                    if (val == 0) rts2[nr2++] = x;
                }
                if (nr2 != d) continue;

                /* Compute S_j = prod(a_i - b_i*r_j) mod q2 */
                unsigned long Sj[MAX_DEG];
                int sq_ok = 1;
                unsigned long Tj[MAX_DEG];
                for (int j = 0; j < d && sq_ok; j++) {
                    unsigned long long pv = 1;
                    for (int i = 0; i < dlen[di]; i++) {
                        long long t2 = ((da[i] % (long long)q2) + q2) % q2;
                        t2 = (t2 + q2 - ((unsigned long long)db[i] % q2 * rts2[j]) % q2) % q2;
                        pv = pv * (unsigned long long)t2 % q2;
                    }
                    Sj[j] = (unsigned long)pv;
                    /* sqrt mod q2 */
                    if (pv == 0) { Tj[j] = 0; continue; }
                    unsigned long long r2=1,b2=pv,e2=(q2-1)/2; unsigned long long m2=q2;
                    while(e2){if(e2&1)r2=r2*b2%m2;b2=b2*b2%m2;e2>>=1;}
                    if (r2 != 1) { sq_ok = 0; break; }
                    /* Tonelli-Shanks */
                    if (q2%4==3){r2=1;b2=pv;e2=(q2+1)/4;while(e2){if(e2&1)r2=r2*b2%m2;b2=b2*b2%m2;e2>>=1;}Tj[j]=(unsigned long)r2;}
                    else{unsigned long Q3=q2-1,S3=0;while(Q3%2==0){Q3/=2;S3++;}
                        unsigned long z2=2;for(;;){r2=1;b2=z2;e2=(q2-1)/2;while(e2){if(e2&1)r2=r2*b2%m2;b2=b2*b2%m2;e2>>=1;}if(r2==m2-1)break;z2++;}
                        unsigned long long M3=S3,c3,t3,R3;
                        r2=1;b2=z2;e2=Q3;while(e2){if(e2&1)r2=r2*b2%m2;b2=b2*b2%m2;e2>>=1;}c3=r2;
                        r2=1;b2=pv;e2=Q3;while(e2){if(e2&1)r2=r2*b2%m2;b2=b2*b2%m2;e2>>=1;}t3=r2;
                        r2=1;b2=pv;e2=(Q3+1)/2;while(e2){if(e2&1)r2=r2*b2%m2;b2=b2*b2%m2;e2>>=1;}R3=r2;
                        for(;;){if(t3==1){Tj[j]=(unsigned long)R3;break;}int i3=0;unsigned long long tt=t3;
                        while(tt!=1){tt=tt*tt%q2;i3++;}unsigned long long bb=c3;
                        for(int jj=0;jj<(int)M3-i3-1;jj++)bb=bb*bb%q2;M3=i3;c3=bb*bb%q2;t3=t3*c3%q2;R3=R3*bb%q2;}}
                    /* Sign will be determined below */
                }
                if (!sq_ok) continue;

                /* Compute prod(a-b*m) mod q2 for verification */
                unsigned long long prod_mod_q = 1;
                unsigned long mmod = mpz_fdiv_ui(poly.m, q2);
                for (int i2 = 0; i2 < dlen[di]; i2++) {
                    long long t2 = ((da[i2] % (long long)q2) + q2) % q2;
                    t2 = (t2 + q2 - ((unsigned long long)db[i2] % q2 * mmod) % q2) % q2;
                    prod_mod_q = prod_mod_q * (unsigned long long)t2 % q2;
                }

                /* Try all 2^d sign combos, find the one where T(m)^2 = prod mod q2 */
                int max_sc = 1 << d;
                unsigned long best_tm = 0;
                int found_sc = 0;
                /* Optimization: if we already know the correct sign combo, use it */
                /* Reset per dependency */
                int known_sc = -1;
                int sc_start = (known_sc >= 0) ? known_sc : 0;
                int sc_end = (known_sc >= 0) ? known_sc + 1 : max_sc;
                for (int sc = sc_start; sc < sc_end && !found_sc; sc++) {
                    unsigned long Tj_sc[MAX_DEG];
                    for (int j = 0; j < d; j++)
                        Tj_sc[j] = (sc & (1<<j)) ? (q2 - Tj[j]) % q2 : Tj[j];

                    /* Lagrange: T(m) mod q2 */
                    unsigned long long result = 0;
                    for (int j = 0; j < d; j++) {
                        unsigned long long num = Tj_sc[j], den = 1;
                        for (int k = 0; k < d; k++) {
                            if (k == j) continue;
                            num = num * ((mmod + q2 - rts2[k]) % q2) % q2;
                            den = den * ((rts2[j] + q2 - rts2[k]) % q2) % q2;
                        }
                        unsigned long long inv2 = 1, bb = den, ee = q2 - 2;
                        while (ee) { if (ee&1) inv2=inv2*bb%q2; bb=bb*bb%q2; ee>>=1; }
                        result = (result + num * inv2 % q2) % q2;
                    }

                    /* Check: T(m)^2 = prod(a-b*m) mod q2 */
                    unsigned long long tm2 = result * result % q2;
                    if (tm2 == prod_mod_q) {
                        best_tm = (unsigned long)result;
                        found_sc = 1;
                        if (known_sc < 0) known_sc = sc;
                    }
                }

                if (found_sc) {
                    crt2_q[ncrt2] = q2;
                    crt2_tm[ncrt2] = best_tm;
                    ncrt2++;
                }
            }

            fprintf(stderr,"  CRT: %d primes for dep %d\n", ncrt2, di);

            /* Simple CRT: use the sign-verified T(m) values.
             * Each crt2_tm[i] was verified: T(m)^2 = prod(a-b*m) mod q_i.
             * CRT them directly and test at the end. */
            if (ncrt2 >= 8) {
                mpz_t crt_v, crt_m;
                mpz_inits(crt_v, crt_m, NULL);
                mpz_set_ui(crt_v, crt2_tm[0]);
                mpz_set_ui(crt_m, crt2_q[0]);

                for (int qi = 1; qi < ncrt2; qi++) {
                    unsigned long q2 = crt2_q[qi];
                    unsigned long tm = crt2_tm[qi];
                    unsigned long a_mod = mpz_fdiv_ui(crt_v, q2);
                    long diff = (long)tm - (long)a_mod;
                    if (diff < 0) diff += q2;
                    unsigned long inv_m2;
                    { unsigned long long iv=1,bb=mpz_fdiv_ui(crt_m,q2),ee=q2-2;
                      while(ee){if(ee&1)iv=iv*bb%q2;bb=bb*bb%q2;ee>>=1;} inv_m2=(unsigned long)iv; }
                    unsigned long t = (unsigned long long)((unsigned long)diff) * inv_m2 % q2;
                    mpz_addmul_ui(crt_v, crt_m, t);
                    mpz_mul_ui(crt_m, crt_m, q2);

                    if (mpz_cmp(crt_m, N) > 0) {
                        mpz_t Y_crt; mpz_init(Y_crt);
                        mpz_mod(Y_crt, crt_v, N);
                        /* Try gcd(X ± Y, N) */
                        mpz_sub(g,X,Y_crt);mpz_mod(g,g,N);mpz_gcd(g,g,N);
                        if(mpz_cmp_ui(g,1)>0&&mpz_cmp(g,N)<0){
                            mpz_t co;mpz_init(co);mpz_divexact(co,N,g);
                            gmp_printf("%Zd\n%Zd\n",g,co);mpz_clear(co);found=1;}
                        if(!found){mpz_add(g,X,Y_crt);mpz_mod(g,g,N);mpz_gcd(g,g,N);
                            if(mpz_cmp_ui(g,1)>0&&mpz_cmp(g,N)<0){
                                mpz_t co;mpz_init(co);mpz_divexact(co,N,g);
                                gmp_printf("%Zd\n%Zd\n",g,co);mpz_clear(co);found=1;}}
                        mpz_clear(Y_crt);
                        break;
                    }
                }
                mpz_clears(crt_v, crt_m, NULL);
            }

            if (0 && hp > 0) { /* DISABLED: old Hensel approach */
                fprintf(stderr,"  p=%lu, roots:", hp);
                for (int j = 0; j < d; j++) fprintf(stderr," %lu", hr[j]);
                fprintf(stderr,"\n");

                /* Number of Hensel lifts needed */
                int nbits = mpz_sizeinbase(N, 2) + 20;
                int pbits = 0; { unsigned long pp = hp; while (pp) { pbits++; pp >>= 1; } }
                int nlifts = 0, lb = pbits;
                while (lb < nbits) { lb *= 2; nlifts++; }

                /* Try all 2^d sign combos */
                int maxsc = 1 << d;
                for (int sc = 0; sc < maxsc && !found; sc++) {
                    mpz_t lr2[MAX_DEG], ls2[MAX_DEG], sv2[MAX_DEG], hmod, htmp;
                    for (int j=0;j<d;j++){mpz_init_set_ui(lr2[j],hr[j]);mpz_init(ls2[j]);mpz_init(sv2[j]);}
                    mpz_init_set_ui(hmod, hp); mpz_init(htmp);

                    /* Initial S_j and sqrt */
                    int ok = 1;
                    for (int j = 0; j < d && ok; j++) {
                        unsigned long long pv = 1;
                        for (int i = 0; i < dlen[di]; i++) {
                            long long t2 = ((da[i] % (long long)hp) + hp) % hp;
                            t2 = (t2 + hp - ((unsigned long long)db[i] % hp * hr[j]) % hp) % hp;
                            pv = pv * (unsigned long long)t2 % hp;
                        }
                        mpz_set_ui(sv2[j], (unsigned long)pv);
                        /* Sqrt mod hp */
                        unsigned long sq = 0;
                        if (pv == 0) { sq = 0; }
                        else {
                            unsigned long long r2=1,b2=pv,e2=(hp-1)/2; unsigned long long m2=hp;
                            while(e2){if(e2&1)r2=r2*b2%m2;b2=b2*b2%m2;e2>>=1;}
                            if (r2 != 1) { ok = 0; break; }
                            if (hp%4==3){r2=1;b2=pv;e2=(hp+1)/4;while(e2){if(e2&1)r2=r2*b2%m2;b2=b2*b2%m2;e2>>=1;}sq=(unsigned long)r2;}
                            else{unsigned long Q3=hp-1,S3=0;while(Q3%2==0){Q3/=2;S3++;}
                                unsigned long z2=2;for(;;){r2=1;b2=z2;e2=(hp-1)/2;while(e2){if(e2&1)r2=r2*b2%m2;b2=b2*b2%m2;e2>>=1;}if(r2==m2-1)break;z2++;}
                                unsigned long long M3=S3,c3,t3,R3;
                                r2=1;b2=z2;e2=Q3;while(e2){if(e2&1)r2=r2*b2%m2;b2=b2*b2%m2;e2>>=1;}c3=r2;
                                r2=1;b2=pv;e2=Q3;while(e2){if(e2&1)r2=r2*b2%m2;b2=b2*b2%m2;e2>>=1;}t3=r2;
                                r2=1;b2=pv;e2=(Q3+1)/2;while(e2){if(e2&1)r2=r2*b2%m2;b2=b2*b2%m2;e2>>=1;}R3=r2;
                                for(;;){if(t3==1){sq=(unsigned long)R3;break;}int i3=0;unsigned long long tt=t3;
                                while(tt!=1){tt=tt*tt%hp;i3++;}unsigned long long bb=c3;
                                for(int jj=0;jj<(int)M3-i3-1;jj++)bb=bb*bb%hp;M3=i3;c3=bb*bb%hp;t3=t3*c3%hp;R3=R3*bb%hp;}}
                        }
                        mpz_set_ui(ls2[j], (sc & (1<<j)) ? hp - sq : sq);
                    }
                    if (!ok) { for(int j=0;j<d;j++){mpz_clear(lr2[j]);mpz_clear(ls2[j]);mpz_clear(sv2[j]);}mpz_clear(hmod);mpz_clear(htmp);continue; }

                    /* Hensel lift nlifts times */
                    for (int lift = 0; lift < nlifts; lift++) {
                        mpz_t nm; mpz_init(nm);
                        mpz_mul(nm, hmod, hmod);
                        /* Lift roots */
                        for (int j = 0; j < d; j++) {
                            mpz_t fv,fpv,iv;mpz_inits(fv,fpv,iv,NULL);
                            mpz_set_ui(fv,0);for(int k=d;k>=0;k--){mpz_mul(fv,fv,lr2[j]);mpz_add(fv,fv,fc[k]);mpz_mod(fv,fv,nm);}
                            mpz_set_ui(fpv,0);for(int k=d;k>=1;k--){mpz_mul(fpv,fpv,lr2[j]);mpz_mul_ui(htmp,fc[k],k);mpz_add(fpv,fpv,htmp);mpz_mod(fpv,fpv,nm);}
                            if(mpz_invert(iv,fpv,nm)){mpz_mul(htmp,fv,iv);mpz_sub(lr2[j],lr2[j],htmp);mpz_mod(lr2[j],lr2[j],nm);}
                            mpz_clears(fv,fpv,iv,NULL);
                        }
                        /* Recompute S_j mod nm */
                        for (int j = 0; j < d; j++) {
                            mpz_set_ui(sv2[j],1);
                            for(int i=0;i<dlen[di];i++){mpz_set_si(htmp,da[i]);mpz_submul_ui(htmp,lr2[j],db[i]);mpz_mod(htmp,htmp,nm);mpz_mul(sv2[j],sv2[j],htmp);mpz_mod(sv2[j],sv2[j],nm);}
                        }
                        /* Lift sqrt */
                        for (int j = 0; j < d; j++) {
                            mpz_t t22,df,tt,iv;mpz_inits(t22,df,tt,iv,NULL);
                            mpz_mul(t22,ls2[j],ls2[j]);mpz_sub(df,t22,sv2[j]);mpz_mod(df,df,nm);
                            mpz_mul_ui(tt,ls2[j],2);mpz_mod(tt,tt,nm);
                            if(mpz_invert(iv,tt,nm)){mpz_mul(htmp,df,iv);mpz_sub(ls2[j],ls2[j],htmp);mpz_mod(ls2[j],ls2[j],nm);}
                            mpz_clears(t22,df,tt,iv,NULL);
                        }
                        mpz_set(hmod, nm); mpz_clear(nm);
                    }

                    /* Verify: T_j^2 = S_j mod hmod for each j */
                    if (sc == 0) {
                        for (int j = 0; j < d; j++) {
                            mpz_t vt2; mpz_init(vt2);
                            mpz_mul(vt2, ls2[j], ls2[j]);
                            mpz_mod(vt2, vt2, hmod);
                            mpz_t vs; mpz_init(vs);
                            mpz_mod(vs, sv2[j], hmod);
                            int match = (mpz_cmp(vt2, vs) == 0);
                            if (!match && j == 0)
                                fprintf(stderr,"  VERIFY: T_%d^2 != S_%d mod hmod (%zu bits)!\n",
                                        j, j, mpz_sizeinbase(hmod, 2));
                            mpz_clears(vt2, vs, NULL);
                        }
                    }

                    /* Verify Lagrange at r_0: should give ls2[0] */
                    if (sc == 0) {
                        mpz_t vl; mpz_init_set_ui(vl, 0);
                        mpz_t vr0; mpz_init(vr0); mpz_mod(vr0, lr2[0], hmod);
                        for (int j = 0; j < d; j++) {
                            mpz_t vn, vd, vi; mpz_inits(vn, vd, vi, NULL);
                            mpz_set(vn, ls2[j]); mpz_set_ui(vd, 1);
                            for (int k = 0; k < d; k++) {
                                if (k == j) continue;
                                mpz_sub(htmp, vr0, lr2[k]); mpz_mod(htmp, htmp, hmod);
                                mpz_mul(vn, vn, htmp); mpz_mod(vn, vn, hmod);
                                mpz_sub(htmp, lr2[j], lr2[k]); mpz_mod(htmp, htmp, hmod);
                                mpz_mul(vd, vd, htmp); mpz_mod(vd, vd, hmod);
                            }
                            if (mpz_invert(vi, vd, hmod)) {
                                mpz_mul(vn, vn, vi); mpz_add(vl, vl, vn); mpz_mod(vl, vl, hmod);
                            } else {
                                fprintf(stderr,"  LAGRANGE: denominator not invertible!\n");
                            }
                            mpz_clears(vn, vd, vi, NULL);
                        }
                        int vmatch = (mpz_cmp(vl, ls2[0]) == 0);
                        fprintf(stderr,"  Lagrange(r_0) = T_0? %d\n", vmatch);
                        mpz_clears(vl, vr0, NULL);
                    }

                    /* Lagrange: T(m) mod hmod */
                    mpz_t mm,Y3;mpz_init(mm);mpz_init_set_ui(Y3,0);
                    mpz_mod(mm, poly.m, hmod);
                    for (int j = 0; j < d; j++) {
                        mpz_t n3,d3,iv;mpz_inits(n3,d3,iv,NULL);
                        mpz_set(n3,ls2[j]);mpz_set_ui(d3,1);
                        for(int k=0;k<d;k++){if(k==j)continue;
                            mpz_sub(htmp,mm,lr2[k]);mpz_mod(htmp,htmp,hmod);mpz_mul(n3,n3,htmp);mpz_mod(n3,n3,hmod);
                            mpz_sub(htmp,lr2[j],lr2[k]);mpz_mod(htmp,htmp,hmod);mpz_mul(d3,d3,htmp);mpz_mod(d3,d3,hmod);}
                        if(mpz_invert(iv,d3,hmod)){mpz_mul(n3,n3,iv);mpz_add(Y3,Y3,n3);mpz_mod(Y3,Y3,hmod);}
                        mpz_clears(n3,d3,iv,NULL);
                    }
                    mpz_mod(Y3, Y3, N);

                    /* Verify T(m)^2 = S(m) mod hmod where S(m) = prod(a-b*m) mod hmod */
                    if (sc == 0) {
                        mpz_t sm, tm2;
                        mpz_inits(sm, tm2, NULL);
                        /* S(m) mod hmod */
                        mpz_set_ui(sm, 1);
                        mpz_t mm2; mpz_init(mm2); mpz_mod(mm2, poly.m, hmod);
                        for (int i = 0; i < dlen[di]; i++) {
                            mpz_t t2; mpz_init(t2);
                            mpz_set_si(t2, da[i]);
                            mpz_submul_ui(t2, mm2, db[i]);
                            mpz_mod(t2, t2, hmod);
                            mpz_mul(sm, sm, t2);
                            mpz_mod(sm, sm, hmod);
                            mpz_clear(t2);
                        }
                        mpz_clear(mm2);
                        /* T(m)^2 mod hmod */
                        mpz_mul(tm2, Y3, Y3); /* Y3 = T(m) mod hmod before reducing to mod N */
                        /* Wait, Y3 was already reduced mod N. Need pre-reduction value. */
                        /* Actually Y3 was computed as mod hmod then reduced mod N later. Let me recompute. */
                        /* The Lagrange step computes Y3 mod hmod, then we reduce mod N. */
                        /* At this point Y3 is mod hmod, not yet reduced. */
                        mpz_mod(tm2, tm2, hmod);
                        /* Also check: does Y3 equal one of the T_j values at the corresponding root? */
                        /* Check T(r_0) = T_0 */
                        mpz_t tr0; mpz_init_set_ui(tr0, 0);
                        mpz_t mm3; mpz_init(mm3); mpz_mod(mm3, lr2[0], hmod);
                        for (int j = d-1; j >= 0; j--) {
                            /* Lagrange at lr2[0] should give ls2[0] */
                        }
                        /* Actually just check Y3 size */
                        fprintf(stderr,"  Y3=%zu bits, hmod=%zu bits, sm=%zu bits\n",
                                mpz_sizeinbase(Y3, 2), mpz_sizeinbase(hmod, 2), mpz_sizeinbase(sm, 2));
                        int match = (mpz_cmp(tm2, sm) == 0);
                        fprintf(stderr,"  T(m)^2 = S(m) mod p^e? %d\n", match);
                        mpz_clears(tr0, mm3, NULL);
                        mpz_clears(sm, tm2, NULL);
                    }

                    /* Verify lifted roots are still roots of f mod hmod */
                    if (sc == 0) {
                        for (int j = 0; j < d; j++) {
                            mpz_t fv2; mpz_init(fv2); mpz_set_ui(fv2, 0);
                            for (int k = d; k >= 0; k--) {
                                mpz_mul(fv2, fv2, lr2[j]);
                                mpz_add(fv2, fv2, fc[k]);
                                mpz_mod(fv2, fv2, hmod);
                            }
                            if (mpz_sgn(fv2) != 0)
                                fprintf(stderr,"  ROOT LIFT FAIL: f(r_%d) != 0 mod p^e!\n", j);
                            mpz_clear(fv2);
                        }
                    }

                    /* Verify: compute prod(a-b*m) mod N and compare with X^2 and Y3^2 */
                    if (sc == 0) {
                        mpz_t vy2, vx2, vprod;
                        mpz_inits(vy2, vx2, vprod, NULL);
                        mpz_mul(vy2, Y3, Y3); mpz_mod(vy2, vy2, N);
                        mpz_mul(vx2, X, X); mpz_mod(vx2, vx2, N);

                        /* Compute prod(a-b*m) mod N directly */
                        mpz_set_ui(vprod, 1);
                        for (int i = 0; i < dlen[di]; i++) {
                            mpz_t term; mpz_init(term);
                            mpz_set_si(term, da[i]);
                            mpz_submul_ui(term, poly.m, db[i]);
                            mpz_mod(term, term, N);
                            mpz_mul(vprod, vprod, term);
                            mpz_mod(vprod, vprod, N);
                            mpz_clear(term);
                        }

                        int y2_eq_prod = (mpz_cmp(vy2, vprod) == 0);
                        int x2_eq_prod = (mpz_cmp(vx2, vprod) == 0);
                        fprintf(stderr,"  dep %d sc %d: Y^2=prod? %d, X^2=prod? %d, Y^2=X^2? %d\n",
                                di, sc, y2_eq_prod, x2_eq_prod, mpz_cmp(vy2,vx2)==0);

                        mpz_clears(vy2, vx2, vprod, NULL);
                    }

                    /* Check */
                    mpz_sub(g,X,Y3);mpz_mod(g,g,N);mpz_gcd(g,g,N);
                    if(mpz_cmp_ui(g,1)>0&&mpz_cmp(g,N)<0){mpz_t co;mpz_init(co);mpz_divexact(co,N,g);gmp_printf("%Zd\n%Zd\n",g,co);mpz_clear(co);found=1;}
                    if(!found){mpz_add(g,X,Y3);mpz_mod(g,g,N);mpz_gcd(g,g,N);
                        if(mpz_cmp_ui(g,1)>0&&mpz_cmp(g,N)<0){mpz_t co;mpz_init(co);mpz_divexact(co,N,g);gmp_printf("%Zd\n%Zd\n",g,co);mpz_clear(co);found=1;}}

                    mpz_clears(Y3,mm,NULL);
                    for(int j=0;j<d;j++){mpz_clear(lr2[j]);mpz_clear(ls2[j]);mpz_clear(sv2[j]);}
                    mpz_clear(hmod);mpz_clear(htmp);
                }

                if (!found) fprintf(stderr,"  Hensel: no factor from %d sign combos\n", maxsc);
            } else {
                fprintf(stderr,"  No suitable Hensel prime found\n");
            }
        }

        /* Skip CRT approach (too slow, exponential blowup) */
        if (0 && !found) {
            fprintf(stderr,"  Trying CRT algebraic sqrt (d=%d)...\n", d);
            /* Find primes q where f(x) mod q has exactly d distinct roots */
            /* Use small primes for brute-force root finding */
            #define MAX_CRT_PRIMES 30
            unsigned long crt_q[MAX_CRT_PRIMES];
            unsigned long crt_roots[MAX_CRT_PRIMES][MAX_DEG];
            int crt_nroots[MAX_CRT_PRIMES];
            int ncrt = 0;

            /* Collect CRT primes */
            for (unsigned long q = 1000; q < 50000 && ncrt < MAX_CRT_PRIMES; q++) {
                mpz_set_ui(tmp, q);
                if (!mpz_probab_prime_p(tmp, 3)) continue;
                /* Find roots of f mod q by brute force */
                int nr = 0;
                unsigned long roots_q[MAX_DEG+1];
                for (unsigned long x = 0; x < q && nr <= d; x++) {
                    unsigned long long val = 0;
                    for (int j = d; j >= 0; j--) {
                        unsigned long cj = mpz_fdiv_ui(fc[j], q);
                        val = (val * x + cj) % q;
                    }
                    if (val == 0) roots_q[nr++] = x;
                }
                if (nr != d) continue; /* need exactly d distinct roots */
                crt_q[ncrt] = q;
                for (int j = 0; j < d; j++) crt_roots[ncrt][j] = roots_q[j];
                crt_nroots[ncrt] = d;
                ncrt++;
            }
            fprintf(stderr,"  Found %d CRT primes with %d roots each\n", ncrt, d);

            if (ncrt >= 3) {
                /* For each CRT prime: compute S_j = prod(a_i - b_i*r_j) mod q
                 * Then T_j = sqrt(S_j) mod q (2 choices each) */

                /* Enumerate sign combinations using first few primes */
                int max_sign_combos = 1;
                for (int i = 0; i < d; i++) max_sign_combos *= 2;

                /* Compute T(m) mod q for each CRT prime, for each sign combo */
                unsigned long *tm_per_q = malloc(ncrt * max_sign_combos * sizeof(unsigned long));

                for (int qi = 0; qi < ncrt; qi++) {
                    unsigned long q = crt_q[qi];
                    unsigned long mmod = mpz_fdiv_ui(poly.m, q);

                    /* S_j = prod(a_i - b_i*r_j) mod q */
                    unsigned long S_j[MAX_DEG];
                    for (int j = 0; j < d; j++) {
                        unsigned long long prod_val = 1;
                        for (int i = 0; i < dlen[di]; i++) {
                            long long aval = da[i];
                            unsigned long long bval = db[i];
                            long long term = ((aval % (long long)q) + q) % q;
                            term = (term + q - (bval % q * crt_roots[qi][j]) % q) % q;
                            prod_val = prod_val * (unsigned long long)term % q;
                        }
                        S_j[j] = (unsigned long)prod_val;
                    }

                    /* For each sign combo, compute T(m) mod q via Lagrange */
                    for (int sc = 0; sc < max_sign_combos; sc++) {
                        unsigned long T_j[MAX_DEG];
                        for (int j = 0; j < d; j++) {
                            /* sqrt(S_j) mod q */
                            unsigned long s = S_j[j];
                            /* Tonelli-Shanks */
                            unsigned long long r2 = 1, b2 = s, e2 = (q+1)/4;
                            if (q % 4 != 3) {
                                /* Full Tonelli-Shanks for q != 3 mod 4 */
                                unsigned long Q2 = q-1, S2 = 0;
                                while (Q2%2==0) { Q2/=2; S2++; }
                                unsigned long z = 2;
                                for (;;) { unsigned long long rr=1,bb=z,ee=(q-1)/2;
                                    while (ee) { if (ee&1) rr=rr*bb%q; bb=bb*bb%q; ee>>=1; }
                                    if (rr==q-1) break; z++; }
                                unsigned long long c2,t2,R2,M3=S2;
                                r2=1;b2=z;e2=Q2; while (e2) { if (e2&1) r2=r2*b2%q; b2=b2*b2%q; e2>>=1; } c2=r2;
                                r2=1;b2=s;e2=Q2; while (e2) { if (e2&1) r2=r2*b2%q; b2=b2*b2%q; e2>>=1; } t2=r2;
                                r2=1;b2=s;e2=(Q2+1)/2; while (e2) { if (e2&1) r2=r2*b2%q; b2=b2*b2%q; e2>>=1; } R2=r2;
                                while (t2!=1) { int ii=0; unsigned long long tt=t2;
                                    while (tt!=1) { tt=tt*tt%q; ii++; }
                                    unsigned long long bb2=c2; for (int jj=0;jj<(int)M3-ii-1;jj++) bb2=bb2*bb2%q;
                                    M3=ii; c2=bb2*bb2%q; t2=t2*c2%q; R2=R2*bb2%q; }
                                r2 = R2;
                            } else {
                                r2=1;b2=s;e2=(q+1)/4;
                                while (e2) { if (e2&1) r2=r2*b2%q; b2=b2*b2%q; e2>>=1; }
                            }
                            /* Apply sign */
                            T_j[j] = (sc & (1<<j)) ? (unsigned long)(q - r2) : (unsigned long)r2;
                        }

                        /* Lagrange interpolation: T(m) mod q */
                        unsigned long long result = 0;
                        for (int j = 0; j < d; j++) {
                            unsigned long long num = T_j[j];
                            unsigned long long den = 1;
                            for (int k = 0; k < d; k++) {
                                if (k == j) continue;
                                num = num * ((mmod + q - crt_roots[qi][k]) % q) % q;
                                den = den * ((crt_roots[qi][j] + q - crt_roots[qi][k]) % q) % q;
                            }
                            /* Invert den */
                            unsigned long long inv = 1, bb = den, ee = q - 2;
                            while (ee) { if (ee&1) inv=inv*bb%q; bb=bb*bb%q; ee>>=1; }
                            result = (result + num * inv % q) % q;
                        }
                        tm_per_q[qi * max_sign_combos + sc] = (unsigned long)result;
                    }
                }

                /* Now try CRT combinations.
                 * Use first 3 primes: 8^3 = 512 combos for d=3.
                 * For each combo, compute T(m) mod (q1*q2*q3) via CRT,
                 * then check gcd(X - T(m), N). */
                fprintf(stderr,"  Trying %d^%d = %d CRT combinations...\n",
                        max_sign_combos, (ncrt < 3 ? ncrt : 3),
                        ncrt < 3 ? 0 : max_sign_combos * max_sign_combos * max_sign_combos);

                /* Need product of primes > N. For 30d, ~10 primes of ~5000 needed.
                 * To limit exponential blowup (8^k), check after each prime addition
                 * and only continue with candidates that give Y^2 = X^2 mod N */
                int nprimes_to_use = ncrt;
                /* Use iterative CRT with pruning */

                /* Start with all sign combos for first prime.
                 * For subsequent primes: for each existing candidate,
                 * check which of the 8 new values it's consistent with (mod q).
                 * This keeps the candidate count bounded. */
                mpz_t *candidates = malloc(max_sign_combos * sizeof(mpz_t));
                int ncand = max_sign_combos;
                for (int i = 0; i < max_sign_combos; i++)
                    mpz_init_set_ui(candidates[i], tm_per_q[0 * max_sign_combos + i]);
                mpz_t mod_acc;
                mpz_init_set_ui(mod_acc, crt_q[0]);

                for (int qi = 1; qi < nprimes_to_use && !found; qi++) {
                    unsigned long q = crt_q[qi];

                    /* For each candidate, find which of the 8 sign combos it matches */
                    int max_new = ncand * max_sign_combos;
                    if (max_new > 50000) { /* too many - just keep best candidates */
                        max_new = 50000;
                    }
                    mpz_t *new_cand = malloc(max_new * sizeof(mpz_t));
                    int nc2 = 0;

                    unsigned long ma_inv_q;
                    { unsigned long long inv=1, bb=mpz_fdiv_ui(mod_acc,q), ee=q-2;
                      while (ee) { if (ee&1) inv=inv*bb%q; bb=bb*bb%q; ee>>=1; }
                      ma_inv_q = (unsigned long)inv; }

                    for (int ci = 0; ci < ncand && nc2 < max_new; ci++) {
                        unsigned long a_mod_q = mpz_fdiv_ui(candidates[ci], q);
                        /* Check all 8 sign combos */
                        for (int sc = 0; sc < max_sign_combos; sc++) {
                            unsigned long tm_q = tm_per_q[qi * max_sign_combos + sc];
                            /* CRT: x = candidates[ci] + mod_acc * ((tm_q - a_mod_q) * inv mod q) */
                            long diff = (long)tm_q - (long)a_mod_q;
                            if (diff < 0) diff += q;
                            unsigned long t_val = (unsigned long long)((unsigned long)diff) * ma_inv_q % q;
                            if (nc2 >= max_new) break;
                            mpz_init(new_cand[nc2]);
                            mpz_mul_ui(new_cand[nc2], mod_acc, t_val);
                            mpz_add(new_cand[nc2], new_cand[nc2], candidates[ci]);
                            nc2++;
                        }
                    }

                    for (int ci = 0; ci < ncand; ci++) mpz_clear(candidates[ci]);
                    free(candidates);
                    candidates = new_cand;
                    ncand = nc2;
                    mpz_mul_ui(mod_acc, mod_acc, q);

                    /* If mod_acc > N, try each candidate */
                    if (mpz_cmp(mod_acc, N) > 0) {
                        fprintf(stderr,"  CRT > N after %d primes, %d candidates\n", qi+1, ncand);
                        for (int ci = 0; ci < ncand && !found; ci++) {
                            mpz_mod(tmp, candidates[ci], N);
                            mpz_sub(g, X, tmp); mpz_mod(g,g,N); mpz_gcd(g,g,N);
                            if (mpz_cmp_ui(g,1)>0 && mpz_cmp(g,N)<0) {
                                mpz_t co; mpz_init(co); mpz_divexact(co,N,g);
                                gmp_printf("%Zd\n%Zd\n",g,co); mpz_clear(co); found=1;
                            }
                            if (!found) {
                                mpz_add(g,X,tmp); mpz_mod(g,g,N); mpz_gcd(g,g,N);
                                if (mpz_cmp_ui(g,1)>0 && mpz_cmp(g,N)<0) {
                                    mpz_t co; mpz_init(co); mpz_divexact(co,N,g);
                                    gmp_printf("%Zd\n%Zd\n",g,co); mpz_clear(co); found=1;
                                }
                            }
                        }
                        if (!found) {
                            fprintf(stderr,"  No factor from %d candidates, trying next dep\n", ncand);
                        }
                        break; /* move to next dependency */
                    }
                }

                for (int ci = 0; ci < ncand; ci++) mpz_clear(candidates[ci]);
                free(candidates);
                mpz_clear(mod_acc);
                free(tm_per_q);
            }
        }

        mpz_clear(g);
        free(da);free(db);free(rex);free(aex);
    }

    if(!found){fprintf(stderr,"Failed\n");return 1;}
    fprintf(stderr,"Done: %.1fs\n",elapsed());
    for(int i=0;i<=d;i++) mpz_clear(fc[i]);
    mpz_clears(X,tmp,rn,an,N,NULL);poly_clear(&poly);fb_free(rfb);fb_free(afb);
    gmp_randclear(rstate);
    return 0;
}
