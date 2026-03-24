/* Minimal QS to debug the extraction issue */
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

static double wt(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9;}

static int*sieve_p(int B,int*c){char*s=(char*)calloc(B+1,1);for(int i=2;i<=B;i++)s[i]=1;
for(int i=2;(long long)i*i<=B;i++)if(s[i])for(int j=i*i;j<=B;j+=i)s[j]=0;
int n=0;for(int i=2;i<=B;i++)if(s[i])n++;int*p=(int*)malloc(n*sizeof(int));int k=0;
for(int i=2;i<=B;i++)if(s[i])p[k++]=i;*c=n;free(s);return p;}

static long msqrt(long n,long p){if(p==2)return n&1;n=((n%p)+p)%p;if(n==0)return 0;
mpz_t a,m,r;mpz_inits(a,m,r,NULL);mpz_set_si(a,n);mpz_set_si(m,p);
mpz_powm_ui(r,a,(p-1)/2,m);if(mpz_cmp_ui(r,1)!=0){mpz_clears(a,m,r,NULL);return-1;}
if((p&3)==3){mpz_powm_ui(r,a,(p+1)/4,m);long R=mpz_get_si(r);mpz_clears(a,m,r,NULL);return R;}
long Q=p-1,S=0;while((Q&1)==0){Q>>=1;S++;}long z;
for(z=2;z<p;z++){mpz_set_si(a,z);mpz_powm_ui(r,a,(p-1)/2,m);if(mpz_cmp_si(r,p-1)==0)break;}
mpz_t c,t,R,b;mpz_inits(c,t,R,b,NULL);
mpz_set_si(a,z);mpz_powm_ui(c,a,Q,m);mpz_set_si(a,n);mpz_powm_ui(t,a,Q,m);
mpz_powm_ui(R,a,(Q+1)/2,m);long Mv=S;
while(1){if(mpz_cmp_ui(t,1)==0){long res=mpz_get_si(R);mpz_clears(a,m,r,c,t,R,b,NULL);return res;}
long i=0;mpz_set(b,t);while(mpz_cmp_ui(b,1)!=0){mpz_mul(b,b,b);mpz_mod(b,b,m);i++;
if(i>=Mv){mpz_clears(a,m,r,c,t,R,b,NULL);return-1;}}
mpz_set(b,c);for(long j=0;j<Mv-i-1;j++){mpz_mul(b,b,b);mpz_mod(b,b,m);}
Mv=i;mpz_mul(c,b,b);mpz_mod(c,c,m);mpz_mul(t,t,c);mpz_mod(t,t,m);
mpz_mul(R,R,b);mpz_mod(R,R,m);}}

typedef struct{int*p;long*sq;int n;}fb_t;
static fb_t bfb(const mpz_t N,int B){int np;int*all=sieve_p(B,&np);
fb_t f;f.p=(int*)malloc(np*sizeof(int));f.sq=(long*)malloc(np*sizeof(long));f.n=0;
for(int i=0;i<np;i++){long nm=mpz_fdiv_ui(N,all[i]);long s=msqrt(nm,all[i]);
if(s>=0){f.p[f.n]=all[i];f.sq[f.n]=s;f.n++;}}free(all);return f;}

typedef struct{mpz_t sv;int*e;}rel_t;
typedef struct{rel_t*d;int n,cap,vl;}rs_t;
static void ri(rs_t*s,int vl){s->vl=vl;s->n=0;s->cap=8192;s->d=(rel_t*)malloc(s->cap*sizeof(rel_t));}
static void ra(rs_t*s,const mpz_t sv,const int*e){
if(s->n>=s->cap){s->cap*=2;s->d=(rel_t*)realloc(s->d,s->cap*sizeof(rel_t));}
mpz_init_set(s->d[s->n].sv,sv);s->d[s->n].e=(int*)malloc(s->vl*sizeof(int));
memcpy(s->d[s->n].e,e,s->vl*sizeof(int));s->n++;}

int main(int argc,char**argv){
    if(argc<2)return 1;
    mpz_t N,factor;mpz_inits(N,factor,NULL);
    mpz_set_str(N,argv[1],10);
    int nd=strlen(argv[1]);double t0=wt();
    fprintf(stderr,"SimpleQS: %d digits\n",nd);

    double ln_N=nd*log(10);double Le=sqrt(ln_N*log(ln_N));
    int B=(int)exp(0.55*Le);if(B<500)B=500;if(B>5000000)B=5000000;
    double lm=1.05*Le;int M=(lm>16)?10000000:(int)exp(lm);
    if(M<20000)M=20000;if(M>10000000)M=10000000;

    fb_t fb=bfb(N,B);
    fprintf(stderr,"FB: %d primes, B=%d, M=%d\n",fb.n,B,M);
    int target=fb.n+30;

    /* QS sieve for FULLY SMOOTH relations only */
    mpz_t m,tmp,val,res;mpz_inits(m,tmp,val,res,NULL);
    mpz_sqrt(m,N);mpz_add_ui(m,m,1);
    int slen=2*M+1;int nb=mpz_sizeinbase(N,2);
    double l2m=1.0+log2(M)+nb/2.0;
    int thr=(int)((l2m)*1024*0.73); /* strict threshold for smooth only */
    int*sv=(int*)calloc(slen,sizeof(int));
    for(int i=0;i<fb.n;i++){
        int p=fb.p[i];long sq=fb.sq[i];long mp=mpz_fdiv_ui(m,p);
        int lp=(int)(log2(p)*1024);
        if(p==2){long r=((sq-mp)%2+2)%2;long s=((r+M)%2+2)%2;
            for(long j=s;j<slen;j+=2)sv[j]+=lp;continue;}
        long r1=((sq-mp)%p+p)%p;long r2=((-sq-mp)%p+p)%p;
        long s1=((r1+M)%p+p)%p;long s2=((r2+M)%p+p)%p;
        for(long j=s1;j<slen;j+=p)sv[j]+=lp;
        if(r1!=r2)for(long j=s2;j<slen;j+=p)sv[j]+=lp;
    }

    rs_t rels;ri(&rels,fb.n+1);
    int*exps=(int*)calloc(fb.n+1,sizeof(int));
    for(int j=0;j<slen;j++){
        if(sv[j]<thr)continue;
        long x=(long)j-M;
        mpz_set_si(tmp,x);mpz_add(val,tmp,m);
        mpz_mul(res,val,val);mpz_sub(res,res,N);
        memset(exps,0,(fb.n+1)*sizeof(int));
        int sign=0;if(mpz_sgn(res)<0){sign=1;mpz_neg(res,res);}exps[0]=sign;
        mpz_set(tmp,res);
        for(int i=0;i<fb.n;i++){unsigned long p=fb.p[i];
            while(mpz_divisible_ui_p(tmp,p)){mpz_divexact_ui(tmp,tmp,p);exps[i+1]++;}}
        if(mpz_cmp_ui(tmp,1)==0){
            mpz_set_si(tmp,x);mpz_add(val,tmp,m);mpz_mod(val,val,N);
            ra(&rels,val,exps);
        }
    }
    free(sv);free(exps);
    fprintf(stderr,"%d fully smooth relations (need %d) %.1fs\n",rels.n,target,wt()-t0);

    if(rels.n<target){fprintf(stderr,"FAIL: not enough\n");return 1;}

    /* Matrix */
    typedef unsigned long wt_t;
    #define WB (sizeof(wt_t)*8)
    int nr=fb.n+1,nc=rels.n;
    int nw=(nc+WB-1)/WB;
    wt_t**mat=(wt_t**)calloc(nr,sizeof(wt_t*));
    for(int i=0;i<nr;i++){mat[i]=(wt_t*)calloc(nw,sizeof(wt_t));
        for(int j=0;j<nc;j++)if(rels.d[j].e[i]&1)mat[i][j/WB]|=(1UL<<(j%WB));}
    int*pc=(int*)malloc(nr*sizeof(int));char*isp=(char*)calloc(nc,1);
    for(int i=0;i<nr;i++){pc[i]=-1;
        for(int j=0;j<nc;j++){if(!isp[j]&&(mat[i][j/WB]&(1UL<<(j%WB)))){
            pc[i]=j;isp[j]=1;
            for(int k=0;k<nr;k++)if(k!=i&&(mat[k][j/WB]&(1UL<<(j%WB))))
                for(int w=0;w<nw;w++)mat[k][w]^=mat[i][w];
            break;}}}

    /* Try deps */
    int ok=0;
    for(int jj=0;jj<nc&&!ok;jj++){
        if(isp[jj])continue;
        int*dep=(int*)malloc((nr+1)*sizeof(int));int cnt=0;dep[cnt++]=jj;
        for(int i=0;i<nr;i++)if(mat[i][jj/WB]&(1UL<<(jj%WB)))if(pc[i]>=0)dep[cnt++]=pc[i];

        /* Extract */
        mpz_t x,y,t2;mpz_inits(x,y,t2,NULL);mpz_set_ui(x,1);
        int*tot=(int*)calloc(fb.n+1,sizeof(int));
        for(int i=0;i<cnt;i++){
            mpz_mul(x,x,rels.d[dep[i]].sv);mpz_mod(x,x,N);
            for(int k=0;k<fb.n+1;k++)tot[k]+=rels.d[dep[i]].e[k];}
        int alleven=1;for(int k=0;k<fb.n+1;k++)if(tot[k]&1)alleven=0;
        if(!alleven){free(tot);free(dep);mpz_clears(x,y,t2,NULL);continue;}

        mpz_set_ui(y,1);
        for(int k=0;k<fb.n;k++){int e=tot[k+1]/2;if(e>0){
            mpz_set_ui(t2,fb.p[k]);mpz_powm_ui(t2,t2,e,N);
            mpz_mul(y,y,t2);mpz_mod(y,y,N);}}
        if((tot[0]/2)&1)mpz_sub(y,N,y);

        mpz_sub(t2,x,y);mpz_gcd(factor,t2,N);
        if(mpz_cmp_ui(factor,1)>0&&mpz_cmp(factor,N)<0)ok=1;
        if(!ok){mpz_add(t2,x,y);mpz_gcd(factor,t2,N);
            if(mpz_cmp_ui(factor,1)>0&&mpz_cmp(factor,N)<0)ok=1;}

        free(tot);free(dep);mpz_clears(x,y,t2,NULL);
    }

    if(ok){mpz_t cf;mpz_init(cf);mpz_tdiv_q(cf,N,factor);
        gmp_printf("FACTOR:%Zd\n",factor);
        gmp_fprintf(stderr,"%Zd x %Zd (%.3fs)\n",factor,cf,wt()-t0);mpz_clear(cf);}
    else fprintf(stderr,"FAIL (%.1fs)\n",wt()-t0);

    for(int i=0;i<nr;i++)free(mat[i]);free(mat);free(pc);free(isp);
    mpz_clears(N,factor,m,tmp,val,res,NULL);
    return ok?0:1;
}
