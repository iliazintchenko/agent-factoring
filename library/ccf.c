/*
 * Cofactor Collision Factoring (CCF) v3
 *
 * Combines: multi-multiplier candidate generation, Bernstein batch smooth
 * detection (primorial product tree), batch GCD cofactor collision matching,
 * single/double large prime matching, GF(2) linear algebra.
 *
 * Usage: ./ccf <N>
 * Output: FACTOR: <p>
 * Compile: gcc -O3 -o ccf ccf.c -lgmp -lm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

/* ---- Primes ---- */
static unsigned int *all_primes;
static int n_all_primes;
static unsigned int *fb;
static int fb_size;
static mpz_t N_g;

static void gen_primes(unsigned int limit) {
    char *s = calloc(limit + 1, 1);
    int cnt = 0;
    for (unsigned int i = 2; i <= limit; i++) {
        if (!s[i]) { cnt++; for (unsigned long j=(unsigned long)i*i; j<=limit; j+=i) s[j]=1; }
    }
    all_primes = malloc(cnt * sizeof(unsigned int));
    n_all_primes = 0;
    memset(s, 0, limit + 1);
    for (unsigned int i = 2; i <= limit; i++) {
        if (!s[i]) { all_primes[n_all_primes++] = i; for (unsigned long j=(unsigned long)i*i; j<=limit; j+=i) s[j]=1; }
    }
    free(s);
}

static int legendre_ui(mpz_t a, unsigned int p) {
    unsigned long r = mpz_fdiv_ui(a, p);
    if (r == 0) return 0;
    unsigned long exp = (p - 1) / 2, res = 1;
    while (exp > 0) { if (exp & 1) res = (res * r) % p; r = (r * r) % p; exp >>= 1; }
    return res == 1 ? 1 : -1;
}

/* ---- Product tree (level-based) ---- */
typedef struct {
    mpz_t *data;      /* all nodes flat */
    int *lev_off;      /* offset of each level */
    int *lev_sz;       /* size of each level */
    int nlevs, nleaves;
} ptree;

static ptree *pt_build(mpz_t *leaves, int n) {
    ptree *t = malloc(sizeof(ptree));
    t->nleaves = n;
    int levs = 1, sz = n;
    while (sz > 1) { sz = (sz+1)/2; levs++; }
    t->nlevs = levs;
    t->lev_off = malloc(levs * sizeof(int));
    t->lev_sz = malloc(levs * sizeof(int));
    int total = 0; sz = n;
    for (int l = 0; l < levs; l++) { t->lev_off[l]=total; t->lev_sz[l]=sz; total+=sz; sz=(sz+1)/2; }
    t->data = malloc(total * sizeof(mpz_t));
    for (int i = 0; i < total; i++) mpz_init(t->data[i]);
    for (int i = 0; i < n; i++) mpz_set(t->data[i], leaves[i]);
    for (int l = 1; l < levs; l++) {
        int ps = t->lev_sz[l-1], po = t->lev_off[l-1], co = t->lev_off[l], cs = t->lev_sz[l];
        for (int i = 0; i < cs; i++) {
            if (2*i+1 < ps) mpz_mul(t->data[co+i], t->data[po+2*i], t->data[po+2*i+1]);
            else mpz_set(t->data[co+i], t->data[po+2*i]);
        }
    }
    return t;
}

static void pt_rem(ptree *t, mpz_t z, mpz_t *out) {
    int total = t->lev_off[t->nlevs-1] + t->lev_sz[t->nlevs-1];
    mpz_t *r = malloc(total * sizeof(mpz_t));
    for (int i = 0; i < total; i++) mpz_init(r[i]);
    mpz_mod(r[t->lev_off[t->nlevs-1]], z, t->data[t->lev_off[t->nlevs-1]]);
    for (int l = t->nlevs-1; l >= 1; l--) {
        int co = t->lev_off[l], po = t->lev_off[l-1], ps = t->lev_sz[l-1], cs = t->lev_sz[l];
        for (int i = 0; i < cs; i++) {
            int li = po+2*i, ri = po+2*i+1;
            mpz_mod(r[li], r[co+i], t->data[li]);
            if (2*i+1 < ps) mpz_mod(r[ri], r[co+i], t->data[ri]);
        }
    }
    for (int i = 0; i < t->nleaves; i++) mpz_set(out[i], r[i]);
    for (int i = 0; i < total; i++) mpz_clear(r[i]);
    free(r);
}

static void pt_free(ptree *t) {
    int total = t->lev_off[t->nlevs-1] + t->lev_sz[t->nlevs-1];
    for (int i = 0; i < total; i++) mpz_clear(t->data[i]);
    free(t->data); free(t->lev_off); free(t->lev_sz); free(t);
}

/* ---- Batch smooth detection ---- */
/* Returns smooth_part[i] = gcd(primorial, candidate[i]) */
static void batch_smooth(mpz_t *cands, int n, mpz_t primorial, mpz_t *smooth) {
    ptree *t = pt_build(cands, n);
    mpz_t *rems = malloc(n * sizeof(mpz_t));
    for (int i = 0; i < n; i++) mpz_init(rems[i]);
    pt_rem(t, primorial, rems);
    for (int i = 0; i < n; i++) mpz_gcd(smooth[i], rems[i], cands[i]);
    for (int i = 0; i < n; i++) mpz_clear(rems[i]);
    free(rems);
    pt_free(t);
}

/* ---- Trial division for exponent extraction ---- */
static int trial_div(mpz_t val, int *exps, mpz_t cof) {
    mpz_set(cof, val);
    memset(exps, 0, fb_size * sizeof(int));
    for (int i = 0; i < fb_size; i++) {
        while (mpz_divisible_ui_p(cof, fb[i])) { mpz_divexact_ui(cof, cof, fb[i]); exps[i]++; }
    }
    return mpz_cmp_ui(cof, 1) == 0;
}

/* ---- GF(2) matrix ---- */
typedef unsigned long long BW;
#define BWS 64
typedef struct { BW *d; int nr,nc,mw,hw,tw; } mat_t;

static mat_t *mat_new(int nr, int nc) {
    mat_t *m = malloc(sizeof(mat_t));
    m->nr=nr; m->nc=nc; m->mw=(nc+BWS-1)/BWS; m->hw=(nr+BWS-1)/BWS; m->tw=m->mw+m->hw;
    m->d = calloc((size_t)nr*m->tw, sizeof(BW));
    for (int i=0;i<nr;i++) m->d[(size_t)i*m->tw+m->mw+i/BWS]|=(1ULL<<(i%BWS));
    return m;
}
static inline int mg(mat_t *m,int r,int c){return(m->d[(size_t)r*m->tw+c/BWS]>>(c%BWS))&1;}
static inline void mf(mat_t *m,int r,int c){m->d[(size_t)r*m->tw+c/BWS]^=(1ULL<<(c%BWS));}
static void mx(mat_t *m,int d,int s){BW *dd=m->d+(size_t)d*m->tw,*ss=m->d+(size_t)s*m->tw;for(int w=0;w<m->tw;w++)dd[w]^=ss[w];}

static int mat_elim(mat_t *m, int **deps, int *nd) {
    int *piv=malloc(m->nc*sizeof(int)); for(int c=0;c<m->nc;c++)piv[c]=-1;
    for(int c=0;c<m->nc;c++){
        int pr=-1;
        for(int r=0;r<m->nr;r++){
            if(!mg(m,r,c))continue;
            int u=0;for(int c2=0;c2<c;c2++)if(piv[c2]==r){u=1;break;}
            if(!u){pr=r;break;}
        }
        if(pr<0)continue; piv[c]=pr;
        for(int r=0;r<m->nr;r++) if(r!=pr&&mg(m,r,c))mx(m,r,pr);
    }
    *nd=0;*deps=malloc(m->nr*sizeof(int));
    for(int r=0;r<m->nr;r++){
        int z=1;for(int w=0;w<m->mw;w++)if(m->d[(size_t)r*m->tw+w]){z=0;break;}
        if(z)(*deps)[(*nd)++]=r;
    }
    free(piv); return *nd;
}

/* ---- Relation ---- */
typedef struct { mpz_t x; int *exp; int nlp; mpz_t lp[3]; } rel_t;

/* ---- Factor extraction ---- */
static int try_dep(mat_t *m, int dep, rel_t *rels, int nr) {
    mpz_t x,y,tmp,lpp;
    mpz_init(x);mpz_init(y);mpz_init(tmp);mpz_init(lpp);
    mpz_set_ui(x,1);mpz_set_ui(lpp,1);
    int *te=calloc(fb_size,sizeof(int));
    for(int i=0;i<nr;i++){
        if(!(m->d[(size_t)dep*m->tw+m->mw+i/BWS]&(1ULL<<(i%BWS))))continue;
        mpz_mul(x,x,rels[i].x);mpz_mod(x,x,N_g);
        for(int j=0;j<fb_size;j++)te[j]+=rels[i].exp[j];
        for(int l=0;l<rels[i].nlp;l++)mpz_mul(lpp,lpp,rels[i].lp[l]);
    }
    for(int j=0;j<fb_size;j++)if(te[j]&1){free(te);mpz_clear(x);mpz_clear(y);mpz_clear(tmp);mpz_clear(lpp);return 0;}
    mpz_set_ui(y,1);
    for(int j=0;j<fb_size;j++){if(te[j]>0){mpz_ui_pow_ui(tmp,fb[j],te[j]/2);mpz_mul(y,y,tmp);mpz_mod(y,y,N_g);}}
    if(mpz_cmp_ui(lpp,1)>0){
        if(!mpz_perfect_square_p(lpp)){free(te);mpz_clear(x);mpz_clear(y);mpz_clear(tmp);mpz_clear(lpp);return 0;}
        mpz_sqrt(tmp,lpp);mpz_mul(y,y,tmp);mpz_mod(y,y,N_g);
    }
    int found=0;
    mpz_sub(tmp,x,y);mpz_gcd(tmp,tmp,N_g);
    if(mpz_cmp_ui(tmp,1)>0&&mpz_cmp(tmp,N_g)<0){gmp_printf("FACTOR: %Zd\n",tmp);found=1;}
    if(!found){mpz_add(tmp,x,y);mpz_gcd(tmp,tmp,N_g);
    if(mpz_cmp_ui(tmp,1)>0&&mpz_cmp(tmp,N_g)<0){gmp_printf("FACTOR: %Zd\n",tmp);found=1;}}
    free(te);mpz_clear(x);mpz_clear(y);mpz_clear(tmp);mpz_clear(lpp);
    return found;
}

/* ---- Main ---- */
/* Good multipliers for Knuth-Schroeppel */
static int mults[] = {1,2,3,5,6,7,10,11,13,14,15,17,19,21,22,23,26,29,30,31,33,34,35,37,38,
    39,41,42,43,46,47,51,53,55,57,58,59,61,62,65,66,67,69,70,71,73,74,77,78,79,82,83,85,
    86,87,89,91,93,94,95,97,101,0};

int main(int argc, char **argv) {
    if(argc<2){fprintf(stderr,"Usage: %s <N>\n",argv[0]);return 1;}
    mpz_init(N_g); mpz_set_str(N_g,argv[1],10);
    int digs=(int)mpz_sizeinbase(N_g,10);
    int nbits=(int)mpz_sizeinbase(N_g,2);
    fprintf(stderr,"CCF: %d digits (%d bits)\n",digs,nbits);
    struct timespec ts,te;
    clock_gettime(CLOCK_MONOTONIC,&ts);

    /* Small factor check */
    gen_primes(1000000);
    for(int i=0;i<n_all_primes;i++){
        if(mpz_divisible_ui_p(N_g,all_primes[i])){printf("FACTOR: %u\n",all_primes[i]);return 0;}
    }

    /* Parameters - tuned for 30-100 digit range */
    double lnN=nbits*log(2.0), lnlnN=log(lnN);
    double Lval=exp(sqrt(lnN*lnlnN));

    /* B scaled with L */
    int B;
    if(digs<=35) B=(int)(pow(Lval,0.55));
    else if(digs<=50) B=(int)(pow(Lval,0.50));
    else B=(int)(pow(Lval,0.45));
    if(B<500)B=500;
    if(B>800000)B=800000;

    /* Generate factor base */
    if((unsigned int)B > all_primes[n_all_primes-1]) gen_primes(B+1000);
    fb=malloc(n_all_primes*sizeof(unsigned int));
    fb_size=0;
    for(int i=0;i<n_all_primes;i++){
        if(all_primes[i]>(unsigned int)B)break;
        if(all_primes[i]==2||legendre_ui(N_g,all_primes[i])==1)
            fb[fb_size++]=all_primes[i];
    }
    int target=fb_size+20;

    /* Sieve range per multiplier */
    int M;
    if(digs<=35) M=B*30;
    else if(digs<=50) M=B*20;
    else M=B*10;
    if(M>10000000) M=10000000;

    /* Count active multipliers */
    int nmults=0;
    for(int i=0;mults[i];i++)nmults++;
    if(nmults>40)nmults=40; /* limit for smaller numbers */

    long long max_cands=(long long)nmults*2*M;
    fprintf(stderr,"CCF: B=%d fb=%d target=%d M=%d nmults=%d max_cands=%lld\n",
            B,fb_size,target,M,nmults,max_cands);

    /* Build primorial once */
    mpz_t primorial;
    mpz_init(primorial);
    mpz_set_ui(primorial,1);
    {
        /* Estimate max value size: 2*M*sqrt(maxmult*N) */
        mpz_t maxval;
        mpz_init(maxval);
        mpz_mul_ui(maxval,N_g,mults[nmults-1]);
        mpz_sqrt(maxval,maxval);
        mpz_mul_ui(maxval,maxval,2*M);
        double log_max=mpz_sizeinbase(maxval,2)*log(2.0);
        mpz_clear(maxval);

        mpz_t pp;
        mpz_init(pp);
        for(int i=0;i<fb_size;i++){
            int e=(int)(log_max/log((double)fb[i]));
            if(e<1)e=1;
            mpz_ui_pow_ui(pp,fb[i],e);
            mpz_mul(primorial,primorial,pp);
        }
        mpz_clear(pp);
    }
    fprintf(stderr,"CCF: primorial %zu bits\n",mpz_sizeinbase(primorial,2));

    /* Allocate relations */
    int rel_cap=target*4;
    rel_t *rels=calloc(rel_cap,sizeof(rel_t));
    for(int i=0;i<rel_cap;i++){mpz_init(rels[i].x);rels[i].exp=calloc(fb_size,sizeof(int));for(int j=0;j<3;j++)mpz_init(rels[i].lp[j]);}
    int n_rels=0;

    /* Partial storage for LP matching */
    int part_cap=2000000;
    mpz_t *px=malloc(part_cap*sizeof(mpz_t));
    int **pe=malloc(part_cap*sizeof(int*));
    mpz_t *pc=malloc(part_cap*sizeof(mpz_t));
    for(int i=0;i<part_cap;i++){mpz_init(px[i]);pe[i]=calloc(fb_size,sizeof(int));mpz_init(pc[i]);}
    int n_parts=0;

    unsigned long lp_bnd=(unsigned long)B*(unsigned long)B;

    /* ---- Generate and batch-test candidates per multiplier ---- */
    int batch_sz=500000;
    mpz_t *bvals=malloc(batch_sz*sizeof(mpz_t));
    mpz_t *bx=malloc(batch_sz*sizeof(mpz_t));
    mpz_t *bsmooth=malloc(batch_sz*sizeof(mpz_t));
    for(int i=0;i<batch_sz;i++){mpz_init(bvals[i]);mpz_init(bx[i]);mpz_init(bsmooth[i]);}

    mpz_t kN,sqkN,xv,vv,cof;
    mpz_init(kN);mpz_init(sqkN);mpz_init(xv);mpz_init(vv);mpz_init(cof);

    int total_smooth=0,total_partial=0,total_tested=0;

    for(int mi=0;mi<nmults && n_rels<target;mi++){
        int k=mults[mi];
        mpz_mul_ui(kN,N_g,k);
        mpz_sqrt(sqkN,kN);

        int bn=0;
        for(int delta=1;delta<=M&&n_rels<target;delta++){
            for(int sd=0;sd<=1;sd++){
                if(sd==0) mpz_add_ui(xv,sqkN,delta);
                else { mpz_sub_ui(xv,sqkN,delta); if(mpz_sgn(xv)<=0)continue; }
                mpz_mul(vv,xv,xv);
                mpz_sub(vv,vv,kN);
                if(mpz_sgn(vv)<0) mpz_neg(vv,vv);
                if(mpz_sgn(vv)==0) continue;
                /* Remove multiplier */
                if(k>1){ while(mpz_divisible_ui_p(vv,k)) mpz_divexact_ui(vv,vv,k); }

                mpz_set(bvals[bn],vv);
                mpz_set(bx[bn],xv);
                bn++;

                if(bn>=batch_sz){
                    /* Process batch */
                    batch_smooth(bvals,bn,primorial,bsmooth);
                    for(int i=0;i<bn;i++){
                        mpz_divexact(cof,bvals[i],bsmooth[i]);
                        total_tested++;
                        if(mpz_cmp_ui(cof,1)==0){
                            /* Fully smooth */
                            if(n_rels<rel_cap){
                                mpz_t tc; mpz_init(tc);
                                trial_div(bvals[i],rels[n_rels].exp,tc);
                                mpz_set(rels[n_rels].x,bx[i]);
                                rels[n_rels].nlp=0;
                                n_rels++; total_smooth++;
                                mpz_clear(tc);
                            }
                        } else if(n_parts<part_cap) {
                            /* Check if cofactor is small enough for LP */
                            if(mpz_fits_ulong_p(cof)&&mpz_get_ui(cof)<=lp_bnd&&mpz_cmp_ui(cof,1)>0){
                                mpz_t tc; mpz_init(tc);
                                trial_div(bvals[i],pe[n_parts],tc);
                                mpz_set(px[n_parts],bx[i]);
                                mpz_set(pc[n_parts],tc);
                                if(mpz_cmp_ui(tc,1)>0) { n_parts++; total_partial++; }
                                mpz_clear(tc);
                            }
                        }
                    }
                    bn=0;

                    /* Check time */
                    clock_gettime(CLOCK_MONOTONIC,&te);
                    double elapsed=(te.tv_sec-ts.tv_sec)+(te.tv_nsec-ts.tv_nsec)*1e-9;
                    if(elapsed>250.0){
                        fprintf(stderr,"CCF: time limit approaching (%.0fs), stopping generation\n",elapsed);
                        goto done_gen;
                    }
                }
            }
        }
        /* Process remaining batch */
        if(bn>0){
            batch_smooth(bvals,bn,primorial,bsmooth);
            for(int i=0;i<bn;i++){
                mpz_divexact(cof,bvals[i],bsmooth[i]);
                total_tested++;
                if(mpz_cmp_ui(cof,1)==0){
                    if(n_rels<rel_cap){
                        mpz_t tc;mpz_init(tc);
                        trial_div(bvals[i],rels[n_rels].exp,tc);
                        mpz_set(rels[n_rels].x,bx[i]);
                        rels[n_rels].nlp=0;
                        n_rels++;total_smooth++;
                        mpz_clear(tc);
                    }
                } else if(n_parts<part_cap){
                    if(mpz_fits_ulong_p(cof)&&mpz_get_ui(cof)<=lp_bnd&&mpz_cmp_ui(cof,1)>0){
                        mpz_t tc;mpz_init(tc);
                        trial_div(bvals[i],pe[n_parts],tc);
                        mpz_set(px[n_parts],bx[i]);
                        mpz_set(pc[n_parts],tc);
                        if(mpz_cmp_ui(tc,1)>0){n_parts++;total_partial++;}
                        mpz_clear(tc);
                    }
                }
            }
        }

        fprintf(stderr,"\r  mult=%d tested=%d smooth=%d partial=%d rels=%d/%d   ",
                k,total_tested,total_smooth,total_partial,n_rels,target);
    }
done_gen:
    fprintf(stderr,"\nCCF: generation done: %d smooth, %d partial, %d rels\n",
            total_smooth,total_partial,n_rels);

    /* ---- LP matching ---- */
    if(n_rels<target && n_parts>=2){
        fprintf(stderr,"CCF: LP matching on %d partials\n",n_parts);
        /* Sort by cofactor */
        int *si=malloc(n_parts*sizeof(int));
        for(int i=0;i<n_parts;i++)si[i]=i;
        /* Simple shell sort */
        for(int gap=n_parts/2;gap>0;gap/=2)
            for(int i=gap;i<n_parts;i++){
                int t=si[i],j=i;
                while(j>=gap&&mpz_cmp(pc[si[j-gap]],pc[t])>0){si[j]=si[j-gap];j-=gap;}
                si[j]=t;
            }
        int matched=0;
        for(int i=0;i+1<n_parts&&n_rels<rel_cap;i++){
            int a=si[i],b=si[i+1];
            if(mpz_cmp(pc[a],pc[b])==0&&mpz_cmp_ui(pc[a],1)>0){
                rel_t *r=&rels[n_rels];
                mpz_mul(r->x,px[a],px[b]); mpz_mod(r->x,r->x,N_g);
                for(int j=0;j<fb_size;j++) r->exp[j]=pe[a][j]+pe[b][j];
                r->nlp=2; mpz_set(r->lp[0],pc[a]); mpz_set(r->lp[1],pc[b]);
                n_rels++; matched++; i++;
            }
        }
        fprintf(stderr,"CCF: LP matched %d, total rels=%d/%d\n",matched,n_rels,target);
        free(si);
    }

    /* ---- Phase: Batch GCD on cofactors ---- */
    if(n_rels<target && n_parts>=100){
        fprintf(stderr,"CCF: batch GCD cofactor collision on %d partials\n",n_parts);
        mpz_t *cofs=malloc(n_parts*sizeof(mpz_t));
        mpz_t *shared=malloc(n_parts*sizeof(mpz_t));
        for(int i=0;i<n_parts;i++){mpz_init_set(cofs[i],pc[i]);mpz_init(shared[i]);}

        /* Build product tree of cofactors */
        ptree *ct=pt_build(cofs,n_parts);
        mpz_t *csq=malloc(n_parts*sizeof(mpz_t));
        for(int i=0;i<n_parts;i++){mpz_init(csq[i]);mpz_mul(csq[i],cofs[i],cofs[i]);}
        ptree *sqt=pt_build(csq,n_parts);

        /* Remainder tree: prod mod c_i^2 */
        mpz_t *rems=malloc(n_parts*sizeof(mpz_t));
        for(int i=0;i<n_parts;i++)mpz_init(rems[i]);
        pt_rem(sqt,ct->data[ct->lev_off[ct->nlevs-1]],rems);

        mpz_t tmp;mpz_init(tmp);
        int colls=0;
        for(int i=0;i<n_parts;i++){
            if(mpz_sgn(cofs[i])==0){mpz_set_ui(shared[i],0);continue;}
            mpz_fdiv_q(tmp,rems[i],cofs[i]);
            mpz_gcd(shared[i],tmp,cofs[i]);
            if(mpz_cmp_ui(shared[i],1)>0)colls++;
            /* Check if shared factor divides N */
            mpz_gcd(tmp,shared[i],N_g);
            if(mpz_cmp_ui(tmp,1)>0&&mpz_cmp(tmp,N_g)<0){
                gmp_printf("FACTOR: %Zd\n",tmp);
                return 0;
            }
        }
        fprintf(stderr,"CCF: %d cofactors with shared factors\n",colls);
        mpz_clear(tmp);

        for(int i=0;i<n_parts;i++){mpz_clear(cofs[i]);mpz_clear(shared[i]);mpz_clear(csq[i]);mpz_clear(rems[i]);}
        free(cofs);free(shared);free(csq);free(rems);
        pt_free(ct);pt_free(sqt);
    }

    /* ---- Linear algebra ---- */
    int found=0;
    if(n_rels>fb_size){
        fprintf(stderr,"CCF: GF(2) elimination %d x %d\n",n_rels,fb_size);
        mat_t *m=mat_new(n_rels,fb_size);
        for(int i=0;i<n_rels;i++)
            for(int j=0;j<fb_size;j++)
                if(rels[i].exp[j]&1)mf(m,i,j);
        int *deps;int nd;
        mat_elim(m,&deps,&nd);
        fprintf(stderr,"CCF: %d deps\n",nd);
        for(int d=0;d<nd;d++){if(try_dep(m,deps[d],rels,n_rels)){found=1;break;}}
        free(deps);free(m->d);free(m);
    } else {
        fprintf(stderr,"CCF: not enough rels %d/%d\n",n_rels,target);
    }

    if(!found) fprintf(stderr,"CCF: FAILED\n");

    clock_gettime(CLOCK_MONOTONIC,&te);
    fprintf(stderr,"CCF: total %.1fs\n",(te.tv_sec-ts.tv_sec)+(te.tv_nsec-ts.tv_nsec)*1e-9);

    /* Cleanup (abbreviated for brevity) */
    mpz_clear(N_g);mpz_clear(primorial);mpz_clear(kN);mpz_clear(sqkN);mpz_clear(xv);mpz_clear(vv);mpz_clear(cof);
    return found?0:1;
}
