/* Minimal SSS test - just smooth relations, no LP, with verification */
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <gmp.h>

static long pmod(long b, long e, long m) {
    long r=1; b=((b%m)+m)%m;
    while(e>0){if(e&1)r=(__int128)r*b%m;b=(__int128)b*b%m;e>>=1;}
    return r;
}
static int leg(long a, long p) { long r=pmod(a,(p-1)/2,p); return r==p-1?-1:(int)r; }
static long ts(long n, long p) {
    if(p==2) return n&1; n=((n%p)+p)%p; if(!n)return 0;
    long Q=p-1,S=0; while(Q%2==0){Q/=2;S++;}
    if(S==1) return pmod(n,(p+1)/4,p);
    long z=2; while(pmod(z,(p-1)/2,p)!=p-1)z++;
    long M=S,c=pmod(z,Q,p),t=pmod(n,Q,p),R=pmod(n,(Q+1)/2,p);
    while(t!=1){long i=0,v=t;while(v!=1){v=(__int128)v*v%p;i++;}
    long b=c;for(long j=0;j<M-i-1;j++)b=(__int128)b*b%p;
    M=i;c=(__int128)b*b%p;t=(__int128)t*c%p;R=(__int128)R*b%p;}
    return R;
}

struct BitRow {
    std::vector<uint64_t> bits,hist; int wb,wh;
    void init(int nb,int nh){wb=(nb+63)/64;wh=(nh+63)/64;bits.assign(wb,0);hist.assign(wh,0);}
    void sbit(int i){bits[i/64]|=1ULL<<(i%64);}
    int gbit(int i)const{return(bits[i/64]>>(i%64))&1;}
    void shist(int i){hist[i/64]|=1ULL<<(i%64);}
    int ghist(int i)const{return(hist[i/64]>>(i%64))&1;}
    void xr(const BitRow&o){for(int i=0;i<wb;i++)bits[i]^=o.bits[i];for(int i=0;i<wh;i++)hist[i]^=o.hist[i];}
    bool zb()const{for(int i=0;i<wb;i++)if(bits[i])return false;return true;}
};

int main(int argc, char*argv[]) {
    mpz_t N,tmp,g,bv; mpz_init(N);mpz_init(tmp);mpz_init(g);mpz_init(bv);
    mpz_set_str(N, argv[1], 10);
    struct timespec ts0; clock_gettime(CLOCK_MONOTONIC,&ts0);
    size_t dig=mpz_sizeinbase(N,10);
    int m=400; if(dig>34)m=600; if(dig>38)m=1000; if(dig>42)m=1400; if(dig>48)m=2400;
    if(dig>52)m=4000; if(dig>56)m=8000; if(dig>62)m=12000; if(dig>70)m=20000;

    mpz_sqrt(bv,N); mpz_mul(tmp,bv,bv); if(mpz_cmp(tmp,N)<0)mpz_add_ui(bv,bv,1);

    std::vector<long> fb;
    {long lim=m*20; std::vector<bool> s(lim+1,true);s[0]=s[1]=false;
    for(long i=2;i*i<=lim;i++)if(s[i])for(long j=i*i;j<=lim;j+=i)s[j]=false;
    for(long i=2;i<=lim&&(int)fb.size()<m;i++)if(s[i]){
        mpz_mod_ui(tmp,N,i);long nm=mpz_get_ui(tmp);
        if(i==2||leg(nm,i)==1)fb.push_back(i);
    }}
    int fbs=fb.size();
    int pls=fbs/5; if(pls<10)pls=10;
    int target=fbs+1+60; /* +1 for sign column, +60 extra for more dependencies */

    std::vector<long> r0(fbs),r1(fbs);
    for(int i=0;i<fbs;i++){
        long p=fb[i]; mpz_mod_ui(tmp,N,p); long r=ts(mpz_get_ui(tmp),p);
        mpz_mod_ui(tmp,bv,p); long bm=mpz_get_ui(tmp);
        r0[i]=((r-bm)%p+p)%p; r1[i]=((p-r-bm)%p+p)%p;
        if(r0[i]>r1[i])std::swap(r0[i],r1[i]);
    }

    /* CRT for plist primes */
    mpz_t nval; mpz_init(nval); mpz_set_ui(nval,1);
    for(int i=0;i<pls;i++) mpz_mul_ui(nval,nval,fb[i]);
    std::vector<mpz_t> coeff(pls),diff(pls);
    for(int i=0;i<pls;i++){
        mpz_init(coeff[i]);mpz_init(diff[i]);
        long p=fb[i]; mpz_t np;mpz_init(np);
        mpz_divexact_ui(np,nval,p);
        mpz_mod_ui(tmp,np,p);
        long inv=pmod(mpz_get_ui(tmp),p-2,p);
        mpz_mul_ui(coeff[i],np,inv);
        mpz_mul_ui(diff[i],coeff[i],r1[i]-r0[i]);
        mpz_clear(np);
    }

    gmp_randstate_t rs; gmp_randinit_default(rs); gmp_randseed_ui(rs,42);

    /* Collect smooth relations with VERIFIED exponents */
    struct Rel { mpz_t xpb; std::vector<int> exps; };
    std::vector<Rel> rels;
    int srch=0;

    while((int)rels.size()<target){
        struct timespec tn; clock_gettime(CLOCK_MONOTONIC,&tn);
        double el=(tn.tv_sec-ts0.tv_sec)+(tn.tv_nsec-ts0.tv_nsec)/1e9;
        if(el>280)break;

        std::vector<int> ch;
        {std::vector<int> pool(pls);for(int i=0;i<pls;i++)pool[i]=i;
        for(int i=0;i<6&&!pool.empty();i++){
            int idx=gmp_urandomm_ui(rs,pool.size());
            ch.push_back(pool[idx]);pool.erase(pool.begin()+idx);
        }}
        if((int)ch.size()<6)continue;

        mpz_t M,xv; mpz_init(M);mpz_init(xv);
        mpz_set_ui(M,1); for(int i:ch)mpz_mul_ui(M,M,fb[i]);
        mpz_set_ui(xv,0);
        for(int i:ch){mpz_mul_ui(tmp,coeff[i],r0[i]);mpz_add(xv,xv,tmp);}
        mpz_mod(xv,xv,M);

        for(int si=0;si<(int)ch.size();si++){
            mpz_add(xv,xv,diff[ch[si]]);mpz_mod(xv,xv,M);

            for(int ji=-1;ji<(int)ch.size();ji++){
                if(ji>=0&&ch[ji]==ch[si])continue;
                mpz_t ms; mpz_init(ms);
                if(ji<0)mpz_set(ms,M); else mpz_divexact_ui(ms,M,fb[ch[ji]]);

                /* Use sorted array instead of hash map for collision counting */
                std::vector<long> kvals;
                kvals.reserve((fbs-pls)*4);
                for(int v=pls;v<fbs;v++){
                    long p=fb[v]; mpz_mod_ui(tmp,ms,p); long mp=mpz_get_ui(tmp);
                    if(!mp)continue; long inv=pmod(mp,p-2,p);
                    mpz_mod_ui(tmp,xv,p); long xm=mpz_get_ui(tmp);
                    long k0=(__int128)((r0[v]-xm+p)%p)*inv%p;
                    long k1=(__int128)((r1[v]-xm+p)%p)*inv%p;
                    kvals.push_back(k0);kvals.push_back(k0-p);
                    kvals.push_back(k1);kvals.push_back(k1-p);
                }
                std::sort(kvals.begin(),kvals.end());

                for(size_t ki=0;ki<kvals.size();){
                    long k=kvals[ki];
                    int cnt=0;
                    while(ki<kvals.size()&&kvals[ki]==k){cnt++;ki++;}
                    if(cnt<=2)continue;
                    mpz_t x2,pv; mpz_init(x2);mpz_init(pv);
                    mpz_set_si(tmp,k); mpz_mul(x2,tmp,ms); mpz_add(x2,x2,xv);
                    mpz_add(pv,x2,bv); /* x+b */
                    mpz_t xpb; mpz_init_set(xpb,pv); /* save x+b */
                    mpz_mul(pv,pv,pv); mpz_sub(pv,pv,N); /* (x+b)^2 - N */

                    if(mpz_sgn(pv)==0){mpz_clear(x2);mpz_clear(pv);mpz_clear(xpb);continue;}

                    std::vector<int> exps(fbs+1,0); /* +1 for sign */
                    if(mpz_sgn(pv)<0){exps[fbs]=1;mpz_neg(pv,pv);}
                    for(int j=0;j<fbs;j++){
                        long p=fb[j];
                        while(mpz_divisible_ui_p(pv,p)){mpz_divexact_ui(pv,pv,p);exps[j]++;}
                    }

                    if(mpz_cmp_ui(pv,1)==0){
                        /* VERIFY: (x+b)^2 mod N should equal (-1)^sign * prod(p^e) mod N */
                        mpz_t lhs,rhs; mpz_init(lhs);mpz_init(rhs);
                        mpz_mul(lhs,xpb,xpb);mpz_mod(lhs,lhs,N);
                        mpz_set_ui(rhs,1);
                        for(int j=0;j<fbs;j++){
                            if(exps[j]==0)continue;
                            mpz_set_ui(tmp,fb[j]);
                            mpz_powm_ui(tmp,tmp,exps[j],N);
                            mpz_mul(rhs,rhs,tmp);mpz_mod(rhs,rhs,N);
                        }
                        if(exps[fbs]){mpz_sub(rhs,N,rhs);} /* multiply by -1 */
                        if(mpz_cmp(lhs,rhs)==0){
                            Rel rel; mpz_init_set(rel.xpb,xpb);
                            rel.exps=exps;
                            rels.push_back(rel);
                        } else {
                            fprintf(stderr,"VERIFY FAILED!\n");
                        }
                        mpz_clear(lhs);mpz_clear(rhs);
                    }
                    mpz_clear(x2);mpz_clear(pv);mpz_clear(xpb);
                }
                mpz_clear(ms);
            }
        }
        mpz_clear(M);mpz_clear(xv);
        srch++;
        if(srch%100==0)fprintf(stderr,"  srch %d: %zu/%d rels\n",srch,rels.size(),target);
    }

    int nrels=rels.size(),ncols=fbs+1;
    fprintf(stderr,"%d verified rels, %d cols\n",nrels,ncols);
    if(nrels<=ncols){printf("FAIL\n");return 1;}
    if(nrels>ncols+200)nrels=ncols+200;

    std::vector<BitRow> mat(nrels);
    for(int i=0;i<nrels;i++){
        mat[i].init(ncols,nrels);mat[i].shist(i);
        for(int j=0;j<ncols;j++)if(rels[i].exps[j]&1)mat[i].sbit(j);
    }
    int rank=0;
    for(int col=0;col<ncols&&rank<nrels;col++){
        int piv=-1;for(int r=rank;r<nrels;r++)if(mat[r].gbit(col)){piv=r;break;}
        if(piv<0)continue;
        if(piv!=rank)std::swap(mat[piv],mat[rank]);
        for(int r=0;r<nrels;r++)if(r!=rank&&mat[r].gbit(col))mat[r].xr(mat[rank]);
        rank++;
    }
    fprintf(stderr,"Rank %d, null %d\n",rank,nrels-rank);

    mpz_t px,py;mpz_init(px);mpz_init(py);
    int found=0;
    for(int row=0;row<nrels&&!found;row++){
        if(!mat[row].zb())continue;
        mpz_set_ui(px,1);
        std::vector<long> se(ncols,0);
        for(int i=0;i<nrels;i++){
            if(!mat[row].ghist(i))continue;
            mpz_mul(px,px,rels[i].xpb);mpz_mod(px,px,N);
            for(int j=0;j<ncols;j++)se[j]+=rels[i].exps[j];
        }
        bool ok=true;for(int j=0;j<ncols;j++)if(se[j]&1){ok=false;break;}
        if(!ok)continue;
        mpz_set_ui(py,1);
        for(int j=0;j<fbs;j++){
            if(!se[j])continue;
            mpz_set_ui(tmp,fb[j]);mpz_powm_ui(tmp,tmp,se[j]/2,N);
            mpz_mul(py,py,tmp);mpz_mod(py,py,N);
        }
        for(int s=0;s<2&&!found;s++){
            if(s==0)mpz_sub(tmp,px,py);else mpz_add(tmp,px,py);
            mpz_gcd(g,tmp,N);
            if(mpz_cmp_ui(g,1)>0&&mpz_cmp(g,N)<0){
                struct timespec tn;clock_gettime(CLOCK_MONOTONIC,&tn);
                double el=(tn.tv_sec-ts0.tv_sec)+(tn.tv_nsec-ts0.tv_nsec)/1e9;
                gmp_printf("FACTOR: %Zd\n",g);
                mpz_t c;mpz_init(c);mpz_divexact(c,N,g);
                gmp_printf("COFACTOR: %Zd\n",c);
                fprintf(stderr,"SSS factored in %.3fs\n",el);
                mpz_clear(c);found=1;
            }
        }
    }
    if(!found){fprintf(stderr,"No factor.\n");printf("FAIL\n");}
    return found?0:1;
}
