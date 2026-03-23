/*
 * Self-Initializing Quadratic Sieve (SIQS).
 *
 * Polynomials g(x) = (a*x + b)^2 - N = a * Q(x), where Q(x) = a*x^2 + 2*b*x + c,
 * c = (b^2 - N)/a, and a is a product of s primes from the factor base, chosen so
 * that a ≈ sqrt(2N)/M. This keeps max|Q(x)| ≈ sqrt(N/2) for x in [-M, M].
 *
 * Single large prime variation.
 * Usage: ./mpqs <number>
 * Single-threaded, seed=42.
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <gmp.h>

static long power_mod(long base, long exp, long mod) {
    long r = 1; base = ((base%mod)+mod)%mod;
    while(exp>0){if(exp&1)r=(__int128)r*base%mod;base=(__int128)base*base%mod;exp>>=1;}
    return r;
}
static int legendre_sym(long a, long p) {
    long r = power_mod(a,(p-1)/2,p); return r==p-1?-1:(int)r;
}
static long tonelli(long n, long p) {
    if(p==2)return n&1; n=((n%p)+p)%p; if(!n)return 0;
    long Q=p-1,S=0; while(Q%2==0){Q/=2;S++;}
    if(S==1)return power_mod(n,(p+1)/4,p);
    long z=2; while(power_mod(z,(p-1)/2,p)!=p-1)z++;
    long M2=S,c=power_mod(z,Q,p),t=power_mod(n,Q,p),R=power_mod(n,(Q+1)/2,p);
    while(t!=1){long i=0,v=t;while(v!=1){v=(__int128)v*v%p;i++;}
    long b=c;for(long j=0;j<M2-i-1;j++)b=(__int128)b*b%p;
    M2=i;c=(__int128)b*b%p;t=(__int128)t*c%p;R=(__int128)R*b%p;}
    return R;
}
static long mod_inv(long a, long m) { return power_mod(a,m-2,m); }

static std::vector<long> sieve_primes(long lim) {
    std::vector<bool> s(lim+1,true);s[0]=s[1]=false;
    for(long i=2;i*i<=lim;i++)if(s[i])for(long j=i*i;j<=lim;j+=i)s[j]=false;
    std::vector<long> p;for(long i=2;i<=lim;i++)if(s[i])p.push_back(i);return p;
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

struct Relation {
    mpz_t val;
    std::vector<int> exps;
    long large_prime;
};

int main(int argc,char*argv[]){
    if(argc<2){fprintf(stderr,"Usage: %s <number>\n",argv[0]);return 1;}
    mpz_t N,tmp,g; mpz_init(N);mpz_init(tmp);mpz_init(g);
    if(mpz_set_str(N,argv[1],10)!=0){fprintf(stderr,"Bad number\n");return 1;}

    struct timespec tstart;
    clock_gettime(CLOCK_MONOTONIC,&tstart);

    size_t digits=mpz_sizeinbase(N,10);
    double ln_n=log(2.0)*mpz_sizeinbase(N,2);
    double ln_ln_n=log(ln_n);
    double opt=exp(0.5*sqrt(ln_n*ln_ln_n));

    long B=(long)(opt*0.9); if(B<2000)B=2000; if(B>8000000)B=8000000;
    long M=65536; if(digits>50)M=131072; if(digits>70)M=262144;
    long lp_mult=80;

    std::vector<long> all_primes=sieve_primes(B);
    std::vector<long> fb; fb.push_back(-1);
    {mpz_t nm;mpz_init(nm);
    for(long p:all_primes){mpz_mod_ui(nm,N,p);long n2=mpz_get_ui(nm);
    if(p==2||legendre_sym(n2,p)==1)fb.push_back(p);}
    mpz_clear(nm);}

    int fb_size=fb.size();
    int target=fb_size+40;
    long lp_bound=B*lp_mult;

    std::vector<long> sqN(fb_size);
    for(int i=1;i<fb_size;i++){long p=fb[i];mpz_mod_ui(tmp,N,p);sqN[i]=tonelli(mpz_get_ui(tmp),p);}

    std::vector<unsigned char> fb_log(fb_size);
    for(int i=1;i<fb_size;i++)fb_log[i]=(unsigned char)(log2((double)fb[i])+0.5);

    /* Compute target_a = sqrt(2N)/M */
    mpz_t target_a;mpz_init(target_a);
    {mpz_t s2n;mpz_init(s2n);mpz_mul_ui(s2n,N,2);mpz_sqrt(s2n,s2n);
    mpz_tdiv_q_ui(target_a,s2n,M);mpz_clear(s2n);}

    fprintf(stderr,"SIQS: %zu dig, B=%ld, fb=%d, M=%ld, target=%d, lp=%ld, target_a=",
            digits,B,fb_size,M,target,lp_bound);
    mpz_out_str(stderr,10,target_a);fprintf(stderr,"\n");

    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);gmp_randseed_ui(rstate,42);

    std::vector<Relation> rels;
    std::unordered_map<long,int> lp_map;
    std::vector<unsigned char> sieve(2*M+1);

    /* Determine how many primes s to multiply for 'a'.
     * We want product ≈ target_a.
     * Use primes from the upper portion of the factor base. */
    double log_target_a=mpz_sizeinbase(target_a,2)*log(2.0);

    /* Pick primes from the top 30% of fb for 'a' factors */
    int a_prime_start=fb_size*2/3; if(a_prime_start<5)a_prime_start=5;
    double avg_log_prime=0;
    int a_prime_count=0;
    for(int i=a_prime_start;i<fb_size;i++){avg_log_prime+=log((double)fb[i]);a_prime_count++;}
    if(a_prime_count>0)avg_log_prime/=a_prime_count;
    int s_primes=(int)(log_target_a/avg_log_prime+0.5);
    if(s_primes<1)s_primes=1;
    if(s_primes>8)s_primes=8;

    fprintf(stderr,"Using %d primes per 'a' (from idx %d-%d, avg logp=%.1f)\n",
            s_primes,a_prime_start,fb_size-1,avg_log_prime);

    int poly_count=0,max_polys=500000;

    int usable_count=0;
    while(usable_count<target && poly_count<max_polys){
        struct timespec tnow;
        clock_gettime(CLOCK_MONOTONIC,&tnow);
        double elapsed=(tnow.tv_sec-tstart.tv_sec)+(tnow.tv_nsec-tstart.tv_nsec)/1e9;
        if(elapsed>280){fprintf(stderr,"Timeout %.1fs, %zu rels\n",elapsed,rels.size());break;}

        /* Generate 'a' = product of s_primes distinct fb primes */
        std::vector<int> a_indices; /* indices into fb[] */
        mpz_t a_val;mpz_init(a_val);mpz_set_ui(a_val,1);

        /* Randomly select s_primes from [a_prime_start, fb_size) */
        {
            std::vector<int> pool;
            for(int i=a_prime_start;i<fb_size;i++)pool.push_back(i);
            for(int i=0;i<s_primes && !pool.empty();i++){
                int idx=gmp_urandomm_ui(rstate,pool.size());
                a_indices.push_back(pool[idx]);
                mpz_mul_ui(a_val,a_val,fb[pool[idx]]);
                pool.erase(pool.begin()+idx);
            }
        }

        if((int)a_indices.size()<s_primes){mpz_clear(a_val);continue;}

        /* Compute b such that b^2 ≡ N mod a, using CRT.
         * For each prime q_i in 'a': b ≡ sqrt(N) mod q_i
         * Use CRT to combine. There are 2^s_primes choices of sign. */

        /* For simplicity, just use one choice: all positive sqroots */
        /* TODO: enumerate sign combinations for more polynomials from same 'a' */

        mpz_t b_val,c_val;mpz_init(b_val);mpz_init(c_val);

        /* CRT: for each prime q, we have b ≡ sqN[q_idx] mod q.
         * b = sum_i (t_i * (a/q_i) * ((a/q_i)^{-1} mod q_i)) where t_i = sqN[q_idx_i]
         * Then reduce b mod a. */
        mpz_set_ui(b_val,0);
        for(int k=0;k<s_primes;k++){
            int qi=a_indices[k];
            long q=fb[qi];
            /* a_over_q = a / q */
            mpz_t a_over_q;mpz_init(a_over_q);
            mpz_divexact_ui(a_over_q,a_val,q);
            /* inv = a_over_q^{-1} mod q */
            mpz_mod_ui(tmp,a_over_q,q);
            long aq_mod_q=mpz_get_ui(tmp);
            long inv=mod_inv(aq_mod_q,q);
            /* contribution = sqN[qi] * inv * a_over_q */
            mpz_mul_ui(tmp,a_over_q,sqN[qi]);
            mpz_mul_ui(tmp,tmp,inv);
            mpz_add(b_val,b_val,tmp);
            mpz_clear(a_over_q);
        }
        mpz_mod(b_val,b_val,a_val);

        /* Ensure b^2 ≡ N mod a */
        mpz_mul(tmp,b_val,b_val);
        mpz_sub(tmp,tmp,N);
        if(!mpz_divisible_p(tmp,a_val)){
            /* CRT failed, skip */
            mpz_clear(a_val);mpz_clear(b_val);mpz_clear(c_val);
            continue;
        }

        /* c = (b^2 - N)/a */
        mpz_mul(c_val,b_val,b_val);
        mpz_sub(c_val,c_val,N);
        mpz_divexact(c_val,c_val,a_val);

        /* Sieve: Q(x) = a*x^2 + 2*b*x + c
         * For fb prime p (not dividing a): x ≡ a^{-1}*(-b ± sqrt(N)) mod p */
        memset(sieve.data(),0,sieve.size());

        long a_long=0;
        if(mpz_fits_slong_p(a_val))a_long=mpz_get_si(a_val);

        for(int i=1;i<fb_size;i++){
            long p=fb[i];
            if(p>2*M)break;

            /* Check if p divides a */
            bool p_divides_a=false;
            for(int k=0;k<s_primes;k++)if(fb[a_indices[k]]==p){p_divides_a=true;break;}

            if(p==2){
                /* Q(x) mod 2: a*x^2 + c mod 2 (2b*x is even) */
                mpz_mod_ui(tmp,c_val,2);
                long cm=mpz_get_ui(tmp);
                if(a_long){
                    long am=a_long&1;
                    /* am*x^2 + cm ≡ 0 mod 2 */
                    /* if am=0: cm≡0 mod 2 (all x work) or cm≡1 (no x works) */
                    /* if am=1: x^2 ≡ -cm mod 2 → x ≡ cm mod 2 (since -1≡1 mod 2) */
                    if(am==0){
                        if(cm==0){for(long j=0;j<=2*M;j++)sieve[j]+=1;}
                    }else{
                        long start=cm;
                        for(long j=M+start;j<=2*M;j+=2)sieve[j]+=1;
                        for(long j=M+start-2;j>=0;j-=2)sieve[j]+=1;
                    }
                }
                continue;
            }

            unsigned char lp2=fb_log[i];

            if(p_divides_a){
                /* When p|a: Q(x) = a*x^2 + 2*b*x + c, and p|a, so Q(x) ≡ 2*b*x + c mod p
                 * x ≡ -c/(2b) mod p */
                mpz_mod_ui(tmp,b_val,p);
                long bm=mpz_get_ui(tmp);
                long twob=(2*bm)%p;
                if(twob==0)continue;
                mpz_mod_ui(tmp,c_val,p);
                long cm=mpz_get_ui(tmp);
                long r=(long)((__int128)(p-cm)*mod_inv(twob,p)%p);
                for(long j=M+r;j<=2*M;j+=p)sieve[j]+=lp2;
                for(long j=(long)M+r-p;j>=0;j-=p)sieve[j]+=lp2;
                continue;
            }

            /* General case: x ≡ a^{-1}(-b ± sqrt(N)) mod p */
            mpz_mod_ui(tmp,a_val,p);
            long am=mpz_get_ui(tmp);
            long ai=mod_inv(am,p);
            mpz_mod_ui(tmp,b_val,p);
            long bm=mpz_get_ui(tmp);
            long neg_b=(p-bm)%p;

            long r1=(__int128)ai*((neg_b+sqN[i])%p)%p;
            long r2=(__int128)ai*((neg_b+p-sqN[i])%p)%p;

            for(long j=M+r1;j<=2*M;j+=p)sieve[j]+=lp2;
            for(long j=(long)M+r1-p;j>=0;j-=p)sieve[j]+=lp2;
            if(r2!=r1){
                for(long j=M+r2;j<=2*M;j+=p)sieve[j]+=lp2;
                for(long j=(long)M+r2-p;j>=0;j-=p)sieve[j]+=lp2;
            }
        }

        /* Threshold: max|Q(x)| ≈ sqrt(N/2) for optimal a.
         * log2(sqrt(N/2)) = bits(N)/2 - 0.5 */
        float log2_maxQ=(float)mpz_sizeinbase(N,2)/2.0f;
        /* But actual Q(x) varies: near center Q(0)=c which could be large or small,
         * and Q(M) ≈ a*M^2 ≈ sqrt(2N)*M/M*M^2... no, a*M^2 with a≈sqrt(2N)/M gives
         * Q(M) ≈ sqrt(2N)*M ≈ sqrt(2N)*M. Hmm, that's not right either.
         *
         * Actually: a ≈ sqrt(2N)/M, so Q(M) ≈ a*M^2 ≈ sqrt(2N)*M. For M=65536:
         * Q(M) ≈ sqrt(2*10^30)*65536 ≈ 10^15 * 65536 ≈ 10^20.
         * Q(0) = c = (b^2-N)/a ≈ N/a ≈ N*M/sqrt(2N) = sqrt(N/2)*M ≈ 10^15*65536 ≈ 10^20.
         * Hmm, those are the same order. The SIQS trick is that |Q(x)| is bounded by
         * roughly M*sqrt(N/(2a)) + ... let me just use a generous threshold. */
        /* Use a rough estimate and be generous with slack */
        float half_bits=(float)mpz_sizeinbase(N,2)/2.0f;
        float log2_M=log2f((float)M);
        float est_max=half_bits+log2_M; /* rough upper bound */
        float slack=28.0f; /* bits of slack */

        for(long idx=0;idx<=2*M;idx++){
            if(sieve[idx]<est_max-slack-10)continue; /* quick reject */

            long xv=idx-M;

            /* More precise estimate for this x */
            float est;
            long ax=(xv<0)?-xv:xv;
            if(ax<2)est=half_bits; /* Q(0)≈c≈sqrt(N/2)*M */
            else{
                /* Q(x) ≈ a*x^2 + ... ≈ max(a*x^2, 2b*x, |c|) */
                /* Rough: log2(a) + 2*log2(|x|) vs half_bits+log2_M */
                float la=log2f((float)mpz_sizeinbase(a_val,2));
                float lx=log2f((float)ax);
                est=la+2*lx; /* a*x^2 term */
                if(est<half_bits)est=half_bits; /* c term dominates for small x */
            }

            if(sieve[idx]<est-slack)continue;

            /* Trial divide */
            mpz_t ax_b,qv;
            mpz_init(ax_b);mpz_init(qv);

            /* ax+b */
            if(a_long && xv>=-1000000000LL && xv<=1000000000LL){
                mpz_set_si(ax_b,xv);
                mpz_mul(ax_b,ax_b,a_val);
            }else{
                mpz_set_si(tmp,xv);
                mpz_mul(ax_b,tmp,a_val);
            }
            mpz_add(ax_b,ax_b,b_val);

            /* (ax+b)^2 - N */
            mpz_mul(qv,ax_b,ax_b);
            mpz_sub(qv,qv,N);

            /* This equals a * Q(x). We'll factor the whole thing and track exponents. */
            std::vector<int> exps(fb_size,0);
            int neg=(mpz_sgn(qv)<0);
            if(neg){exps[0]=1;mpz_neg(qv,qv);}

            /* Trial divide by all fb primes */
            for(int j=1;j<fb_size;j++){
                long p=fb[j];
                while(mpz_divisible_ui_p(qv,p)){
                    mpz_divexact_ui(qv,qv,p);
                    exps[j]++;
                }
            }

            bool smooth=(mpz_cmp_ui(qv,1)==0);
            long rem=0;
            bool lp_ok=false;

            if(!smooth && mpz_fits_ulong_p(qv)){
                rem=mpz_get_ui(qv);
                if(rem>1 && rem<=(unsigned long)lp_bound){
                    if(mpz_probab_prime_p(qv,3))lp_ok=true;
                }
            }

            if(smooth){
                Relation rel;
                mpz_init_set(rel.val,ax_b);
                rel.exps=exps;
                rel.large_prime=0;
                rels.push_back(rel);
                usable_count++;
            }else if(lp_ok){
                auto it=lp_map.find(rem);
                if(it!=lp_map.end()){
                    int oi=it->second;
                    Relation&orel=rels[oi];
                    Relation rel;
                    mpz_init(rel.val);
                    mpz_mul(rel.val,ax_b,orel.val);
                    mpz_mod(rel.val,rel.val,N);
                    rel.exps.resize(fb_size);
                    for(int j=0;j<fb_size;j++)
                        rel.exps[j]=exps[j]+orel.exps[j];
                    rel.large_prime=-rem;
                    rels.push_back(rel);
                    usable_count++;
                }else{
                    Relation rel;
                    mpz_init_set(rel.val,ax_b);
                    rel.exps=exps;
                    rel.large_prime=rem;
                    rels.push_back(rel);
                    lp_map[rem]=rels.size()-1;
                }
            }

            mpz_clear(ax_b);mpz_clear(qv);
        }

        mpz_clear(a_val);mpz_clear(b_val);mpz_clear(c_val);
        poly_count++;

        if(poly_count%200==0){
            clock_gettime(CLOCK_MONOTONIC,&tnow);
            elapsed=(tnow.tv_sec-tstart.tv_sec)+(tnow.tv_nsec-tstart.tv_nsec)/1e9;
            int full=0;
            for(auto&r:rels)if(r.large_prime==0||r.large_prime<0)full++;
            fprintf(stderr,"  poly %d: %d usable/%zu total (%.1fs)\n",poly_count,full,rels.size(),elapsed);
        }
    }

    /* Collect usable relations */
    std::vector<int> usable;
    for(int i=0;i<(int)rels.size();i++){
        if(rels[i].large_prime==0||rels[i].large_prime<0)
            usable.push_back(i);
    }

    fprintf(stderr,"%zu usable rels (need %d), %d polys\n",usable.size(),target,poly_count);

    if((int)usable.size()<=fb_size){
        fprintf(stderr,"Not enough.\n");
        printf("FAIL\n");goto cleanup;
    }

    if((int)usable.size()>target+200)usable.resize(target+200);

    {
        int nrels=usable.size();
        fprintf(stderr,"LinAlg: %dx%d...\n",nrels,fb_size);

        std::vector<BitRow> mat(nrels);
        for(int i=0;i<nrels;i++){
            mat[i].init(fb_size,nrels);mat[i].shist(i);
            auto&exps=rels[usable[i]].exps;
            for(int j=0;j<fb_size;j++)if(exps[j]&1)mat[i].sbit(j);
        }

        int rank=0;
        for(int col=0;col<fb_size&&rank<nrels;col++){
            int pivot=-1;
            for(int row=rank;row<nrels;row++)if(mat[row].gbit(col)){pivot=row;break;}
            if(pivot==-1)continue;
            if(pivot!=rank)std::swap(mat[pivot],mat[rank]);
            for(int row=0;row<nrels;row++)
                if(row!=rank&&mat[row].gbit(col))mat[row].xr(mat[rank]);
            rank++;
        }
        fprintf(stderr,"Rank %d, null %d\n",rank,nrels-rank);

        int found=0;
        mpz_t px,py;mpz_init(px);mpz_init(py);

        for(int row=0;row<nrels&&!found;row++){
            if(!mat[row].zb())continue;

            mpz_set_ui(px,1);
            std::vector<long> se(fb_size,0);
            std::vector<long> lps; /* large primes from combined rels */

            for(int i=0;i<nrels;i++){
                if(!mat[row].ghist(i))continue;
                int ri=usable[i];
                mpz_mul(px,px,rels[ri].val);
                mpz_mod(px,px,N);
                for(int j=0;j<fb_size;j++)se[j]+=rels[ri].exps[j];
                if(rels[ri].large_prime<0)lps.push_back(-rels[ri].large_prime);
            }

            bool ok=true;
            for(int j=0;j<fb_size;j++)if(se[j]&1){ok=false;break;}
            if(!ok)continue;

            mpz_set_ui(py,1);
            for(int j=1;j<fb_size;j++){
                if(se[j]==0)continue;
                mpz_set_ui(tmp,fb[j]);
                mpz_powm_ui(tmp,tmp,se[j]/2,N);
                mpz_mul(py,py,tmp);
                mpz_mod(py,py,N);
            }
            /* Include large primes: each LP appears twice (from the two LP rels that were combined),
             * so LP^2 in the product. We need LP^1 in y. */
            for(long lp:lps){
                mpz_set_ui(tmp,lp);
                mpz_mul(py,py,tmp);
                mpz_mod(py,py,N);
            }

            for(int sign=0;sign<2&&!found;sign++){
                if(sign==0)mpz_sub(tmp,px,py);else mpz_add(tmp,px,py);
                mpz_gcd(g,tmp,N);
                if(mpz_cmp_ui(g,1)>0&&mpz_cmp(g,N)<0){
                    struct timespec tnow;clock_gettime(CLOCK_MONOTONIC,&tnow);
                    double elapsed=(tnow.tv_sec-tstart.tv_sec)+(tnow.tv_nsec-tstart.tv_nsec)/1e9;
                    gmp_printf("FACTOR: %Zd\n",g);
                    mpz_t cof;mpz_init(cof);mpz_divexact(cof,N,g);
                    gmp_printf("COFACTOR: %Zd\n",cof);
                    fprintf(stderr,"SIQS factored in %.3fs (%d polys, %zu rels)\n",elapsed,poly_count,usable.size());
                    mpz_clear(cof);found=1;
                }
            }
        }
        mpz_clear(px);mpz_clear(py);
        if(!found){fprintf(stderr,"No factor from deps.\n");printf("FAIL\n");}
    }

cleanup:
    for(auto&r:rels)mpz_clear(r.val);
    mpz_clear(N);mpz_clear(tmp);mpz_clear(g);mpz_clear(target_a);
    gmp_randclear(rstate);
    return 0;
}
