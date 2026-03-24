/* Verify a single MPQS relation */
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static long msqrt(long n,long p){
    if(p==2)return n&1; n=((n%p)+p)%p; if(n==0)return 0;
    mpz_t a,m,r; mpz_inits(a,m,r,NULL);
    mpz_set_si(a,n); mpz_set_si(m,p);
    mpz_powm_ui(r,a,(p-1)/2,m);
    if(mpz_cmp_ui(r,1)!=0){mpz_clears(a,m,r,NULL); return -1;}
    if((p&3)==3){mpz_powm_ui(r,a,(p+1)/4,m); long R=mpz_get_si(r); mpz_clears(a,m,r,NULL); return R;}
    /* Full Tonelli-Shanks - skip for now, most primes are 3 mod 4 at small sizes */
    long Q=p-1,S=0; while((Q&1)==0){Q>>=1;S++;}
    long z; for(z=2;z<p;z++){mpz_set_si(a,z);mpz_powm_ui(r,a,(p-1)/2,m);if(mpz_cmp_si(r,p-1)==0)break;}
    mpz_t c,t,R,b; mpz_inits(c,t,R,b,NULL);
    mpz_set_si(a,z);mpz_powm_ui(c,a,Q,m);mpz_set_si(a,n);mpz_powm_ui(t,a,Q,m);mpz_powm_ui(R,a,(Q+1)/2,m);
    long Mv=S;
    while(1){if(mpz_cmp_ui(t,1)==0){long res=mpz_get_si(R);mpz_clears(a,m,r,c,t,R,b,NULL);return res;}
    long i=0;mpz_set(b,t);while(mpz_cmp_ui(b,1)!=0){mpz_mul(b,b,b);mpz_mod(b,b,m);i++;if(i>=Mv){mpz_clears(a,m,r,c,t,R,b,NULL);return-1;}}
    mpz_set(b,c);for(long j=0;j<Mv-i-1;j++){mpz_mul(b,b,b);mpz_mod(b,b,m);}
    Mv=i;mpz_mul(c,b,b);mpz_mod(c,c,m);mpz_mul(t,t,c);mpz_mod(t,t,m);mpz_mul(R,R,b);mpz_mod(R,R,m);}
}

int main() {
    mpz_t N, A, B, C, q, tmp, q_inv, gx, sv, sv2, gx_mod;
    mpz_inits(N, A, B, C, q, tmp, q_inv, gx, sv, sv2, gx_mod, NULL);
    mpz_set_str(N, "202930257038936745993812575957", 10);

    /* Generate one MPQS polynomial */
    mpz_mul_ui(tmp, N, 2); mpz_sqrt(tmp, tmp); mpz_sqrt(tmp, tmp); /* (2N)^{1/4} */
    int M = 10000000;
    mpz_tdiv_q_ui(tmp, tmp, (unsigned long)sqrt(M));
    mpz_add_ui(tmp, tmp, 42);
    mpz_nextprime(q, tmp);

    /* Find a q where N is QR */
    for (int att = 0; att < 100; att++) {
        if (!mpz_fits_ulong_p(q)) { mpz_nextprime(q, q); continue; }
        unsigned long qv = mpz_get_ui(q);
        long nq = mpz_fdiv_ui(N, qv);
        long sq = msqrt(nq, qv);
        if (sq <= 0) { mpz_nextprime(q, q); continue; }

        gmp_printf("q = %Zd, sq = %ld\n", q, sq);

        mpz_mul(A, q, q);
        mpz_set_ui(B, sq);

        /* Hensel lift */
        mpz_set_ui(tmp, (unsigned long)sq * sq);
        mpz_sub(tmp, N, tmp);
        mpz_tdiv_q_ui(tmp, tmp, qv);
        mpz_t inv, pz; mpz_inits(inv, pz, NULL);
        mpz_set_ui(inv, 2 * sq % qv); mpz_set_ui(pz, qv);
        mpz_invert(inv, inv, pz);
        mpz_mul(tmp, tmp, inv); mpz_mod_ui(tmp, tmp, qv);
        mpz_mul_ui(tmp, tmp, qv); mpz_add(B, B, tmp);
        mpz_clears(inv, pz, NULL);

        mpz_mul(C, B, B); mpz_sub(C, C, N); mpz_tdiv_q(C, C, A);

        gmp_printf("A = %Zd\nB = %Zd\nC = %Zd\n", A, B, C);

        /* Verify B^2 ≡ N (mod A) */
        mpz_mul(tmp, B, B); mpz_sub(tmp, tmp, N);
        printf("B^2 - N divisible by A? %d\n", mpz_divisible_p(tmp, A) != 0);

        /* Compute q^{-1} mod N */
        mpz_set_ui(tmp, qv);
        mpz_invert(q_inv, tmp, N);

        /* Test x = 0: g(0) = C */
        long x = 0;
        gmp_printf("\nTesting x = %ld:\n", x);
        mpz_set(gx, C);
        gmp_printf("g(x) = %Zd\n", gx);

        /* sv = (Ax + B) * q^{-1} mod N */
        mpz_set_si(tmp, x); mpz_mul(sv, A, tmp); mpz_add(sv, sv, B);
        gmp_printf("Ax+B = %Zd\n", sv);
        mpz_mul(sv, sv, q_inv); mpz_mod(sv, sv, N);
        gmp_printf("sv = (Ax+B)*q^{-1} mod N = %Zd\n", sv);

        /* Verify: sv^2 ≡ g(x) (mod N) */
        mpz_mul(sv2, sv, sv); mpz_mod(sv2, sv2, N);
        mpz_mod(gx_mod, gx, N);
        if (mpz_sgn(gx_mod) < 0) mpz_add(gx_mod, gx_mod, N);
        gmp_printf("sv^2 mod N = %Zd\ng(x) mod N = %Zd\n", sv2, gx_mod);
        printf("Match? %d\n", mpz_cmp(sv2, gx_mod) == 0);

        /* Also test: (Ax+B)^2 ≡ A*g(x) + N (over Z) */
        mpz_set_si(tmp, x); mpz_mul(tmp, A, tmp); mpz_add(tmp, tmp, B);
        mpz_mul(tmp, tmp, tmp); /* (Ax+B)^2 */
        mpz_t agn; mpz_init(agn);
        mpz_mul(agn, A, gx); mpz_add(agn, agn, N);
        gmp_printf("(Ax+B)^2 = %Zd\nA*g(x)+N = %Zd\nEqual? %d\n", tmp, agn, mpz_cmp(tmp, agn) == 0);
        mpz_clear(agn);

        break;
    }

    mpz_clears(N, A, B, C, q, tmp, q_inv, gx, sv, sv2, gx_mod, NULL);
    return 0;
}
