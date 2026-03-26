#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>

static gmp_randstate_t rstate;

/* ---- Trial Division ---- */
int trial_division(mpz_t n, mpz_t factor, unsigned long limit) {
    if (mpz_divisible_ui_p(n, 2)) { mpz_set_ui(factor, 2); return 1; }
    if (mpz_divisible_ui_p(n, 3)) { mpz_set_ui(factor, 3); return 1; }
    for (unsigned long i = 5, w = 2; i <= limit; i += w, w = 6 - w)
        if (mpz_divisible_ui_p(n, i)) { mpz_set_ui(factor, i); return 1; }
    return 0;
}

/* ---- Pollard Rho (Brent, batched GCD) ---- */
int pollard_rho(mpz_t n, mpz_t factor, unsigned long max_iters) {
    mpz_t y, c, g, q, x, ys, tmp;
    mpz_inits(y, c, g, q, x, ys, tmp, NULL);
    int found = 0;
    for (int attempt = 0; attempt < 5 && !found; attempt++) {
        mpz_urandomm(y, rstate, n);
        mpz_urandomm(c, rstate, n);
        if (mpz_sgn(c) == 0) mpz_set_ui(c, 1);
        mpz_set_ui(q, 1);
        unsigned long r = 1, total = 0;
        while (total < max_iters && !found) {
            mpz_set(x, y);
            for (unsigned long i = 0; i < r; i++) {
                mpz_mul(y, y, y); mpz_add(y, y, c); mpz_mod(y, y, n);
            }
            for (unsigned long k = 0; k < r && !found; ) {
                mpz_set(ys, y);
                unsigned long batch = r - k; if (batch > 256) batch = 256;
                for (unsigned long i = 0; i < batch; i++) {
                    mpz_mul(y, y, y); mpz_add(y, y, c); mpz_mod(y, y, n);
                    mpz_sub(tmp, x, y); mpz_abs(tmp, tmp);
                    mpz_mul(q, q, tmp); mpz_mod(q, q, n);
                }
                mpz_gcd(g, q, n);
                if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, n) < 0) {
                    mpz_set(factor, g); found = 1;
                } else if (mpz_cmp(g, n) == 0) {
                    mpz_set(y, ys);
                    for (unsigned long i = 0; i < batch; i++) {
                        mpz_mul(y, y, y); mpz_add(y, y, c); mpz_mod(y, y, n);
                        mpz_sub(tmp, x, y); mpz_abs(tmp, tmp);
                        mpz_gcd(g, tmp, n);
                        if (mpz_cmp_ui(g, 1) > 0) {
                            if (mpz_cmp(g, n) < 0) { mpz_set(factor, g); found = 1; }
                            break;
                        }
                    }
                }
                k += batch; total += batch;
                if (total >= max_iters) break;
            }
            r *= 2;
        }
    }
    mpz_clears(y, c, g, q, x, ys, tmp, NULL);
    return found;
}

/* ---- Primes ---- */
static unsigned int *primes; static int nprimes;
void sieve(unsigned long lim) {
    char *s = calloc(lim+1, 1); s[0]=s[1]=1;
    for (unsigned long i=2; i*i<=lim; i++) if(!s[i]) for(unsigned long j=i*i;j<=lim;j+=i)s[j]=1;
    nprimes=0; for(unsigned long i=2;i<=lim;i++) if(!s[i]) nprimes++;
    primes=malloc(nprimes*sizeof(unsigned int)); int k=0;
    for(unsigned long i=2;i<=lim;i++) if(!s[i]) primes[k++]=(unsigned int)i;
    free(s);
}

/* Precompute product of prime powers <= B as mpz */
void compute_exponent(mpz_t e, unsigned long B) {
    mpz_set_ui(e, 1);
    for (int i = 0; i < nprimes && primes[i] <= B; i++) {
        unsigned long p = primes[i], pp = p;
        while (pp <= B) { mpz_mul_ui(e, e, p); if (pp > B/p) break; pp *= p; }
    }
}

/* ---- Pollard p-1 ---- */
int pollard_pm1(mpz_t n, mpz_t factor, mpz_t e_pm1) {
    mpz_t a, g, tmp;
    mpz_inits(a, g, tmp, NULL);
    mpz_set_ui(a, 2);
    mpz_powm(a, a, e_pm1, n);
    mpz_sub_ui(tmp, a, 1);
    mpz_gcd(g, tmp, n);
    int found = (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, n) < 0);
    if (found) mpz_set(factor, g);
    mpz_clears(a, g, tmp, NULL);
    return found;
}

/* ---- Williams p+1 ---- */
int williams_pp1(mpz_t n, mpz_t factor, mpz_t e_pp1) {
    mpz_t V, V1, tmp, g;
    mpz_inits(V, V1, tmp, g, NULL);
    for (int seed = 3; seed <= 13; seed++) {
        mpz_set_ui(V, 2); mpz_set_ui(V1, seed);
        int bits = mpz_sizeinbase(e_pp1, 2);
        for (int i = bits - 1; i >= 0; i--) {
            if (mpz_tstbit(e_pp1, i)) {
                mpz_mul(V, V, V1); mpz_sub_ui(V, V, seed); mpz_mod(V, V, n);
                mpz_mul(V1, V1, V1); mpz_sub_ui(V1, V1, 2); mpz_mod(V1, V1, n);
            } else {
                mpz_mul(V1, V, V1); mpz_sub_ui(V1, V1, seed); mpz_mod(V1, V1, n);
                mpz_mul(V, V, V); mpz_sub_ui(V, V, 2); mpz_mod(V, V, n);
            }
        }
        mpz_sub_ui(tmp, V, 2);
        mpz_gcd(g, tmp, n);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, n) < 0) {
            mpz_set(factor, g);
            mpz_clears(V, V1, tmp, g, NULL);
            return 1;
        }
    }
    mpz_clears(V, V1, tmp, g, NULL);
    return 0;
}

/* ---- ECM (Montgomery form, stage 1 only) ---- */
typedef struct { mpz_t X, Z; } pt;
void pt_init(pt *P) { mpz_init(P->X); mpz_init(P->Z); }
void pt_clear(pt *P) { mpz_clear(P->X); mpz_clear(P->Z); }
void pt_set(pt *R, pt *P) { mpz_set(R->X, P->X); mpz_set(R->Z, P->Z); }

static mpz_t _t1, _t2, _t3, _t4; /* shared temps */

static inline void ecm_dbl(pt *R, pt *P, mpz_t a24, mpz_t n) {
    mpz_add(_t1, P->X, P->Z); mpz_mul(_t1, _t1, _t1); mpz_mod(_t1, _t1, n);
    mpz_sub(_t2, P->X, P->Z); mpz_mul(_t2, _t2, _t2); mpz_mod(_t2, _t2, n);
    mpz_sub(_t3, _t1, _t2);
    mpz_mul(R->X, _t1, _t2); mpz_mod(R->X, R->X, n);
    mpz_mul(_t4, a24, _t3); mpz_mod(_t4, _t4, n);
    mpz_add(_t4, _t4, _t2);
    mpz_mul(R->Z, _t3, _t4); mpz_mod(R->Z, R->Z, n);
}

static inline void ecm_add(pt *R, pt *P, pt *Q, pt *D, mpz_t n) {
    mpz_add(_t1, P->X, P->Z); mpz_sub(_t2, Q->X, Q->Z);
    mpz_mul(_t3, _t1, _t2); mpz_mod(_t3, _t3, n);
    mpz_sub(_t1, P->X, P->Z); mpz_add(_t2, Q->X, Q->Z);
    mpz_mul(_t4, _t1, _t2); mpz_mod(_t4, _t4, n);
    mpz_add(_t1, _t3, _t4); mpz_mul(_t1, _t1, _t1); mpz_mod(_t1, _t1, n);
    mpz_mul(R->X, _t1, D->Z); mpz_mod(R->X, R->X, n);
    mpz_sub(_t1, _t3, _t4); mpz_mul(_t1, _t1, _t1); mpz_mod(_t1, _t1, n);
    mpz_mul(R->Z, _t1, D->X); mpz_mod(R->Z, R->Z, n);
}

/* Scalar multiply by mpz using Montgomery ladder */
void ecm_scalar(pt *R, pt *P, mpz_t k, mpz_t a24, mpz_t n) {
    pt R0, R1, tmp;
    pt_init(&R0); pt_init(&R1); pt_init(&tmp);
    pt_set(&R0, P);
    ecm_dbl(&R1, P, a24, n);
    for (int i = mpz_sizeinbase(k, 2) - 2; i >= 0; i--) {
        if (mpz_tstbit(k, i)) {
            ecm_add(&tmp, &R1, &R0, P, n); pt_set(&R0, &tmp);
            ecm_dbl(&tmp, &R1, a24, n); pt_set(&R1, &tmp);
        } else {
            ecm_add(&tmp, &R0, &R1, P, n); pt_set(&R1, &tmp);
            ecm_dbl(&tmp, &R0, a24, n); pt_set(&R0, &tmp);
        }
    }
    pt_set(R, &R0);
    pt_clear(&R0); pt_clear(&R1); pt_clear(&tmp);
}

int ecm_one_curve(mpz_t n, mpz_t factor, mpz_t e_ecm) {
    mpz_t sigma, u, v, A, a24, g, t;
    mpz_inits(sigma, u, v, A, a24, g, t, NULL);
    pt P;
    pt_init(&P);

    mpz_urandomm(sigma, rstate, n);
    mpz_add_ui(sigma, sigma, 6);
    mpz_mul(u, sigma, sigma); mpz_sub_ui(u, u, 5); mpz_mod(u, u, n);
    mpz_mul_ui(v, sigma, 4); mpz_mod(v, v, n);
    mpz_powm_ui(P.X, u, 3, n);
    mpz_powm_ui(P.Z, v, 3, n);

    mpz_sub(A, v, u); mpz_powm_ui(A, A, 3, n);
    mpz_mul_ui(t, u, 3); mpz_add(t, t, v); mpz_mod(t, t, n);
    mpz_mul(A, A, t); mpz_mod(A, A, n);
    mpz_powm_ui(t, u, 3, n); mpz_mul(t, t, v); mpz_mul_ui(t, t, 4); mpz_mod(t, t, n);
    if (!mpz_invert(t, t, n)) {
        mpz_gcd(g, t, n);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, n) < 0) {
            mpz_set(factor, g); pt_clear(&P); mpz_clears(sigma,u,v,A,a24,g,t,NULL); return 1;
        }
        pt_clear(&P); mpz_clears(sigma,u,v,A,a24,g,t,NULL); return 0;
    }
    mpz_mul(A, A, t); mpz_sub_ui(A, A, 2); mpz_mod(A, A, n);
    mpz_add_ui(a24, A, 2);
    mpz_set_ui(t, 4);
    if (!mpz_invert(t, t, n)) { pt_clear(&P); mpz_clears(sigma,u,v,A,a24,g,t,NULL); return 0; }
    mpz_mul(a24, a24, t); mpz_mod(a24, a24, n);

    ecm_scalar(&P, &P, e_ecm, a24, n);

    mpz_gcd(g, P.Z, n);
    int found = (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, n) < 0);
    if (found) mpz_set(factor, g);
    pt_clear(&P); mpz_clears(sigma,u,v,A,a24,g,t,NULL);
    return found;
}

int ecm_factor(mpz_t n, mpz_t factor, mpz_t e_ecm, int ncurves) {
    for (int c = 0; c < ncurves; c++)
        if (ecm_one_curve(n, factor, e_ecm)) return 1;
    return 0;
}

/* ---- Main ---- */
int main(void) {
    gmp_randinit_mt(rstate);
    gmp_randseed_ui(rstate, (unsigned long)time(NULL) ^ 0xCAFE);
    mpz_inits(_t1, _t2, _t3, _t4, NULL);

    fprintf(stderr, "Sieving...\n");
    sieve(1000000UL);
    fprintf(stderr, "%d primes up to %u\n", nprimes, primes[nprimes-1]);

    /* Precompute exponents for various B1 values */
    mpz_t e_pm1, e_pp1, e_ecm_2k, e_ecm_11k, e_ecm_50k, e_ecm_250k, e_ecm_1m;
    mpz_inits(e_pm1, e_pp1, e_ecm_2k, e_ecm_11k, e_ecm_50k, e_ecm_250k, e_ecm_1m, NULL);

    fprintf(stderr, "Precomputing exponents...\n");
    compute_exponent(e_pm1, 1000000UL);
    compute_exponent(e_pp1, 100000UL);
    compute_exponent(e_ecm_2k, 2000UL);
    compute_exponent(e_ecm_11k, 11000UL);
    compute_exponent(e_ecm_50k, 50000UL);
    compute_exponent(e_ecm_250k, 250000UL);
    compute_exponent(e_ecm_1m, 1000000UL);
    fprintf(stderr, "Ready.\n");

    char line[512];
    while (fgets(line, sizeof(line), stdin)) {
        int digits; char Nstr[256];
        if (sscanf(line, "%d %255s", &digits, Nstr) != 2) continue;

        mpz_t n, factor;
        mpz_init(n); mpz_init(factor);
        mpz_set_str(n, Nstr, 10);

        struct timespec t0, t1;
        clock_gettime(CLOCK_MONOTONIC, &t0);
        char method[32] = ""; int factored = 0;

        fprintf(stderr, "[%d] %s\n", digits, Nstr);

        if (!factored && trial_division(n, factor, 1000000UL))
            { strcpy(method, "trial"); factored = 1; }

        /* Scale effort by digit count */
        int rho_iters = (digits <= 40) ? 3000000 : 1000000;
        int ecm2k_n   = (digits <= 40) ? 25 : 10;
        int ecm11k_n  = (digits <= 50) ? 50 : 20;
        int ecm50k_n  = (digits <= 55) ? 30 : 15;
        int ecm250k_n = (digits <= 60) ? 15 : 8;
        int ecm1m_n   = (digits <= 65) ? 8 : 0;

        if (!factored && pollard_rho(n, factor, rho_iters))
            { strcpy(method, "rho"); factored = 1; }

        if (!factored && pollard_pm1(n, factor, e_pm1))
            { strcpy(method, "p-1"); factored = 1; }

        if (!factored && williams_pp1(n, factor, e_pp1))
            { strcpy(method, "p+1"); factored = 1; }

        if (!factored && ecm2k_n > 0 && ecm_factor(n, factor, e_ecm_2k, ecm2k_n))
            { strcpy(method, "ECM-2K"); factored = 1; }

        if (!factored && ecm11k_n > 0 && ecm_factor(n, factor, e_ecm_11k, ecm11k_n))
            { strcpy(method, "ECM-11K"); factored = 1; }

        if (!factored && ecm50k_n > 0 && ecm_factor(n, factor, e_ecm_50k, ecm50k_n))
            { strcpy(method, "ECM-50K"); factored = 1; }

        if (!factored && ecm250k_n > 0 && ecm_factor(n, factor, e_ecm_250k, ecm250k_n))
            { strcpy(method, "ECM-250K"); factored = 1; }

        if (!factored && ecm1m_n > 0 && ecm_factor(n, factor, e_ecm_1m, ecm1m_n))
            { strcpy(method, "ECM-1M"); factored = 1; }

        clock_gettime(CLOCK_MONOTONIC, &t1);
        double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;

        if (factored) {
            mpz_t q; mpz_init(q);
            mpz_divexact(q, n, factor);
            if (mpz_cmp(factor, q) > 0) mpz_swap(factor, q);
            gmp_printf("%d\t%s\t%s\t%Zd\t%Zd\t%.6f\n", digits, Nstr, method, factor, q, elapsed);
            fprintf(stderr, "  [%s] %.3fs\n", method, elapsed);
            mpz_clear(q);
        } else {
            printf("%d\t%s\tFAILED\t\t\t%.6f\n", digits, Nstr, elapsed);
            fprintf(stderr, "  FAILED %.3fs\n", elapsed);
        }
        fflush(stdout);
        mpz_clear(n); mpz_clear(factor);
    }

    mpz_clears(e_pm1, e_pp1, e_ecm_2k, e_ecm_11k, e_ecm_50k, e_ecm_250k, e_ecm_1m, NULL);
    mpz_clears(_t1, _t2, _t3, _t4, NULL);
    free(primes); gmp_randclear(rstate);
    return 0;
}
