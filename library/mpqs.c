/*
 * Multiple Polynomial Quadratic Sieve (MPQS)
 *
 * Uses polynomials f(x) = Ax^2 + 2Bx + C where:
 *   A = prime q from factor base, q ≈ sqrt(2N)/M
 *   B = sqrt(N) mod q
 *   C = (B^2 - N) / q
 *   (qx + B)^2 ≡ q * f(x) (mod N)
 *
 * With large prime variation and double-large-prime (2LP).
 *
 * Usage: ./mpqs <number>
 * Single-threaded. Seed=42.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

#define MAX_FB 12000
#define MAX_RELS 50000

/* ---- Factor base ---- */
typedef struct {
    unsigned int p;
    int logp;
    unsigned long sqrt_n; /* sqrt(N) mod p */
} fb_t;

static fb_t FB[MAX_FB];
static int fb_count;

/* ---- Relations ---- */
typedef struct {
    mpz_t sqrt_val;   /* qx + B (mod N) */
    int *expo;         /* exponent vector over FB (fb_count entries) */
    int neg;
    unsigned long lp1, lp2; /* 0 = none; one or two large primes */
} rel_t;

static rel_t RELS[MAX_RELS];
static int rel_count = 0;

/* ---- Primes ---- */
static int *PRIMES;
static int NPRIMES;

static void gen_primes(int lim) {
    char *s = calloc(lim+1,1);
    int c=0;
    for(int i=2;i<=lim;i++) if(!s[i]){c++;for(long j=(long)i*i;j<=lim;j+=i)s[j]=1;}
    PRIMES=malloc(c*sizeof(int)); NPRIMES=0;
    memset(s,0,lim+1);
    for(int i=2;i<=lim;i++) if(!s[i]){PRIMES[NPRIMES++]=i;for(long j=(long)i*i;j<=lim;j+=i)s[j]=1;}
    free(s);
}

/* Tonelli-Shanks */
static unsigned long tsqrt(unsigned long a, unsigned long p) {
    if(p==2) return a&1;
    if(a==0) return 0;
    mpz_t az,pz,rz;
    mpz_init_set_ui(az,a); mpz_init_set_ui(pz,p); mpz_init(rz);
    mpz_powm_ui(rz,az,(p-1)/2,pz);
    if(mpz_cmp_ui(rz,1)!=0){mpz_clear(az);mpz_clear(pz);mpz_clear(rz);return(unsigned long)-1;}
    if((p&3)==3){mpz_powm_ui(rz,az,(p+1)/4,pz);unsigned long r=mpz_get_ui(rz);mpz_clear(az);mpz_clear(pz);mpz_clear(rz);return r;}
    unsigned long Q=p-1,S=0;while(!(Q&1)){Q>>=1;S++;}
    unsigned long z=2;
    mpz_t zz; mpz_init(zz);
    while(1){mpz_set_ui(zz,z);mpz_powm_ui(zz,zz,(p-1)/2,pz);if(mpz_cmp_ui(zz,p-1)==0)break;z++;}
    mpz_t M,cc,t,R,b;
    mpz_init_set_ui(M,S);mpz_init(cc);mpz_init(t);mpz_init(R);mpz_init(b);
    mpz_set_ui(cc,z);mpz_powm_ui(cc,cc,Q,pz);
    mpz_set_ui(t,a);mpz_powm_ui(t,t,Q,pz);
    mpz_set_ui(R,a);mpz_powm_ui(R,R,(Q+1)/2,pz);
    while(1){
        if(mpz_cmp_ui(t,1)==0){unsigned long r=mpz_get_ui(R);mpz_clear(M);mpz_clear(cc);mpz_clear(t);mpz_clear(R);mpz_clear(b);mpz_clear(zz);mpz_clear(az);mpz_clear(pz);mpz_clear(rz);return r;}
        mpz_set(zz,t);unsigned long i=0;
        while(mpz_cmp_ui(zz,1)!=0){mpz_mul(zz,zz,zz);mpz_mod(zz,zz,pz);i++;}
        unsigned long m=mpz_get_ui(M);
        mpz_set(b,cc);for(unsigned long j=0;j<m-i-1;j++){mpz_mul(b,b,b);mpz_mod(b,b,pz);}
        mpz_set_ui(M,i);
        mpz_mul(cc,b,b);mpz_mod(cc,cc,pz);
        mpz_mul(t,t,cc);mpz_mod(t,t,pz);
        mpz_mul(R,R,b);mpz_mod(R,R,pz);
    }
}

/* ---- Parameters ---- */
typedef struct { int fb_sz; int M; double thresh_sub; int lp_mult; } par_t;

static par_t get_par(int d) {
    if(d<=30)  return(par_t){100,  32768,  10, 40};
    if(d<=35)  return(par_t){200,  49152,  12, 50};
    if(d<=40)  return(par_t){350,  65536,  13, 50};
    if(d<=45)  return(par_t){600,  65536,  14, 50};
    if(d<=50)  return(par_t){1000, 98304,  15, 60};
    if(d<=55)  return(par_t){1600, 131072, 16, 60};
    if(d<=60)  return(par_t){2500, 131072, 17, 60};
    if(d<=65)  return(par_t){3500, 196608, 18, 60};
    if(d<=70)  return(par_t){4500, 262144, 18, 60};
    if(d<=75)  return(par_t){6000, 327680, 19, 60};
    if(d<=80)  return(par_t){7500, 393216, 20, 60};
    if(d<=85)  return(par_t){9000, 393216, 20, 60};
    if(d<=90)  return(par_t){10500,524288, 21, 60};
    if(d<=95)  return(par_t){11500,524288, 21, 60};
    return(par_t){12000,524288,21,60};
}

static mpz_t N;

/* Build factor base */
static int build_fb(int target) {
    int cnt = 0;
    for (int i = 0; i < NPRIMES && cnt < target; i++) {
        unsigned int p = PRIMES[i];
        unsigned long nm = mpz_fdiv_ui(N, p);
        if (p == 2) {
            FB[cnt].p = 2; FB[cnt].logp = 1; FB[cnt].sqrt_n = nm & 1; cnt++;
            continue;
        }
        unsigned long sr = tsqrt(nm, p);
        if (sr == (unsigned long)-1) continue;
        FB[cnt].p = p;
        FB[cnt].logp = (int)(log2((double)p) + 0.5);
        if (FB[cnt].logp < 1) FB[cnt].logp = 1;
        FB[cnt].sqrt_n = sr;
        cnt++;
    }
    return cnt;
}

/* Find index of prime q in FB, or -1 */
static int fb_index(unsigned int q) {
    for (int i = 0; i < fb_count; i++)
        if (FB[i].p == q) return i;
    return -1;
}

/* Store a relation */
static void store_rel(mpz_t sv, int *expo, int neg, unsigned long lp1, unsigned long lp2) {
    if (rel_count >= MAX_RELS) return;
    rel_t *r = &RELS[rel_count];
    mpz_init_set(r->sqrt_val, sv);
    r->expo = malloc(fb_count * sizeof(int));
    memcpy(r->expo, expo, fb_count * sizeof(int));
    r->neg = neg;
    r->lp1 = lp1;
    r->lp2 = lp2;
    rel_count++;
}

/* ---- Large prime graph for 2LP ---- */
#define HT_SIZE (1 << 20)
typedef struct ht_node { unsigned long key; int rel_idx; struct ht_node *next; } ht_node;
static ht_node *HT[HT_SIZE];

static int ht_hash(unsigned long k) { return (int)((k * 2654435761UL) >> 12) & (HT_SIZE-1); }

static int ht_find(unsigned long k) {
    for (ht_node *n = HT[ht_hash(k)]; n; n = n->next)
        if (n->key == k) return n->rel_idx;
    return -1;
}
static void ht_add(unsigned long k, int idx) {
    int h = ht_hash(k);
    ht_node *n = malloc(sizeof(ht_node));
    n->key = k; n->rel_idx = idx; n->next = HT[h]; HT[h] = n;
}

/* ---- Gaussian elimination mod 2 ---- */
/* Returns rank; matrix has identity augmented for tracking */
static int gauss(unsigned char **mat, int nr, int nc) {
    int rank = 0;
    for (int c = 0; c < nc && rank < nr; c++) {
        int pr = -1;
        for (int r = rank; r < nr; r++) if (mat[r][c]) { pr = r; break; }
        if (pr < 0) continue;
        if (pr != rank) { unsigned char *t = mat[pr]; mat[pr] = mat[rank]; mat[rank] = t; }
        for (int r = 0; r < nr; r++)
            if (r != rank && mat[r][c])
                for (int k = 0; k < nc + nr; k++) mat[r][k] ^= mat[rank][k];
        rank++;
    }
    return rank;
}

int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }
    mpz_init_set_str(N, argv[1], 10);
    int digits = strlen(argv[1]);
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    /* Trial division */
    for (unsigned long p = 2; p < 1000000; p++) {
        if (mpz_divisible_ui_p(N, p)) {
            mpz_t q; mpz_init(q); mpz_divexact_ui(q, N, p);
            clock_gettime(CLOCK_MONOTONIC, &t1);
            gmp_printf("%lu %Zd\n", p, q);
            fprintf(stderr, "Trial: time=%.3fs\n", (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)/1e9);
            return 0;
        }
    }

    par_t par = get_par(digits);
    gen_primes(par.fb_sz * 50);
    fb_count = build_fb(par.fb_sz);
    int M = par.M;
    int sieve_len = 2 * M;

    unsigned long lp_bound = (unsigned long)FB[fb_count-1].p * par.lp_mult;
    int needed = fb_count + 50;

    fprintf(stderr, "MPQS: %d digits, FB=%d (max=%u), M=%d, LP=%lu\n",
            digits, fb_count, FB[fb_count-1].p, M, lp_bound);

    /* Target A: A ≈ sqrt(2N) / M */
    mpz_t target_a, tmp1, tmp2;
    mpz_init(target_a); mpz_init(tmp1); mpz_init(tmp2);
    mpz_mul_ui(target_a, N, 2);
    mpz_sqrt(target_a, target_a);
    mpz_tdiv_q_ui(target_a, target_a, M);

    /* For MPQS, max|f(x)| ≈ M * sqrt(N/2) when A ≈ sqrt(2N)/M */
    double log2fmax = log2((double)M) + (mpz_sizeinbase(N, 2) - 1) / 2.0;
    double log2_lpb = log2((double)lp_bound);
    int thresh = (int)(log2fmax - log2_lpb - 3.0);
    if (thresh < 15) thresh = 15;

    fprintf(stderr, "Threshold: %d (log2fmax=%.1f, log2_lpb=%.1f)\n", thresh, log2fmax, log2_lpb);

    unsigned char *sieve = malloc(sieve_len);
    int *soln1 = malloc(fb_count * sizeof(int));
    int *soln2 = malloc(fb_count * sizeof(int));
    int *expo = calloc(fb_count, sizeof(int));
    mpz_t A, B, C, fval, sqv, cofactor;
    mpz_init(A); mpz_init(B); mpz_init(C);
    mpz_init(fval); mpz_init(sqv); mpz_init(cofactor);

    gmp_randstate_t rs;
    gmp_randinit_default(rs);
    gmp_randseed_ui(rs, 42);

    int smooth_cnt = 0, partial_cnt = 0, combined_cnt = 0;
    int polys = 0;
    memset(HT, 0, sizeof(HT));

    /* Choose A primes from FB that are near sqrt(target_a) - actually A should be near target_a directly
     * We'll try A = single FB prime near target_a. If target_a is too large for a single prime,
     * we'll use a random product of 2 primes. */

    /* Actually for single-prime MPQS, A should be near target_a, and A must be a prime
     * where N is a QR. We pick from among our FB primes that are closest to target_a.
     * But FB primes are typically much smaller than target_a (which can be huge).
     *
     * Better approach: pick a prime q (not necessarily in FB) near target_a, check N is QR mod q.
     * Then in the exponent vector, we include q's contribution to f(x)*q.
     * Actually: (qx+B)^2 = q*f(x) + N, so (qx+B)^2 ≡ q*f(x) (mod N).
     * For the congruence of squares, we need q*f(x) to be smooth.
     * q is a large prime (near target_a), so q*f(x) has one large factor q.
     * If f(x) itself is smooth, then q*f(x) has q as a large prime factor.
     * We can treat q as a "polynomial large prime" that's the same for all relations
     * from this polynomial.
     *
     * But this means relations from different polynomials have different q's,
     * and we can only combine relations from the same polynomial.
     * That's bad - it means we can't mix relations across polynomials.
     *
     * The standard MPQS approach: A should be the SQUARE of a prime, A = q^2.
     * Then (qx + B)^2 = A*x^2 + 2qBx + B^2 = ... actually this doesn't simplify.
     *
     * Let me re-derive. Standard MPQS:
     * Choose A (could be prime or composite, must be a QR mod N).
     * Find B with B^2 ≡ N (mod A) and |B| ≤ A/2.
     * Set C = (B^2 - N)/A (this is an integer).
     * Polynomial g(x) = Ax^2 + 2Bx + C.
     * Then A*g(x) = (Ax + B)^2 - N.
     * So (Ax + B)^2 ≡ A*g(x) (mod N).
     *
     * Now: if A is itself a perfect square, A = a^2, then:
     * (Ax + B)^2 ≡ a^2 * g(x) (mod N)
     * (Ax + B / a)... no, B is not divisible by a in general.
     *
     * Hmm. The standard approach for combining across polynomials:
     *
     * From one polynomial: (Ax+B)^2 ≡ A*g(x) (mod N)
     * If A is the same for all polynomials (SIQS approach), then A factors out
     * and we need g(x) to be smooth.
     * When combining relations from different polynomials with the same A:
     * ∏(Ax_i + B_i)^2 ≡ A^k * ∏g(x_i) (mod N) where k = number of relations
     * If k is even and ∏g(x_i) is a perfect square, then:
     * ∏(Ax_i + B_i)^2 ≡ (A^{k/2} * Y)^2 (mod N)
     *
     * For relations from DIFFERENT A values:
     * (A1*x1 + B1)^2 ≡ A1*g1(x1) (mod N)
     * (A2*x2 + B2)^2 ≡ A2*g2(x2) (mod N)
     * Product: ∏(Ai*xi + Bi)^2 ≡ ∏(Ai*gi(xi)) (mod N)
     * We need ∏(Ai*gi(xi)) to be a perfect square.
     * So the exponent vector should include the prime factorization of Ai as well.
     *
     * For MPQS with A = single prime q (different for each poly):
     * The exponent vector for each relation includes:
     * - Sign of g(x)
     * - Exponents of FB primes in g(x) (from trial division)
     * - Exponent of q (always 1, since we have q * g(x))
     *
     * But q is different for each polynomial, so we'd need columns for every q.
     * That's too many columns!
     *
     * The solution: choose A to be a PRODUCT of FB primes.
     * If A = q1 * q2 * ... * qk where qi are FB primes, then q1..qk are already
     * in the exponent vector. We just add 1 to each qi's exponent.
     *
     * This is the standard SIQS approach: A = product of ~k FB primes where
     * k is chosen so that A ≈ target_a.
     *
     * For simplicity, I'll use A = product of 2-3 FB primes from the upper range.
     */

    /* Select A primes from upper portion of FB */
    /* target_a has about digits/2 digits; we need FB primes that multiply to near target_a */
    int target_a_bits = mpz_sizeinbase(target_a, 2);
    fprintf(stderr, "target_a has %d bits\n", target_a_bits);

    /* We need k primes of size ~target_a_bits/k bits each */
    /* FB primes go up to FB[fb_count-1].p */
    int max_fb_bits = (int)(log2((double)FB[fb_count-1].p));
    int num_a_primes;
    if (target_a_bits <= max_fb_bits + 2) num_a_primes = 1;
    else if (target_a_bits <= 2 * max_fb_bits + 2) num_a_primes = 2;
    else if (target_a_bits <= 3 * max_fb_bits + 2) num_a_primes = 3;
    else if (target_a_bits <= 4 * max_fb_bits + 2) num_a_primes = 4;
    else if (target_a_bits <= 5 * max_fb_bits + 2) num_a_primes = 5;
    else if (target_a_bits <= 6 * max_fb_bits + 2) num_a_primes = 6;
    else if (target_a_bits <= 7 * max_fb_bits + 2) num_a_primes = 7;
    else if (target_a_bits <= 8 * max_fb_bits + 2) num_a_primes = 8;
    else if (target_a_bits <= 9 * max_fb_bits + 2) num_a_primes = 9;
    else num_a_primes = 10;

    /* Pick A primes from the upper half of FB */
    int a_start = fb_count / 2;  /* start of range for A primes */
    int a_range = fb_count - a_start;

    fprintf(stderr, "Using %d primes for A (a_range=%d, max_fb_bits=%d)\n",
            num_a_primes, a_range, max_fb_bits);

    /* A-prime indices (into FB) */
    int *a_indices = malloc(num_a_primes * sizeof(int));

    while (smooth_cnt + combined_cnt < needed && polys < 2000000) {
        /* Pick num_a_primes distinct FB indices from [a_start, fb_count) */
        for (int i = 0; i < num_a_primes; i++) {
            int idx;
            int retry = 0;
            do {
                idx = a_start + (int)(gmp_urandomm_ui(rs, a_range));
                retry = 0;
                for (int j = 0; j < i; j++)
                    if (a_indices[j] == idx) { retry = 1; break; }
            } while (retry);
            a_indices[i] = idx;
        }

        /* Compute A = product of selected primes */
        mpz_set_ui(A, 1);
        for (int i = 0; i < num_a_primes; i++)
            mpz_mul_ui(A, A, FB[a_indices[i]].p);

        /* Find B with B^2 ≡ N (mod A) using CRT */
        /* For each prime q = FB[a_indices[i]].p:
         * B ≡ sqrt(N) mod q (we have this already as FB[].sqrt_n)
         * Use CRT to combine */

        /* Start with B mod q0, then CRT with each subsequent qi */
        unsigned long b_mod_q0 = FB[a_indices[0]].sqrt_n;
        mpz_set_ui(B, b_mod_q0);
        mpz_set_ui(tmp1, FB[a_indices[0]].p); /* modulus so far */

        int crt_ok = 1;
        for (int i = 1; i < num_a_primes && crt_ok; i++) {
            unsigned long qi = FB[a_indices[i]].p;
            unsigned long si = FB[a_indices[i]].sqrt_n;

            /* B ≡ b_curr (mod mod_so_far)
             * B ≡ si (mod qi)
             * CRT: B = b_curr + mod_so_far * t where t ≡ (si - b_curr) * mod_so_far^{-1} (mod qi) */

            unsigned long b_curr_mod_qi = mpz_fdiv_ui(B, qi);
            long diff = (long)si - (long)b_curr_mod_qi;
            if (diff < 0) diff += qi;

            /* Compute mod_so_far^{-1} mod qi */
            unsigned long msf_mod_qi = mpz_fdiv_ui(tmp1, qi);
            mpz_set_ui(tmp2, msf_mod_qi);
            mpz_t qi_mpz;
            mpz_init_set_ui(qi_mpz, qi);
            if (!mpz_invert(tmp2, tmp2, qi_mpz)) {
                crt_ok = 0;
                mpz_clear(qi_mpz);
                break;
            }
            unsigned long inv = mpz_get_ui(tmp2);
            mpz_clear(qi_mpz);

            unsigned long t = ((unsigned long)diff * inv) % qi;

            /* B = B + tmp1 * t */
            mpz_mul_ui(tmp2, tmp1, t);
            mpz_add(B, B, tmp2);

            /* tmp1 = tmp1 * qi */
            mpz_mul_ui(tmp1, tmp1, qi);

            /* Keep B in [0, tmp1) */
            mpz_mod(B, B, tmp1);
        }

        if (!crt_ok) continue;

        /* Verify B^2 ≡ N (mod A) */
        mpz_mul(tmp1, B, B);
        mpz_sub(tmp1, tmp1, N);
        if (!mpz_divisible_p(tmp1, A)) {
            /* Try negating B for some primes - simplest: try B = A - B */
            mpz_sub(B, A, B);
            mpz_mul(tmp1, B, B);
            mpz_sub(tmp1, tmp1, N);
            if (!mpz_divisible_p(tmp1, A)) {
                /* Can happen due to CRT sign choices. Skip this A. */
                polys++;
                continue;
            }
        }

        /* C = (B^2 - N) / A */
        mpz_mul(C, B, B);
        mpz_sub(C, C, N);
        mpz_divexact(C, C, A);

        /* Compute sieve starting positions:
         * g(x) = Ax^2 + 2Bx + C, solve g(x) ≡ 0 (mod p) for each FB prime p
         * x ≡ (-B ± sqrt(N)) * A^{-1} (mod p)
         * (since B^2 - AC = N, discriminant = 4(B^2 - AC) = 4N)
         */
        for (int i = 0; i < fb_count; i++) {
            unsigned int p = FB[i].p;

            /* Check if p divides A */
            int p_in_A = 0;
            for (int j = 0; j < num_a_primes; j++)
                if (FB[a_indices[j]].p == p) { p_in_A = 1; break; }

            if (p_in_A) {
                /* g(x) ≡ 2Bx + C (mod p) since A ≡ 0 (mod p)
                 * x ≡ -C * (2B)^{-1} (mod p) - but need to handle carefully
                 * Actually if p | A, then g(x) = (A/p)*px^2 + 2Bx + C
                 * Since p | A, p | (B^2-N)/... this gets complicated
                 * Simplest: just skip this prime in sieving but add 1 to its exponent */
                soln1[i] = -1;
                soln2[i] = -1;
                continue;
            }

            if (p == 2) {
                /* g(x) mod 2 = C mod 2 (since A even => not possible, A is product of odd FB primes unless 2 is one) */
                /* Since FB[0] = 2 and we skip if p_in_A, if 2 is not in A: */
                unsigned long c2 = mpz_fdiv_ui(C, 2);
                unsigned long a2 = mpz_fdiv_ui(A, 2);
                if (a2 == 0) {
                    /* A even (shouldn't happen for product of odd primes) */
                    soln1[i] = (c2 == 0) ? 0 : -1;
                    soln2[i] = -1;
                } else {
                    /* Ax^2 + C ≡ 0 mod 2 => x^2 + C ≡ 0 mod 2 */
                    soln1[i] = (c2 == 0) ? 0 : 1;
                    soln2[i] = -1; /* only one solution mod 2 */
                }
                continue;
            }

            /* Compute A^{-1} mod p */
            unsigned long a_mod_p = mpz_fdiv_ui(A, p);
            mpz_set_ui(tmp1, a_mod_p);
            mpz_set_ui(tmp2, p);
            mpz_invert(tmp1, tmp1, tmp2);
            unsigned long ainv = mpz_get_ui(tmp1);

            unsigned long b_mod_p = mpz_fdiv_ui(B, p);
            unsigned long sr = FB[i].sqrt_n;

            /* x1 = (-B + sr) * A^{-1} mod p */
            long x1 = (long)((p - b_mod_p + sr) % p);
            x1 = (long)(((unsigned long)x1 * ainv) % p);

            /* x2 = (-B - sr) * A^{-1} mod p */
            long x2 = (long)((2*p - b_mod_p - sr) % p);
            x2 = (long)(((unsigned long)x2 * ainv) % p);

            /* Convert to sieve coordinates: sieve[j] for j in [0, sieve_len)
             * represents x = j - M. So we want j = x + M.
             * x ≡ x1 (mod p) => j = ((x1 + M) % p) */
            soln1[i] = (int)(((long)x1 + M) % (long)p);
            if (soln1[i] < 0) soln1[i] += p;
            soln2[i] = (int)(((long)x2 + M) % (long)p);
            if (soln2[i] < 0) soln2[i] += p;
        }

        /* Sieve */
        memset(sieve, 0, sieve_len);
        for (int i = 0; i < fb_count; i++) {
            if (soln1[i] < 0) continue;
            unsigned int p = FB[i].p;
            int lp = FB[i].logp;
            for (int j = soln1[i]; j < sieve_len; j += p) sieve[j] += lp;
            if (soln2[i] >= 0 && soln2[i] != soln1[i])
                for (int j = soln2[i]; j < sieve_len; j += p) sieve[j] += lp;
        }

        /* Scan candidates */
        for (int j = 0; j < sieve_len; j++) {
            if (sieve[j] < thresh) continue;
            if (rel_count >= MAX_RELS - 10) break;

            long x = (long)j - M;

            /* f(x) = Ax^2 + 2Bx + C */
            mpz_set_si(fval, x);
            mpz_mul_si(tmp1, fval, x);   /* x^2 */
            mpz_mul(tmp1, tmp1, A);       /* Ax^2 */
            mpz_mul(tmp2, B, fval);
            mpz_mul_ui(tmp2, tmp2, 2);    /* 2Bx */
            mpz_add(fval, tmp1, tmp2);
            mpz_add(fval, fval, C);       /* + C */

            /* sqrt_val = Ax + B */
            mpz_set_si(sqv, x);
            mpz_mul(sqv, sqv, A);
            mpz_add(sqv, sqv, B);

            int neg = (mpz_sgn(fval) < 0);
            if (neg) mpz_neg(fval, fval);

            /* Trial divide fval by FB */
            mpz_set(cofactor, fval);
            memset(expo, 0, fb_count * sizeof(int));

            /* Include A's prime factors with exponent +1 each
             * (since (Ax+B)^2 ≡ A*f(x) (mod N), the actual smooth value is A*f(x)) */
            for (int k = 0; k < num_a_primes; k++)
                expo[a_indices[k]] += 1;

            for (int k = 0; k < fb_count; k++) {
                unsigned int p = FB[k].p;
                while (mpz_divisible_ui_p(cofactor, p)) {
                    mpz_divexact_ui(cofactor, cofactor, p);
                    expo[k]++;
                }
            }

            /* Check smoothness */
            if (mpz_cmp_ui(cofactor, 1) == 0) {
                /* Fully smooth */
                store_rel(sqv, expo, neg, 0, 0);
                smooth_cnt++;
            } else if (mpz_fits_ulong_p(cofactor)) {
                unsigned long cf = mpz_get_ui(cofactor);
                if (cf <= lp_bound) {
                    /* Single large prime */
                    int match = ht_find(cf);
                    if (match >= 0) {
                        /* Combine partials */
                        store_rel(sqv, expo, neg, cf, 0);
                        int i1 = match, i2 = rel_count - 1;

                        /* Create combined relation */
                        rel_t *r = &RELS[rel_count];
                        mpz_init(r->sqrt_val);
                        mpz_mul(r->sqrt_val, RELS[i1].sqrt_val, RELS[i2].sqrt_val);
                        mpz_mod(r->sqrt_val, r->sqrt_val, N);
                        r->expo = malloc(fb_count * sizeof(int));
                        for (int k = 0; k < fb_count; k++)
                            r->expo[k] = RELS[i1].expo[k] + RELS[i2].expo[k];
                        r->neg = RELS[i1].neg ^ RELS[i2].neg;
                        r->lp1 = 0;
                        r->lp2 = cf; /* store for Y computation */
                        rel_count++;
                        combined_cnt++;
                    } else {
                        store_rel(sqv, expo, neg, cf, 0);
                        ht_add(cf, rel_count - 1);
                        partial_cnt++;
                    }
                }
            }
        }

        polys++;
        if (polys % 200 == 0) {
            clock_gettime(CLOCK_MONOTONIC, &t1);
            double el = (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)/1e9;
            fprintf(stderr, "Poly %d: smooth=%d combined=%d partial=%d need=%d (%.1fs)\n",
                    polys, smooth_cnt, combined_cnt, partial_cnt, needed, el);
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double sieve_t = (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)/1e9;
    fprintf(stderr, "Sieve done: smooth=%d combined=%d (%.1fs, %d polys)\n",
            smooth_cnt, combined_cnt, sieve_t, polys);

    /* Collect usable relations (lp1 == 0: either smooth or combined) */
    int *usable = malloc(rel_count * sizeof(int));
    int nusable = 0;
    for (int i = 0; i < rel_count; i++)
        if (RELS[i].lp1 == 0) usable[nusable++] = i;

    fprintf(stderr, "Usable: %d (need %d)\n", nusable, fb_count + 1);

    if (nusable < fb_count + 1) {
        printf("FAIL\n");
        fprintf(stderr, "Not enough relations\n");
        return 1;
    }

    /* Gaussian elimination */
    int nrows = nusable;
    int ncols = fb_count + 1; /* +1 for sign */

    unsigned char **mat = malloc(nrows * sizeof(unsigned char *));
    for (int i = 0; i < nrows; i++) {
        mat[i] = calloc(ncols + nrows, 1);
        int ri = usable[i];
        mat[i][0] = RELS[ri].neg & 1;
        for (int j = 0; j < fb_count; j++)
            mat[i][j+1] = RELS[ri].expo[j] & 1;
        mat[i][ncols + i] = 1;
    }

    int rank = gauss(mat, nrows, ncols);
    fprintf(stderr, "Matrix %dx%d, rank=%d, deps=%d\n", nrows, ncols, rank, nrows-rank);

    /* Try null-space vectors */
    int found = 0;
    mpz_t X, Y, g, pm;
    mpz_init(X); mpz_init(Y); mpz_init(g); mpz_init(pm);

    for (int row = rank; row < nrows && !found; row++) {
        mpz_set_ui(X, 1);
        int *total = calloc(fb_count, sizeof(int));

        for (int i = 0; i < nrows; i++) {
            if (!mat[row][ncols + i]) continue;
            int ri = usable[i];
            mpz_mul(X, X, RELS[ri].sqrt_val);
            mpz_mod(X, X, N);
            for (int j = 0; j < fb_count; j++)
                total[j] += RELS[ri].expo[j];
        }

        /* Y = product of FB[j].p ^ (total[j]/2) mod N */
        mpz_set_ui(Y, 1);
        for (int j = 0; j < fb_count; j++) {
            if (total[j] <= 0) continue;
            mpz_set_ui(pm, FB[j].p);
            mpz_powm_ui(pm, pm, total[j] / 2, N);
            mpz_mul(Y, Y, pm);
            mpz_mod(Y, Y, N);
        }

        /* Include large primes from combined relations */
        for (int i = 0; i < nrows; i++) {
            if (!mat[row][ncols + i]) continue;
            int ri = usable[i];
            if (RELS[ri].lp2 > 0) {
                mpz_set_ui(pm, RELS[ri].lp2);
                mpz_mul(Y, Y, pm);
                mpz_mod(Y, Y, N);
            }
        }

        mpz_sub(g, X, Y); mpz_gcd(g, g, N);
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) found = 1;
        else {
            mpz_add(g, X, Y); mpz_gcd(g, g, N);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) found = 1;
        }
        free(total);
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double elapsed = (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)/1e9;

    if (found) {
        mpz_t cof; mpz_init(cof); mpz_divexact(cof, N, g);
        gmp_printf("%Zd %Zd\n", g, cof);
        fprintf(stderr, "MPQS time=%.3fs\n", elapsed);
        return 0;
    }
    fprintf(stderr, "MPQS linalg failed, time=%.3fs\n", elapsed);
    printf("FAIL\n");
    return 1;
}
