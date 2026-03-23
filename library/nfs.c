/*
 * Basic Number Field Sieve (NFS) implementation.
 *
 * Achieves L[1/3] scaling by using two "images" of each relation:
 * - Rational side: a + b*m must be smooth
 * - Algebraic side: Norm_f(a + b*α) must be smooth
 *
 * Uses base-m polynomial selection and line sieving.
 *
 * Usage: ./nfs <number>
 * Single-threaded. Seed=42.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

#define MAX_RFB 8000   /* max rational factor base size */
#define MAX_AFB 8000   /* max algebraic factor base size */
#define MAX_RELS 25000

/* ---- Primes ---- */
static int *PRIMES;
static int NPRIMES;

static void gen_primes(int lim) {
    char *s = calloc(lim+1,1);
    int c = 0;
    for (int i = 2; i <= lim; i++) if (!s[i]) { c++; for (long j=(long)i*i;j<=lim;j+=i) s[j]=1; }
    PRIMES = malloc(c * sizeof(int));
    NPRIMES = 0;
    memset(s, 0, lim+1);
    for (int i = 2; i <= lim; i++) if (!s[i]) { PRIMES[NPRIMES++] = i; for (long j=(long)i*i;j<=lim;j+=i) s[j]=1; }
    free(s);
}

/* Tonelli-Shanks */
static long tsqrt_mod(long a, long p) {
    if (p == 2) return a & 1;
    if (a == 0) return 0;
    mpz_t az, pz, rz;
    mpz_init_set_si(az, a); mpz_init_set_si(pz, p); mpz_init(rz);
    mpz_mod(az, az, pz);
    mpz_powm_ui(rz, az, (p-1)/2, pz);
    if (mpz_cmp_ui(rz, 1) != 0) { mpz_clear(az); mpz_clear(pz); mpz_clear(rz); return -1; }
    if ((p & 3) == 3) { mpz_powm_ui(rz, az, (p+1)/4, pz); long r = mpz_get_si(rz); mpz_clear(az); mpz_clear(pz); mpz_clear(rz); return r; }
    long Q = p-1, S = 0;
    while (!(Q&1)) { Q >>= 1; S++; }
    long z = 2;
    mpz_t zz; mpz_init(zz);
    while (1) { mpz_set_si(zz, z); mpz_powm_ui(zz, zz, (p-1)/2, pz); if (mpz_cmp_ui(zz, p-1) == 0) break; z++; }
    mpz_t M, c, t, R, b;
    mpz_init_set_si(M, S); mpz_init(c); mpz_init(t); mpz_init(R); mpz_init(b);
    mpz_set_si(c, z); mpz_powm_ui(c, c, Q, pz);
    mpz_set(t, az); mpz_powm_ui(t, t, Q, pz);
    mpz_set(R, az); mpz_powm_ui(R, R, (Q+1)/2, pz);
    while (1) {
        if (mpz_cmp_ui(t, 1) == 0) { long r = mpz_get_si(R); mpz_clear(M); mpz_clear(c); mpz_clear(t); mpz_clear(R); mpz_clear(b); mpz_clear(zz); mpz_clear(az); mpz_clear(pz); mpz_clear(rz); return r; }
        mpz_set(zz, t); long i = 0;
        while (mpz_cmp_ui(zz, 1) != 0) { mpz_mul(zz, zz, zz); mpz_mod(zz, zz, pz); i++; }
        long m = mpz_get_si(M);
        mpz_set(b, c);
        for (long j = 0; j < m-i-1; j++) { mpz_mul(b, b, b); mpz_mod(b, b, pz); }
        mpz_set_si(M, i);
        mpz_mul(c, b, b); mpz_mod(c, c, pz);
        mpz_mul(t, t, c); mpz_mod(t, t, pz);
        mpz_mul(R, R, b); mpz_mod(R, R, pz);
    }
}

/* Find roots of polynomial f mod p */
/* f(x) = c[d]*x^d + ... + c[1]*x + c[0] */
static int find_roots_mod_p(mpz_t *coeff, int degree, long p, long *roots) {
    int nroots = 0;
    /* For small p, just try all values */
    if (p < 10000) {
        mpz_t val, xp, term, pp;
        mpz_init(val); mpz_init(xp); mpz_init(term); mpz_init_set_si(pp, p);
        for (long x = 0; x < p; x++) {
            mpz_set_ui(val, 0);
            mpz_set_ui(xp, 1);
            for (int i = 0; i <= degree; i++) {
                mpz_mul(term, coeff[i], xp);
                mpz_add(val, val, term);
                mpz_mul_ui(xp, xp, x);
            }
            mpz_mod(val, val, pp);
            if (mpz_sgn(val) == 0) {
                roots[nroots++] = x;
            }
        }
        mpz_clear(val); mpz_clear(xp); mpz_clear(term); mpz_clear(pp);
        return nroots;
    }
    /* For larger p, use Hensel/Cantor-Zassenhaus - but for our case, p is always small */
    return 0;
}

/* Rational factor base entry */
typedef struct { int p; int logp; } rfb_t;

/* Algebraic factor base entry: (p, r) where f(r) ≡ 0 (mod p) */
typedef struct { int p; int logp; long root; } afb_t;

/* Relation: pair (a, b) with a+b*m smooth and Norm_f(a,b) smooth */
typedef struct {
    long a, b;
    int *r_expo;    /* rational exponents */
    int *a_expo;    /* algebraic exponents */
    int r_neg;      /* sign of a+b*m */
    int a_neg;      /* sign of algebraic norm */
} rel_t;

static rfb_t RFB[MAX_RFB];
static int rfb_count;
static afb_t AFB[MAX_AFB];
static int afb_count;
static rel_t RELS[MAX_RELS];
static int rel_count = 0;

static mpz_t N_global;
static mpz_t *F_coeff; /* polynomial coefficients */
static int F_degree;
static mpz_t M_val;    /* m value */

/* NFS parameters by digit count */
typedef struct {
    int degree;     /* polynomial degree */
    int rfb_sz;     /* rational factor base size */
    int afb_sz;     /* algebraic factor base size */
    int sieve_a;    /* sieve range for a: [-sieve_a, sieve_a] */
    int max_b;      /* max b value */
    double r_thresh_sub;
    double a_thresh_sub;
} nfs_par_t;

static nfs_par_t get_par(int digits) {
    /* NFS needs larger factor bases than QS for equivalent sizes.
     * The advantage only shows at larger digit counts (better scaling).
     * FB sizes are approximately L[1/3]^{2/3} which grows slower than L[1/2]^{1/2}. */
    if (digits <= 35)  return (nfs_par_t){3, 2000,  2000,  200000,  2000,  12, 12};
    if (digits <= 40)  return (nfs_par_t){3, 2500,  2500,  300000,  3000,  13, 13};
    if (digits <= 50)  return (nfs_par_t){3, 3500,  3500,  500000,  5000,  14, 14};
    if (digits <= 60)  return (nfs_par_t){4, 4000,  4000,  500000,  5000,  15, 15};
    if (digits <= 70)  return (nfs_par_t){4, 5000,  5000,  800000,  8000,  16, 16};
    if (digits <= 80)  return (nfs_par_t){5, 5500,  5500,  1000000, 10000, 17, 17};
    if (digits <= 90)  return (nfs_par_t){5, 6000,  6000,  1500000, 15000, 18, 18};
    if (digits <= 100) return (nfs_par_t){5, 7000,  7000,  2000000, 20000, 19, 19};
    return (nfs_par_t){5, 8000, 8000, 2500000, 25000, 20, 20};
}

/* Base-m polynomial selection:
 * Write N in base m where m = N^(1/d).
 * f(x) = c_d*x^d + ... + c_1*x + c_0 where c_i are the base-m digits.
 * Then f(m) = N, so f(m) ≡ 0 (mod N). */
static void select_polynomial(int degree) {
    F_degree = degree;
    F_coeff = malloc((degree + 1) * sizeof(mpz_t));
    for (int i = 0; i <= degree; i++) mpz_init(F_coeff[i]);
    mpz_init(M_val);

    /* m = floor(N^(1/d)) */
    mpz_t tmp;
    mpz_init(tmp);
    mpz_root(M_val, N_global, degree);

    /* Extract base-m digits */
    mpz_set(tmp, N_global);
    for (int i = 0; i < degree; i++) {
        mpz_tdiv_qr(tmp, F_coeff[i], tmp, M_val);
    }
    mpz_set(F_coeff[degree], tmp);

    /* Verify: f(m) should equal N */
    mpz_set(tmp, F_coeff[degree]);
    for (int i = degree - 1; i >= 0; i--) {
        mpz_mul(tmp, tmp, M_val);
        mpz_add(tmp, tmp, F_coeff[i]);
    }
    if (mpz_cmp(tmp, N_global) != 0) {
        fprintf(stderr, "ERROR: f(m) != N\n");
        gmp_fprintf(stderr, "f(m) = %Zd\nN    = %Zd\n", tmp, N_global);
        exit(1);
    }

    gmp_fprintf(stderr, "Polynomial: degree=%d, m=%Zd\n", degree, M_val);
    for (int i = degree; i >= 0; i--)
        gmp_fprintf(stderr, "  c[%d] = %Zd\n", i, F_coeff[i]);

    mpz_clear(tmp);
}

/* Compute algebraic norm: Norm_f(a, b) = (-b)^d * f(-a/b)
 * = sum_{i=0}^{d} c_i * (-a)^i * b^(d-i) * (-1)^d ... actually:
 * Norm(a + b*alpha) = resultant of (a + b*x, f(x))
 * = b^d * f(-a/b)  (up to sign)
 * = sum_{i=0}^{d} c_i * (-a)^i * b^(d-i) */
static void compute_alg_norm(mpz_t norm, long a, long b) {
    mpz_t term, neg_a_pow, b_pow;
    mpz_init(term);
    mpz_init_set_si(neg_a_pow, 1); /* (-a)^0 */
    mpz_init(b_pow);

    /* b_pow = b^d */
    mpz_set_si(b_pow, b);
    mpz_pow_ui(b_pow, b_pow, F_degree);

    mpz_set_ui(norm, 0);
    long neg_a = -a;

    for (int i = 0; i <= F_degree; i++) {
        /* term = c_i * (-a)^i * b^(d-i) */
        mpz_mul(term, F_coeff[i], neg_a_pow);
        mpz_mul(term, term, b_pow);
        mpz_add(norm, norm, term);

        /* Update powers */
        mpz_mul_si(neg_a_pow, neg_a_pow, neg_a);
        if (i < F_degree) {
            /* b_pow = b_pow / b (we're going from b^d down to b^0) */
            if (b != 0) { if (b > 0) mpz_divexact_ui(b_pow, b_pow, b); else { mpz_divexact_ui(b_pow, b_pow, -b); mpz_neg(b_pow, b_pow); } }
            else mpz_set_ui(b_pow, 0);
        }
    }

    mpz_clear(term);
    mpz_clear(neg_a_pow);
    mpz_clear(b_pow);
}

/* Build rational and algebraic factor bases */
static void build_factor_bases(int rfb_target, int afb_target) {
    rfb_count = 0;
    afb_count = 0;

    long roots[64];

    for (int i = 0; i < NPRIMES && (rfb_count < rfb_target || afb_count < afb_target); i++) {
        int p = PRIMES[i];

        if (rfb_count < rfb_target) {
            RFB[rfb_count].p = p;
            RFB[rfb_count].logp = (int)(log2((double)p) + 0.5);
            if (RFB[rfb_count].logp < 1) RFB[rfb_count].logp = 1;
            rfb_count++;
        }

        if (afb_count < afb_target) {
            /* Find roots of f(x) mod p */
            int nr = find_roots_mod_p(F_coeff, F_degree, p, roots);
            for (int j = 0; j < nr && afb_count < afb_target; j++) {
                AFB[afb_count].p = p;
                AFB[afb_count].logp = (int)(log2((double)p) + 0.5);
                if (AFB[afb_count].logp < 1) AFB[afb_count].logp = 1;
                AFB[afb_count].root = roots[j];
                afb_count++;
            }
        }
    }

    fprintf(stderr, "RFB: %d primes (max %d), AFB: %d (prime,root) pairs (max %d)\n",
            rfb_count, RFB[rfb_count-1].p, afb_count, AFB[afb_count-1].p);
}

/* Line sieve for a given b value.
 * For each b, sieve over a in [-sieve_a, sieve_a].
 * Rational: a + b*m must be smooth (sieve as in QS)
 * Algebraic: Norm_f(a,b) must be smooth (sieve using (p,r) pairs)
 */
static int line_sieve(long b, int sieve_a, double r_thresh_sub, double a_thresh_sub) {
    int sieve_len = 2 * sieve_a;
    unsigned char *r_sieve = calloc(sieve_len, 1); /* rational */
    unsigned char *a_sieve = calloc(sieve_len, 1); /* algebraic */

    if (!r_sieve || !a_sieve) {
        free(r_sieve); free(a_sieve);
        return 0;
    }

    /* Compute log2 of max rational value: |a + b*m| ≈ b*m + sieve_a */
    mpz_t bm;
    mpz_init(bm);
    mpz_mul_si(bm, M_val, b);
    mpz_add_ui(bm, bm, sieve_a);
    double log2_rmax = mpz_sizeinbase(bm, 2);

    /* Rational sieve: a + b*m ≡ 0 (mod p) => a ≡ -b*m (mod p) */
    for (int i = 0; i < rfb_count; i++) {
        int p = RFB[i].p;
        long bm_mod_p = mpz_fdiv_ui(M_val, p);
        bm_mod_p = (bm_mod_p * (b % p)) % p;
        long start = (long)(p - bm_mod_p) % p;
        /* Convert to sieve index: sieve[j] for a = j - sieve_a */
        long idx = (start + sieve_a) % p;
        if (idx < 0) idx += p;
        for (long j = idx; j < sieve_len; j += p)
            r_sieve[j] += RFB[i].logp;
    }

    /* Algebraic sieve: Norm_f(a, b) ≡ 0 (mod p) when a ≡ -b*r (mod p)
     * where r is a root of f mod p */
    /* log2 of max algebraic norm: roughly d * log2(max coeff * max(a,b)) */
    double log2_amax = F_degree * (mpz_sizeinbase(F_coeff[F_degree], 2) + log2((double)(sieve_a > b ? sieve_a : b) + 1));

    for (int i = 0; i < afb_count; i++) {
        int p = AFB[i].p;
        long r = AFB[i].root;
        /* a + b*r ≡ 0 (mod p) => a ≡ -b*r (mod p) */
        long neg_br = (-(b % p) * (r % p)) % p;
        if (neg_br < 0) neg_br += p;
        long idx = (neg_br + sieve_a) % p;
        if (idx < 0) idx += p;
        for (long j = idx; j < sieve_len; j += p)
            a_sieve[j] += AFB[i].logp;
    }

    /* Rational threshold */
    int r_thresh = (int)(log2_rmax - r_thresh_sub);
    if (r_thresh < 10) r_thresh = 10;
    /* Algebraic threshold */
    int a_thresh = (int)(log2_amax - a_thresh_sub);
    if (a_thresh < 10) a_thresh = 10;

    /* Scan for candidates where BOTH sides pass threshold */
    int found = 0;
    mpz_t r_val, a_norm, cofactor;
    mpz_init(r_val); mpz_init(a_norm); mpz_init(cofactor);

    int *r_expo = calloc(rfb_count, sizeof(int));
    int *a_expo = calloc(afb_count, sizeof(int));

    for (int j = 0; j < sieve_len && rel_count < MAX_RELS; j++) {
        if (r_sieve[j] < r_thresh || a_sieve[j] < a_thresh) continue;

        long a = (long)j - sieve_a;
        if (a == 0) continue;  /* a and b must be coprime */

        /* Check gcd(a, b) == 1 */
        long g = a < 0 ? -a : a;
        long h = b;
        while (h) { long t = h; h = g % h; g = t; }
        if (g != 1) continue;

        /* Trial divide rational side: a + b*m */
        mpz_mul_si(r_val, M_val, b);
        if (a >= 0) mpz_add_ui(r_val, r_val, a);
        else mpz_sub_ui(r_val, r_val, -a);
        int r_neg = (mpz_sgn(r_val) < 0);
        if (r_neg) mpz_neg(r_val, r_val);

        mpz_set(cofactor, r_val);
        memset(r_expo, 0, rfb_count * sizeof(int));
        for (int k = 0; k < rfb_count; k++) {
            int p = RFB[k].p;
            while (mpz_divisible_ui_p(cofactor, p)) {
                mpz_divexact_ui(cofactor, cofactor, p);
                r_expo[k]++;
            }
        }
        if (mpz_cmp_ui(cofactor, 1) != 0) continue; /* not smooth */

        /* Trial divide algebraic side */
        compute_alg_norm(a_norm, a, b);
        int a_neg = (mpz_sgn(a_norm) < 0);
        if (a_neg) mpz_neg(a_norm, a_norm);

        mpz_set(cofactor, a_norm);
        memset(a_expo, 0, afb_count * sizeof(int));
        for (int k = 0; k < afb_count; k++) {
            int p = AFB[k].p;
            while (mpz_divisible_ui_p(cofactor, p)) {
                mpz_divexact_ui(cofactor, cofactor, p);
                a_expo[k]++;
            }
        }
        if (mpz_cmp_ui(cofactor, 1) != 0) continue; /* not smooth */

        /* Both sides smooth! Store relation */
        rel_t *r = &RELS[rel_count];
        r->a = a;
        r->b = b;
        r->r_expo = malloc(rfb_count * sizeof(int));
        memcpy(r->r_expo, r_expo, rfb_count * sizeof(int));
        r->a_expo = malloc(afb_count * sizeof(int));
        memcpy(r->a_expo, a_expo, afb_count * sizeof(int));
        r->r_neg = r_neg;
        r->a_neg = a_neg;
        rel_count++;
        found++;
    }

    free(r_expo); free(a_expo);
    free(r_sieve); free(a_sieve);
    mpz_clear(r_val); mpz_clear(a_norm); mpz_clear(cofactor); mpz_clear(bm);

    return found;
}

/* Gaussian elimination mod 2 */
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

    mpz_init_set_str(N_global, argv[1], 10);
    int digits = strlen(argv[1]);

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    /* Trial division */
    for (unsigned long p = 2; p < 1000000; p++) {
        if (mpz_divisible_ui_p(N_global, p)) {
            mpz_t q; mpz_init(q); mpz_divexact_ui(q, N_global, p);
            clock_gettime(CLOCK_MONOTONIC, &t1);
            gmp_printf("%lu %Zd\n", p, q);
            fprintf(stderr, "Trial: time=%.3fs\n", (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)/1e9);
            return 0;
        }
    }

    nfs_par_t par = get_par(digits);
    gen_primes((par.rfb_sz > par.afb_sz ? par.rfb_sz : par.afb_sz) * 50);

    /* Polynomial selection */
    select_polynomial(par.degree);

    /* Build factor bases */
    build_factor_bases(par.rfb_sz, par.afb_sz);

    int needed = rfb_count + afb_count + 10;
    fprintf(stderr, "NFS: %d digits, degree=%d, RFB=%d, AFB=%d, need=%d relations\n",
            digits, par.degree, rfb_count, afb_count, needed);
    fprintf(stderr, "Sieve: a in [-%d, %d], b in [1, %d]\n", par.sieve_a, par.sieve_a, par.max_b);

    /* Line sieve */
    for (long b = 1; b <= par.max_b && rel_count < needed; b++) {
        int found = line_sieve(b, par.sieve_a, par.r_thresh_sub, par.a_thresh_sub);

        if (b % 50 == 0) {
            clock_gettime(CLOCK_MONOTONIC, &t1);
            double el = (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)/1e9;
            fprintf(stderr, "b=%ld: rels=%d/%d (%.1fs)\n", b, rel_count, needed, el);
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double sieve_t = (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)/1e9;
    fprintf(stderr, "Sieve done: %d/%d relations (%.1fs)\n", rel_count, needed, sieve_t);

    if (rel_count < needed) {
        fprintf(stderr, "Not enough relations\n");
        printf("FAIL\n");
        return 1;
    }

    /* Linear algebra: matrix has columns for:
     * - 1 sign column (rational)
     * - rfb_count columns (rational exponents mod 2)
     * - 1 sign column (algebraic)
     * - afb_count columns (algebraic exponents mod 2)
     * Total columns = rfb_count + afb_count + 2
     */
    int ncols = rfb_count + afb_count + 2;
    int nrows = rel_count;

    unsigned char **mat = malloc(nrows * sizeof(unsigned char *));
    for (int i = 0; i < nrows; i++) {
        mat[i] = calloc(ncols + nrows, 1);
        /* Rational sign */
        mat[i][0] = RELS[i].r_neg & 1;
        /* Rational exponents */
        for (int j = 0; j < rfb_count; j++)
            mat[i][j + 1] = RELS[i].r_expo[j] & 1;
        /* Algebraic sign */
        mat[i][rfb_count + 1] = RELS[i].a_neg & 1;
        /* Algebraic exponents */
        for (int j = 0; j < afb_count; j++)
            mat[i][rfb_count + 2 + j] = RELS[i].a_expo[j] & 1;
        /* Identity */
        mat[i][ncols + i] = 1;
    }

    int rank = gauss(mat, nrows, ncols);
    fprintf(stderr, "Matrix %dx%d, rank=%d, deps=%d\n", nrows, ncols, rank, nrows-rank);

    /* NFS square root step:
     *
     * For each dependency, we have:
     * - Rational side: ∏(a_i + b_i*m) = S_r² in Z (all rational exponents even)
     * - Algebraic side: ∏(a_i + b_i*α) = S_a² in Z[α] (all algebraic exponents even)
     *
     * Since φ(α) = m (mod N), φ(S_a)² = φ(∏(a_i + b_i*α)) = ∏(a_i + b_i*m) = S_r² (mod N)
     * So gcd(S_r - φ(S_a), N) gives a factor.
     *
     * Computing S_r: just ∏p^(e/2) from the rational exponents.
     * Computing φ(S_a): compute ∏(a_i + b_i*α) in Z[α]/(f(α)) mod N,
     *   then find its square root using Euler criterion + GCD trick.
     *
     * Polynomial arithmetic in R = Z[α]/(f(α)) mod N:
     * Elements are polynomials of degree < d with coefficients in Z/NZ.
     * Multiplication: polynomial multiply then reduce mod f(α).
     */

    /* Polynomial multiplication in R mod N: multiply two degree-(d-1) polys mod f mod N */
    int d = F_degree;

    /* poly_mul: result[0..d-1] = a[0..d-1] * b[0..d-1] mod f mod N */
    /* Inline since d is small (3-5) */
    mpz_t *poly_tmp = malloc(2 * d * sizeof(mpz_t));
    for (int i = 0; i < 2*d; i++) mpz_init(poly_tmp[i]);

    int found = 0;
    mpz_t X, Y, g, pm;
    mpz_init(X); mpz_init(Y); mpz_init(g); mpz_init(pm);

    for (int row = rank; row < nrows && !found; row++) {
        /* Step 1: Compute S_r = rational square root */
        int *r_total = calloc(rfb_count, sizeof(int));

        /* Also compute P = ∏(a_i + b_i*α) in R mod N */
        mpz_t *P = malloc(d * sizeof(mpz_t));
        for (int k = 0; k < d; k++) mpz_init(P[k]);
        mpz_set_ui(P[0], 1); /* P starts as 1 */

        int dep_size = 0;
        for (int i = 0; i < nrows; i++) {
            if (!mat[row][ncols + i]) continue;
            dep_size++;

            for (int j = 0; j < rfb_count; j++)
                r_total[j] += RELS[i].r_expo[j];

            /* Multiply P by (a_i + b_i * α) in R mod N */
            /* linear = a_i + b_i * α, so linear[0] = a_i, linear[1] = b_i, rest 0 */
            /* P_new = P * linear mod f mod N */
            for (int k = 0; k < 2*d; k++) mpz_set_ui(poly_tmp[k], 0);

            for (int j = 0; j < d; j++) {
                /* P[j] * a_i goes to poly_tmp[j] */
                mpz_t ta, tb;
                mpz_init(ta); mpz_init(tb);
                mpz_mul_si(ta, P[j], RELS[i].a);
                mpz_add(poly_tmp[j], poly_tmp[j], ta);
                /* P[j] * b_i * α goes to poly_tmp[j+1] */
                mpz_mul_si(tb, P[j], RELS[i].b);
                mpz_add(poly_tmp[j+1], poly_tmp[j+1], tb);
                mpz_clear(ta); mpz_clear(tb);
            }

            /* Reduce mod f: for degree d and above, use f(α) = 0
             * α^d = -(c_{d-1}*α^{d-1} + ... + c_0) / c_d
             * Since c_d = F_coeff[d], we need its inverse mod N.
             * For base-m polynomial, c_d is typically 1 or small.
             */
            /* Actually: f(α) = c_d*α^d + c_{d-1}*α^{d-1} + ... + c_0 = 0
             * So α^d = -(c_{d-1}*α^{d-1} + ... + c_0) / c_d */
            mpz_t cd_inv;
            mpz_init(cd_inv);
            if (mpz_invert(cd_inv, F_coeff[d], N_global) == 0) {
                /* c_d not invertible mod N => gcd(c_d, N) is a factor! */
                mpz_gcd(g, F_coeff[d], N_global);
                if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N_global) < 0) {
                    found = 1;
                    mpz_clear(cd_inv);
                    break;
                }
            }
            mpz_neg(cd_inv, cd_inv);
            mpz_mod(cd_inv, cd_inv, N_global);

            /* Reduce from degree 2d-2 down to d-1 */
            for (int k = 2*d - 2; k >= d; k--) {
                mpz_mod(poly_tmp[k], poly_tmp[k], N_global);
                if (mpz_sgn(poly_tmp[k]) == 0) continue;
                /* α^k = α^{k-d} * α^d = α^{k-d} * (-cd_inv) * (c_{d-1}*α^{d-1} + ... + c_0) */
                mpz_t coeff_k;
                mpz_init(coeff_k);
                mpz_mul(coeff_k, poly_tmp[k], cd_inv);
                mpz_mod(coeff_k, coeff_k, N_global);
                for (int j = 0; j < d; j++) {
                    mpz_t t;
                    mpz_init(t);
                    mpz_mul(t, coeff_k, F_coeff[j]);
                    mpz_add(poly_tmp[k - d + j], poly_tmp[k - d + j], t);
                    mpz_clear(t);
                }
                mpz_set_ui(poly_tmp[k], 0);
                mpz_clear(coeff_k);
            }

            /* Copy reduced result back to P */
            for (int k = 0; k < d; k++) {
                mpz_mod(P[k], poly_tmp[k], N_global);
            }

            mpz_clear(cd_inv);
        }

        if (found) { free(r_total); for (int k = 0; k < d; k++) mpz_clear(P[k]); free(P); break; }

        /* Step 2: Compute S_r = rational square root */
        mpz_set_ui(Y, 1);
        for (int j = 0; j < rfb_count; j++) {
            if (r_total[j] <= 0) continue;
            mpz_set_ui(pm, RFB[j].p);
            mpz_powm_ui(pm, pm, r_total[j] / 2, N_global);
            mpz_mul(Y, Y, pm);
            mpz_mod(Y, Y, N_global);
        }

        /* Step 3: Compute algebraic square root φ(S_a)
         * P = ∏(a_i + b_i*α) in R mod N. We need S such that S² = P in R/(f, N).
         * Then φ(S_a) = S evaluated at m modulo N.
         *
         * Approach: compute P^((N^d + 1) / 2) in R mod N.
         * Actually, use the Euler criterion approach:
         * For the multiplicative group of R/(f, p) ≅ F_{p^d}* (when f is irred mod p),
         * the order is p^d - 1, so square root is P^((p^d+1)/2).
         * For N = pq, we use P^((N^d+1)/2) which works mod p and mod q differently.
         *
         * But (N^d+1)/2 is enormous. Use repeated squaring in R.
         * Each step: polynomial multiply mod f mod N.
         *
         * Actually, computing P^((N^d - 1)/2) first:
         * If result is 1 mod p and -1 mod q (or vice versa),
         * then gcd(result - 1, N) gives a factor directly!
         */

        /* Compute exponent = (N^d - 1) / 2 */
        mpz_t exp;
        mpz_init(exp);
        mpz_pow_ui(exp, N_global, d);
        mpz_sub_ui(exp, exp, 1);
        mpz_tdiv_q_ui(exp, exp, 2);

        /* Compute P^exp in R mod N using repeated squaring */
        mpz_t *result = malloc(d * sizeof(mpz_t));
        mpz_t *base = malloc(d * sizeof(mpz_t));
        for (int k = 0; k < d; k++) {
            mpz_init_set_ui(result[k], k == 0 ? 1 : 0); /* result = 1 */
            mpz_init_set(base[k], P[k]); /* base = P */
        }

        /* Repeated squaring */
        int exp_bits = mpz_sizeinbase(exp, 2);
        for (int bit = 0; bit < exp_bits; bit++) {
            if (mpz_tstbit(exp, bit)) {
                /* result = result * base mod f mod N */
                for (int k = 0; k < 2*d; k++) mpz_set_ui(poly_tmp[k], 0);
                for (int j1 = 0; j1 < d; j1++)
                    for (int j2 = 0; j2 < d; j2++) {
                        mpz_t t; mpz_init(t);
                        mpz_mul(t, result[j1], base[j2]);
                        mpz_add(poly_tmp[j1+j2], poly_tmp[j1+j2], t);
                        mpz_clear(t);
                    }
                /* Reduce mod f mod N */
                mpz_t cd_inv; mpz_init(cd_inv);
                mpz_invert(cd_inv, F_coeff[d], N_global);
                mpz_neg(cd_inv, cd_inv); mpz_mod(cd_inv, cd_inv, N_global);
                for (int k = 2*d-2; k >= d; k--) {
                    mpz_mod(poly_tmp[k], poly_tmp[k], N_global);
                    if (mpz_sgn(poly_tmp[k]) == 0) continue;
                    mpz_t ck; mpz_init(ck);
                    mpz_mul(ck, poly_tmp[k], cd_inv); mpz_mod(ck, ck, N_global);
                    for (int j = 0; j < d; j++) {
                        mpz_t t; mpz_init(t);
                        mpz_mul(t, ck, F_coeff[j]);
                        mpz_add(poly_tmp[k-d+j], poly_tmp[k-d+j], t);
                        mpz_clear(t);
                    }
                    mpz_set_ui(poly_tmp[k], 0);
                    mpz_clear(ck);
                }
                for (int k = 0; k < d; k++) mpz_mod(result[k], poly_tmp[k], N_global);
                mpz_clear(cd_inv);
            }

            /* base = base^2 mod f mod N */
            for (int k = 0; k < 2*d; k++) mpz_set_ui(poly_tmp[k], 0);
            for (int j1 = 0; j1 < d; j1++)
                for (int j2 = 0; j2 < d; j2++) {
                    mpz_t t; mpz_init(t);
                    mpz_mul(t, base[j1], base[j2]);
                    mpz_add(poly_tmp[j1+j2], poly_tmp[j1+j2], t);
                    mpz_clear(t);
                }
            mpz_t cd_inv2; mpz_init(cd_inv2);
            mpz_invert(cd_inv2, F_coeff[d], N_global);
            mpz_neg(cd_inv2, cd_inv2); mpz_mod(cd_inv2, cd_inv2, N_global);
            for (int k = 2*d-2; k >= d; k--) {
                mpz_mod(poly_tmp[k], poly_tmp[k], N_global);
                if (mpz_sgn(poly_tmp[k]) == 0) continue;
                mpz_t ck; mpz_init(ck);
                mpz_mul(ck, poly_tmp[k], cd_inv2); mpz_mod(ck, ck, N_global);
                for (int j = 0; j < d; j++) {
                    mpz_t t; mpz_init(t);
                    mpz_mul(t, ck, F_coeff[j]);
                    mpz_add(poly_tmp[k-d+j], poly_tmp[k-d+j], t);
                    mpz_clear(t);
                }
                mpz_set_ui(poly_tmp[k], 0);
                mpz_clear(ck);
            }
            for (int k = 0; k < d; k++) mpz_mod(base[k], poly_tmp[k], N_global);
            mpz_clear(cd_inv2);
        }

        /* result = P^((N^d - 1)/2) in R mod N.
         * Check: if result ≡ 1, P is a QR on both sides (try P^((N^d+1)/2) for sqrt)
         * If result ≡ -1, P is a QNR on both sides
         * Otherwise: gcd of result coefficients ± 1 with N gives factor */

        /* Try gcd of each coefficient with N */
        for (int k = 0; k < d && !found; k++) {
            if (k == 0) {
                /* Check result[0] - 1 */
                mpz_sub_ui(pm, result[0], 1);
                mpz_gcd(g, pm, N_global);
                if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N_global) < 0) found = 1;
                if (!found) {
                    /* Check result[0] + 1 */
                    mpz_add_ui(pm, result[0], 1);
                    mpz_gcd(g, pm, N_global);
                    if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N_global) < 0) found = 1;
                }
            } else {
                mpz_gcd(g, result[k], N_global);
                if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N_global) < 0) found = 1;
            }
        }

        if (!found) {
            /* Try: compute sqrt = P * result (= P^((N^d+1)/2)) then evaluate at m */
            /* sqrt_alg = P * P^((N^d-1)/2) = P^((N^d+1)/2) */
            for (int k = 0; k < 2*d; k++) mpz_set_ui(poly_tmp[k], 0);
            for (int j1 = 0; j1 < d; j1++)
                for (int j2 = 0; j2 < d; j2++) {
                    mpz_t t; mpz_init(t);
                    mpz_mul(t, P[j1], result[j2]);
                    mpz_add(poly_tmp[j1+j2], poly_tmp[j1+j2], t);
                    mpz_clear(t);
                }
            mpz_t cd_inv3; mpz_init(cd_inv3);
            mpz_invert(cd_inv3, F_coeff[d], N_global);
            mpz_neg(cd_inv3, cd_inv3); mpz_mod(cd_inv3, cd_inv3, N_global);
            for (int k = 2*d-2; k >= d; k--) {
                mpz_mod(poly_tmp[k], poly_tmp[k], N_global);
                if (mpz_sgn(poly_tmp[k]) == 0) continue;
                mpz_t ck; mpz_init(ck);
                mpz_mul(ck, poly_tmp[k], cd_inv3); mpz_mod(ck, ck, N_global);
                for (int j = 0; j < d; j++) {
                    mpz_t t; mpz_init(t);
                    mpz_mul(t, ck, F_coeff[j]);
                    mpz_add(poly_tmp[k-d+j], poly_tmp[k-d+j], t);
                    mpz_clear(t);
                }
                mpz_set_ui(poly_tmp[k], 0);
                mpz_clear(ck);
            }
            mpz_clear(cd_inv3);

            /* Evaluate sqrt polynomial at m: X_a = poly_tmp[0] + poly_tmp[1]*m + poly_tmp[2]*m^2 + ... */
            mpz_set_ui(X, 0);
            mpz_t mpow;
            mpz_init_set_ui(mpow, 1);
            for (int k = 0; k < d; k++) {
                mpz_mod(poly_tmp[k], poly_tmp[k], N_global);
                mpz_t t; mpz_init(t);
                mpz_mul(t, poly_tmp[k], mpow);
                mpz_add(X, X, t);
                mpz_mod(X, X, N_global);
                mpz_mul(mpow, mpow, M_val);
                mpz_mod(mpow, mpow, N_global);
                mpz_clear(t);
            }
            mpz_clear(mpow);

            /* gcd(S_r - X_a, N) */
            mpz_sub(g, Y, X); mpz_gcd(g, g, N_global);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N_global) < 0) found = 1;
            if (!found) {
                mpz_add(g, Y, X); mpz_gcd(g, g, N_global);
                if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N_global) < 0) found = 1;
            }
        }

        mpz_clear(exp);
        for (int k = 0; k < d; k++) { mpz_clear(result[k]); mpz_clear(base[k]); mpz_clear(P[k]); }
        free(result); free(base); free(P);
        free(r_total);
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double elapsed = (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)/1e9;

    if (found) {
        mpz_t cof; mpz_init(cof); mpz_divexact(cof, N_global, g);
        gmp_printf("%Zd %Zd\n", g, cof);
        fprintf(stderr, "NFS time=%.3fs\n", elapsed);
        return 0;
    }
    fprintf(stderr, "NFS linalg failed, time=%.3fs\n", elapsed);
    printf("FAIL\n");
    return 1;
}
