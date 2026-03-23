/*
 * lattice_factor.c - Lattice-based factoring approach
 *
 * Uses lattice reduction (LLL) to find smooth representations of
 * numbers related to N. The key idea from Schnorr-Lenstra:
 *
 * Construct a lattice where short vectors correspond to relations
 * of the form: a_1*log(p_1) + a_2*log(p_2) + ... ≈ log(N)
 * which implies p_1^a_1 * p_2^a_2 * ... ≈ N
 *
 * More precisely, for factor base {p_1,...,p_k}, construct the lattice
 * with basis vectors:
 *   b_i = (0,...,0, C, 0,...,0, round(C*ln(p_i)/ln(N)))
 *                   ^i-th position
 * for i=1..k, plus:
 *   b_{k+1} = (0,...,0, round(C))
 *
 * Short vectors in this lattice correspond to:
 *   sum a_i * ln(p_i) ≈ 0 (mod ln(N))
 * i.e., product of p_i^a_i ≈ N^m for small m
 * i.e., product of p_i^a_i mod N is small and potentially smooth
 *
 * This is a novel approach to finding smooth congruences without sieving.
 *
 * Compile: gcc -O3 -march=native -o lattice_factor library/lattice_factor.c -lgmp -lm
 * Usage: ./lattice_factor <N>
 *
 * NOTE: Uses GMP's built-in functions. For real LLL, we'd need
 * fplll or a custom implementation. This uses a simplified approach.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

#define SEED 42

/*
 * Simplified lattice reduction approach:
 *
 * Instead of full Schnorr-Lenstra, use a practical variant:
 *
 * 1. Choose a factor base B = {p_1, ..., p_k}
 * 2. For random vectors e = (e_1, ..., e_k) with small entries,
 *    compute v = p_1^e_1 * ... * p_k^e_k mod N
 * 3. If v (or N-v) is B-smooth, we have a congruence
 * 4. Use LLL to find short vectors e that make v small
 *
 * The lattice: columns are log representations of primes mod N.
 * A short vector in the lattice corresponds to a small product mod N.
 *
 * Practical variant using "smooth number cannon":
 * Choose random exponents, compute product mod N, check smoothness.
 * Use LLL to bias the search toward vectors that give small products.
 */

/* Matrix operations for LLL (using long long for simplicity) */
typedef struct {
    int rows, cols;
    long long **data;
} matrix_t;

static matrix_t *mat_alloc(int rows, int cols) {
    matrix_t *m = malloc(sizeof(matrix_t));
    m->rows = rows;
    m->cols = cols;
    m->data = malloc(sizeof(long long*) * rows);
    for (int i = 0; i < rows; i++) {
        m->data[i] = calloc(cols, sizeof(long long));
    }
    return m;
}

static void mat_free(matrix_t *m) {
    for (int i = 0; i < m->rows; i++) free(m->data[i]);
    free(m->data);
    free(m);
}

/* Gram-Schmidt orthogonalization (using doubles for speed) */
static void gram_schmidt(matrix_t *B, double **mu, double *Bstar_norm) {
    int n = B->rows;
    int d = B->cols;

    double **bstar = malloc(sizeof(double*) * n);
    for (int i = 0; i < n; i++) {
        bstar[i] = malloc(sizeof(double) * d);
        for (int j = 0; j < d; j++)
            bstar[i][j] = B->data[i][j];

        for (int k = 0; k < i; k++) {
            double dot_bi_bk = 0, dot_bk_bk = Bstar_norm[k];
            for (int j = 0; j < d; j++)
                dot_bi_bk += B->data[i][j] * bstar[k][j];

            mu[i][k] = (dot_bk_bk > 1e-10) ? dot_bi_bk / dot_bk_bk : 0;

            for (int j = 0; j < d; j++)
                bstar[i][j] -= mu[i][k] * bstar[k][j];
        }

        Bstar_norm[i] = 0;
        for (int j = 0; j < d; j++)
            Bstar_norm[i] += bstar[i][j] * bstar[i][j];
    }

    for (int i = 0; i < n; i++) free(bstar[i]);
    free(bstar);
}

/* LLL reduction with delta = 3/4 */
static void lll_reduce(matrix_t *B) {
    int n = B->rows;
    int d = B->cols;
    double delta = 0.75;

    double **mu = malloc(sizeof(double*) * n);
    for (int i = 0; i < n; i++)
        mu[i] = calloc(n, sizeof(double));
    double *Bstar_norm = calloc(n, sizeof(double));

    gram_schmidt(B, mu, Bstar_norm);

    int k = 1;
    while (k < n) {
        /* Size reduction */
        for (int j = k - 1; j >= 0; j--) {
            if (fabs(mu[k][j]) > 0.5) {
                long long r = (long long)round(mu[k][j]);
                for (int i = 0; i < d; i++)
                    B->data[k][i] -= r * B->data[j][i];
                /* Update mu */
                gram_schmidt(B, mu, Bstar_norm);
            }
        }

        /* Lovász condition */
        if (Bstar_norm[k] >= (delta - mu[k][k-1] * mu[k][k-1]) * Bstar_norm[k-1]) {
            k++;
        } else {
            /* Swap k and k-1 */
            long long *tmp = B->data[k];
            B->data[k] = B->data[k-1];
            B->data[k-1] = tmp;

            gram_schmidt(B, mu, Bstar_norm);
            if (k > 1) k--;
        }
    }

    for (int i = 0; i < n; i++) free(mu[i]);
    free(mu);
    free(Bstar_norm);
}

/* Check if v is B-smooth by trial division */
static int is_smooth(mpz_t v, int *primes, int nprimes,
                    unsigned char *exponents, unsigned int *lp,
                    unsigned long lp_bound) {
    mpz_t tmp;
    mpz_init(tmp);
    mpz_abs(tmp, v);

    memset(exponents, 0, nprimes + 1);
    *lp = 0;

    if (mpz_sgn(v) < 0) exponents[0] = 1; /* sign bit */

    for (int i = 0; i < nprimes; i++) {
        while (mpz_divisible_ui_p(tmp, primes[i])) {
            exponents[i + 1] ^= 1;
            mpz_divexact_ui(tmp, tmp, primes[i]);
        }
    }

    int result = 0;
    if (mpz_cmp_ui(tmp, 1) == 0) {
        result = 1; /* Fully smooth */
    } else if (mpz_fits_ulong_p(tmp) && mpz_get_ui(tmp) <= lp_bound) {
        *lp = mpz_get_ui(tmp);
        result = 2; /* SLP */
    }

    mpz_clear(tmp);
    return result;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    mpz_t N;
    mpz_init(N);
    mpz_set_str(N, argv[1], 10);

    int digits = mpz_sizeinbase(N, 10);
    fprintf(stderr, "Factoring %d-digit number using lattice approach\n", digits);

    /* Check small factors */
    for (int p = 2; p < 10000; p++) {
        if (mpz_divisible_ui_p(N, p)) {
            mpz_t cof;
            mpz_init(cof);
            mpz_divexact_ui(cof, N, p);
            gmp_printf("%Zd = %d * %Zd\n", N, p, cof);
            mpz_clear(cof); mpz_clear(N);
            return 0;
        }
    }

    /* Factor base */
    double ln_n = digits * log(10);
    double ln_ln_n = log(ln_n);
    int smooth_bound = (int)exp(0.5 * sqrt(ln_n * ln_ln_n));
    if (smooth_bound < 300) smooth_bound = 300;
    if (smooth_bound > 500000) smooth_bound = 500000;

    char *sieve = calloc(smooth_bound + 1, 1);
    int *primes = malloc(sizeof(int) * (smooth_bound + 100));
    int nprimes = 0;
    for (int i = 2; i <= smooth_bound; i++) {
        if (!sieve[i]) {
            /* Only include primes where N is QR mod p */
            if (i == 2 || mpz_kronecker_ui(N, i) >= 0) {
                primes[nprimes++] = i;
            }
            for (long j = (long)i*i; j <= smooth_bound; j += i)
                sieve[j] = 1;
        }
    }
    free(sieve);

    fprintf(stderr, "Factor base: %d primes up to %d\n", nprimes, smooth_bound);

    unsigned long lp_bound = (unsigned long)primes[nprimes-1] * 50;
    int target = nprimes + 30;

    /*
     * APPROACH 1: Random smooth number search with lattice guidance
     *
     * The lattice approach: construct vectors (e_1,...,e_k) such that
     * p_1^e_1 * ... * p_k^e_k is close to N (or a power of N).
     * Then the residue mod N is small and more likely smooth.
     *
     * Practical implementation:
     * 1. Use a subset of FB primes
     * 2. Construct lattice: each row represents a prime's contribution
     * 3. LLL reduce the lattice
     * 4. Short vectors in reduced basis give small residues mod N
     * 5. Test each for smoothness
     */

    /* Use a smaller subset for lattice (LLL is O(d^6)) */
    int lattice_dim = nprimes;
    if (lattice_dim > 40) lattice_dim = 40; /* LLL practical limit */

    fprintf(stderr, "Lattice dimension: %d\n", lattice_dim);

    /* Scale factor for log representation */
    double log_N = mpz_sizeinbase(N, 2) * log(2);
    long long C = 1LL << 30; /* Large constant for precision */

    /*
     * Lattice basis: (lattice_dim+1) x (lattice_dim+1) matrix
     *
     * Row i (for prime p_i): e_i with 1 in position i, and
     *   round(C * ln(p_i)) in the last column
     *
     * Row k+1 (for N): 0s except round(C * ln(N)) in last column
     *
     * Short vectors: sum(a_i * ln(p_i)) ≈ m * ln(N)
     * => product(p_i^a_i) ≈ N^m
     * => product(p_i^a_i) mod N is small
     */

    int dim = lattice_dim + 1;
    matrix_t *basis = mat_alloc(dim, dim);

    /* Fill lattice */
    for (int i = 0; i < lattice_dim; i++) {
        basis->data[i][i] = 1; /* Identity in first lattice_dim columns */
        basis->data[i][lattice_dim] = (long long)round(C * log(primes[i]));
    }
    basis->data[lattice_dim][lattice_dim] = (long long)round(C * log_N);

    fprintf(stderr, "Running LLL reduction...\n");
    lll_reduce(basis);
    fprintf(stderr, "LLL done.\n");

    /* Each row of the reduced basis gives exponents (a_1,...,a_k)
     * such that product(p_i^a_i) ≈ N^m for some small m.
     * We test: v = product(p_i^{a_i}) mod N for smoothness.
     *
     * Also: perturb reduced basis vectors to generate more candidates.
     */

    int rel_count = 0;
    int total_tested = 0;

    /* Storage for relations */
    typedef struct {
        mpz_t v;       /* value mod N */
        mpz_t prod;    /* product of p_i^a_i mod N */
        unsigned char *exp; /* exponent vector */
        unsigned int lp;
    } lattice_rel_t;

    lattice_rel_t *rels = malloc(sizeof(lattice_rel_t) * (target + 100));
    unsigned char *exp_buf = calloc(nprimes + 1, 1);

    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, SEED);

    mpz_t prod, tmp;
    mpz_init(prod);
    mpz_init(tmp);

    /* Test each reduced basis vector */
    for (int row = 0; row < dim && rel_count < target; row++) {
        /* Extract exponents from first lattice_dim columns */
        mpz_set_ui(prod, 1);
        int valid = 1;

        for (int i = 0; i < lattice_dim; i++) {
            long long e = basis->data[row][i];
            if (e == 0) continue;

            mpz_set_ui(tmp, primes[i]);
            if (e > 0) {
                if (e > 1000) { valid = 0; break; } /* Exponent too large */
                mpz_powm_ui(tmp, tmp, e, N);
            } else {
                /* Negative exponent: use modular inverse */
                if (-e > 1000) { valid = 0; break; }
                mpz_invert(tmp, tmp, N);
                if (mpz_sgn(tmp) == 0) { valid = 0; break; }
                mpz_powm_ui(tmp, tmp, -e, N);
            }
            mpz_mul(prod, prod, tmp);
            mpz_mod(prod, prod, N);
        }

        if (!valid) continue;
        total_tested++;

        /* Check if prod or N-prod is smooth */
        unsigned int lp = 0;
        int smooth = is_smooth(prod, primes, nprimes, exp_buf, &lp, lp_bound);

        if (!smooth) {
            /* Try N - prod */
            mpz_sub(tmp, N, prod);
            smooth = is_smooth(tmp, primes, nprimes, exp_buf, &lp, lp_bound);
            if (smooth) mpz_set(prod, tmp);
        }

        if (smooth) {
            lattice_rel_t *r = &rels[rel_count];
            mpz_init_set(r->v, prod);
            mpz_init_set(r->prod, prod);
            r->exp = malloc(nprimes + 1);
            memcpy(r->exp, exp_buf, nprimes + 1);
            r->lp = lp;
            rel_count++;
            fprintf(stderr, "  Row %d: found relation (type %d, lp=%u)\n", row, smooth, lp);
        }
    }

    fprintf(stderr, "Basis vectors: %d tested, %d relations\n", total_tested, rel_count);

    /* Generate more candidates by random combinations of basis vectors */
    fprintf(stderr, "Generating random lattice combinations...\n");
    int max_iters = 1000000;

    for (int iter = 0; iter < max_iters && rel_count < target; iter++) {
        clock_gettime(CLOCK_MONOTONIC, &t1);
        double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
        if (elapsed > 290.0) break;

        /* Random small combination of basis vectors */
        long long combo[50] = {0};
        for (int j = 0; j < dim; j++) combo[j] = 0;

        /* Sparse random combination: pick 2-4 basis vectors with small coefficients */
        int n_terms = 2 + gmp_urandomm_ui(rstate, 3);
        for (int t = 0; t < n_terms; t++) {
            int r = gmp_urandomm_ui(rstate, dim);
            int c = (int)gmp_urandomm_ui(rstate, 5) - 2; /* -2 to 2 */
            for (int j = 0; j < dim; j++)
                combo[j] += c * basis->data[r][j];
        }

        /* Compute product of p_i^{combo_i} mod N */
        mpz_set_ui(prod, 1);
        int valid = 1;
        for (int i = 0; i < lattice_dim; i++) {
            long long e = combo[i];
            if (e == 0) continue;
            if (e > 500 || e < -500) { valid = 0; break; }

            mpz_set_ui(tmp, primes[i]);
            if (e > 0) {
                mpz_powm_ui(tmp, tmp, e, N);
            } else {
                mpz_invert(tmp, tmp, N);
                if (mpz_sgn(tmp) == 0) { valid = 0; break; }
                mpz_powm_ui(tmp, tmp, -e, N);
            }
            mpz_mul(prod, prod, tmp);
            mpz_mod(prod, prod, N);
        }
        if (!valid) continue;
        total_tested++;

        /* Check smoothness */
        unsigned int lp = 0;
        int smooth = is_smooth(prod, primes, nprimes, exp_buf, &lp, lp_bound);
        if (!smooth) {
            mpz_sub(tmp, N, prod);
            smooth = is_smooth(tmp, primes, nprimes, exp_buf, &lp, lp_bound);
            if (smooth) mpz_set(prod, tmp);
        }

        if (smooth) {
            lattice_rel_t *r = &rels[rel_count];
            mpz_init_set(r->v, prod);
            mpz_init_set(r->prod, prod);
            r->exp = malloc(nprimes + 1);
            memcpy(r->exp, exp_buf, nprimes + 1);
            r->lp = lp;
            rel_count++;
        }

        if ((iter + 1) % 100000 == 0) {
            clock_gettime(CLOCK_MONOTONIC, &t1);
            elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
            fprintf(stderr, "[%.1fs] iter=%d tested=%d rels=%d/%d (%.0f/s)\n",
                    elapsed, iter+1, total_tested, rel_count, target,
                    rel_count / elapsed);
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
    fprintf(stderr, "Total: %d tested, %d relations in %.1fs\n",
            total_tested, rel_count, elapsed);

    if (rel_count < nprimes + 1) {
        fprintf(stderr, "Not enough relations (%d/%d)\n", rel_count, nprimes + 1);
        gmp_printf("FAIL %Zd\n", N);
    } else {
        fprintf(stderr, "TODO: linear algebra step not yet implemented\n");
        /* Would need GF(2) Gaussian elimination on the exponent matrix,
         * then compute x = product of (p_i^a_i mod N) for dependency,
         * y = sqrt of product of v values, gcd(x-y, N) */
        gmp_printf("FAIL %Zd\n", N);
    }

    /* Cleanup */
    free(exp_buf);
    mpz_clear(prod); mpz_clear(tmp);
    for (int i = 0; i < rel_count; i++) {
        mpz_clear(rels[i].v);
        mpz_clear(rels[i].prod);
        free(rels[i].exp);
    }
    free(rels);
    mat_free(basis);
    free(primes);
    gmp_randclear(rstate);
    mpz_clear(N);

    return 0;
}
