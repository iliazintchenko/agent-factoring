/*
 * lattice_factor.c - Lattice-based integer factoring (Schnorr-Seysen style)
 *
 * Compile: gcc -O3 -march=native -o lattice_factor library/lattice_factor.c -lgmp -lm
 *
 * Algorithm:
 *   1. Select factor base of small primes where N is a QR mod p
 *   2. Compute sqrt(N) mod p for each factor base prime (Tonelli-Shanks)
 *   3. Sieve for x values near sqrt(N) where x^2-N is smooth (basic QS-like)
 *   4. Use lattice (LLL) to find additional x values via CRT where
 *      x^2-N is divisible by many factor base primes
 *   5. Gaussian elimination mod 2 on the exponent matrix
 *   6. Combine relations: x^2 = y^2 (mod N), gcd(x-y, N)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <math.h>
#include <time.h>

/* ================================================================== */
/* Tonelli-Shanks: compute r such that r^2 = n (mod p), p odd prime   */
/* Returns 1 on success, 0 if n is not a QR mod p                     */
/* ================================================================== */
static int tonelli_shanks(mpz_t r, const mpz_t n, const mpz_t p) {
    if (mpz_cmp_ui(p, 2) == 0) {
        mpz_set(r, n);
        mpz_mod(r, r, p);
        return 1;
    }

    int leg = mpz_legendre(n, p);
    if (leg == -1) return 0;
    if (leg == 0) { mpz_set_ui(r, 0); return 1; }

    mpz_t Q, pm1, z, c, t, R_val, tmp, b;
    mpz_inits(Q, pm1, z, c, t, R_val, tmp, b, NULL);

    mpz_sub_ui(pm1, p, 1);
    mpz_set(Q, pm1);
    unsigned long S = 0;
    while (mpz_even_p(Q)) { mpz_fdiv_q_2exp(Q, Q, 1); S++; }

    if (S == 1) {
        mpz_add_ui(tmp, p, 1);
        mpz_fdiv_q_2exp(tmp, tmp, 2);
        mpz_powm(r, n, tmp, p);
        mpz_clears(Q, pm1, z, c, t, R_val, tmp, b, NULL);
        return 1;
    }

    mpz_set_ui(z, 2);
    while (mpz_legendre(z, p) != -1) mpz_add_ui(z, z, 1);

    unsigned long M = S;
    mpz_powm(c, z, Q, p);
    mpz_add_ui(tmp, Q, 1);
    mpz_fdiv_q_2exp(tmp, tmp, 1);
    mpz_powm(R_val, n, tmp, p);
    mpz_powm(t, n, Q, p);

    while (1) {
        if (mpz_cmp_ui(t, 1) == 0) { mpz_set(r, R_val); break; }
        if (mpz_cmp_ui(t, 0) == 0) { mpz_set_ui(r, 0); break; }

        unsigned long i = 0;
        mpz_set(tmp, t);
        while (mpz_cmp_ui(tmp, 1) != 0) { mpz_powm_ui(tmp, tmp, 2, p); i++; }

        mpz_set(b, c);
        for (unsigned long j = 0; j < M - i - 1; j++) mpz_powm_ui(b, b, 2, p);

        M = i;
        mpz_powm_ui(c, b, 2, p);
        mpz_mul(R_val, R_val, b); mpz_mod(R_val, R_val, p);
        mpz_mul(t, t, c); mpz_mod(t, t, p);
    }

    mpz_clears(Q, pm1, z, c, t, R_val, tmp, b, NULL);
    return 1;
}

/* ================================================================== */
/* Small primes sieve                                                  */
/* ================================================================== */
static int *sieve_primes(int limit, int *count) {
    char *is_comp = calloc(limit + 1, 1);
    int *primes = malloc(sizeof(int) * (limit / 2 + 10));
    int cnt = 0;
    for (int i = 2; i <= limit; i++) {
        if (!is_comp[i]) {
            primes[cnt++] = i;
            if ((long long)i * i <= limit)
                for (int j = i * i; j <= limit; j += i) is_comp[j] = 1;
        }
    }
    free(is_comp);
    *count = cnt;
    return primes;
}

/* ================================================================== */
/* LLL lattice reduction - GMP rationals for exact Gram-Schmidt       */
/* basis: n x m matrix stored row-major as mpz_t[n*m]                 */
/* Uses delta = 3/4                                                    */
/* ================================================================== */
#define MAT(B, i, j, m) ((B)[(i)*(m) + (j)])

static void lll_reduce(mpz_t *basis, int n, int m) {
    mpq_t *mu = malloc(sizeof(mpq_t) * n * n);
    mpq_t *Bsq = malloc(sizeof(mpq_t) * n);
    mpq_t delta_q, tmp_q, tmp_q2, half_q, neg_half_q;
    mpz_t tmp_z, rmu;

    for (int i = 0; i < n * n; i++) mpq_init(mu[i]);
    for (int i = 0; i < n; i++) mpq_init(Bsq[i]);
    mpq_init(delta_q); mpq_init(tmp_q); mpq_init(tmp_q2);
    mpq_init(half_q); mpq_init(neg_half_q);
    mpz_init(tmp_z); mpz_init(rmu);

    mpq_set_ui(delta_q, 3, 4);
    mpq_set_ui(half_q, 1, 2);
    mpq_set_si(neg_half_q, -1, 2);

    /* Dot product of basis rows ri and rj -> result as mpq (integer-valued) */
    #define DOT(res, ri, rj) do { \
        mpz_set_ui(mpq_numref(res), 0); \
        for (int _d = 0; _d < m; _d++) \
            mpz_addmul(mpq_numref(res), MAT(basis,ri,_d,m), MAT(basis,rj,_d,m)); \
        mpz_set_ui(mpq_denref(res), 1); \
    } while(0)

    /* Compute full Gram-Schmidt data from scratch starting at `from` */
    auto void compute_gs(int from);
    void compute_gs(int from) {
        for (int i = from; i < n; i++) {
            DOT(Bsq[i], i, i);
            for (int j = 0; j < i; j++) {
                DOT(mu[i*n+j], i, j);
                for (int l = 0; l < j; l++) {
                    mpq_mul(tmp_q, mu[j*n+l], mu[i*n+l]);
                    mpq_mul(tmp_q, tmp_q, Bsq[l]);
                    mpq_sub(mu[i*n+j], mu[i*n+j], tmp_q);
                }
                if (mpq_sgn(Bsq[j]) != 0)
                    mpq_div(mu[i*n+j], mu[i*n+j], Bsq[j]);
            }
            for (int j = 0; j < i; j++) {
                mpq_mul(tmp_q, mu[i*n+j], mu[i*n+j]);
                mpq_mul(tmp_q, tmp_q, Bsq[j]);
                mpq_sub(Bsq[i], Bsq[i], tmp_q);
            }
        }
    }

    compute_gs(0);

    int k = 1;
    while (k < n) {
        /* Size-reduce b_k */
        for (int j = k - 1; j >= 0; j--) {
            if (mpq_cmp(mu[k*n+j], half_q) > 0 || mpq_cmp(mu[k*n+j], neg_half_q) < 0) {
                /* round(mu) */
                mpz_mul_2exp(tmp_z, mpq_numref(mu[k*n+j]), 1);
                mpz_add(tmp_z, tmp_z, mpq_denref(mu[k*n+j]));
                mpz_mul_2exp(rmu, mpq_denref(mu[k*n+j]), 1);
                mpz_fdiv_q(rmu, tmp_z, rmu);
                if (mpz_sgn(rmu) == 0) continue;

                for (int d = 0; d < m; d++)
                    mpz_submul(MAT(basis,k,d,m), rmu, MAT(basis,j,d,m));

                mpq_set_z(tmp_q, rmu);
                mpq_sub(mu[k*n+j], mu[k*n+j], tmp_q);
                for (int l = 0; l < j; l++) {
                    mpq_mul(tmp_q2, tmp_q, mu[j*n+l]);
                    mpq_sub(mu[k*n+l], mu[k*n+l], tmp_q2);
                }
            }
        }

        /* Lovász condition */
        mpq_mul(tmp_q, mu[k*n+(k-1)], mu[k*n+(k-1)]);
        mpq_sub(tmp_q, delta_q, tmp_q);
        mpq_mul(tmp_q, tmp_q, Bsq[k-1]);

        if (mpq_cmp(Bsq[k], tmp_q) >= 0) {
            k++;
        } else {
            /* Swap rows k and k-1 */
            for (int d = 0; d < m; d++)
                mpz_swap(MAT(basis,k,d,m), MAT(basis,k-1,d,m));

            /* Update Gram-Schmidt data using swap formulas */
            mpq_set(tmp_q, mu[k*n+(k-1)]);

            mpq_t Bnew_km1, mu_new, Bnew_k;
            mpq_init(Bnew_km1); mpq_init(mu_new); mpq_init(Bnew_k);

            /* Bnew[k-1] = B[k] + mu^2*B[k-1] */
            mpq_mul(tmp_q2, tmp_q, tmp_q);
            mpq_mul(tmp_q2, tmp_q2, Bsq[k-1]);
            mpq_add(Bnew_km1, Bsq[k], tmp_q2);

            if (mpq_sgn(Bnew_km1) == 0) {
                mpq_clear(Bnew_km1); mpq_clear(mu_new); mpq_clear(Bnew_k);
                compute_gs(0);
                k = (k > 1) ? k - 1 : 1;
                continue;
            }

            /* mu_new[k][k-1] */
            mpq_mul(mu_new, tmp_q, Bsq[k-1]);
            mpq_div(mu_new, mu_new, Bnew_km1);

            /* Bnew[k] */
            mpq_mul(Bnew_k, Bsq[k-1], Bsq[k]);
            mpq_div(Bnew_k, Bnew_k, Bnew_km1);

            /* Update mu for rows i > k */
            for (int i = k+1; i < n; i++) {
                mpq_t oikm1, oik;
                mpq_init(oikm1); mpq_init(oik);
                mpq_set(oikm1, mu[i*n+(k-1)]);
                mpq_set(oik, mu[i*n+k]);

                mpq_mul(mu[i*n+(k-1)], oikm1, mu_new);
                mpq_mul(tmp_q2, oik, Bsq[k]);
                mpq_div(tmp_q2, tmp_q2, Bnew_km1);
                mpq_add(mu[i*n+(k-1)], mu[i*n+(k-1)], tmp_q2);

                mpq_mul(mu[i*n+k], tmp_q, oik);
                mpq_sub(mu[i*n+k], oikm1, mu[i*n+k]);

                mpq_clear(oikm1); mpq_clear(oik);
            }

            for (int j = 0; j < k-1; j++)
                mpq_swap(mu[(k-1)*n+j], mu[k*n+j]);

            mpq_set(mu[k*n+(k-1)], mu_new);
            mpq_set(Bsq[k-1], Bnew_km1);
            mpq_set(Bsq[k], Bnew_k);

            mpq_clear(Bnew_km1); mpq_clear(mu_new); mpq_clear(Bnew_k);
            k = (k > 1) ? k-1 : 1;
        }
    }

    for (int i = 0; i < n*n; i++) mpq_clear(mu[i]);
    for (int i = 0; i < n; i++) mpq_clear(Bsq[i]);
    mpq_clear(delta_q); mpq_clear(tmp_q); mpq_clear(tmp_q2);
    mpq_clear(half_q); mpq_clear(neg_half_q);
    mpz_clear(tmp_z); mpz_clear(rmu);
    free(mu); free(Bsq);
}

/* ================================================================== */
/* Factor base selection                                               */
/* ================================================================== */
static int select_factor_base(const mpz_t N, int *fb, int fb_max, int prime_limit) {
    int nprimes;
    int *primes = sieve_primes(prime_limit, &nprimes);
    int cnt = 0;
    mpz_t p_mpz;
    mpz_init(p_mpz);
    for (int i = 0; i < nprimes && cnt < fb_max; i++) {
        if (primes[i] == 2) { fb[cnt++] = 2; continue; }
        mpz_set_ui(p_mpz, primes[i]);
        if (mpz_legendre(N, p_mpz) >= 0) fb[cnt++] = primes[i];
    }
    mpz_clear(p_mpz);
    free(primes);
    return cnt;
}

/* ================================================================== */
/* Trial-divide val over factor base, fill exponent vector             */
/* Returns 1 if fully smooth                                           */
/* ================================================================== */
static int trial_factor(const mpz_t val, const int *fb, int fb_size, int *exponents) {
    memset(exponents, 0, sizeof(int) * fb_size);
    if (mpz_sgn(val) == 0) return 0;

    mpz_t tmp;
    mpz_init(tmp);
    mpz_abs(tmp, val);

    for (int i = 0; i < fb_size; i++) {
        while (mpz_divisible_ui_p(tmp, fb[i])) {
            mpz_divexact_ui(tmp, tmp, fb[i]);
            exponents[i]++;
        }
    }
    int smooth = (mpz_cmp_ui(tmp, 1) == 0);
    mpz_clear(tmp);
    return smooth;
}

/* ================================================================== */
/* Relation storage                                                    */
/* ================================================================== */
typedef struct {
    mpz_t x;          /* x value */
    int *exponents;    /* exponent vector (fb_size entries) */
    int sign;          /* 1 if x^2 - N < 0 */
} relation_t;

/* ================================================================== */
/* Gaussian elimination mod 2                                          */
/* Returns list of ALL dependencies found (array of arrays)            */
/* ================================================================== */
typedef unsigned long bitmask_t;
#define BPW (sizeof(bitmask_t)*8)

typedef struct {
    int *indices;
    int count;
} dependency_t;

static dependency_t *gauss_find_deps(relation_t *rels, int nrels, int fb_size, int *ndeps_out) {
    int ncols = fb_size + 1;  /* +1 for sign */
    int nw_exp = (ncols + BPW - 1) / BPW;
    int nw_hist = (nrels + BPW - 1) / BPW;
    int nw = nw_exp + nw_hist;

    bitmask_t **rows = malloc(sizeof(bitmask_t*) * nrels);
    for (int i = 0; i < nrels; i++) {
        rows[i] = calloc(nw, sizeof(bitmask_t));
        for (int j = 0; j < fb_size; j++)
            if (rels[i].exponents[j] & 1)
                rows[i][j/BPW] |= (1UL << (j%BPW));
        if (rels[i].sign & 1)
            rows[i][fb_size/BPW] |= (1UL << (fb_size%BPW));
        rows[i][nw_exp + i/BPW] |= (1UL << (i%BPW));
    }

    /* Row echelon form */
    int *used = calloc(nrels, sizeof(int));
    for (int col = 0; col < ncols; col++) {
        int pr = -1;
        for (int r = 0; r < nrels; r++) {
            if (!used[r] && (rows[r][col/BPW] & (1UL << (col%BPW)))) {
                pr = r; break;
            }
        }
        if (pr == -1) continue;
        used[pr] = 1;
        for (int r = 0; r < nrels; r++) {
            if (r != pr && (rows[r][col/BPW] & (1UL << (col%BPW)))) {
                for (int w = 0; w < nw; w++) rows[r][w] ^= rows[pr][w];
            }
        }
    }

    /* Collect all zero rows (dependencies) */
    int max_deps = nrels;
    dependency_t *deps = malloc(sizeof(dependency_t) * max_deps);
    int ndeps = 0;

    for (int r = 0; r < nrels; r++) {
        int all_zero = 1;
        for (int w = 0; w < nw_exp; w++)
            if (rows[r][w]) { all_zero = 0; break; }
        if (!all_zero) continue;

        int cnt = 0;
        for (int b = 0; b < nrels; b++)
            if (rows[r][nw_exp + b/BPW] & (1UL << (b%BPW))) cnt++;
        if (cnt < 1) continue;

        deps[ndeps].indices = malloc(sizeof(int) * cnt);
        deps[ndeps].count = cnt;
        int idx = 0;
        for (int b = 0; b < nrels; b++)
            if (rows[r][nw_exp + b/BPW] & (1UL << (b%BPW)))
                deps[ndeps].indices[idx++] = b;
        ndeps++;
    }

    for (int i = 0; i < nrels; i++) free(rows[i]);
    free(rows);
    free(used);
    *ndeps_out = ndeps;
    return deps;
}

/* ================================================================== */
/* Try to extract a factor from a dependency                           */
/* Returns 1 if nontrivial factor found                                */
/* ================================================================== */
static int try_dependency(const mpz_t N, relation_t *rels, const int *fb,
                          int fb_size, const dependency_t *dep, mpz_t factor_out) {
    /* x_prod = product of x_i mod N */
    /* y = product of p_j^(sum_exp_j/2) mod N  (with sign adjustment) */
    mpz_t x_prod, y, g;
    mpz_inits(x_prod, y, g, NULL);
    mpz_set_ui(x_prod, 1);

    int *sum_exp = calloc(fb_size, sizeof(int));
    int sign_sum = 0;

    for (int d = 0; d < dep->count; d++) {
        int ri = dep->indices[d];
        mpz_mul(x_prod, x_prod, rels[ri].x);
        mpz_mod(x_prod, x_prod, N);
        for (int j = 0; j < fb_size; j++)
            sum_exp[j] += rels[ri].exponents[j];
        sign_sum += rels[ri].sign;
    }

    /* Verify all exponents are even */
    int ok = 1;
    for (int j = 0; j < fb_size; j++) {
        if (sum_exp[j] % 2 != 0) { ok = 0; break; }
    }
    if (sign_sum % 2 != 0) ok = 0;

    if (!ok) {
        free(sum_exp);
        mpz_clears(x_prod, y, g, NULL);
        return 0;
    }

    /* Compute y = product of p_j^(sum_exp_j/2) mod N */
    mpz_set_ui(y, 1);
    for (int j = 0; j < fb_size; j++) {
        int he = sum_exp[j] / 2;
        if (he > 0) {
            mpz_t pe;
            mpz_init(pe);
            mpz_ui_pow_ui(pe, fb[j], he);
            mpz_mul(y, y, pe);
            mpz_mod(y, y, N);
            mpz_clear(pe);
        }
    }

    /* Try gcd(x_prod - y, N) */
    mpz_sub(g, x_prod, y);
    mpz_mod(g, g, N);
    mpz_gcd(g, g, N);
    if (mpz_cmp_ui(g, 1) != 0 && mpz_cmp(g, N) != 0) {
        mpz_set(factor_out, g);
        free(sum_exp);
        mpz_clears(x_prod, y, g, NULL);
        return 1;
    }

    /* Try gcd(x_prod + y, N) */
    mpz_add(g, x_prod, y);
    mpz_mod(g, g, N);
    mpz_gcd(g, g, N);
    if (mpz_cmp_ui(g, 1) != 0 && mpz_cmp(g, N) != 0) {
        mpz_set(factor_out, g);
        free(sum_exp);
        mpz_clears(x_prod, y, g, NULL);
        return 1;
    }

    free(sum_exp);
    mpz_clears(x_prod, y, g, NULL);
    return 0;
}

/* ================================================================== */
/* Main lattice factoring routine                                      */
/* ================================================================== */
static int lattice_factor(const mpz_t N, mpz_t factor_out) {
    size_t digits = mpz_sizeinbase(N, 10);

    /* Configure factor base size based on digit count.
     * Use the L(N) = exp(sqrt(ln N * ln ln N)) heuristic.
     * The optimal smoothness bound B ≈ L(N)^{1/2√2}.
     * fb_max limits how many primes we select; prime_limit is the sieve
     * bound for generating candidate primes. */
    /* B = L(N)^{1/sqrt(2)} is the standard QS smoothness bound.
     * We scale it down slightly since lattice methods with trial
     * factoring are slower per candidate than proper QS sieving. */
    double lnN = digits * log(10.0);
    double L = exp(sqrt(lnN * log(lnN)));
    /* Balance factor base size vs smoothness probability.
     * Larger B = more smooth numbers but more relations needed.
     * With a proper sieve, higher B is better. */
    double B = pow(L, 0.55);
    if (B < 100) B = 100;
    if (B > 200000) B = 200000;

    int prime_limit = (int)(B * 1.2);
    int fb_max = (int)(B / log(B)) + 20;
    if (fb_max > 5000) fb_max = 5000;

    int *fb = malloc(sizeof(int) * fb_max);
    int fb_size = select_factor_base(N, fb, fb_max, prime_limit);
    if (fb_size < 3) {
        fprintf(stderr, "Factor base too small (%d)\n", fb_size);
        free(fb);
        return 0;
    }

    printf("  Factor base: %d primes (2 .. %d)\n", fb_size, fb[fb_size-1]);

    /* Compute sqrt(N) mod p for each FB prime */
    mpz_t *roots = malloc(sizeof(mpz_t) * fb_size);
    mpz_t p_mpz, n_mod_p;
    mpz_init(p_mpz); mpz_init(n_mod_p);

    for (int i = 0; i < fb_size; i++) {
        mpz_init(roots[i]);
        mpz_set_ui(p_mpz, fb[i]);
        mpz_mod(n_mod_p, N, p_mpz);
        if (fb[i] == 2) {
            mpz_set_ui(roots[i], mpz_odd_p(N) ? 1 : 0);
        } else {
            tonelli_shanks(roots[i], n_mod_p, p_mpz);
        }
    }

    mpz_t sqrt_N;
    mpz_init(sqrt_N);
    mpz_sqrt(sqrt_N, N);

    /* Allocate relations */
    int needed = fb_size + 5;
    int max_rels = needed * 4;
    relation_t *rels = malloc(sizeof(relation_t) * max_rels);
    for (int i = 0; i < max_rels; i++) {
        mpz_init(rels[i].x);
        rels[i].exponents = calloc(fb_size, sizeof(int));
        rels[i].sign = 0;
    }
    int nrels = 0;

    /* ------------------------------------------------------------ */
    /* Phase 1: Logarithmic sieve near sqrt(N) for smooth x^2-N    */
    /* Uses sieve-by-primes: accumulate log contributions, then     */
    /* only trial-factor candidates whose log sum is close to the   */
    /* expected log of x^2-N.                                       */
    /* ------------------------------------------------------------ */
    {
        long sieve_size = 500000;
        if (digits > 15) sieve_size = 2000000;
        if (digits > 25) sieve_size = 8000000;
        if (digits > 35) sieve_size = 30000000;
        if (digits > 45) sieve_size = 60000000;

        /* sieve_start = sqrt(N) - sieve_size/2 */
        mpz_t sieve_start;
        mpz_init(sieve_start);
        mpz_sub_ui(sieve_start, sqrt_N, sieve_size / 2);
        if (mpz_sgn(sieve_start) <= 0) mpz_set_ui(sieve_start, 2);

        /* Allocate log sieve array */
        float *logsum = calloc(sieve_size, sizeof(float));

        /* For the sieve, we accumulate log(p) for each prime dividing x^2-N.
         * A position is a candidate if the accumulated log is close to
         * log(|x^2-N|). Since |x^2-N| varies with position, we use a
         * per-position approach: store expected_log[j] and compare.
         *
         * expected_log[j] = log(|sieve_start+j|^2 - N|)
         * We precompute this approximately. For positions near sqrt(N):
         *   |x^2-N| ≈ |2*sqrt(N)*(x - sqrt(N)) + (x-sqrt(N))^2|
         *
         * Threshold: accept if logsum[j] > expected_log[j] - T
         * where T accounts for one large prime we might miss. */
        double log_largest_fb = log((double)fb[fb_size-1]);
        /* Large tolerance: we expect smooth numbers to have most of their
         * log accounted for by sieved primes. Allow missing up to half
         * the expected log since prime powers and near-misses are common. */
        double T = log_largest_fb * 3.0 + 5.0;

        /* Precompute approximate expected log for sieve positions.
         * Position j in the sieve corresponds to x = sieve_start + j.
         * We compute delta = x - sqrt(N), then |x^2-N| ≈ |2*sqrt(N)*delta + delta^2|.
         * For efficiency, compute log(|x^2-N|) approximately using double. */
        double dsqrt = mpz_get_d(sqrt_N);
        double dstart = mpz_get_d(sieve_start);
        float *expected_log = malloc(sizeof(float) * sieve_size);
        for (long j = 0; j < sieve_size; j++) {
            double dx = dstart + j;
            double delta = dx - dsqrt;
            double residual = fabs(2.0 * dsqrt * delta + delta * delta);
            if (residual < 1.0) residual = 1.0;
            expected_log[j] = (float)log(residual);
        }

        /* For each FB prime, find starting offsets and sieve */
        mpz_t tmp_off;
        mpz_init(tmp_off);

        for (int fi = 0; fi < fb_size; fi++) {
            int p = fb[fi];
            float logp = logf((float)p);

            /* Find offset where (sieve_start + offset)^2 ≡ N (mod p)
             * i.e., sieve_start + offset ≡ ±root (mod p)
             * offset ≡ root - sieve_start (mod p)  and  offset ≡ -root - sieve_start (mod p) */

            /* Compute sieve_start mod p */
            unsigned long s_mod_p = mpz_fdiv_ui(sieve_start, p);
            unsigned long r = mpz_get_ui(roots[fi]);

            /* For p=2, handle specially */
            if (p == 2) {
                /* x^2 - N mod 2: if N is odd, x^2-N is even when x is odd */
                unsigned long start = (s_mod_p % 2 == 1) ? 0 : 1;
                for (long j = start; j < sieve_size; j += 2)
                    logsum[j] += logp;
                continue;
            }

            /* Two roots: r and p-r */
            for (int sign = 0; sign < 2; sign++) {
                unsigned long root_val = (sign == 0) ? r : (p - r);
                if (root_val >= (unsigned long)p) root_val %= p;

                long off = (long)root_val - (long)s_mod_p;
                if (off < 0) off += p;

                /* Sieve with this root */
                for (long j = off; j < sieve_size; j += p)
                    logsum[j] += logp;

                /* Also sieve prime powers */
                long pp = (long)p * p;
                while (pp < sieve_size * 2 && pp > 0) {
                    /* Find offset for p^k */
                    /* We need (sieve_start + j)^2 ≡ N (mod p^k) */
                    /* This requires lifting the root, which is complex.
                     * For simplicity, just check divisibility during trial factor. */
                    break;
                }
            }
        }

        /* Now scan for candidates where logsum ≈ expected_log and trial-factor */
        mpz_t sx, sv;
        mpz_init(sx); mpz_init(sv);
        int *ev = calloc(fb_size, sizeof(int));
        int candidates = 0;

        for (long j = 0; j < sieve_size && nrels < max_rels; j++) {
            if (logsum[j] < expected_log[j] - T) continue;
            candidates++;

            mpz_add_ui(sx, sieve_start, j);
            mpz_mul(sv, sx, sx);
            mpz_sub(sv, sv, N);
            if (mpz_sgn(sv) == 0) continue;

            int sign = (mpz_sgn(sv) < 0) ? 1 : 0;

            if (trial_factor(sv, fb, fb_size, ev)) {
                mpz_set(rels[nrels].x, sx);
                memcpy(rels[nrels].exponents, ev, sizeof(int)*fb_size);
                rels[nrels].sign = sign;
                nrels++;
            }
        }

        free(ev);
        free(logsum);
        free(expected_log);
        mpz_clears(sx, sv, sieve_start, tmp_off, NULL);
        printf("  Sieve phase: %d smooth relations found (checked %d candidates)\n", nrels, candidates);
    }

    /* ------------------------------------------------------------ */
    /* Phase 2: Lattice-based CRT approach for more relations        */
    /* Use subsets of FB primes, CRT to find x ≡ r_i (mod p_i),     */
    /* then search near sqrt(N) for smooth x^2-N.                   */
    /* ------------------------------------------------------------ */
    if (nrels < needed) {
        printf("  Lattice-CRT phase (have %d, need %d)...\n", nrels, needed);

        mpz_t crt_x, crt_mod, tmp_crt, inv_z;
        mpz_inits(crt_x, crt_mod, tmp_crt, inv_z, NULL);

        int *ev = calloc(fb_size, sizeof(int));
        mpz_t base_x, test_x, test_val;
        mpz_inits(base_x, test_x, test_val, NULL);

        /* Try subsets of increasing size */
        for (int ss = 3; ss <= fb_size && ss <= 12 && nrels < needed; ss++) {
            for (int attempt = 0; attempt < 100 && nrels < needed; attempt++) {
                /* Select ss primes from FB (skip p=2 for CRT simplicity) */
                int *sidx = malloc(sizeof(int) * ss);
                mpz_set_ui(crt_x, 0);
                mpz_set_ui(crt_mod, 1);
                int ok = 1;

                for (int i = 0; i < ss && ok; i++) {
                    int idx = 1 + ((attempt * 7 + i * 13 + ss * 3) % (fb_size - 1));
                    sidx[i] = idx;
                    mpz_set_ui(p_mpz, fb[idx]);

                    /* Check coprimality with current modulus */
                    mpz_t g;
                    mpz_init(g);
                    mpz_gcd(g, crt_mod, p_mpz);
                    if (mpz_cmp_ui(g, 1) != 0) { mpz_clear(g); ok = 0; break; }
                    mpz_clear(g);

                    /* CRT step */
                    mpz_sub(tmp_crt, roots[idx], crt_x);
                    mpz_mod(tmp_crt, tmp_crt, p_mpz);
                    if (!mpz_invert(inv_z, crt_mod, p_mpz)) { ok = 0; break; }
                    mpz_mul(tmp_crt, tmp_crt, inv_z);
                    mpz_mod(tmp_crt, tmp_crt, p_mpz);
                    mpz_addmul(crt_x, crt_mod, tmp_crt);
                    mpz_mul(crt_mod, crt_mod, p_mpz);
                    mpz_mod(crt_x, crt_x, crt_mod);
                }

                if (!ok) { free(sidx); continue; }

                /* Find k so that crt_x + k*crt_mod ≈ sqrt(N) */
                mpz_sub(tmp_crt, sqrt_N, crt_x);
                mpz_fdiv_q(tmp_crt, tmp_crt, crt_mod);

                mpz_set(base_x, crt_x);
                mpz_addmul(base_x, tmp_crt, crt_mod);

                /* Search near base_x */
                mpz_set(test_x, base_x);
                mpz_submul_ui(test_x, crt_mod, 200);

                for (int delta = 0; delta < 400 && nrels < max_rels; delta++) {
                    mpz_mul(test_val, test_x, test_x);
                    mpz_sub(test_val, test_val, N);
                    if (mpz_sgn(test_val) != 0) {
                        int sign = (mpz_sgn(test_val) < 0) ? 1 : 0;
                        if (trial_factor(test_val, fb, fb_size, ev)) {
                            /* Check for duplicates */
                            int dup = 0;
                            for (int r = 0; r < nrels; r++)
                                if (mpz_cmp(rels[r].x, test_x) == 0) { dup = 1; break; }
                            if (!dup) {
                                mpz_set(rels[nrels].x, test_x);
                                memcpy(rels[nrels].exponents, ev, sizeof(int)*fb_size);
                                rels[nrels].sign = sign;
                                nrels++;
                            }
                        }
                    }
                    mpz_add(test_x, test_x, crt_mod);
                }

                free(sidx);
            }
        }

        free(ev);
        mpz_clears(crt_x, crt_mod, tmp_crt, inv_z, base_x, test_x, test_val, NULL);
        printf("  After lattice-CRT: %d relations total\n", nrels);
    }

    /* ------------------------------------------------------------ */
    /* Phase 3: LLL lattice for candidate generation                 */
    /* Build lattice where short vectors give x with x^2-N divisible */
    /* by many FB primes. We use the approach:                       */
    /*   Row 0: (M, 0, 0, ..., 0)  where M = prod of subset primes  */
    /*   Row i: (r_i, 0, ..., p_i, ..., 0) with r_i in col 0       */
    /* Short vectors give small linear combinations.                 */
    /* ------------------------------------------------------------ */
    if (nrels < needed && fb_size >= 4) {
        printf("  LLL lattice phase...\n");

        int ldim = (fb_size < 12) ? fb_size : 12;
        int mdim = ldim + 1;  /* matrix is mdim x mdim */

        mpz_t *lat = malloc(sizeof(mpz_t) * mdim * mdim);
        for (int i = 0; i < mdim * mdim; i++) mpz_init_set_ui(lat[i], 0);

        /* Compute M = product of first ldim FB primes (skipping 2) */
        mpz_t M;
        mpz_init_set_ui(M, 1);
        for (int i = 0; i < ldim && i + 1 < fb_size; i++)
            mpz_mul_ui(M, M, fb[i + 1]);

        /* Row 0: (M, 0, ..., 0) */
        mpz_set(MAT(lat, 0, 0, mdim), M);

        /* Row i (1..ldim): (r_i, 0, ..., p_i, ..., 0) */
        for (int i = 0; i < ldim; i++) {
            int fi = i + 1;  /* skip p=2 */
            if (fi >= fb_size) break;
            mpz_set(MAT(lat, i+1, 0, mdim), roots[fi]);
            mpz_set_ui(MAT(lat, i+1, i+1, mdim), fb[fi]);
        }

        printf("  Running LLL on %d x %d lattice...\n", mdim, mdim);
        lll_reduce(lat, mdim, mdim);
        printf("  LLL complete.\n");

        /* Extract x candidates from short vectors (column 0) */
        int *ev = calloc(fb_size, sizeof(int));
        mpz_t sx, sv, adj;
        mpz_inits(sx, sv, adj, NULL);

        for (int row = 0; row < mdim && nrels < max_rels; row++) {
            mpz_abs(sx, MAT(lat, row, 0, mdim));
            if (mpz_sgn(sx) == 0) continue;

            /* Adjust to be near sqrt(N): find k such that sx + k*M ≈ sqrt(N) */
            mpz_sub(adj, sqrt_N, sx);
            mpz_fdiv_q(adj, adj, M);
            mpz_addmul(sx, adj, M);

            /* Search near this x */
            mpz_t tx, tv;
            mpz_init(tx); mpz_init(tv);
            mpz_set(tx, sx);
            mpz_submul_ui(tx, M, 20);

            for (int d = 0; d < 40 && nrels < max_rels; d++) {
                mpz_mul(tv, tx, tx);
                mpz_sub(tv, tv, N);
                if (mpz_sgn(tv) != 0) {
                    int sign = (mpz_sgn(tv) < 0) ? 1 : 0;
                    if (trial_factor(tv, fb, fb_size, ev)) {
                        int dup = 0;
                        for (int r = 0; r < nrels; r++)
                            if (mpz_cmp(rels[r].x, tx) == 0) { dup = 1; break; }
                        if (!dup) {
                            mpz_set(rels[nrels].x, tx);
                            memcpy(rels[nrels].exponents, ev, sizeof(int)*fb_size);
                            rels[nrels].sign = sign;
                            nrels++;
                        }
                    }
                }
                mpz_add(tx, tx, M);
            }
            mpz_clear(tx); mpz_clear(tv);
        }

        free(ev);
        mpz_clears(sx, sv, adj, M, NULL);
        for (int i = 0; i < mdim * mdim; i++) mpz_clear(lat[i]);
        free(lat);
        printf("  After LLL phase: %d relations total\n", nrels);
    }

    printf("  Collected %d relations (needed %d)\n", nrels, needed);

    if (nrels < needed) {
        fprintf(stderr, "  Not enough relations (%d < %d)\n", nrels, needed);
        goto fail;
    }

    /* ------------------------------------------------------------ */
    /* Phase 4: Gaussian elimination + factor extraction              */
    /* ------------------------------------------------------------ */
    {
        int ndeps = 0;
        dependency_t *deps = gauss_find_deps(rels, nrels, fb_size, &ndeps);
        printf("  Gaussian elimination found %d dependencies\n", ndeps);

        for (int d = 0; d < ndeps; d++) {
            if (try_dependency(N, rels, fb, fb_size, &deps[d], factor_out)) {
                printf("  Dependency %d yielded a factor!\n", d);
                for (int i = 0; i < ndeps; i++) free(deps[i].indices);
                free(deps);
                goto success;
            }
        }

        for (int i = 0; i < ndeps; i++) free(deps[i].indices);
        free(deps);
        printf("  No dependency yielded a nontrivial factor\n");
    }

fail:
    for (int i = 0; i < fb_size; i++) mpz_clear(roots[i]);
    free(roots);
    mpz_clears(p_mpz, n_mod_p, sqrt_N, NULL);
    for (int i = 0; i < max_rels; i++) { mpz_clear(rels[i].x); free(rels[i].exponents); }
    free(rels);
    free(fb);
    return 0;

success:
    for (int i = 0; i < fb_size; i++) mpz_clear(roots[i]);
    free(roots);
    mpz_clears(p_mpz, n_mod_p, sqrt_N, NULL);
    for (int i = 0; i < max_rels; i++) { mpz_clear(rels[i].x); free(rels[i].exponents); }
    free(rels);
    free(fb);
    return 1;
}

/* ================================================================== */
/* Pollard's rho (for comparison / fallback)                           */
/* ================================================================== */
static int pollard_rho(const mpz_t N, mpz_t factor_out) {
    mpz_t x, y, d, c, tmp;
    mpz_inits(x, y, d, c, tmp, NULL);
    gmp_randstate_t st;
    gmp_randinit_default(st);
    gmp_randseed_ui(st, 42);

    for (int a = 0; a < 50; a++) {
        mpz_urandomm(x, st, N);
        mpz_set(y, x);
        mpz_urandomm(c, st, N);
        if (mpz_sgn(c) == 0) mpz_set_ui(c, 1);
        mpz_set_ui(d, 1);

        while (mpz_cmp_ui(d, 1) == 0) {
            mpz_mul(x, x, x); mpz_add(x, x, c); mpz_mod(x, x, N);
            mpz_mul(y, y, y); mpz_add(y, y, c); mpz_mod(y, y, N);
            mpz_mul(y, y, y); mpz_add(y, y, c); mpz_mod(y, y, N);
            mpz_sub(tmp, x, y); mpz_abs(tmp, tmp);
            mpz_gcd(d, tmp, N);
        }
        if (mpz_cmp(d, N) != 0) {
            mpz_set(factor_out, d);
            mpz_clears(x, y, d, c, tmp, NULL);
            gmp_randclear(st);
            return 1;
        }
    }
    mpz_clears(x, y, d, c, tmp, NULL);
    gmp_randclear(st);
    return 0;
}

/* ================================================================== */
/* Main                                                                */
/* ================================================================== */
int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N> [rho]\n", argv[0]);
        return 1;
    }

    mpz_t N, factor;
    mpz_inits(N, factor, NULL);
    if (mpz_set_str(N, argv[1], 10) != 0) {
        fprintf(stderr, "Invalid number: %s\n", argv[1]);
        return 1;
    }

    if (mpz_cmp_ui(N, 4) < 0) {
        gmp_printf("%Zd is too small\n", N);
        return 1;
    }
    if (mpz_even_p(N)) { printf("2\n"); mpz_clears(N, factor, NULL); return 0; }
    if (mpz_perfect_power_p(N)) {
        mpz_t root;
        mpz_init(root);
        for (unsigned long e = 2; e <= 64; e++)
            if (mpz_root(root, N, e)) {
                gmp_printf("%Zd (perfect %lu-th power)\n", root, e);
                mpz_clear(root); mpz_clears(N, factor, NULL); return 0;
            }
        mpz_clear(root);
    }

    int use_rho = (argc >= 3 && strcmp(argv[2], "rho") == 0);

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    int found;
    if (use_rho) {
        printf("Pollard's rho:\n");
        found = pollard_rho(N, factor);
    } else {
        gmp_printf("Lattice factoring: N = %Zd (%zd digits)\n", N, mpz_sizeinbase(N, 10));
        found = lattice_factor(N, factor);
        if (!found) {
            printf("Lattice failed, falling back to Pollard's rho...\n");
            found = pollard_rho(N, factor);
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)/1e9;

    if (found) {
        mpz_t cofactor;
        mpz_init(cofactor);
        mpz_divexact(cofactor, N, factor);
        gmp_printf("Factor: %Zd\n", factor);
        gmp_printf("Cofactor: %Zd\n", cofactor);
        mpz_clear(cofactor);
    } else {
        printf("Failed to find a factor.\n");
    }
    printf("Time: %.4f seconds\n", elapsed);

    mpz_clears(N, factor, NULL);
    return found ? 0 : 1;
}
