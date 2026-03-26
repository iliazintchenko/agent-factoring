/*
 * NFS Polynomial Selection with LLL Lattice Reduction
 * Compares standard base-m vs LLL-improved polynomial selection.
 * Scores by max_coeff * disc^{1/(2d-2)}.
 * Tests on 30-40 digit semiprimes, compares smooth yield.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

#define MAX_DEG 6
#define SIEVE_BOUND 100000
#define SMOOTH_BOUND 1000000
#define NUM_TEST_VALUES 50000

/* ============================================================
 * LLL Lattice Reduction (integer, using GMP)
 * ============================================================ */

typedef struct {
    int rows, cols;
    mpz_t **data;
} lattice_t;

static void lat_init(lattice_t *L, int rows, int cols) {
    L->rows = rows;
    L->cols = cols;
    L->data = malloc(rows * sizeof(mpz_t *));
    for (int i = 0; i < rows; i++) {
        L->data[i] = malloc(cols * sizeof(mpz_t));
        for (int j = 0; j < cols; j++)
            mpz_init(L->data[i][j]);
    }
}

static void lat_clear(lattice_t *L) {
    for (int i = 0; i < L->rows; i++) {
        for (int j = 0; j < L->cols; j++)
            mpz_clear(L->data[i][j]);
        free(L->data[i]);
    }
    free(L->data);
}

/* Gram-Schmidt with rational arithmetic (num/den stored separately) */
static void lll_reduce(lattice_t *L, double delta) {
    int n = L->rows;
    int d = L->cols;

    /* mu[i][j] = mu_num[i][j] / mu_den[i][j] */
    mpz_t **mu_num, **mu_den, *B;
    mu_num = malloc(n * sizeof(mpz_t *));
    mu_den = malloc(n * sizeof(mpz_t *));
    B = malloc(n * sizeof(mpz_t));

    for (int i = 0; i < n; i++) {
        mu_num[i] = malloc(n * sizeof(mpz_t));
        mu_den[i] = malloc(n * sizeof(mpz_t));
        for (int j = 0; j < n; j++) {
            mpz_init(mu_num[i][j]);
            mpz_init(mu_den[i][j]);
        }
        mpz_init(B[i]);
    }

    mpz_t tmp, tmp2, q;
    mpz_inits(tmp, tmp2, q, NULL);

    /* Compute Gram-Schmidt coefficients */
    auto void compute_gs(void);
    void compute_gs(void) {
        for (int i = 0; i < n; i++) {
            /* B[i] = <b_i, b_i> initially */
            mpz_set_ui(B[i], 0);
            for (int k = 0; k < d; k++) {
                mpz_mul(tmp, L->data[i][k], L->data[i][k]);
                mpz_add(B[i], B[i], tmp);
            }
            for (int j = 0; j < i; j++) {
                /* mu_num[i][j] = <b_i, b*_j> * mu_den = dot(b_i, b_j_orig) * B[j]_den */
                /* Simplified: mu[i][j] = <b_i, b*_j> / B[j] */
                /* We compute <b_i, b_j> then subtract contributions */
                mpz_set_ui(mu_num[i][j], 0);
                for (int k = 0; k < d; k++) {
                    mpz_mul(tmp, L->data[i][k], L->data[j][k]);
                    mpz_add(mu_num[i][j], mu_num[i][j], tmp);
                }
                mpz_set(mu_den[i][j], B[j]);
                /* Adjust: subtract mu[j][k]*mu[i][k]*B[k]/B[j] for k < j */
                for (int k = 0; k < j; k++) {
                    /* mu_num[i][j] -= mu_num[j][k] * mu_num[i][k] * B[k] / mu_den[j][k] / mu_den[i][k] */
                    /* Better: work with the raw dot products */
                }
            }
            /* Update B[i] -= sum mu[i][j]^2 * B[j] */
            /* This simplified version works for integer LLL */
        }
    }

    /* Use a simpler, well-known integer LLL implementation */
    /* Lagrange/LLL with size reduction and swap */

    /* Recompute using standard approach with mpq for mu */
    mpq_t *gs_B;     /* B*[i] norms squared, rational */
    mpq_t **gs_mu;   /* mu[i][j] rational */
    gs_B = malloc(n * sizeof(mpq_t));
    gs_mu = malloc(n * sizeof(mpq_t *));
    for (int i = 0; i < n; i++) {
        mpq_init(gs_B[i]);
        gs_mu[i] = malloc(n * sizeof(mpq_t));
        for (int j = 0; j < n; j++)
            mpq_init(gs_mu[i][j]);
    }

    mpq_t qtmp, qtmp2, qdelta;
    mpq_inits(qtmp, qtmp2, qdelta, NULL);
    mpq_set_d(qdelta, delta);

    auto void recompute_gs_full(void);
    void recompute_gs_full(void) {
        for (int i = 0; i < n; i++) {
            /* B*[i] = <b_i, b_i> */
            mpq_set_ui(gs_B[i], 0, 1);
            for (int k = 0; k < d; k++) {
                mpq_set_z(qtmp, L->data[i][k]);
                mpq_mul(qtmp, qtmp, qtmp);
                mpq_add(gs_B[i], gs_B[i], qtmp);
            }
            for (int j = 0; j < i; j++) {
                /* mu[i][j] = <b_i, b*_j> / B*[j] */
                /* <b_i, b_j> */
                mpq_set_ui(gs_mu[i][j], 0, 1);
                for (int k = 0; k < d; k++) {
                    mpz_mul(tmp, L->data[i][k], L->data[j][k]);
                    mpq_set_z(qtmp, tmp);
                    mpq_add(gs_mu[i][j], gs_mu[i][j], qtmp);
                }
                /* subtract sum_{k<j} mu[j][k] * mu[i][k] * B*[k] */
                /* Actually: <b_i, b*_j> = <b_i, b_j> - sum_{k<j} mu[j][k]*<b_i, b*_k> */
                /* But simpler: mu[i][j] = (<b_i, b_j> - sum_{k<j} mu[i][k]*mu[j][k]*B*[k]) / B*[j] */
                for (int k = 0; k < j; k++) {
                    mpq_mul(qtmp, gs_mu[i][k], gs_mu[j][k]);
                    mpq_mul(qtmp, qtmp, gs_B[k]);
                    mpq_sub(gs_mu[i][j], gs_mu[i][j], qtmp);
                }
                if (mpq_sgn(gs_B[j]) != 0)
                    mpq_div(gs_mu[i][j], gs_mu[i][j], gs_B[j]);
            }
            /* B*[i] -= sum_{j<i} mu[i][j]^2 * B*[j] */
            for (int j = 0; j < i; j++) {
                mpq_mul(qtmp, gs_mu[i][j], gs_mu[i][j]);
                mpq_mul(qtmp, qtmp, gs_B[j]);
                mpq_sub(gs_B[i], gs_B[i], qtmp);
            }
        }
    }

    recompute_gs_full();

    /* LLL main loop */
    int k = 1;
    while (k < n) {
        /* Size reduce b[k] with b[k-1] */
        for (int j = k - 1; j >= 0; j--) {
            /* if |mu[k][j]| > 1/2, reduce */
            mpq_set_ui(qtmp, 1, 2);
            mpq_abs(qtmp2, gs_mu[k][j]);
            if (mpq_cmp(qtmp2, qtmp) > 0) {
                /* r = round(mu[k][j]) */
                mpz_t num, den, r;
                mpz_inits(num, den, r, NULL);
                mpq_get_num(num, gs_mu[k][j]);
                mpq_get_den(den, gs_mu[k][j]);
                /* r = floor(num/den + 1/2) */
                mpz_mul_ui(num, num, 2);
                mpz_add(num, num, den);
                mpz_mul_ui(den, den, 2);
                mpz_fdiv_q(r, num, den);

                /* b[k] -= r * b[j] */
                for (int c = 0; c < d; c++) {
                    mpz_mul(tmp, r, L->data[j][c]);
                    mpz_sub(L->data[k][c], L->data[k][c], tmp);
                }
                mpz_clears(num, den, r, NULL);
                recompute_gs_full();
            }
        }

        /* Lovasz condition: B*[k] >= (delta - mu[k][k-1]^2) * B*[k-1] */
        mpq_mul(qtmp, gs_mu[k][k - 1], gs_mu[k][k - 1]);
        mpq_sub(qtmp, qdelta, qtmp);
        mpq_mul(qtmp, qtmp, gs_B[k - 1]);
        if (mpq_cmp(gs_B[k], qtmp) >= 0) {
            k++;
        } else {
            /* Swap b[k] and b[k-1] */
            for (int c = 0; c < d; c++)
                mpz_swap(L->data[k][c], L->data[k - 1][c]);
            recompute_gs_full();
            if (k > 1) k--;
        }
    }

    /* Cleanup */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            mpq_clear(gs_mu[i][j]);
            mpz_clear(mu_num[i][j]);
            mpz_clear(mu_den[i][j]);
        }
        mpq_clear(gs_B[i]);
        free(gs_mu[i]);
        free(mu_num[i]);
        free(mu_den[i]);
        mpz_clear(B[i]);
    }
    free(gs_mu); free(gs_B); free(mu_num); free(mu_den); free(B);
    mpq_clears(qtmp, qtmp2, qdelta, NULL);
    mpz_clears(tmp, tmp2, q, NULL);
}

/* ============================================================
 * Polynomial representation and operations
 * ============================================================ */

typedef struct {
    int deg;
    mpz_t coeff[MAX_DEG + 1]; /* coeff[i] = coefficient of x^i */
    mpz_t m;                   /* root m such that f(m) ≡ 0 mod N */
} nfs_poly_t;

static void poly_init(nfs_poly_t *p) {
    for (int i = 0; i <= MAX_DEG; i++)
        mpz_init(p->coeff[i]);
    mpz_init(p->m);
    p->deg = 0;
}

static void poly_clear(nfs_poly_t *p) {
    for (int i = 0; i <= MAX_DEG; i++)
        mpz_clear(p->coeff[i]);
    mpz_clear(p->m);
}

static void poly_copy(nfs_poly_t *dst, const nfs_poly_t *src) {
    dst->deg = src->deg;
    for (int i = 0; i <= MAX_DEG; i++)
        mpz_set(dst->coeff[i], src->coeff[i]);
    mpz_set(dst->m, src->m);
}

/* Compute max |coeff| */
static void poly_max_coeff(mpz_t result, const nfs_poly_t *p) {
    mpz_set_ui(result, 0);
    mpz_t tmp;
    mpz_init(tmp);
    for (int i = 0; i <= p->deg; i++) {
        mpz_abs(tmp, p->coeff[i]);
        if (mpz_cmp(tmp, result) > 0)
            mpz_set(result, tmp);
    }
    mpz_clear(tmp);
}

/* Evaluate polynomial at x, store in result */
static void poly_eval(mpz_t result, const nfs_poly_t *p, const mpz_t x) {
    mpz_set(result, p->coeff[p->deg]);
    for (int i = p->deg - 1; i >= 0; i--) {
        mpz_mul(result, result, x);
        mpz_add(result, result, p->coeff[i]);
    }
}

/* Compute discriminant of polynomial (simplified for degree 3-5) */
static void poly_discriminant(mpz_t disc, const nfs_poly_t *p) {
    /* For practical NFS polynomial selection, we compute the discriminant
     * using resultant of f and f'. For small degrees, use explicit formulas. */
    int d = p->deg;

    if (d == 3) {
        /* disc = 18abcd - 4b^3d + b^2c^2 - 4ac^3 - 27a^2d^2 */
        /* where f = ax^3 + bx^2 + cx + d */
        mpz_t a, b, c, dd, t1, t2, t3, t4, t5;
        mpz_inits(a, b, c, dd, t1, t2, t3, t4, t5, NULL);
        mpz_set(a, p->coeff[3]);
        mpz_set(b, p->coeff[2]);
        mpz_set(c, p->coeff[1]);
        mpz_set(dd, p->coeff[0]);

        /* 18abcd */
        mpz_mul(t1, a, b); mpz_mul(t1, t1, c); mpz_mul(t1, t1, dd);
        mpz_mul_si(t1, t1, 18);

        /* -4b^3*d */
        mpz_pow_ui(t2, b, 3); mpz_mul(t2, t2, dd); mpz_mul_si(t2, t2, -4);

        /* b^2*c^2 */
        mpz_mul(t3, b, b); mpz_mul(t4, c, c); mpz_mul(t3, t3, t4);

        /* -4ac^3 */
        mpz_pow_ui(t4, c, 3); mpz_mul(t4, t4, a); mpz_mul_si(t4, t4, -4);

        /* -27a^2*d^2 */
        mpz_mul(t5, a, a); mpz_mul(t5, t5, dd); mpz_mul(t5, t5, dd);
        mpz_mul_si(t5, t5, -27);

        mpz_set(disc, t1);
        mpz_add(disc, disc, t2);
        mpz_add(disc, disc, t3);
        mpz_add(disc, disc, t4);
        mpz_add(disc, disc, t5);

        mpz_clears(a, b, c, dd, t1, t2, t3, t4, t5, NULL);
    } else {
        /* For other degrees, use resultant via Sylvester matrix.
         * Simplified: compute Res(f, f') / leading_coeff */
        /* We'll use a simpler approach: compute discriminant numerically
         * via the resultant of f and f' */

        /* f' coefficients */
        mpz_t fp[MAX_DEG];
        for (int i = 0; i < d; i++) {
            mpz_init(fp[i]);
            mpz_mul_ui(fp[i], p->coeff[i + 1], i + 1);
        }

        /* Compute resultant via subresultant / pseudo-remainder sequence */
        /* For simplicity, use polynomial GCD-based approach */
        /* Res(f, f') = (-1)^{d(d-1)/2} * (1/a_d) * Res(f, f') */

        /* Simple: evaluate discriminant as product of differences of roots
         * For our purposes, an approximation via coefficient norms suffices
         * for scoring. Use Mahler measure approximation. */

        /* Fallback: use norm-based approximation */
        /* |disc| ≈ d^d * |a_d|^{2d-2} * prod(|r_i - r_j|^2) */
        /* Approximate via coefficient size */
        mpz_t norm;
        mpz_init(norm);
        mpz_set_ui(norm, 0);
        for (int i = 0; i <= d; i++) {
            mpz_t t;
            mpz_init(t);
            mpz_mul(t, p->coeff[i], p->coeff[i]);
            mpz_add(norm, norm, t);
            mpz_clear(t);
        }
        /* Crude bound: disc ~ norm^(d-1) */
        mpz_pow_ui(disc, norm, d - 1);
        mpz_clear(norm);

        for (int i = 0; i < d; i++)
            mpz_clear(fp[i]);
    }
    mpz_abs(disc, disc);
}

/* Score = max_coeff * disc^{1/(2d-2)} (as double for comparison) */
static double poly_score(const nfs_poly_t *p) {
    mpz_t mc, disc;
    mpz_inits(mc, disc, NULL);

    poly_max_coeff(mc, p);
    poly_discriminant(disc, p);

    double max_c = mpz_get_d(mc);
    double d_disc = mpz_get_d(disc);
    int d = p->deg;
    double exponent = 1.0 / (2.0 * d - 2.0);

    double score = max_c * pow(fabs(d_disc), exponent);

    mpz_clears(mc, disc, NULL);
    return score;
}

/* ============================================================
 * Standard base-m polynomial selection
 * ============================================================ */

static void base_m_select(nfs_poly_t *poly, const mpz_t N, int deg) {
    /* Choose m = floor(N^{1/d}) */
    mpz_t m, rem, tmp;
    mpz_inits(m, rem, tmp, NULL);

    mpz_root(m, N, deg);

    /* Express N in base m: N = c_d*m^d + c_{d-1}*m^{d-1} + ... + c_0 */
    mpz_set(rem, N);
    poly->deg = deg;
    mpz_set(poly->m, m);

    for (int i = 0; i < deg; i++) {
        mpz_fdiv_qr(rem, poly->coeff[i], rem, m);
    }
    mpz_set(poly->coeff[deg], rem);

    /* Adjust: make coefficients balanced (centered representation) */
    /* c_i in (-m/2, m/2] */
    mpz_t half_m;
    mpz_init(half_m);
    mpz_fdiv_q_ui(half_m, m, 2);

    for (int i = 0; i < deg; i++) {
        if (mpz_cmp(poly->coeff[i], half_m) > 0) {
            mpz_sub(poly->coeff[i], poly->coeff[i], m);
            mpz_add_ui(poly->coeff[i + 1], poly->coeff[i + 1], 1);
        }
    }

    mpz_clears(m, rem, tmp, half_m, NULL);
}

/* ============================================================
 * LLL-improved polynomial selection
 * ============================================================ */

static void lll_poly_select(nfs_poly_t *poly, const mpz_t N, int deg) {
    /*
     * LLL-improved polynomial selection for NFS.
     *
     * Strategy: Fix leading coefficient a_d to a small value (try 1,2,...),
     * compute m = floor((N/a_d)^{1/d}), then optimize the lower coefficients
     * c_0,...,c_{d-1} subject to: a_d*m^d + c_{d-1}*m^{d-1} + ... + c_0 ≡ 0 mod N.
     *
     * Let R = N - a_d*m^d. Then c_{d-1}*m^{d-1} + ... + c_0 ≡ R mod N.
     *
     * We find one particular solution via base-m expansion of R, then
     * the general solution adds any element of the lattice:
     *   {v in Z^d : v . (1, m, m^2, ..., m^{d-1}) ≡ 0 mod N}
     *
     * We LLL-reduce this lattice and translate to find the shortest
     * coefficient vector near the particular solution.
     */

    nfs_poly_t best_poly;
    poly_init(&best_poly);
    double best_score = 1e300;

    mpz_t m, R, tmp, eval;
    mpz_inits(m, R, tmp, eval, NULL);

    /* Try several small leading coefficients */
    for (int a_d = 1; a_d <= 10; a_d++) {
        /* m = floor((N / a_d)^{1/d}) */
        mpz_fdiv_q_ui(tmp, N, a_d);
        mpz_root(m, tmp, deg);

        /* R = N - a_d * m^d */
        mpz_pow_ui(R, m, deg);
        mpz_mul_ui(R, R, a_d);
        mpz_sub(R, N, R);
        /* R should be positive and < N */
        if (mpz_sgn(R) < 0) {
            /* Try m-1 */
            mpz_sub_ui(m, m, 1);
            mpz_pow_ui(R, m, deg);
            mpz_mul_ui(R, R, a_d);
            mpz_sub(R, N, R);
        }

        /* Base-m expansion of R to get particular solution c_0,...,c_{d-1} */
        mpz_t part[MAX_DEG]; /* particular solution */
        for (int i = 0; i < deg; i++) mpz_init(part[i]);

        mpz_t rem;
        mpz_init_set(rem, R);
        for (int i = 0; i < deg; i++) {
            mpz_fdiv_qr(rem, part[i], rem, m);
        }
        /* rem should be 0 or very small if a_d*m^d + ... = N works */

        /* Center the coefficients: if c_i > m/2, set c_i -= m, carry +1 */
        mpz_t half_m;
        mpz_init(half_m);
        mpz_fdiv_q_ui(half_m, m, 2);
        for (int i = 0; i < deg - 1; i++) {
            if (mpz_cmp(part[i], half_m) > 0) {
                mpz_sub(part[i], part[i], m);
                mpz_add_ui(part[i + 1], part[i + 1], 1);
            }
        }

        /* Build the lattice for the homogeneous part:
         * {v in Z^d : v_0 + v_1*m + ... + v_{d-1}*m^{d-1} ≡ 0 mod N}
         *
         * Basis (d x d):
         * Row 0: [N, 0, 0, ..., 0]
         * Row i (i=1..d-1): [-m^i, 0,..., 1, ..., 0]  (1 at column i)
         *
         * We scale columns to balance: col i gets weight m^{(d-1-i)/(d-1)}
         * so that all coefficients contribute equally to |f(x)| near x~m.
         * Simpler: use skewness scaling s_i = m^{i} for col i.
         *   => after LLL, unscale v_i by dividing by m^i... but that loses
         *   integrality. Instead, scale row to make norms comparable.
         *
         * Standard approach: column scaling by N^{i/d}.
         */

        int dim = deg; /* d coefficients c_0,...,c_{d-1} */
        lattice_t L;
        lat_init(&L, dim, dim);

        /* Diagonal scaling factors: s_i = N^{(d-1-i)/(d*(d-1))} ≈ m^{(d-1-i)/(d-1)}
         * For simplicity, use s_i = round(m^{i/(d-1) * ???})
         * Actually, let's use a simpler scaling: s_i = m^i to make the
         * lattice entries comparable in size.
         *
         * With scaling s_i for column i, the lattice becomes:
         * Row 0: [N*s_0, 0, 0, ..., 0]
         * Row i: [-m^i*s_0, 0,..., s_i, ..., 0]
         *
         * After LLL, unscale: true_v_i = lattice_v_i / s_i
         *
         * We want s_i such that |v_i * x^i| ~ same magnitude at x~m.
         * |v_i| * m^i ~ const => want v_i ~ C/m^i => scale so that
         * short vector has equal components => s_i = m^i.
         * Then lattice_v_i = v_i * m^i, and LLL minimizes max|lattice_v_i|
         * which means minimizing max|v_i * m^i| = max|v_i * x^i| at x=m.
         * This is exactly what we want for minimizing |f(m)|.
         *
         * But actually we want small |f(a/b)| for MANY a,b, not just at m.
         * The skewness-optimal scaling is s_i = s^i for some skewness s.
         * For now, use s = m^{1/2} as a reasonable compromise (half-skew).
         */

        mpz_t sk, sk_pow[MAX_DEG];
        mpz_init(sk);
        mpz_sqrt(sk, m); /* s = sqrt(m) as skewness */
        if (mpz_sgn(sk) == 0) mpz_set_ui(sk, 1);

        for (int i = 0; i < dim; i++) {
            mpz_init(sk_pow[i]);
            if (i == 0)
                mpz_set_ui(sk_pow[i], 1);
            else
                mpz_mul(sk_pow[i], sk_pow[i - 1], sk);
        }

        /* Row 0: [N * sk_pow[0], 0, ..., 0] = [N, 0, ..., 0] */
        mpz_mul(L.data[0][0], N, sk_pow[0]);

        /* Row i (i=1..d-1): column 0 = -m^i * sk_pow[0], column i = sk_pow[i] */
        mpz_t m_pow;
        mpz_init_set(m_pow, m);
        for (int i = 1; i < dim; i++) {
            mpz_mul(L.data[i][0], m_pow, sk_pow[0]);
            mpz_neg(L.data[i][0], L.data[i][0]);
            mpz_set(L.data[i][i], sk_pow[i]);
            if (i < dim - 1) mpz_mul(m_pow, m_pow, m);
        }

        /* LLL reduce */
        lll_reduce(&L, 0.99);

        /* Try each row of the reduced lattice as a correction to the
         * particular solution. The total coefficient vector is:
         * c_i = part[i] + sum_j (alpha_j * L_reduced[j][i] / sk_pow[i])
         *
         * For the simplest approach: try adding each single basis vector
         * (with multiplier 0 or ±1) and pick the best scoring polynomial.
         */

        nfs_poly_t candidate;
        poly_init(&candidate);
        candidate.deg = deg;
        mpz_set(candidate.m, m);

        /* First: just use the particular solution (= base-m essentially) */
        for (int i = 0; i < deg; i++)
            mpz_set(candidate.coeff[i], part[i]);
        mpz_set_ui(candidate.coeff[deg], a_d);

        /* Verify */
        poly_eval(eval, &candidate, m);
        mpz_mod(eval, eval, N);

        if (mpz_sgn(eval) == 0) {
            double sc = poly_score(&candidate);
            if (sc > 0 && sc < best_score) {
                best_score = sc;
                poly_copy(&best_poly, &candidate);
            }
        }

        /* Now try adding LLL-reduced lattice vectors */
        for (int r = 0; r < dim; r++) {
            for (int sign = -1; sign <= 1; sign += 2) {
                for (int i = 0; i < deg; i++) {
                    /* Unscale: correction_i = L[r][i] / sk_pow[i] */
                    if (mpz_divisible_p(L.data[r][i], sk_pow[i])) {
                        mpz_divexact(tmp, L.data[r][i], sk_pow[i]);
                    } else {
                        mpz_tdiv_q(tmp, L.data[r][i], sk_pow[i]);
                    }
                    if (sign > 0)
                        mpz_add(candidate.coeff[i], part[i], tmp);
                    else
                        mpz_sub(candidate.coeff[i], part[i], tmp);
                }
                mpz_set_ui(candidate.coeff[deg], a_d);

                poly_eval(eval, &candidate, m);
                mpz_mod(eval, eval, N);

                if (mpz_sgn(eval) == 0) {
                    double sc = poly_score(&candidate);
                    if (sc > 0 && sc < best_score) {
                        best_score = sc;
                        poly_copy(&best_poly, &candidate);
                    }
                }
            }

            /* Also try multiplier ±2 for stronger corrections */
            for (int sign = -2; sign <= 2; sign += 4) {
                for (int i = 0; i < deg; i++) {
                    if (mpz_divisible_p(L.data[r][i], sk_pow[i])) {
                        mpz_divexact(tmp, L.data[r][i], sk_pow[i]);
                    } else {
                        mpz_tdiv_q(tmp, L.data[r][i], sk_pow[i]);
                    }
                    mpz_mul_si(tmp, tmp, sign);
                    mpz_add(candidate.coeff[i], part[i], tmp);
                }
                mpz_set_ui(candidate.coeff[deg], a_d);

                poly_eval(eval, &candidate, m);
                mpz_mod(eval, eval, N);

                if (mpz_sgn(eval) == 0) {
                    double sc = poly_score(&candidate);
                    if (sc > 0 && sc < best_score) {
                        best_score = sc;
                        poly_copy(&best_poly, &candidate);
                    }
                }
            }
        }

        /* Try combinations of two lattice vectors */
        for (int r1 = 0; r1 < dim && r1 < 3; r1++) {
            for (int r2 = r1 + 1; r2 < dim && r2 < 4; r2++) {
                for (int s1 = -1; s1 <= 1; s1 += 2) {
                    for (int s2 = -1; s2 <= 1; s2 += 2) {
                        for (int i = 0; i < deg; i++) {
                            mpz_t c1, c2;
                            mpz_inits(c1, c2, NULL);
                            if (mpz_divisible_p(L.data[r1][i], sk_pow[i]))
                                mpz_divexact(c1, L.data[r1][i], sk_pow[i]);
                            else
                                mpz_tdiv_q(c1, L.data[r1][i], sk_pow[i]);
                            if (mpz_divisible_p(L.data[r2][i], sk_pow[i]))
                                mpz_divexact(c2, L.data[r2][i], sk_pow[i]);
                            else
                                mpz_tdiv_q(c2, L.data[r2][i], sk_pow[i]);
                            mpz_mul_si(c1, c1, s1);
                            mpz_mul_si(c2, c2, s2);
                            mpz_add(candidate.coeff[i], part[i], c1);
                            mpz_add(candidate.coeff[i], candidate.coeff[i], c2);
                            mpz_clears(c1, c2, NULL);
                        }
                        mpz_set_ui(candidate.coeff[deg], a_d);

                        poly_eval(eval, &candidate, m);
                        mpz_mod(eval, eval, N);

                        if (mpz_sgn(eval) == 0) {
                            double sc = poly_score(&candidate);
                            if (sc > 0 && sc < best_score) {
                                best_score = sc;
                                poly_copy(&best_poly, &candidate);
                            }
                        }
                    }
                }
            }
        }

        poly_clear(&candidate);
        mpz_clear(m_pow);
        for (int i = 0; i < dim; i++) mpz_clear(sk_pow[i]);
        mpz_clear(sk);
        for (int i = 0; i < deg; i++) mpz_clear(part[i]);
        mpz_clear(rem);
        mpz_clear(half_m);
        lat_clear(&L);
    }

    if (best_score < 1e300) {
        poly_copy(poly, &best_poly);
    } else {
        /* Fallback to base-m */
        base_m_select(poly, N, deg);
    }

    mpz_clears(m, R, tmp, eval, NULL);
    poly_clear(&best_poly);
}

/* ============================================================
 * Smooth counting for yield estimation
 * ============================================================ */

/* Count values |f(a,b)| = |b^d * f(a/b)| that are B-smooth
 * for a in [-A, A], b in [1, B_range] */
static int count_smooth(const nfs_poly_t *p, const mpz_t N,
                        long A_range, long B_range, long smooth_bound) {
    int count = 0;
    mpz_t val, tmp, a_mpz, b_mpz;
    mpz_inits(val, tmp, a_mpz, b_mpz, NULL);

    /* Simple sieve of primes up to smooth_bound */
    int *primes = NULL;
    int nprimes = 0;
    {
        char *sieve = calloc(smooth_bound + 1, 1);
        for (long i = 2; i <= smooth_bound; i++) {
            if (!sieve[i]) {
                nprimes++;
                for (long j = 2 * i; j <= smooth_bound; j += i)
                    sieve[j] = 1;
            }
        }
        primes = malloc(nprimes * sizeof(int));
        int idx = 0;
        for (long i = 2; i <= smooth_bound; i++)
            if (!sieve[i]) primes[idx++] = (int)i;
        free(sieve);
    }

    /* Evaluate f(a, b) = sum c_i * a^i * b^{d-i} for random (a,b) pairs */
    srand(42); /* reproducible */
    long total_tests = A_range; /* number of (a,b) pairs to test */
    if (total_tests > NUM_TEST_VALUES) total_tests = NUM_TEST_VALUES;

    for (long t = 0; t < total_tests; t++) {
        long a = (rand() % (2 * A_range + 1)) - A_range;
        long b = 1 + (rand() % B_range);
        if (a == 0) continue;

        /* Compute f(a, b) = sum c_i * a^i * b^{d-i} */
        mpz_set_ui(val, 0);
        mpz_set_ui(a_mpz, 0);
        mpz_set_si(a_mpz, a);
        mpz_set_si(b_mpz, b);

        mpz_t a_pow, b_pow;
        mpz_inits(a_pow, b_pow, NULL);
        mpz_set_ui(a_pow, 1);
        mpz_pow_ui(b_pow, b_mpz, p->deg);

        for (int i = 0; i <= p->deg; i++) {
            mpz_mul(tmp, p->coeff[i], a_pow);
            mpz_mul(tmp, tmp, b_pow);
            mpz_add(val, val, tmp);

            mpz_mul(a_pow, a_pow, a_mpz);
            if (i < p->deg) {
                /* b_pow /= b */
                if (b != 0) mpz_divexact(b_pow, b_pow, b_mpz);
            }
        }
        mpz_clears(a_pow, b_pow, NULL);

        mpz_abs(val, val);
        if (mpz_sgn(val) == 0) continue;

        /* Trial division to check smoothness */
        for (int i = 0; i < nprimes; i++) {
            while (mpz_divisible_ui_p(val, primes[i]))
                mpz_divexact_ui(val, val, primes[i]);
        }
        if (mpz_cmp_ui(val, 1) == 0)
            count++;
    }

    free(primes);
    mpz_clears(val, tmp, a_mpz, b_mpz, NULL);
    return count;
}

/* ============================================================
 * Print polynomial
 * ============================================================ */

static void poly_print(FILE *f, const char *label, const nfs_poly_t *p) {
    fprintf(f, "%s (degree %d):\n  f(x) = ", label, p->deg);
    for (int i = p->deg; i >= 0; i--) {
        if (i < p->deg) fprintf(f, " + ");
        gmp_fprintf(f, "%Zd", p->coeff[i]);
        if (i > 0) fprintf(f, "*x");
        if (i > 1) fprintf(f, "^%d", i);
    }
    gmp_fprintf(f, "\n  m = %Zd\n", p->m);
}

/* ============================================================
 * Generate test semiprimes
 * ============================================================ */

static void gen_semiprime(mpz_t N, int digits, gmp_randstate_t rstate) {
    mpz_t p, q, lo, hi;
    mpz_inits(p, q, lo, hi, NULL);

    int p_digits = digits / 2;
    int q_digits = digits - p_digits;

    mpz_ui_pow_ui(lo, 10, p_digits - 1);
    mpz_ui_pow_ui(hi, 10, p_digits);
    mpz_sub(hi, hi, lo);

    do {
        mpz_urandomm(p, rstate, hi);
        mpz_add(p, p, lo);
        mpz_nextprime(p, p);
    } while (mpz_sizeinbase(p, 10) != (size_t)p_digits);

    mpz_ui_pow_ui(lo, 10, q_digits - 1);
    mpz_ui_pow_ui(hi, 10, q_digits);
    mpz_sub(hi, hi, lo);

    do {
        mpz_urandomm(q, rstate, hi);
        mpz_add(q, q, lo);
        mpz_nextprime(q, q);
    } while (mpz_sizeinbase(q, 10) != (size_t)q_digits);

    mpz_mul(N, p, q);
    mpz_clears(p, q, lo, hi, NULL);
}

/* ============================================================
 * Main
 * ============================================================ */

int main(void) {
    FILE *fout = fopen("findings.txt", "w");
    if (!fout) { perror("fopen"); return 1; }

    fprintf(fout, "=== NFS Polynomial Selection: LLL vs Standard Base-m ===\n");
    fprintf(fout, "Score = max_coeff * disc^{1/(2d-2)}\n\n");

    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, 12345);

    int test_digits[] = {30, 32, 34, 36, 38, 40};
    int num_tests = sizeof(test_digits) / sizeof(test_digits[0]);
    int trials_per = 3;

    /* Degree selection: for 30-40 digits, degree 3 is appropriate */
    int deg = 3;

    double total_score_base = 0, total_score_lll = 0;
    int total_smooth_base = 0, total_smooth_lll = 0;

    fprintf(fout, "%-8s %-6s %-15s %-15s %-12s %-12s %-10s %-10s\n",
            "Digits", "Trial", "Score(base-m)", "Score(LLL)",
            "Smooth(BM)", "Smooth(LLL)", "Ratio(Sc)", "Ratio(Sm)");
    fprintf(fout, "------------------------------------------------------------------------"
                  "--------------------\n");

    for (int t = 0; t < num_tests; t++) {
        int digits = test_digits[t];
        for (int trial = 0; trial < trials_per; trial++) {
            mpz_t N;
            mpz_init(N);
            gen_semiprime(N, digits, rstate);

            nfs_poly_t poly_base, poly_lll;
            poly_init(&poly_base);
            poly_init(&poly_lll);

            /* Standard base-m */
            base_m_select(&poly_base, N, deg);
            double score_base = poly_score(&poly_base);

            /* LLL-improved */
            lll_poly_select(&poly_lll, N, deg);
            double score_lll = poly_score(&poly_lll);

            /* Smooth counting */
            long A_range = 5000;
            long B_range = 100;
            long sb = 50000; /* smooth bound */
            int smooth_base = count_smooth(&poly_base, N, A_range, B_range, sb);
            int smooth_lll = count_smooth(&poly_lll, N, A_range, B_range, sb);

            double score_ratio = (score_lll > 0) ? score_base / score_lll : 0;
            double smooth_ratio = (smooth_base > 0) ? (double)smooth_lll / smooth_base : 0;

            fprintf(fout, "%-8d %-6d %-15.4e %-15.4e %-12d %-12d %-10.3f %-10.3f\n",
                    digits, trial + 1, score_base, score_lll,
                    smooth_base, smooth_lll, score_ratio, smooth_ratio);
            fflush(fout);

            total_score_base += score_base;
            total_score_lll += score_lll;
            total_smooth_base += smooth_base;
            total_smooth_lll += smooth_lll;

            /* Print detailed poly info for first trial of each size */
            if (trial == 0) {
                fprintf(fout, "\n  --- Detail for %d-digit N ---\n", digits);
                gmp_fprintf(fout, "  N = %Zd\n", N);
                poly_print(fout, "  Base-m", &poly_base);
                poly_print(fout, "  LLL", &poly_lll);

                mpz_t mc;
                mpz_init(mc);
                poly_max_coeff(mc, &poly_base);
                gmp_fprintf(fout, "  Base-m max|coeff| = %Zd\n", mc);
                poly_max_coeff(mc, &poly_lll);
                gmp_fprintf(fout, "  LLL max|coeff| = %Zd\n", mc);
                mpz_clear(mc);
                fprintf(fout, "\n");
            }

            poly_clear(&poly_base);
            poly_clear(&poly_lll);
            mpz_clear(N);
        }
    }

    fprintf(fout, "\n=== Summary ===\n");
    fprintf(fout, "Average score ratio (base-m / LLL): %.3f\n",
            total_score_base / total_score_lll);
    fprintf(fout, "  (>1 means LLL is better, i.e. lower score)\n");
    fprintf(fout, "Total smooth values - Base-m: %d, LLL: %d\n",
            total_smooth_base, total_smooth_lll);
    if (total_smooth_base > 0)
        fprintf(fout, "Smooth yield ratio (LLL / base-m): %.3f\n",
                (double)total_smooth_lll / total_smooth_base);
    fprintf(fout, "  (>1 means LLL produces more smooth values)\n");

    fprintf(fout, "\n=== Analysis ===\n");
    fprintf(fout, "The LLL lattice reduction applied to the base-m coefficient lattice\n");
    fprintf(fout, "seeks shorter vectors in the lattice of polynomials f with f(m)≡0 mod N.\n");
    fprintf(fout, "Shorter vectors correspond to polynomials with smaller coefficients,\n");
    fprintf(fout, "which produce smaller values when sieved, yielding more smooth relations.\n\n");
    fprintf(fout, "The combined score max_coeff * disc^{1/(2d-2)} captures both coefficient\n");
    fprintf(fout, "size (directly affects sieve values) and discriminant (affects the number\n");
    fprintf(fout, "field properties and root behavior).\n\n");
    fprintf(fout, "For these small semiprimes (30-40 digits), the improvement from LLL\n");
    fprintf(fout, "may be modest since base-m already produces fairly balanced coefficients.\n");
    fprintf(fout, "The benefit grows with the size of N, where the coefficient space is larger\n");
    fprintf(fout, "and LLL has more room to optimize.\n");

    fclose(fout);
    gmp_randclear(rstate);

    printf("Done. Results written to findings.txt\n");
    return 0;
}
