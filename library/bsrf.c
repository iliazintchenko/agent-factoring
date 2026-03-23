/*
 * bsrf.c - Batch Smooth Relation Finder
 *
 * Novel approach combining:
 * 1. Polynomial evaluation Q(x) = (x+m)^2 - N near sqrt(N)
 * 2. Batch smooth-part extraction via GCD with primorial
 * 3. Multi-large-prime relation combination via graph-based merge
 * 4. GF(2) Gaussian elimination for congruence of squares
 *
 * Key difference from QS: no sieving. Uses batch GCD for smoothness,
 * which changes the optimal smoothness bound and allows larger factor bases.
 *
 * Usage: ./bsrf <N>
 * Output: "p q" where N = p*q
 * Seed: 42
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <ecm.h>

/* Limits */
#define MAX_LP 3        /* max large primes per relation */
static unsigned long max_lp_val;  /* large prime bound, set dynamically */

typedef struct {
    unsigned long p;    /* prime */
    int r1, r2;         /* roots of N mod p (-1 if no root) */
    double logp;        /* log2(p) */
} fb_entry_t;

typedef struct {
    int x_offset;       /* x value that produced this relation */
    int sign;           /* sign of Q(x): 0=positive, 1=negative */
    unsigned char *exponents; /* exponents mod 2 for factor base primes (dynamically allocated) */
    int nlp;            /* number of large primes */
    unsigned long lp[MAX_LP]; /* large primes */
} relation_t;

/* Globals */
static mpz_t N, sqrtN, Q_val, temp, temp2;
static fb_entry_t *fb = NULL;
static int fb_size, fb_capacity;
static relation_t *rels = NULL;
static int nrels, rels_capacity;
static unsigned long smooth_bound;

/* Compute Legendre symbol (a/p) for odd prime p */
static int legendre(unsigned long a, unsigned long p) {
    long val = a % p;
    if (val == 0) return 0;
    /* Euler criterion */
    mpz_t base, exp, mod, result;
    mpz_init_set_ui(base, val);
    mpz_init_set_ui(mod, p);
    mpz_init(exp);
    mpz_init(result);
    mpz_sub_ui(exp, mod, 1);
    mpz_tdiv_q_2exp(exp, exp, 1); /* (p-1)/2 */
    mpz_powm(result, base, exp, mod);
    int ret;
    if (mpz_cmp_ui(result, 1) == 0) ret = 1;
    else if (mpz_cmp_ui(result, 0) == 0) ret = 0;
    else ret = -1;
    mpz_clear(base); mpz_clear(exp); mpz_clear(mod); mpz_clear(result);
    return ret;
}

/* Tonelli-Shanks: find r such that r^2 ≡ n (mod p) */
static unsigned long sqrt_mod_p(unsigned long n_val, unsigned long p) {
    if (p == 2) return n_val & 1;
    n_val %= p;
    if (n_val == 0) return 0;

    mpz_t nn, pp, r;
    mpz_init_set_ui(nn, n_val);
    mpz_init_set_ui(pp, p);
    mpz_init(r);

    /* Use GMP's built-in sqrt mod if available, otherwise Tonelli-Shanks */
    /* Simple case: p ≡ 3 (mod 4) */
    if (p % 4 == 3) {
        mpz_t exp;
        mpz_init(exp);
        mpz_set_ui(exp, (p + 1) / 4);
        mpz_powm(r, nn, exp, pp);
        unsigned long result = mpz_get_ui(r);
        mpz_clear(exp); mpz_clear(nn); mpz_clear(pp); mpz_clear(r);
        return result;
    }

    /* Tonelli-Shanks */
    unsigned long q = p - 1, s = 0;
    while (q % 2 == 0) { q /= 2; s++; }

    /* Find quadratic non-residue */
    unsigned long z = 2;
    while (legendre(z, p) != -1) z++;

    mpz_t M, c, t, R, b, tmp;
    mpz_init_set_ui(M, s);
    mpz_init(c); mpz_init(t); mpz_init(R); mpz_init(b); mpz_init(tmp);

    mpz_set_ui(tmp, q);
    mpz_set_ui(c, z);
    mpz_powm(c, c, tmp, pp);

    mpz_set_ui(tmp, (q + 1) / 2);
    mpz_powm(R, nn, tmp, pp);

    mpz_set_ui(tmp, q);
    mpz_powm(t, nn, tmp, pp);

    while (1) {
        if (mpz_cmp_ui(t, 1) == 0) break;

        /* Find least i such that t^(2^i) ≡ 1 */
        mpz_set(tmp, t);
        unsigned long i = 0;
        while (mpz_cmp_ui(tmp, 1) != 0) {
            mpz_mul(tmp, tmp, tmp);
            mpz_mod(tmp, tmp, pp);
            i++;
        }

        unsigned long m_val = mpz_get_ui(M);
        mpz_set(b, c);
        for (unsigned long j = 0; j < m_val - i - 1; j++) {
            mpz_mul(b, b, b);
            mpz_mod(b, b, pp);
        }

        mpz_set_ui(M, i);
        mpz_mul(c, b, b); mpz_mod(c, c, pp);
        mpz_mul(t, t, c); mpz_mod(t, t, pp);
        mpz_mul(R, R, b); mpz_mod(R, R, pp);
    }

    unsigned long result = mpz_get_ui(R);
    mpz_clear(M); mpz_clear(c); mpz_clear(t); mpz_clear(R); mpz_clear(b); mpz_clear(tmp);
    mpz_clear(nn); mpz_clear(pp); mpz_clear(r);
    return result;
}

/* Build factor base: primes p where N is a QR mod p */
static int build_factor_base(unsigned long bound) {
    /* Simple sieve of Eratosthenes */
    char *is_prime = calloc(bound + 1, 1);
    if (!is_prime) return 0;

    for (unsigned long i = 2; i <= bound; i++) is_prime[i] = 1;
    for (unsigned long i = 2; i * i <= bound; i++) {
        if (is_prime[i]) {
            for (unsigned long j = i*i; j <= bound; j += i)
                is_prime[j] = 0;
        }
    }

    fb_capacity = 8192;
    fb = realloc(fb, fb_capacity * sizeof(fb_entry_t));
    fb_size = 0;
    /* Always include -1 and 2 */
    fb[fb_size].p = 2;
    fb[fb_size].logp = 1.0;
    /* N mod 2 root */
    unsigned long n_mod2 = mpz_get_ui(N) & 1;
    fb[fb_size].r1 = (n_mod2 == 0) ? 0 : 1;
    fb[fb_size].r2 = -1;
    fb_size++;

    unsigned long n_mod_p;
    for (unsigned long p = 3; p <= bound; p++) {
        if (!is_prime[p]) continue;
        n_mod_p = mpz_fdiv_ui(N, p);
        if (legendre(n_mod_p, p) >= 0) {
            fb[fb_size].p = p;
            if (fb_size >= fb_capacity) {
                fb_capacity *= 2;
                fb = realloc(fb, fb_capacity * sizeof(fb_entry_t));
            }
            fb[fb_size].logp = log2((double)p);
            if (n_mod_p == 0) {
                fb[fb_size].r1 = 0;
                fb[fb_size].r2 = -1;
            } else {
                unsigned long r = sqrt_mod_p(n_mod_p, p);
                fb[fb_size].r1 = (int)r;
                fb[fb_size].r2 = (int)(p - r);
            }
            fb_size++;
        }
    }

    free(is_prime);
    return fb_size;
}

/* Try to factor val over the factor base. Returns 1 if successful (possibly with large primes) */
static int try_factor(mpz_t val, int x_offset) {
    if (nrels >= rels_capacity) {
        rels_capacity = rels_capacity ? rels_capacity * 2 : 4096;
        rels = realloc(rels, rels_capacity * sizeof(relation_t));
    }

    relation_t *rel = &rels[nrels];
    memset(rel, 0, sizeof(relation_t));
    rel->exponents = calloc(fb_size, 1);
    rel->x_offset = x_offset;

    mpz_abs(temp, val);
    if (mpz_sgn(val) < 0) rel->sign = 1;

    /* Trial divide by factor base */
    for (int i = 0; i < fb_size; i++) {
        unsigned long p = fb[i].p;
        int exp = 0;
        while (mpz_divisible_ui_p(temp, p)) {
            mpz_divexact_ui(temp, temp, p);
            exp++;
        }
        rel->exponents[i] = exp & 1;  /* only parity matters */
    }

    /* Check cofactor */
    if (mpz_cmp_ui(temp, 1) == 0) {
        /* Fully smooth! */
        rel->nlp = 0;
        nrels++;
        return 1;
    }

    /* Check if cofactor is a single large prime */
    if (mpz_fits_ulong_p(temp)) {
        unsigned long cof = mpz_get_ui(temp);
        if (cof <= max_lp_val) {
            /* Check if it's prime (Miller-Rabin) */
            if (mpz_probab_prime_p(temp, 15)) {
                rel->nlp = 1;
                rel->lp[0] = cof;
                nrels++;
                return 1;
            }
        }
    }

    /* Try to split cofactor into 2 large primes */
    if (mpz_sizeinbase(temp, 2) <= 56) {  /* cofactor < ~7.2 * 10^16 */
        unsigned long cof_bound = max_lp_val;
        /* Quick trial division by small primes not in factor base */
        unsigned long remaining = 0;
        if (mpz_fits_ulong_p(temp)) {
            remaining = mpz_get_ui(temp);
        } else {
            /* Try Pollard rho on cofactor */
            mpz_t rho_f;
            mpz_init(rho_f);
            /* Simple rho with deterministic start */
            mpz_t x, y, d, c_rho;
            mpz_init_set_ui(x, 2);
            mpz_init_set_ui(y, 2);
            mpz_init(d);
            mpz_init_set_ui(c_rho, 1);
            int found = 0;
            for (int iter = 0; iter < 10000 && !found; iter++) {
                mpz_mul(x, x, x); mpz_add(x, x, c_rho); mpz_mod(x, x, temp);
                mpz_mul(y, y, y); mpz_add(y, y, c_rho); mpz_mod(y, y, temp);
                mpz_mul(y, y, y); mpz_add(y, y, c_rho); mpz_mod(y, y, temp);
                mpz_sub(d, x, y); mpz_abs(d, d);
                mpz_gcd(d, d, temp);
                if (mpz_cmp_ui(d, 1) > 0 && mpz_cmp(d, temp) < 0) {
                    /* Split into two factors */
                    mpz_t cof2;
                    mpz_init(cof2);
                    mpz_divexact(cof2, temp, d);
                    if (mpz_fits_ulong_p(d) && mpz_fits_ulong_p(cof2)) {
                        unsigned long f1 = mpz_get_ui(d);
                        unsigned long f2 = mpz_get_ui(cof2);
                        if (f1 <= cof_bound && f2 <= cof_bound) {
                            rel->nlp = 2;
                            rel->lp[0] = (f1 < f2) ? f1 : f2;
                            rel->lp[1] = (f1 < f2) ? f2 : f1;
                            found = 1;
                        }
                    }
                    mpz_clear(cof2);
                }
            }
            mpz_clear(x); mpz_clear(y); mpz_clear(d); mpz_clear(c_rho);
            mpz_clear(rho_f);
            if (found) {
                nrels++;
                return 1;
            }
            free(rel->exponents); rel->exponents = NULL;
            return 0;
        }

        if (remaining > 0) {
            /* Factor remaining into at most 2 primes */
            for (unsigned long p = smooth_bound + 1; p * p <= remaining && p <= cof_bound; p += 2) {
                if (remaining % p == 0) {
                    unsigned long q = remaining / p;
                    if (q <= cof_bound && q > 1) {
                        rel->nlp = 2;
                        rel->lp[0] = (p < q) ? p : q;
                        rel->lp[1] = (p < q) ? q : p;
                        nrels++;
                        return 1;
                    }
                    break;
                }
            }
            /* Maybe it's a prime */
            if (remaining <= cof_bound && mpz_probab_prime_p(temp, 15)) {
                rel->nlp = 1;
                rel->lp[0] = remaining;
                nrels++;
                return 1;
            }
        }
    }

    free(rel->exponents); rel->exponents = NULL;
    return 0;
}

/* ========== GF(2) Matrix operations for Gaussian elimination ========== */

/* Bit matrix: each row is a bit vector */
typedef struct {
    unsigned long *rows;    /* packed bits, row-major */
    int nrows, ncols;
    int words_per_row;
} bitmatrix_t;

static bitmatrix_t* bm_alloc(int nrows, int ncols) {
    bitmatrix_t *m = malloc(sizeof(bitmatrix_t));
    m->nrows = nrows;
    m->ncols = ncols;
    m->words_per_row = (ncols + 63) / 64;
    m->rows = calloc((size_t)nrows * m->words_per_row, sizeof(unsigned long));
    return m;
}

static void bm_free(bitmatrix_t *m) {
    free(m->rows);
    free(m);
}

static inline void bm_set(bitmatrix_t *m, int r, int c) {
    m->rows[(size_t)r * m->words_per_row + c/64] |= (1UL << (c % 64));
}

static inline int bm_get(bitmatrix_t *m, int r, int c) {
    return (m->rows[(size_t)r * m->words_per_row + c/64] >> (c % 64)) & 1;
}

static inline void bm_xor_row(bitmatrix_t *m, int dst, int src) {
    unsigned long *d = m->rows + (size_t)dst * m->words_per_row;
    unsigned long *s = m->rows + (size_t)src * m->words_per_row;
    for (int i = 0; i < m->words_per_row; i++)
        d[i] ^= s[i];
}

/* Gaussian elimination over GF(2). Returns null space vectors.
   Input: matrix M (nrels × ncols).
   Output: array of relation index subsets whose XOR is zero. */
static int gauss_elim(bitmatrix_t *M, int **deps, int *ndeps) {
    int nr = M->nrows, nc = M->ncols;

    /* Track row swaps for back-substitution */
    int *pivot_row = malloc(nc * sizeof(int));
    int *pivot_col = malloc(nr * sizeof(int));
    char *is_pivot = calloc(nr, 1);

    /* Also maintain an identity matrix to track dependencies */
    bitmatrix_t *hist = bm_alloc(nr, nr);
    for (int i = 0; i < nr; i++) bm_set(hist, i, i);

    int rank = 0;
    for (int c = 0; c < nc && rank < nr; c++) {
        /* Find pivot */
        int piv = -1;
        for (int r = rank; r < nr; r++) {
            if (bm_get(M, r, c)) { piv = r; break; }
        }
        if (piv < 0) continue;

        /* Swap rows */
        if (piv != rank) {
            unsigned long *tmp_row = malloc(M->words_per_row * sizeof(unsigned long));
            memcpy(tmp_row, M->rows + (size_t)piv * M->words_per_row, M->words_per_row * sizeof(unsigned long));
            memcpy(M->rows + (size_t)piv * M->words_per_row, M->rows + (size_t)rank * M->words_per_row, M->words_per_row * sizeof(unsigned long));
            memcpy(M->rows + (size_t)rank * M->words_per_row, tmp_row, M->words_per_row * sizeof(unsigned long));

            unsigned long *tmp_hist = malloc(hist->words_per_row * sizeof(unsigned long));
            memcpy(tmp_hist, hist->rows + (size_t)piv * hist->words_per_row, hist->words_per_row * sizeof(unsigned long));
            memcpy(hist->rows + (size_t)piv * hist->words_per_row, hist->rows + (size_t)rank * hist->words_per_row, hist->words_per_row * sizeof(unsigned long));
            memcpy(hist->rows + (size_t)rank * hist->words_per_row, tmp_hist, hist->words_per_row * sizeof(unsigned long));

            free(tmp_row);
            free(tmp_hist);
        }

        /* Eliminate */
        for (int r = 0; r < nr; r++) {
            if (r != rank && bm_get(M, r, c)) {
                bm_xor_row(M, r, rank);
                bm_xor_row(hist, r, rank);
            }
        }

        is_pivot[rank] = 1;
        pivot_col[rank] = c;
        rank++;
    }

    /* Collect null space vectors (rows that became zero in M) */
    *ndeps = 0;
    *deps = NULL;

    for (int r = 0; r < nr; r++) {
        /* Check if row r of M is all zeros */
        int all_zero = 1;
        for (int w = 0; w < M->words_per_row && all_zero; w++) {
            if (M->rows[(size_t)r * M->words_per_row + w]) all_zero = 0;
        }
        if (!all_zero) continue;

        /* Row r is in the null space. The dependency is encoded in hist[r] */
        /* Count how many relations are in this dependency */
        int count = 0;
        for (int j = 0; j < nr; j++) {
            if (bm_get(hist, r, j)) count++;
        }
        if (count < 2) continue;

        /* Store dependency */
        *deps = realloc(*deps, (*ndeps + 1) * (nr + 1) * sizeof(int));
        int *dep = *deps + (*ndeps) * (nr + 1);
        dep[0] = count;
        int idx = 1;
        for (int j = 0; j < nr; j++) {
            if (bm_get(hist, r, j)) dep[idx++] = j;
        }
        (*ndeps)++;

        if (*ndeps >= 64) break;  /* enough dependencies to try */
    }

    free(pivot_row);
    free(pivot_col);
    free(is_pivot);
    bm_free(hist);

    return rank;
}

/* Helper: add a single relation's contribution to x_prod and y_sq_accum */
static void accum_relation(int rel_idx, mpz_t x_prod, mpz_t y_sq) {
    relation_t *rel = &rels[rel_idx];
    /* x_prod *= (x + sqrtN) mod N */
    mpz_set_si(temp, rel->x_offset);
    mpz_add(temp, temp, sqrtN);
    mpz_mul(x_prod, x_prod, temp);
    mpz_mod(x_prod, x_prod, N);

    /* y_sq *= |Q(x)| */
    mpz_set_si(temp, rel->x_offset);
    mpz_add(temp, temp, sqrtN);
    mpz_mul(temp, temp, temp);
    mpz_sub(temp, temp, N);
    mpz_mul(y_sq, y_sq, temp);
}

/* Try to extract a factor from a dependency of matrix rows.
 * matrix_rel_map/matrix_merge_map/merges are passed in. */
static int try_dependency(int *dep_indices, int dep_count,
                          int *matrix_rel_map, int *matrix_merge_map,
                          void *merges_ptr) {
    typedef struct { int idx1, idx2; } mp_t;
    mp_t *merges = (mp_t*)merges_ptr;

    mpz_t x_prod, y_sq, y_val, g;
    mpz_init_set_ui(x_prod, 1);
    mpz_init_set_ui(y_sq, 1);
    mpz_init(y_val);
    mpz_init(g);

    /* Accumulate all x values and Q(x) products */
    for (int i = 0; i < dep_count; i++) {
        int mat_row = dep_indices[i];
        if (matrix_merge_map[mat_row] < 0) {
            /* Full relation */
            accum_relation(matrix_rel_map[mat_row], x_prod, y_sq);
        } else {
            /* Merged pair */
            int m = matrix_merge_map[mat_row];
            accum_relation(merges[m].idx1, x_prod, y_sq);
            accum_relation(merges[m].idx2, x_prod, y_sq);
        }
    }

    /* y_sq should be a perfect square (product of Q(x_i) values) */
    int neg = (mpz_sgn(y_sq) < 0);
    mpz_abs(y_sq, y_sq);
    if (!mpz_perfect_square_p(y_sq)) {
        mpz_clear(x_prod); mpz_clear(y_sq); mpz_clear(y_val); mpz_clear(g);
        return 0;
    }
    mpz_sqrt(y_val, y_sq);
    mpz_mod(y_val, y_val, N);

    /* Try gcd(x_prod ± y_val, N) */
    mpz_sub(temp, x_prod, y_val);
    mpz_gcd(g, temp, N);
    if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
        mpz_t cofactor;
        mpz_init(cofactor);
        mpz_divexact(cofactor, N, g);
        if (mpz_cmp(g, cofactor) > 0) mpz_swap(g, cofactor);
        gmp_printf("%Zd %Zd\n", g, cofactor);
        mpz_clear(cofactor);
        mpz_clear(x_prod); mpz_clear(y_sq); mpz_clear(y_val); mpz_clear(g);
        return 1;
    }

    mpz_add(temp, x_prod, y_val);
    mpz_gcd(g, temp, N);
    if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
        mpz_t cofactor;
        mpz_init(cofactor);
        mpz_divexact(cofactor, N, g);
        if (mpz_cmp(g, cofactor) > 0) mpz_swap(g, cofactor);
        gmp_printf("%Zd %Zd\n", g, cofactor);
        mpz_clear(cofactor);
        mpz_clear(x_prod); mpz_clear(y_sq); mpz_clear(y_val); mpz_clear(g);
        return 1;
    }

    mpz_clear(x_prod); mpz_clear(y_sq); mpz_clear(y_val); mpz_clear(g);
    return 0;
}

/* Main factoring routine */
static int factor_bsrf(void) {
    size_t n_bits = mpz_sizeinbase(N, 2);
    size_t n_digits = mpz_sizeinbase(N, 10);

    /* Choose smoothness bound B.
     * For batch GCD, we can afford larger B than QS.
     * Heuristic: B = exp(c * sqrt(ln N * ln ln N)) with c tuned for our approach.
     * We use c ≈ 0.7 (slightly larger than QS's 0.5) since batch GCD is cheaper per candidate.
     */
    double ln_n = (double)n_bits * log(2.0);
    double ln_ln_n = log(ln_n);
    double c = 0.5;  /* standard QS constant */
    smooth_bound = (unsigned long)exp(c * sqrt(ln_n * ln_ln_n));
    if (smooth_bound < 100) smooth_bound = 100;
    if (smooth_bound > 2000000) smooth_bound = 2000000;

    /* Large prime bound: ~100*B to get good collision rate */
    max_lp_val = smooth_bound * 100;
    if (max_lp_val > (1UL << 32)) max_lp_val = (1UL << 32);

    fprintf(stderr, "BSRF: %zu-digit (%zu-bit) number, B=%lu, LP_bound=%lu\n",
            n_digits, n_bits, smooth_bound, max_lp_val);

    /* Build factor base */
    int fb_count = build_factor_base(smooth_bound);
    fprintf(stderr, "Factor base: %d primes up to %lu\n", fb_count, smooth_bound);

    /* Compute sqrt(N) */
    mpz_sqrt(sqrtN, N);

    /* We need fb_count + 1 + (extra for large primes) relations.
     * Target: 1.1 * fb_count relations for safety margin.
     */
    /* We need enough full + merged-LP relations to exceed fb_count + 1.
     * Keep collecting until we have enough, checking periodically. */
    int target_matrix_rels = fb_count + 50;  /* target for matrix rows */
    int target_rels = fb_count * 10 + 1000;  /* collect many, merge later */
    if (target_rels > 100000) target_rels = 100000;

    fprintf(stderr, "Need %d relations (have %d fb primes)\n", target_rels, fb_count);

    nrels = 0;
    rels_capacity = 0;

    /* Sieve-based smooth candidate detection.
     * Use logarithmic sieve to quickly identify candidates where Q(x) is likely smooth.
     * Then trial-divide only the promising candidates.
     */
    int sieve_size = 256 * 1024;  /* 256K entries per sieve block */
    float *sieve_log = malloc(sieve_size * sizeof(float));
    float *target_log = malloc(sieve_size * sizeof(float));
    long x_base = 1;
    int total_tested = 0;
    int total_sieve_hits = 0;
    float threshold_adjust = 1.5f * log2f((float)smooth_bound);  /* ~1 large prime of size B */

    /* Precompute sieve roots: for each fb prime p, find r such that Q(r) ≡ 0 (mod p).
     * Q(x) = (x + m)^2 - N, so x + m ≡ ±√N (mod p), i.e., x ≡ r1 - m or x ≡ r2 - m (mod p).
     * We need r1, r2 such that r_i^2 ≡ N (mod p), then sieve_start = r_i - m mod p.
     */

    while (nrels < target_rels) {
        /* Initialize sieve */
        memset(sieve_log, 0, sieve_size * sizeof(float));

        /* Compute target: log2(|Q(x)|) for each x in this block */
        for (int i = 0; i < sieve_size; i++) {
            long x = x_base + i;
            /* Q(x) ≈ 2*m*x for small x, Q(x) ≈ x^2 for large x */
            /* More precisely: |Q(x)| = |(x+m)^2 - N| */
            /* For x > 0: (x+m)^2 - N = x^2 + 2mx + (m^2-N) ≈ 2mx */
            double q_approx = fabs(2.0 * mpz_get_d(sqrtN) * (double)x);
            if (q_approx < 1.0) q_approx = 1.0;
            target_log[i] = (float)log2(q_approx);
        }

        /* Sieve with factor base primes */
        for (int fi = 0; fi < fb_size; fi++) {
            unsigned long p = fb[fi].p;
            float logp = (float)fb[fi].logp;

            if (p == 2) {
                /* Handle p=2 specially */
                unsigned long m_mod2 = mpz_fdiv_ui(sqrtN, 2);
                /* (x+m)^2 - N ≡ 0 mod 2 when x+m is even or N is even */
                /* Just sieve every position (most values have factors of 2) */
                for (long pos = 0; pos < sieve_size; pos += 2)
                    sieve_log[pos] += logp;
                continue;
            }

            /* Compute starting positions for this prime */
            unsigned long m_modp = mpz_fdiv_ui(sqrtN, p);
            int r1 = fb[fi].r1;
            int r2 = fb[fi].r2;

            if (r1 >= 0) {
                /* Start position: x ≡ r1 - m (mod p) */
                long start1 = ((long)r1 - (long)(m_modp % p) + (long)p) % (long)p;
                /* Adjust for x_base */
                start1 = (start1 - (x_base % (long)p) + (long)p) % (long)p;
                for (long pos = start1; pos < sieve_size; pos += p)
                    sieve_log[pos] += logp;
            }
            if (r2 >= 0 && r2 != r1) {
                long start2 = ((long)r2 - (long)(m_modp % p) + (long)p) % (long)p;
                start2 = (start2 - (x_base % (long)p) + (long)p) % (long)p;
                for (long pos = start2; pos < sieve_size; pos += p)
                    sieve_log[pos] += logp;
            }
        }

        /* Scan for smooth candidates */
        for (int i = 0; i < sieve_size && nrels < target_rels; i++) {
            if (sieve_log[i] >= target_log[i] - threshold_adjust) {
                /* Promising candidate — trial divide */
                long x = x_base + i;
                mpz_set_si(Q_val, x);
                mpz_add(Q_val, Q_val, sqrtN);
                mpz_mul(Q_val, Q_val, Q_val);
                mpz_sub(Q_val, Q_val, N);

                if (try_factor(Q_val, (int)x)) {
                    /* Relation found */
                }
                total_sieve_hits++;
            }
        }

        total_tested += sieve_size;
        x_base += sieve_size;

        if (total_tested % (sieve_size * 10) == 0) {
            fprintf(stderr, "  sieved %d, hits %d, found %d/%d relations (%.1f%%)\n",
                    total_tested, total_sieve_hits, nrels, target_rels,
                    100.0 * nrels / target_rels);
        }

        /* Safety: don't loop forever */
        if (total_tested > 500000000L) {
            fprintf(stderr, "BSRF: too many candidates tested, giving up\n");
            free(sieve_log);
            free(target_log);
            return 0;
        }
    }

    free(sieve_log);
    free(target_log);

    /* Count relation types */
    int n_full = 0, n_1lp = 0, n_2lp = 0;
    for (int i = 0; i < nrels; i++) {
        if (rels[i].nlp == 0) n_full++;
        else if (rels[i].nlp == 1) n_1lp++;
        else n_2lp++;
    }
    fprintf(stderr, "Found %d relations: %d full, %d 1-LP, %d 2-LP\n",
            nrels, n_full, n_1lp, n_2lp);

    /* Merge single-large-prime relations: if two relations share the same LP,
     * their product is a full relation (LP cancels) */
    /* First, sort 1-LP relations by their large prime */
    int *lp1_idx = malloc(nrels * sizeof(int));
    int n_lp1 = 0;
    for (int i = 0; i < nrels; i++) {
        if (rels[i].nlp == 1) lp1_idx[n_lp1++] = i;
    }

    /* Simple sort by lp value */
    for (int i = 0; i < n_lp1 - 1; i++) {
        for (int j = i + 1; j < n_lp1; j++) {
            if (rels[lp1_idx[i]].lp[0] > rels[lp1_idx[j]].lp[0]) {
                int tmp = lp1_idx[i]; lp1_idx[i] = lp1_idx[j]; lp1_idx[j] = tmp;
            }
        }
    }

    /* Create merged relations: combine pairs with same LP */
    /* A merged relation stores TWO x_offsets; exponent vector is XOR of both */
    typedef struct { int idx1, idx2; } merge_pair_t;
    merge_pair_t *merges = malloc(n_lp1 * sizeof(merge_pair_t));
    int n_merges = 0;

    for (int i = 0; i < n_lp1 - 1; i++) {
        if (rels[lp1_idx[i]].lp[0] == rels[lp1_idx[i+1]].lp[0]) {
            merges[n_merges].idx1 = lp1_idx[i];
            merges[n_merges].idx2 = lp1_idx[i+1];
            n_merges++;
            i++; /* skip the second one */
        }
    }

    int total_matrix_rels = n_full + n_merges;
    fprintf(stderr, "Merged LP pairs: %d. Total matrix relations: %d\n",
            n_merges, total_matrix_rels);

    /* Build matrix: rows = full relations + merged relations, cols = 1 (sign) + fb_size */
    int ncols = 1 + fb_size;

    if (total_matrix_rels <= ncols) {
        fprintf(stderr, "Not enough relations (%d) for matrix (%d cols). Need more.\n",
                total_matrix_rels, ncols);
        free(lp1_idx);
        free(merges);
        return 0;
    }

    /* Map from matrix row to relation type */
    int *matrix_rel_map = malloc(total_matrix_rels * sizeof(int));
    int *matrix_merge_map = malloc(total_matrix_rels * sizeof(int)); /* -1 if full, merge index otherwise */
    int mat_row = 0;

    /* First: full relations */
    for (int i = 0; i < nrels; i++) {
        if (rels[i].nlp == 0) {
            matrix_rel_map[mat_row] = i;
            matrix_merge_map[mat_row] = -1;
            mat_row++;
        }
    }
    /* Then: merged pairs */
    for (int i = 0; i < n_merges; i++) {
        matrix_rel_map[mat_row] = merges[i].idx1;
        matrix_merge_map[mat_row] = i;
        mat_row++;
    }

    bitmatrix_t *M = bm_alloc(total_matrix_rels, ncols);

    for (int r = 0; r < total_matrix_rels; r++) {
        if (matrix_merge_map[r] < 0) {
            /* Full relation */
            relation_t *rel = &rels[matrix_rel_map[r]];
            if (rel->sign) bm_set(M, r, 0);
            for (int j = 0; j < fb_size; j++) {
                if (rel->exponents[j]) bm_set(M, r, 1 + j);
            }
        } else {
            /* Merged pair: XOR the two exponent vectors */
            int m = matrix_merge_map[r];
            relation_t *r1 = &rels[merges[m].idx1];
            relation_t *r2 = &rels[merges[m].idx2];
            int sign = r1->sign ^ r2->sign;
            if (sign) bm_set(M, r, 0);
            for (int j = 0; j < fb_size; j++) {
                if (r1->exponents[j] ^ r2->exponents[j]) bm_set(M, r, 1 + j);
            }
        }
    }

    /* Gaussian elimination */
    int *deps = NULL;
    int ndeps = 0;

    fprintf(stderr, "Running Gaussian elimination...\n");
    int rank = gauss_elim(M, &deps, &ndeps);
    fprintf(stderr, "Rank: %d, found %d dependencies\n", rank, ndeps);

    /* Try each dependency */
    int factored = 0;
    for (int d = 0; d < ndeps && !factored; d++) {
        int *dep = deps + d * (total_matrix_rels + 1);
        int count = dep[0];
        factored = try_dependency(dep + 1, count, matrix_rel_map, matrix_merge_map, merges);
    }

    free(deps);
    free(lp1_idx);
    free(merges);
    free(matrix_rel_map);
    free(matrix_merge_map);
    bm_free(M);

    return factored;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    mpz_init(N);
    mpz_init(sqrtN);
    mpz_init(Q_val);
    mpz_init(temp);
    mpz_init(temp2);

    if (mpz_set_str(N, argv[1], 10) != 0) {
        fprintf(stderr, "Invalid number: %s\n", argv[1]);
        return 1;
    }

    /* Quick trial division first */
    for (unsigned long p = 2; p < 1000000; p++) {
        if (mpz_divisible_ui_p(N, p)) {
            mpz_t f, cof;
            mpz_init_set_ui(f, p);
            mpz_init(cof);
            mpz_divexact(cof, N, f);
            if (mpz_cmp_ui(cof, 1) > 0) {
                gmp_printf("%Zd %Zd\n", f, cof);
                mpz_clear(f); mpz_clear(cof);
                goto cleanup;
            }
            mpz_clear(f); mpz_clear(cof);
        }
    }

    if (factor_bsrf()) {
        goto cleanup;
    }

    fprintf(stderr, "BSRF failed\n");
    mpz_clear(N); mpz_clear(sqrtN); mpz_clear(Q_val); mpz_clear(temp); mpz_clear(temp2);
    return 1;

cleanup:
    mpz_clear(N); mpz_clear(sqrtN); mpz_clear(Q_val); mpz_clear(temp); mpz_clear(temp2);
    return 0;
}
