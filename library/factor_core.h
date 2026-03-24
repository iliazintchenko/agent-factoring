/*
 * factor_core.h - Core utilities for factoring experiments
 *
 * Provides: GMP wrappers, smoothness detection, relation management,
 * matrix solving over GF(2), and factor extraction.
 */
#ifndef FACTOR_CORE_H
#define FACTOR_CORE_H

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

/* ===== Timing ===== */
static inline double wall_time(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* ===== Prime sieve (simple Eratosthenes) ===== */
static int *sieve_primes(int B, int *count) {
    char *is_prime = (char *)calloc(B + 1, 1);
    if (!is_prime) return NULL;
    for (int i = 2; i <= B; i++) is_prime[i] = 1;
    for (int i = 2; (long long)i * i <= B; i++)
        if (is_prime[i])
            for (int j = i * i; j <= B; j += i)
                is_prime[j] = 0;
    int cnt = 0;
    for (int i = 2; i <= B; i++) if (is_prime[i]) cnt++;
    int *primes = (int *)malloc(cnt * sizeof(int));
    int idx = 0;
    for (int i = 2; i <= B; i++) if (is_prime[i]) primes[idx++] = i;
    *count = cnt;
    free(is_prime);
    return primes;
}

/* ===== Tonelli-Shanks: compute sqrt(n) mod p ===== */
static int tonelli_shanks(mpz_t result, const mpz_t n, unsigned long p) {
    mpz_t n_mod, tmp, z, M, c, t, R, b;
    mpz_inits(n_mod, tmp, z, M, c, t, R, b, NULL);

    mpz_set_ui(tmp, p);
    mpz_mod(n_mod, n, tmp);

    if (mpz_jacobi(n_mod, tmp) != 1) {
        mpz_clears(n_mod, tmp, z, M, c, t, R, b, NULL);
        return 0; /* not a QR */
    }

    if (p == 2) {
        mpz_set_ui(result, mpz_odd_p(n_mod) ? 1 : 0);
        mpz_clears(n_mod, tmp, z, M, c, t, R, b, NULL);
        return 1;
    }

    /* Factor out powers of 2: p-1 = Q * 2^S */
    unsigned long S = 0;
    unsigned long Q = p - 1;
    while ((Q & 1) == 0) { S++; Q >>= 1; }

    /* Find a non-residue z */
    unsigned long zz;
    for (zz = 2; zz < p; zz++) {
        mpz_set_ui(tmp, zz);
        mpz_set_ui(z, p);
        if (mpz_jacobi(tmp, z) == -1) break;
    }

    mpz_set_ui(M, S);
    mpz_set_ui(tmp, p);
    mpz_set_ui(c, zz);
    mpz_powm_ui(c, c, Q, tmp);          /* c = z^Q mod p */
    mpz_powm_ui(t, n_mod, Q, tmp);      /* t = n^Q mod p */
    mpz_powm_ui(R, n_mod, (Q + 1) / 2, tmp); /* R = n^((Q+1)/2) mod p */

    while (1) {
        if (mpz_cmp_ui(t, 0) == 0) { mpz_set_ui(result, 0); break; }
        if (mpz_cmp_ui(t, 1) == 0) { mpz_set(result, R); break; }

        /* Find least i such that t^(2^i) = 1 */
        unsigned long i = 0;
        mpz_set(b, t);
        while (mpz_cmp_ui(b, 1) != 0 && i < mpz_get_ui(M)) {
            mpz_mul(b, b, b);
            mpz_mod_ui(b, b, p);
            i++;
        }

        if (i == mpz_get_ui(M)) {
            mpz_clears(n_mod, tmp, z, M, c, t, R, b, NULL);
            return 0;
        }

        /* b = c^(2^(M-i-1)) */
        mpz_set(b, c);
        for (unsigned long j = 0; j < mpz_get_ui(M) - i - 1; j++) {
            mpz_mul(b, b, b);
            mpz_mod_ui(b, b, p);
        }

        mpz_set_ui(M, i);
        mpz_mul(c, b, b); mpz_mod_ui(c, c, p);
        mpz_mul(t, t, c); mpz_mod_ui(t, t, p);
        mpz_mul(R, R, b); mpz_mod_ui(R, R, p);
    }

    mpz_clears(n_mod, tmp, z, M, c, t, R, b, NULL);
    return 1;
}

/* ===== Relation storage ===== */
typedef struct {
    mpz_t x;          /* x value: x^2 ≡ product (mod N) */
    int *exponents;    /* exponent vector over factor base */
    int fb_size;       /* factor base size */
    int sign;          /* 1 if positive, -1 if negative */
} relation_t;

typedef struct {
    relation_t *rels;
    int count;
    int capacity;
    int fb_size;
} relation_set_t;

static void relset_init(relation_set_t *rs, int fb_size, int initial_cap) {
    rs->fb_size = fb_size;
    rs->count = 0;
    rs->capacity = initial_cap;
    rs->rels = (relation_t *)malloc(initial_cap * sizeof(relation_t));
}

static void relset_add(relation_set_t *rs, const mpz_t x, const int *exps, int sign) {
    if (rs->count >= rs->capacity) {
        rs->capacity *= 2;
        rs->rels = (relation_t *)realloc(rs->rels, rs->capacity * sizeof(relation_t));
    }
    relation_t *r = &rs->rels[rs->count];
    mpz_init_set(r->x, x);
    r->exponents = (int *)malloc(rs->fb_size * sizeof(int));
    memcpy(r->exponents, exps, rs->fb_size * sizeof(int));
    r->fb_size = rs->fb_size;
    r->sign = sign;
    rs->count++;
}

static void relset_free(relation_set_t *rs) {
    for (int i = 0; i < rs->count; i++) {
        mpz_clear(rs->rels[i].x);
        free(rs->rels[i].exponents);
    }
    free(rs->rels);
}

/* ===== GF(2) matrix solver (Gaussian elimination) ===== */
/* Finds a subset of relations whose exponent vectors sum to 0 mod 2 */
typedef unsigned long word_t;
#define WORD_BITS (sizeof(word_t) * 8)

static int solve_gf2(relation_set_t *rs, int **result_indices, int *result_count) {
    int nrows = rs->fb_size + 1; /* +1 for sign */
    int ncols = rs->count;
    if (ncols <= nrows) return 0; /* need more relations than factor base */

    int nwords = (ncols + WORD_BITS - 1) / WORD_BITS;
    word_t **matrix = (word_t **)calloc(nrows, sizeof(word_t *));
    for (int i = 0; i < nrows; i++)
        matrix[i] = (word_t *)calloc(nwords, sizeof(word_t));

    /* Fill matrix: row i = parity of exponent of prime i across all relations */
    for (int j = 0; j < ncols; j++) {
        relation_t *r = &rs->rels[j];
        /* Row 0: sign bit */
        if (r->sign < 0)
            matrix[0][j / WORD_BITS] |= (1UL << (j % WORD_BITS));
        /* Rows 1..fb_size: exponent parities */
        for (int i = 0; i < rs->fb_size; i++) {
            if (r->exponents[i] & 1)
                matrix[i + 1][j / WORD_BITS] |= (1UL << (j % WORD_BITS));
        }
    }

    /* Gaussian elimination to find null space */
    int *pivot_col = (int *)malloc(nrows * sizeof(int));
    int *is_pivot = (int *)calloc(ncols, sizeof(int));
    int rank = 0;

    for (int i = 0; i < nrows && rank < nrows; i++) {
        /* Find pivot */
        int pc = -1;
        for (int j = 0; j < ncols; j++) {
            if (!is_pivot[j] && (matrix[i][j / WORD_BITS] & (1UL << (j % WORD_BITS)))) {
                pc = j;
                break;
            }
        }
        if (pc == -1) continue;

        pivot_col[rank] = pc;
        is_pivot[pc] = 1;

        /* Eliminate */
        for (int k = 0; k < nrows; k++) {
            if (k != i && (matrix[k][pc / WORD_BITS] & (1UL << (pc % WORD_BITS)))) {
                for (int w = 0; w < nwords; w++)
                    matrix[k][w] ^= matrix[i][w];
            }
        }
        rank++;
    }

    /* Find free variables (non-pivot columns) */
    *result_count = 0;
    *result_indices = NULL;

    for (int j = 0; j < ncols; j++) {
        if (is_pivot[j]) continue;

        /* This free variable gives a null space vector */
        int *indices = (int *)malloc(ncols * sizeof(int));
        int cnt = 0;
        indices[cnt++] = j;

        /* Add pivot columns that have a 1 in column j */
        for (int i = 0; i < nrows; i++) {
            if (matrix[i][j / WORD_BITS] & (1UL << (j % WORD_BITS))) {
                /* Find which pivot column corresponds to this row */
                for (int r = 0; r < rank; r++) {
                    if (pivot_col[r] >= 0) {
                        /* Check if this row's pivot is at pivot_col[r] */
                        int pc = pivot_col[r];
                        if (matrix[i][pc / WORD_BITS] & (1UL << (pc % WORD_BITS))) {
                            indices[cnt++] = pc;
                            break;
                        }
                    }
                }
            }
        }

        *result_indices = indices;
        *result_count = cnt;
        break; /* Just return the first null space vector */
    }

    for (int i = 0; i < nrows; i++) free(matrix[i]);
    free(matrix);
    free(pivot_col);
    free(is_pivot);

    return (*result_count > 0) ? 1 : 0;
}

/* ===== Factor extraction ===== */
/* Given a set of relation indices forming a null space vector,
   compute x = product of all x_i mod N, y = sqrt(product of all f(x_i)) mod N,
   then gcd(x - y, N) */
static int extract_factor(relation_set_t *rs, int *indices, int cnt,
                          const mpz_t N, mpz_t factor) {
    mpz_t x_prod, y_squared, y, tmp;
    mpz_inits(x_prod, y_squared, y, tmp, NULL);

    mpz_set_ui(x_prod, 1);

    /* Accumulate exponents */
    int *total_exp = (int *)calloc(rs->fb_size, sizeof(int));
    int total_sign = 0;

    for (int i = 0; i < cnt; i++) {
        relation_t *r = &rs->rels[indices[i]];
        mpz_mul(x_prod, x_prod, r->x);
        mpz_mod(x_prod, x_prod, N);
        if (r->sign < 0) total_sign++;
        for (int j = 0; j < rs->fb_size; j++)
            total_exp[j] += r->exponents[j];
    }

    /* Compute y = product of primes^(total_exp/2) mod N */
    mpz_set_ui(y, 1);
    /* Note: we need the factor base primes here */
    /* This will be called with the factor base available externally */

    free(total_exp);
    mpz_clears(x_prod, y_squared, y, tmp, NULL);
    return 0; /* placeholder */
}

#endif /* FACTOR_CORE_H */
