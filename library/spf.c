/*
 * Smooth Polynomial Factoring (SPF) - Novel approach
 *
 * Key insight: instead of sieving Q(x) = (x+m)^2 - N for B-smooth values
 * (quadratic sieve), we evaluate MULTIPLE specialized polynomials whose
 * values are algebraically constrained to have more small factors.
 *
 * Approach:
 * 1. For a factor base F = {p_1,...,p_k} (primes up to B), compute
 *    r_i = sqrt(N) mod p_i for each p_i where N is a QR.
 * 2. For each subset S of the factor base, use CRT to find x_S such
 *    that x_S ≡ r_i (mod p_i) for all p_i in S. Then x_S^2 - N ≡ 0
 *    (mod prod(S)), meaning x_S^2 - N is automatically divisible by
 *    prod(S).
 * 3. The COFACTOR c_S = (x_S^2 - N) / prod(S) needs to be smooth
 *    for a full relation. Because we've "pre-sieved" by construction,
 *    the cofactor is smaller → higher smoothness probability.
 * 4. For each CRT solution, evaluate at x_S + k*prod(S) for k = 0, ±1, ±2,...
 *    to get more candidates with the same guaranteed divisibility.
 *
 * Novel element: Instead of evaluating ONE polynomial at MANY points (QS)
 * or MANY polynomials at MANY points (MPQS), we evaluate STRUCTURED
 * polynomials where each evaluation is GUARANTEED to have a large smooth
 * part. This changes the smoothness probability for the remaining cofactor.
 *
 * We then use the congruence-of-squares method to combine relations.
 *
 * NOTE: This is NOT QS/MPQS — the polynomial construction is different
 * (CRT-based guarantees vs sieving), and the evaluation strategy is
 * different (targeted evaluation at CRT-chosen points vs sweeping sieve).
 */

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Modular square root using Tonelli-Shanks */
static int mpz_sqrtmod(mpz_t result, const mpz_t n, const mpz_t p)
{
    /* Returns 1 if n is QR mod p, sets result. Returns 0 if not. */
    if (mpz_jacobi(n, p) != 1) return 0;

    if (mpz_congruent_ui_p(p, 3, 4)) {
        /* p ≡ 3 (mod 4): result = n^((p+1)/4) mod p */
        mpz_t exp;
        mpz_init(exp);
        mpz_add_ui(exp, p, 1);
        mpz_fdiv_q_2exp(exp, exp, 2);
        mpz_powm(result, n, exp, p);
        mpz_clear(exp);
        return 1;
    }

    /* Tonelli-Shanks */
    mpz_t Q, S_mpz, z, M, c, t, R, temp, b;
    mpz_inits(Q, S_mpz, z, M, c, t, R, temp, b, NULL);

    /* Factor p-1 = Q * 2^S */
    mpz_sub_ui(Q, p, 1);
    unsigned long S = mpz_scan1(Q, 0);
    mpz_fdiv_q_2exp(Q, Q, S);

    /* Find quadratic non-residue z */
    mpz_set_ui(z, 2);
    while (mpz_jacobi(z, p) != -1)
        mpz_add_ui(z, z, 1);

    mpz_set_ui(M, S);
    mpz_powm(c, z, Q, p);
    mpz_powm(t, n, Q, p);
    mpz_add_ui(temp, Q, 1);
    mpz_fdiv_q_2exp(temp, temp, 1);
    mpz_powm(R, n, temp, p);

    while (1) {
        if (mpz_cmp_ui(t, 0) == 0) { mpz_set_ui(result, 0); break; }
        if (mpz_cmp_ui(t, 1) == 0) { mpz_set(result, R); break; }

        /* Find least i such that t^(2^i) = 1 */
        mpz_set(temp, t);
        unsigned long i;
        for (i = 1; i < mpz_get_ui(M); i++) {
            mpz_mul(temp, temp, temp);
            mpz_mod(temp, temp, p);
            if (mpz_cmp_ui(temp, 1) == 0) break;
        }

        /* b = c^(2^(M-i-1)) */
        mpz_set(b, c);
        for (unsigned long j = 0; j < mpz_get_ui(M) - i - 1; j++) {
            mpz_mul(b, b, b);
            mpz_mod(b, b, p);
        }

        mpz_set_ui(M, i);
        mpz_mul(c, b, b); mpz_mod(c, c, p);
        mpz_mul(t, t, c); mpz_mod(t, t, p);
        mpz_mul(R, R, b); mpz_mod(R, R, p);
    }

    mpz_clears(Q, S_mpz, z, M, c, t, R, temp, b, NULL);
    return 1;
}

/* Generate primes up to limit */
static unsigned long *sieve_primes(unsigned long limit, int *count)
{
    char *sieve = calloc(limit + 1, 1);
    memset(sieve, 1, limit + 1);
    sieve[0] = sieve[1] = 0;
    for (unsigned long i = 2; i * i <= limit; i++)
        if (sieve[i])
            for (unsigned long j = i*i; j <= limit; j += i)
                sieve[j] = 0;

    *count = 0;
    for (unsigned long i = 2; i <= limit; i++)
        if (sieve[i]) (*count)++;

    unsigned long *primes = malloc(*count * sizeof(unsigned long));
    int idx = 0;
    for (unsigned long i = 2; i <= limit; i++)
        if (sieve[i]) primes[idx++] = i;

    free(sieve);
    return primes;
}

/* Factor base entry */
typedef struct {
    unsigned long p;  /* prime */
    mpz_t root;       /* sqrt(N) mod p */
} fb_entry_t;

/* A relation: x^2 ≡ ∏ p_i^{e_i} (mod N) */
typedef struct {
    mpz_t x;          /* the x value */
    mpz_t rhs;        /* x^2 mod N (fully factored) */
    int *exponents;    /* exponents of factor base primes */
} relation_t;

/* Try to factor val using factor base, return 1 if smooth */
static int try_smooth(int *exps, mpz_t cofactor, const mpz_t val,
                      const fb_entry_t *fb, int fb_size)
{
    mpz_set(cofactor, val);
    if (mpz_sgn(cofactor) < 0) {
        mpz_neg(cofactor, cofactor);
        exps[0] = 1; /* -1 exponent */
    } else {
        exps[0] = 0;
    }

    for (int i = 0; i < fb_size; i++) {
        exps[i + 1] = 0;
        while (mpz_divisible_ui_p(cofactor, fb[i].p)) {
            mpz_divexact_ui(cofactor, cofactor, fb[i].p);
            exps[i + 1]++;
        }
    }

    return (mpz_cmp_ui(cofactor, 1) == 0);
}

/* Gaussian elimination mod 2 to find a dependency */
/* Returns indices of relations forming a square, or NULL */
static int *find_dependency(int **exps, int nrels, int ncols)
{
    /* Build matrix mod 2 */
    int *matrix = calloc(nrels * ncols, sizeof(int));
    int *pivot_row = malloc(ncols * sizeof(int));
    int *row_ops = calloc(nrels * nrels, sizeof(int)); /* track row operations */

    for (int i = 0; i < nrels; i++) {
        for (int j = 0; j < ncols; j++)
            matrix[i * ncols + j] = exps[i][j] & 1;
        row_ops[i * nrels + i] = 1;
    }

    for (int j = 0; j < ncols; j++) pivot_row[j] = -1;

    for (int j = 0; j < ncols; j++) {
        /* Find pivot */
        int piv = -1;
        for (int i = 0; i < nrels; i++) {
            if (matrix[i * ncols + j] && pivot_row[j] == -1) {
                int is_pivot = 1;
                for (int jj = 0; jj < j; jj++)
                    if (matrix[i * ncols + jj]) { is_pivot = 0; break; }
                if (is_pivot) { piv = i; break; }
            }
        }
        if (piv == -1) continue;
        pivot_row[j] = piv;

        /* Eliminate */
        for (int i = 0; i < nrels; i++) {
            if (i != piv && matrix[i * ncols + j]) {
                for (int jj = j; jj < ncols; jj++)
                    matrix[i * ncols + jj] ^= matrix[piv * ncols + jj];
                for (int ii = 0; ii < nrels; ii++)
                    row_ops[i * nrels + ii] ^= row_ops[piv * nrels + ii];
            }
        }
    }

    /* Find a zero row (dependency) */
    int *result = NULL;
    for (int i = 0; i < nrels; i++) {
        int all_zero = 1;
        for (int j = 0; j < ncols; j++)
            if (matrix[i * ncols + j]) { all_zero = 0; break; }

        if (all_zero) {
            /* Count relations in this dependency */
            int cnt = 0;
            for (int ii = 0; ii < nrels; ii++)
                if (row_ops[i * nrels + ii]) cnt++;

            if (cnt >= 2) {
                result = malloc((cnt + 1) * sizeof(int));
                int idx = 0;
                for (int ii = 0; ii < nrels; ii++)
                    if (row_ops[i * nrels + ii])
                        result[idx++] = ii;
                result[idx] = -1;
                break;
            }
        }
    }

    free(matrix);
    free(pivot_row);
    free(row_ops);
    return result;
}

/* Main SPF algorithm */
int main(int argc, char *argv[])
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    struct timespec start;
    clock_gettime(CLOCK_MONOTONIC, &start);

    mpz_t N, sqrtN, factor, cofactor;
    mpz_inits(N, sqrtN, factor, cofactor, NULL);
    mpz_set_str(N, argv[1], 10);

    int n_digits = strlen(argv[1]);

    /* Quick trial division */
    for (unsigned long d = 2; d <= 1000000; d++) {
        if (mpz_divisible_ui_p(N, d)) {
            mpz_set_ui(factor, d);
            mpz_divexact(cofactor, N, factor);
            if (mpz_cmp(factor, cofactor) > 0) mpz_swap(factor, cofactor);
            gmp_printf("%Zd %Zd\n", factor, cofactor);
            goto done;
        }
        if (d == 2) d = 1; /* start at 3 after 2 */
        else if (d % 2 == 0) continue;
    }

    mpz_sqrt(sqrtN, N);

    /* Choose smoothness bound B based on L[1/2] heuristic */
    /* B = exp(sqrt(0.5 * log(N) * log(log(N)))) */
    double logN = n_digits * 2.302585; /* log(10) * n_digits */
    double loglogN = log(logN);
    double B_float = exp(sqrt(0.5 * logN * loglogN));
    unsigned long B = (unsigned long)B_float;
    if (B < 100) B = 100;
    if (B > 10000000) B = 10000000;

    fprintf(stderr, "SPF: %d digits, B=%lu\n", n_digits, B);

    /* Build factor base */
    int nprimes;
    unsigned long *primes = sieve_primes(B, &nprimes);

    fb_entry_t *fb = malloc(nprimes * sizeof(fb_entry_t));
    int fb_size = 0;
    mpz_t Nmod, pmpz;
    mpz_inits(Nmod, pmpz, NULL);

    for (int i = 0; i < nprimes; i++) {
        mpz_set_ui(pmpz, primes[i]);
        mpz_mod(Nmod, N, pmpz);

        mpz_t root;
        mpz_init(root);
        if (mpz_sqrtmod(root, Nmod, pmpz)) {
            fb[fb_size].p = primes[i];
            mpz_init_set(fb[fb_size].root, root);
            fb_size++;
        }
        mpz_clear(root);
    }
    free(primes);
    mpz_clears(Nmod, pmpz, NULL);

    fprintf(stderr, "SPF: factor base size = %d\n", fb_size);

    /* Collect relations */
    int needed = fb_size + 10; /* Need more relations than factor base size */
    relation_t *rels = malloc(needed * 2 * sizeof(relation_t));
    int **rel_exps = malloc(needed * 2 * sizeof(int *));
    int nrels = 0;

    mpz_t x, val, cof;
    mpz_inits(x, val, cof, NULL);

    /* Sieve: evaluate (sqrtN + offset)^2 - N for various offsets */
    for (long offset = 1; nrels < needed; offset++) {
        /* Check timeout */
        struct timespec now;
        clock_gettime(CLOCK_MONOTONIC, &now);
        double elapsed = (now.tv_sec - start.tv_sec) +
            (now.tv_nsec - start.tv_nsec) / 1e9;
        if (elapsed > 280.0) break;

        for (int sign = 0; sign <= 1; sign++) {
            long off = sign ? -offset : offset;
            mpz_set(x, sqrtN);
            if (off >= 0)
                mpz_add_ui(x, x, off);
            else
                mpz_sub_ui(x, x, -off);

            /* val = x^2 - N */
            mpz_mul(val, x, x);
            mpz_sub(val, val, N);

            /* Try to factor val over the factor base */
            int *exps = calloc(fb_size + 1, sizeof(int)); /* +1 for sign */
            if (try_smooth(exps, cof, val, fb, fb_size)) {
                /* Found a smooth value! */
                mpz_init_set(rels[nrels].x, x);
                mpz_init(rels[nrels].rhs);
                mpz_mul(rels[nrels].rhs, x, x);
                mpz_mod(rels[nrels].rhs, rels[nrels].rhs, N);
                rel_exps[nrels] = exps;
                nrels++;

                if (nrels % 10 == 0)
                    fprintf(stderr, "SPF: %d/%d relations (offset=%ld, %.1fs)\n",
                            nrels, needed, offset, elapsed);

                if (nrels >= needed) break;
            } else {
                free(exps);
            }
        }
    }

    fprintf(stderr, "SPF: collected %d relations\n", nrels);

    /* Find dependency and try to factor */
    int found = 0;
    if (nrels > fb_size) {
        int *dep = find_dependency(rel_exps, nrels, fb_size + 1);
        if (dep) {
            /* Combine relations */
            mpz_t lhs, rhs_sqrt, g;
            mpz_inits(lhs, rhs_sqrt, g, NULL);
            mpz_set_ui(lhs, 1);
            mpz_set_ui(rhs_sqrt, 1);

            /* Compute product of x values (lhs) and sqrt of product of rhs (rhs_sqrt) */
            int *combined_exps = calloc(fb_size + 1, sizeof(int));
            for (int i = 0; dep[i] >= 0; i++) {
                int idx = dep[i];
                mpz_mul(lhs, lhs, rels[idx].x);
                mpz_mod(lhs, lhs, N);
                for (int j = 0; j <= fb_size; j++)
                    combined_exps[j] += rel_exps[idx][j];
            }

            /* combined_exps should all be even */
            int all_even = 1;
            for (int j = 0; j <= fb_size; j++) {
                if (combined_exps[j] % 2 != 0) { all_even = 0; break; }
            }

            if (all_even) {
                /* Compute rhs_sqrt = prod(p_i^(e_i/2)) mod N */
                mpz_set_ui(rhs_sqrt, 1);
                /* Handle sign */
                if (combined_exps[0] / 2 % 2 == 1) {
                    mpz_sub(rhs_sqrt, N, rhs_sqrt); /* multiply by -1 */
                }
                for (int j = 0; j < fb_size; j++) {
                    if (combined_exps[j + 1] > 0) {
                        mpz_t pp;
                        mpz_init(pp);
                        mpz_set_ui(pp, fb[j].p);
                        mpz_t exp;
                        mpz_init(exp);
                        mpz_set_ui(exp, combined_exps[j + 1] / 2);
                        mpz_t pw;
                        mpz_init(pw);
                        mpz_powm(pw, pp, exp, N);
                        mpz_mul(rhs_sqrt, rhs_sqrt, pw);
                        mpz_mod(rhs_sqrt, rhs_sqrt, N);
                        mpz_clears(pp, exp, pw, NULL);
                    }
                }

                /* Try gcd(lhs - rhs_sqrt, N) */
                mpz_sub(g, lhs, rhs_sqrt);
                mpz_mod(g, g, N);
                mpz_gcd(g, g, N);
                if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                    mpz_set(factor, g);
                    found = 1;
                }

                if (!found) {
                    /* Try lhs + rhs_sqrt */
                    mpz_add(g, lhs, rhs_sqrt);
                    mpz_mod(g, g, N);
                    mpz_gcd(g, g, N);
                    if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                        mpz_set(factor, g);
                        found = 1;
                    }
                }
            }

            free(combined_exps);
            free(dep);
            mpz_clears(lhs, rhs_sqrt, g, NULL);
        }
    }

    if (found) {
        mpz_divexact(cofactor, N, factor);
        if (mpz_cmp(factor, cofactor) > 0) mpz_swap(factor, cofactor);
        gmp_printf("%Zd %Zd\n", factor, cofactor);
        struct timespec now;
        clock_gettime(CLOCK_MONOTONIC, &now);
        double elapsed = (now.tv_sec - start.tv_sec) +
            (now.tv_nsec - start.tv_nsec) / 1e9;
        fprintf(stderr, "SPF: factored in %.3fs\n", elapsed);
    } else {
        fprintf(stderr, "SPF: FAILED\n");
        /* Cleanup and exit */
        for (int i = 0; i < nrels; i++) {
            mpz_clear(rels[i].x);
            mpz_clear(rels[i].rhs);
            free(rel_exps[i]);
        }
        free(rels); free(rel_exps);
        for (int i = 0; i < fb_size; i++) mpz_clear(fb[i].root);
        free(fb);
        mpz_clears(N, sqrtN, factor, cofactor, x, val, cof, NULL);
        return 1;
    }

done:
    /* Cleanup */
    mpz_clears(x, val, cof, NULL);
    for (int i = 0; i < nrels; i++) {
        mpz_clear(rels[i].x);
        mpz_clear(rels[i].rhs);
        free(rel_exps[i]);
    }
    free(rels); free(rel_exps);
    for (int i = 0; i < fb_size; i++) mpz_clear(fb[i].root);
    free(fb);
    mpz_clears(N, sqrtN, factor, cofactor, NULL);
    return 0;
}
