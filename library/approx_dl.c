/*
 * approx_dl.c — Approximate discrete log relations for factoring
 *
 * NOVEL IDEA: Instead of exact DL (which is as hard as factoring),
 * compute APPROXIMATE DL relations and look for structural patterns.
 *
 * For small primes p_1, ..., p_k and fixed base g = 2:
 * Find exponents e_i such that 2^{e_i} ≡ p_i + ε_i (mod N)
 * where ε_i is "small".
 *
 * Since 2^{e_i} mod N = CRT(2^{e_i} mod p, 2^{e_i} mod q):
 * - 2^{e_i} mod p is well-defined (element of Z/pZ)
 * - 2^{e_i} mod q is well-defined (element of Z/qZ)
 *
 * For 2^{e_i} ≡ p_i (mod N): this holds iff both
 *   2^{e_i} ≡ p_i (mod p) AND 2^{e_i} ≡ p_i (mod q)
 *
 * The first congruence constrains e_i mod ord_p(2).
 * The second constrains e_i mod ord_q(2).
 *
 * If we find MANY such relations, the pattern of e_i values
 * reveals information about ord_p(2) and ord_q(2).
 *
 * For APPROXIMATE relations (ε_i ≠ 0):
 *   2^{e_i} ≡ p_i + ε_i (mod p): constrains e_i mod ord_p(2) given p_i, ε_i
 *   2^{e_i} ≡ p_i + ε_i (mod q): constrains e_i mod ord_q(2)
 *
 * The PATTERN of which e_i give small ε_i encodes information about
 * ord_p(2) and ord_q(2), which encode p and q.
 *
 * EXPERIMENT: Compute 2^e mod N for many e values and record when
 * the result is close to a small prime. Analyze the distribution
 * of successful e values.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N> [max_exponent] [tolerance]\n", argv[0]);
        return 1;
    }

    mpz_t N, g, result, target, diff, gcd_val;
    mpz_inits(N, g, result, target, diff, gcd_val, NULL);
    mpz_set_str(N, argv[1], 10);
    mpz_set_ui(g, 2); /* base */

    long max_exp = argc > 2 ? atol(argv[2]) : 10000000;
    long tolerance = argc > 3 ? atol(argv[3]) : 1000;

    int n_digits = mpz_sizeinbase(N, 10);

    /* Small primes to search for */
    unsigned long small_primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
                                     31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
                                     73, 79, 83, 89, 97};
    int n_primes = sizeof(small_primes) / sizeof(small_primes[0]);

    /* Compute 2^e mod N using repeated squaring */
    /* But we want to do this incrementally: 2^{e+1} = 2 * 2^e mod N */
    mpz_set_ui(result, 1);

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    int hits = 0;
    int factor_found = 0;

    for (long e = 0; e < max_exp && !factor_found; e++) {
        mpz_mul_ui(result, result, 2);
        mpz_mod(result, result, N);

        /* Check if result is close to any small prime */
        for (int i = 0; i < n_primes; i++) {
            mpz_set_ui(target, small_primes[i]);

            /* Check result ≈ target (mod N) */
            mpz_sub(diff, result, target);
            mpz_mod(diff, diff, N);

            /* diff should be small or close to N */
            if (mpz_cmp_ui(diff, tolerance) < 0 ||
                (mpz_sub(diff, N, diff), mpz_cmp_ui(diff, tolerance) < 0)) {

                /* Found approximate relation: 2^e ≈ p_i (mod N) */
                mpz_sub(diff, result, target);
                mpz_gcd(gcd_val, diff, N);

                if (mpz_cmp_ui(gcd_val, 1) > 0 && mpz_cmp(gcd_val, N) < 0) {
                    /* FACTOR FOUND via approximate DL! */
                    mpz_t cofactor;
                    mpz_init(cofactor);
                    mpz_divexact(cofactor, N, gcd_val);
                    gmp_printf("%Zd %Zd\n", gcd_val, cofactor);

                    clock_gettime(CLOCK_MONOTONIC, &t1);
                    double elapsed = (t1.tv_sec - t0.tv_sec) +
                                   (t1.tv_nsec - t0.tv_nsec) / 1e9;
                    fprintf(stderr, "ApproxDL: Factor found at e=%ld, target=%lu, "
                            "diff=", e, small_primes[i]);
                    gmp_fprintf(stderr, "%Zd", diff);
                    fprintf(stderr, " (%.3fs)\n", elapsed);
                    mpz_clear(cofactor);
                    factor_found = 1;
                    break;
                }

                hits++;
                if (hits <= 20) {
                    fprintf(stderr, "Hit: e=%ld, 2^e ≈ %lu (mod N), diff ≈ ",
                            e, small_primes[i]);
                    gmp_fprintf(stderr, "%Zd", diff);
                    fprintf(stderr, "\n");
                }
            }
        }

        if (e % 1000000 == 0 && e > 0) {
            clock_gettime(CLOCK_MONOTONIC, &t1);
            double elapsed = (t1.tv_sec - t0.tv_sec) +
                           (t1.tv_nsec - t0.tv_nsec) / 1e9;
            fprintf(stderr, "ApproxDL: e=%ld, hits=%d, %.1fs\n", e, hits, elapsed);
        }
    }

    if (!factor_found) {
        clock_gettime(CLOCK_MONOTONIC, &t1);
        double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
        fprintf(stderr, "ApproxDL: No factor found after %ld steps (%d hits, %.1fs)\n",
                max_exp, hits, elapsed);
    }

    mpz_clears(N, g, result, target, diff, gcd_val, NULL);
    return factor_found ? 0 : 1;
}
