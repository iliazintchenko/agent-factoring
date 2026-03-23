/*
 * smoothness_profile.c — Measure partial-smoothness distribution of QS polynomial values
 *
 * For N given on command line, evaluates Q(x) = (x + ceil(sqrt(N)))^2 - N
 * for x in [-M, M], trial divides by primes up to B, and reports how many
 * values have 0, 1, 2, ... large prime factors.
 *
 * Usage: ./smoothness_profile <N> [B] [M]
 * Defaults: B = 100000, M = 500000
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

/* Simple prime sieve up to limit */
static int *sieve_primes(int limit, int *count) {
    char *is_composite = calloc(limit + 1, 1);
    int *primes = malloc(sizeof(int) * (limit / 2 + 10));
    int cnt = 0;
    for (int i = 2; i <= limit; i++) {
        if (!is_composite[i]) {
            primes[cnt++] = i;
            if ((long long)i * i <= limit) {
                for (int j = i * i; j <= limit; j += i)
                    is_composite[j] = 1;
            }
        }
    }
    free(is_composite);
    *count = cnt;
    return primes;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N> [B] [M]\n", argv[0]);
        return 1;
    }

    mpz_t N, sqrtN, qx, rem, cofactor;
    mpz_inits(N, sqrtN, qx, rem, cofactor, NULL);
    mpz_set_str(N, argv[1], 10);

    int B = (argc > 2) ? atoi(argv[2]) : 100000;
    long M = (argc > 3) ? atol(argv[3]) : 500000;

    /* Compute ceil(sqrt(N)) */
    mpz_sqrt(sqrtN, N);
    mpz_mul(qx, sqrtN, sqrtN);
    if (mpz_cmp(qx, N) < 0)
        mpz_add_ui(sqrtN, sqrtN, 1);

    /* Generate primes up to B */
    int nprimes;
    int *primes = sieve_primes(B, &nprimes);

    /* Find which primes have N as a QR (for the factor base) */
    int *fb_primes = malloc(sizeof(int) * nprimes);
    int *fb_roots = malloc(sizeof(int) * nprimes);  /* sqrt(N) mod p */
    int fb_size = 0;

    /* Always include 2 */
    fb_primes[fb_size] = 2;
    fb_roots[fb_size] = 1;
    fb_size++;

    mpz_t tmp;
    mpz_init(tmp);
    for (int i = 1; i < nprimes; i++) {  /* skip 2, already added */
        int p = primes[i];
        /* Check if N is a QR mod p */
        mpz_set_ui(tmp, p);
        int legendre = mpz_legendre(N, tmp);
        if (legendre == 1) {
            fb_primes[fb_size] = p;
            /* Compute sqrt(N) mod p using Tonelli-Shanks (via GMP) */
            /* We'll use a simpler approach: brute force for small p, GMP for larger */
            unsigned long n_mod_p = mpz_fdiv_ui(N, p);
            /* Find r such that r^2 = n_mod_p (mod p) */
            unsigned long r = 0;
            for (unsigned long x = 0; (unsigned long)x < (unsigned long)p; x++) {
                if ((x * x) % p == n_mod_p) { r = x; break; }
            }
            fb_roots[fb_size] = (int)r;
            fb_size++;
        }
    }
    mpz_clear(tmp);
    free(primes);

    fprintf(stderr, "N = %s\n", argv[1]);
    fprintf(stderr, "Factor base: %d primes up to %d\n", fb_size, B);
    fprintf(stderr, "Sieve range: [-%ld, %ld]\n", M, M);

    /* Profile: for each x in [-M, M], compute Q(x) and trial divide */
    int max_large = 10;
    long *counts = calloc(max_large + 1, sizeof(long));
    long total = 0;
    long smooth_count = 0;

    mpz_t xval, base;
    mpz_inits(xval, base, NULL);

    for (long x = -M; x <= M; x++) {
        /* Q(x) = (x + sqrtN)^2 - N */
        mpz_set_si(xval, x);
        mpz_add(base, xval, sqrtN);
        mpz_mul(qx, base, base);
        mpz_sub(qx, qx, N);

        /* Take absolute value */
        int sign = mpz_sgn(qx);
        if (sign == 0) continue;
        if (sign < 0) mpz_neg(qx, qx);

        /* Trial divide by factor base */
        mpz_set(cofactor, qx);
        for (int i = 0; i < fb_size; i++) {
            unsigned long p = fb_primes[i];
            while (mpz_divisible_ui_p(cofactor, p)) {
                mpz_divexact_ui(cofactor, cofactor, p);
            }
        }

        /* Count number of "large" prime factors in cofactor */
        /* cofactor = product of primes > B */
        int n_large;
        if (mpz_cmp_ui(cofactor, 1) == 0) {
            n_large = 0;  /* fully smooth */
            smooth_count++;
        } else {
            /* Estimate number of large prime factors from size */
            double log_cofactor = mpz_sizeinbase(cofactor, 2) * 0.693;
            double log_B = log((double)B);
            /* Rough estimate: n_large ≈ log(cofactor) / log(B) */
            n_large = (int)(log_cofactor / log_B + 0.5);
            if (n_large < 1) n_large = 1;
            if (n_large > max_large) n_large = max_large;
        }
        if (n_large <= max_large)
            counts[n_large]++;
        total++;
    }

    printf("PROFILE N=%s digits=%zu B=%d M=%ld\n",
           argv[1], strlen(argv[1]), B, M);
    printf("total_values=%ld factor_base_size=%d\n", total, fb_size);
    printf("fully_smooth=%ld (%.6f%%)\n", smooth_count, 100.0 * smooth_count / total);
    for (int k = 0; k <= max_large; k++) {
        printf("large_primes<=%d: %ld (%.6f%%)\n", k, counts[k], 100.0 * counts[k] / total);
    }

    /* Cleanup */
    free(counts);
    free(fb_primes);
    free(fb_roots);
    mpz_clears(N, sqrtN, qx, rem, cofactor, xval, base, NULL);
    return 0;
}
