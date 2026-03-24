/*
 * cofactor_analysis.c — Analyze cofactor distribution after smooth extraction
 *
 * For QS candidates Q(x) = (x+m)^2 - N, extract the B-smooth part and
 * study the cofactor distribution. This tells us:
 *
 * 1. What fraction of Q(x) values are "almost smooth" (cofactor is prime < B²)?
 * 2. What fraction have cofactors that could be factored by ECM cheaply?
 * 3. How does the cofactor size distribution change with N?
 *
 * If cofactors are typically much smaller than √N, a two-phase approach
 * (smooth extraction + ECM) could have better scaling than pure sieving.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

void analyze_cofactors(const char *N_str, unsigned long B, int M) {
    mpz_t N, sqrtN, m, x_val, qx, abs_qx, smooth, cofactor, primorial, g;
    mpz_inits(N, sqrtN, m, x_val, qx, abs_qx, smooth, cofactor, primorial, g, NULL);
    mpz_set_str(N, N_str, 10);
    int n_digits = mpz_sizeinbase(N, 10);

    mpz_sqrt(sqrtN, N);
    mpz_add_ui(m, sqrtN, 1);

    /* Build primorial with prime powers */
    mpz_t bound_est;
    mpz_init(bound_est);
    mpz_mul_ui(bound_est, sqrtN, 2 * M);

    mpz_set_ui(primorial, 1);
    mpz_t pp, pk;
    mpz_inits(pp, pk, NULL);
    mpz_set_ui(pp, 2);
    while (mpz_cmp_ui(pp, B) <= 0) {
        mpz_set(pk, pp);
        while (1) {
            mpz_t pk2;
            mpz_init(pk2);
            mpz_mul(pk2, pk, pp);
            if (mpz_cmp(pk2, bound_est) > 0) { mpz_clear(pk2); break; }
            mpz_set(pk, pk2);
            mpz_clear(pk2);
        }
        mpz_mul(primorial, primorial, pk);
        mpz_nextprime(pp, pp);
    }
    mpz_clears(pp, pk, bound_est, NULL);

    /* Histogram: cofactor size in bits */
    int max_bits = mpz_sizeinbase(sqrtN, 2) + 30;
    int *hist = calloc(max_bits + 1, sizeof(int));
    int n_fully_smooth = 0;
    int n_one_large_prime = 0;  /* cofactor is prime < B² */
    int n_two_large_prime = 0;  /* cofactor is product of 2 primes, each < B */
    int n_tested = 0;
    double total_cofactor_ratio = 0.0;

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    for (int x = 1; x <= M; x++) {
        mpz_add_ui(x_val, m, x);
        mpz_mul(qx, x_val, x_val);
        mpz_sub(qx, qx, N);
        mpz_abs(abs_qx, qx);

        /* Extract B-smooth part */
        mpz_set(cofactor, abs_qx);
        while (1) {
            mpz_gcd(g, cofactor, primorial);
            if (mpz_cmp_ui(g, 1) == 0) break;
            mpz_divexact(cofactor, cofactor, g);
        }

        n_tested++;
        int cof_bits = mpz_sizeinbase(cofactor, 2);
        int val_bits = mpz_sizeinbase(abs_qx, 2);

        if (cof_bits <= max_bits) hist[cof_bits]++;

        double ratio = (double)cof_bits / val_bits;
        total_cofactor_ratio += ratio;

        if (mpz_cmp_ui(cofactor, 1) == 0) {
            n_fully_smooth++;
        } else if (mpz_probab_prime_p(cofactor, 10)) {
            unsigned long B2 = B * B;
            mpz_t B2_mpz;
            mpz_init_set_ui(B2_mpz, B);
            mpz_mul_ui(B2_mpz, B2_mpz, B);
            if (mpz_cmp(cofactor, B2_mpz) < 0) {
                n_one_large_prime++;
            }
            mpz_clear(B2_mpz);
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;

    printf("=== Cofactor Analysis: %d digits, B=%lu, M=%d ===\n",
           n_digits, B, M);
    printf("Time: %.2fs\n", elapsed);
    printf("Fully B-smooth: %d / %d (%.4f%%)\n",
           n_fully_smooth, n_tested, 100.0 * n_fully_smooth / n_tested);
    printf("One large prime (< B²): %d / %d (%.4f%%)\n",
           n_one_large_prime, n_tested, 100.0 * n_one_large_prime / n_tested);
    printf("Avg cofactor ratio (cofactor_bits / value_bits): %.3f\n",
           total_cofactor_ratio / n_tested);

    /* Print histogram (compressed) */
    printf("\nCofactor size distribution (bits -> count):\n");
    int bucket_size = 5;
    for (int b = 0; b <= max_bits; b += bucket_size) {
        int count = 0;
        for (int i = b; i < b + bucket_size && i <= max_bits; i++)
            count += hist[i];
        if (count > 0)
            printf("  %3d-%3d bits: %d\n", b, b + bucket_size - 1, count);
    }

    free(hist);
    mpz_clears(N, sqrtN, m, x_val, qx, abs_qx, smooth, cofactor, primorial, g, NULL);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N_decimal> [B] [M]\n", argv[0]);
        return 1;
    }

    mpz_t N;
    mpz_init_set_str(N, argv[1], 10);
    int n_bits = mpz_sizeinbase(N, 2);

    /* Default B: L[1/2, 0.5] */
    double ln_N = n_bits * log(2.0);
    double ln_ln_N = log(ln_N);
    unsigned long B = argc > 2 ? atol(argv[2]) :
        (unsigned long)(exp(0.5 * sqrt(ln_N * ln_ln_N)));
    int M = argc > 3 ? atoi(argv[3]) : 500000;

    if (B < 1000) B = 1000;

    mpz_clear(N);

    analyze_cofactors(argv[1], B, M);
    return 0;
}
