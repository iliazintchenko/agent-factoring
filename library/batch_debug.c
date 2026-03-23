/*
 * batch_debug.c - Debug: compare batch smoothness vs trial division
 * Tests whether the product tree / remainder tree approach actually
 * finds the same smooth numbers as brute force trial division.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

#define SEED 42

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    mpz_t N, kN;
    mpz_init(N);
    mpz_init(kN);
    mpz_set_str(N, argv[1], 10);
    mpz_set(kN, N); /* k=1 for simplicity */

    int digits = mpz_sizeinbase(N, 10);
    double ln_n = digits * log(10);
    int smooth_bound = (int)(exp(0.5 * sqrt(ln_n * log(ln_n))) * 0.7);
    if (smooth_bound < 300) smooth_bound = 300;

    /* Generate primes up to smooth_bound */
    char *sieve = calloc(smooth_bound + 1, 1);
    int *primes = malloc(sizeof(int) * smooth_bound);
    int pcount = 0;
    for (int i = 2; i <= smooth_bound; i++) {
        if (!sieve[i]) {
            primes[pcount++] = i;
            for (long j = (long)i*i; j <= smooth_bound; j += i)
                sieve[j] = 1;
        }
    }
    free(sieve);

    fprintf(stderr, "%d primes up to %d\n", pcount, smooth_bound);

    /* Compute primorial = product of prime powers */
    mpz_t primorial;
    mpz_init_set_ui(primorial, 1);
    for (int i = 0; i < pcount; i++) {
        unsigned long pk = primes[i];
        while (pk <= 1000000000000000UL / primes[i]) pk *= primes[i];
        mpz_mul_ui(primorial, primorial, pk);
    }
    fprintf(stderr, "Primorial: %zu bits\n", mpz_sizeinbase(primorial, 2));

    /* Generate some Q(x) = (x + sqrt(N))^2 - N values */
    mpz_t sqrtN, qx, tmp;
    mpz_init(sqrtN);
    mpz_init(qx);
    mpz_init(tmp);
    mpz_sqrt(sqrtN, N);

    unsigned long lp_bound = smooth_bound * 50;

    int total = 1000;
    int td_smooth = 0, td_slp = 0;
    int batch_smooth = 0, batch_slp = 0;
    int both_smooth = 0;

    for (int x = 1; x <= total; x++) {
        /* Q(x) = (x + sqrtN)^2 - N */
        mpz_add_ui(qx, sqrtN, x);
        mpz_mul(qx, qx, qx);
        mpz_sub(qx, qx, N);

        /* Trial division */
        mpz_abs(tmp, qx);
        for (int i = 0; i < pcount; i++) {
            while (mpz_divisible_ui_p(tmp, primes[i]))
                mpz_divexact_ui(tmp, tmp, primes[i]);
        }
        int td_result = 0;
        if (mpz_cmp_ui(tmp, 1) == 0) { td_smooth++; td_result = 1; }
        else if (mpz_fits_ulong_p(tmp) && mpz_get_ui(tmp) <= lp_bound) { td_slp++; td_result = 2; }

        /* Batch method: GCD with primorial */
        mpz_abs(tmp, qx);
        /* Iterated GCD */
        int iters = 0;
        while (iters < 30) {
            mpz_t g;
            mpz_init(g);
            mpz_gcd(g, tmp, primorial);
            if (mpz_cmp_ui(g, 1) <= 0) { mpz_clear(g); break; }
            mpz_divexact(tmp, tmp, g);
            mpz_clear(g);
            iters++;
        }

        int batch_result = 0;
        if (mpz_cmp_ui(tmp, 1) == 0) { batch_smooth++; batch_result = 1; }
        else if (mpz_fits_ulong_p(tmp) && mpz_get_ui(tmp) <= lp_bound) { batch_slp++; batch_result = 2; }

        if (td_result && batch_result) both_smooth++;
        if (td_result && !batch_result) {
            gmp_fprintf(stderr, "x=%d: TD found (type %d) but batch missed. Remaining: %Zd (%zu bits)\n",
                       x, td_result, tmp, mpz_sizeinbase(tmp, 2));
        }
        if (!td_result && batch_result) {
            fprintf(stderr, "x=%d: batch found but TD missed?!\n", x);
        }
    }

    /* Also check QS Q(x) values are of expected size */
    mpz_add_ui(qx, sqrtN, 1);
    mpz_mul(qx, qx, qx);
    mpz_sub(qx, qx, N);
    fprintf(stderr, "\nQ(1) size: %zu bits (%zu digits)\n",
            mpz_sizeinbase(qx, 2), mpz_sizeinbase(qx, 10));

    mpz_add_ui(qx, sqrtN, 1000);
    mpz_mul(qx, qx, qx);
    mpz_sub(qx, qx, N);
    fprintf(stderr, "Q(1000) size: %zu bits (%zu digits)\n",
            mpz_sizeinbase(qx, 2), mpz_sizeinbase(qx, 10));

    mpz_add_ui(qx, sqrtN, 10000);
    mpz_mul(qx, qx, qx);
    mpz_sub(qx, qx, N);
    fprintf(stderr, "Q(10000) size: %zu bits (%zu digits)\n\n",
            mpz_sizeinbase(qx, 2), mpz_sizeinbase(qx, 10));

    fprintf(stderr, "Results for first %d Q(x) values:\n", total);
    fprintf(stderr, "  Trial division: %d smooth, %d SLP\n", td_smooth, td_slp);
    fprintf(stderr, "  Batch GCD:      %d smooth, %d SLP\n", batch_smooth, batch_slp);
    fprintf(stderr, "  Both found:     %d\n", both_smooth);
    fprintf(stderr, "  Smoothness probability: %.2f%%\n", 100.0 * (td_smooth + td_slp) / total);

    free(primes);
    mpz_clear(N); mpz_clear(kN); mpz_clear(primorial);
    mpz_clear(sqrtN); mpz_clear(qx); mpz_clear(tmp);
    return 0;
}
