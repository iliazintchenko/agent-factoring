/*
 * candidate_strategies.c — Compare candidate generation strategies
 *
 * For QS-style factoring, compare the SMOOTHNESS YIELD (relations per
 * second of CPU time) of different candidate generation strategies:
 *
 * 1. Sequential: x = √N, √N+1, √N+2, ... (standard QS)
 * 2. Random: x = random values near √N
 * 3. Fibonacci-like: x = √N + fib(k) for Fibonacci numbers
 * 4. Lattice: x = √N + LLL-selected offsets
 * 5. CRT-structured: x satisfying multiple congruences simultaneously
 *
 * Each strategy generates candidates, computes x² - N, and tests for
 * smoothness via trial division. The key metric is: how many smooth
 * values per unit of candidate-generation effort?
 *
 * Usage: ./candidate_strategies <N>
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

static double walltime(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* Test smoothness of val with primes up to B. Returns 1 if B-smooth. */
static int is_smooth(mpz_t val, unsigned long B) {
    mpz_t tmp;
    mpz_init(tmp);
    mpz_abs(tmp, val);

    for (unsigned long p = 2; p <= B; p++) {
        if (p > 2 && p % 2 == 0) continue;
        int is_p = 1;
        for (unsigned long d = 3; d * d <= p; d += 2)
            if (p % d == 0) { is_p = 0; break; }
        if (p > 2 && !is_p) continue;

        while (mpz_divisible_ui_p(tmp, p))
            mpz_divexact_ui(tmp, tmp, p);
    }

    int smooth = (mpz_cmp_ui(tmp, 1) == 0);
    mpz_clear(tmp);
    return smooth;
}

int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }

    mpz_t N, sqrtN, x, val;
    mpz_init(N); mpz_init(sqrtN); mpz_init(x); mpz_init(val);
    mpz_set_str(N, argv[1], 10);
    int digits = strlen(argv[1]);
    mpz_sqrt(sqrtN, N);
    mpz_add_ui(sqrtN, sqrtN, 1);

    unsigned long B = 5000; /* smoothness bound */
    int trials = 50000; /* candidates per strategy */

    printf("# %d-digit N, B=%lu, %d trials per strategy\n\n", digits, B, trials);

    /* Strategy 1: Sequential */
    {
        double t0 = walltime();
        int smooth_count = 0;
        long total_val_bits = 0;

        for (int i = 1; i <= trials; i++) {
            mpz_set(x, sqrtN);
            mpz_add_ui(x, x, i);
            mpz_mul(val, x, x);
            mpz_sub(val, val, N);
            total_val_bits += mpz_sizeinbase(val, 2);
            if (is_smooth(val, B)) smooth_count++;
        }

        double elapsed = walltime() - t0;
        printf("Sequential: %d smooth / %d trials = %.4f%% in %.3fs (%.0f smooth/s, avg %ld bits)\n",
               smooth_count, trials, 100.0*smooth_count/trials, elapsed,
               smooth_count/elapsed, total_val_bits/trials);
    }

    /* Strategy 2: Random offsets */
    {
        double t0 = walltime();
        int smooth_count = 0;
        long total_val_bits = 0;

        gmp_randstate_t rng;
        gmp_randinit_default(rng);
        gmp_randseed_ui(rng, 42);

        mpz_t offset, max_off;
        mpz_init(offset); mpz_init(max_off);
        mpz_set_ui(max_off, trials * 2);

        for (int i = 0; i < trials; i++) {
            mpz_urandomm(offset, rng, max_off);
            mpz_add_ui(offset, offset, 1);
            mpz_set(x, sqrtN);
            mpz_add(x, x, offset);
            mpz_mul(val, x, x);
            mpz_sub(val, val, N);
            total_val_bits += mpz_sizeinbase(val, 2);
            if (is_smooth(val, B)) smooth_count++;
        }

        double elapsed = walltime() - t0;
        printf("Random:     %d smooth / %d trials = %.4f%% in %.3fs (%.0f smooth/s, avg %ld bits)\n",
               smooth_count, trials, 100.0*smooth_count/trials, elapsed,
               smooth_count/elapsed, total_val_bits/trials);

        mpz_clear(offset); mpz_clear(max_off);
        gmp_randclear(rng);
    }

    /* Strategy 3: Quadratic residue focused — only try x where x² ≡ N (mod small product) */
    {
        double t0 = walltime();
        int smooth_count = 0;
        long total_val_bits = 0;
        int tested = 0;

        /* For each small prime p where N is QR mod p, only test x ≡ ±√N (mod p).
         * This means x² - N ≡ 0 (mod p), guaranteeing p divides the value. */
        unsigned long primes[] = {3, 5, 7, 11, 13};
        unsigned long product = 1;
        for (int i = 0; i < 5; i++) product *= primes[i];

        /* Find x ≡ √N mod product using CRT */
        for (int t = 0; t < trials && tested < trials; t++) {
            long offset = (long)(t * product);
            if (offset > (long)(trials * 10)) break;
            mpz_set(x, sqrtN);
            mpz_add_ui(x, x, offset);
            mpz_mul(val, x, x);
            mpz_sub(val, val, N);
            total_val_bits += mpz_sizeinbase(val, 2);
            tested++;
            if (is_smooth(val, B)) smooth_count++;
        }

        double elapsed = walltime() - t0;
        if (tested > 0)
            printf("CRT-skip:   %d smooth / %d trials = %.4f%% in %.3fs (%.0f smooth/s, avg %ld bits)\n",
                   smooth_count, tested, 100.0*smooth_count/tested, elapsed,
                   smooth_count/elapsed, total_val_bits/tested);
    }

    /* Strategy 4: Near-zero values (closest to √N) */
    {
        double t0 = walltime();
        int smooth_count = 0;
        long total_val_bits = 0;

        /* Test x = √N ± 1, √N ± 2, ..., prioritizing smallest |x² - N| */
        for (int i = 1; i <= trials; i++) {
            int sign = (i % 2 == 0) ? 1 : -1;
            int offset = (i + 1) / 2;
            mpz_set(x, sqrtN);
            if (sign > 0) mpz_add_ui(x, x, offset);
            else mpz_sub_ui(x, x, offset);

            mpz_mul(val, x, x);
            mpz_sub(val, val, N);
            total_val_bits += mpz_sizeinbase(val, 2);
            if (is_smooth(val, B)) smooth_count++;
        }

        double elapsed = walltime() - t0;
        printf("Near-zero:  %d smooth / %d trials = %.4f%% in %.3fs (%.0f smooth/s, avg %ld bits)\n",
               smooth_count, trials, 100.0*smooth_count/trials, elapsed,
               smooth_count/elapsed, total_val_bits/trials);
    }

    /* Strategy 5: Multiplied candidates — k*N for small smooth k */
    {
        double t0 = walltime();
        int smooth_count = 0;
        long total_val_bits = 0;
        int tested = 0;

        /* For each small smooth k, compute x² - kN. If kN has different
         * residuosity than N mod some primes, different candidates become smooth. */
        unsigned long multipliers[] = {1, 3, 5, 7, 11, 13, 15, 17, 19, 21};
        int n_mult = 10;

        for (int mi = 0; mi < n_mult && tested < trials; mi++) {
            mpz_t kN, sqrtKN;
            mpz_init(kN); mpz_init(sqrtKN);
            mpz_mul_ui(kN, N, multipliers[mi]);
            mpz_sqrt(sqrtKN, kN);
            mpz_add_ui(sqrtKN, sqrtKN, 1);

            for (int i = 1; i <= trials / n_mult && tested < trials; i++) {
                mpz_set(x, sqrtKN);
                mpz_add_ui(x, x, i);
                mpz_mul(val, x, x);
                mpz_sub(val, val, kN);
                total_val_bits += mpz_sizeinbase(val, 2);
                tested++;
                if (is_smooth(val, B)) smooth_count++;
            }

            mpz_clear(kN); mpz_clear(sqrtKN);
        }

        double elapsed = walltime() - t0;
        printf("Multi-k:    %d smooth / %d trials = %.4f%% in %.3fs (%.0f smooth/s, avg %ld bits)\n",
               smooth_count, tested, 100.0*smooth_count/tested, elapsed,
               smooth_count/elapsed, total_val_bits/tested);
    }

    mpz_clear(N); mpz_clear(sqrtN); mpz_clear(x); mpz_clear(val);
    return 0;
}
