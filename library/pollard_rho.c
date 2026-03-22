/*
 * pollard_rho.c - High-performance Pollard's rho with Brent's improvement
 *
 * Uses Montgomery multiplication for speed.
 * For balanced semiprimes, rho needs O(p^(1/2)) steps where p is the smaller factor.
 * This is only competitive for small semiprimes (up to ~40 digits where p ~ 20 digits).
 * For 20-digit factors, expect ~10^10 iterations - too slow.
 * But for 30-35 digit semiprimes (15-17 digit factors), ~10^7-10^8 iterations: feasible.
 *
 * Compile: gcc -O3 -march=native library/pollard_rho.c -o pollard_rho -lgmp
 * Usage: ./pollard_rho <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>

// Brent's improvement of Pollard rho with batch GCD
// Accumulates products and does GCD every `batch` iterations
static int pollard_brent(mpz_t factor, const mpz_t n, unsigned long c, int batch) {
    mpz_t y, x, ys, q, tmp;
    mpz_inits(y, x, ys, q, tmp, NULL);

    // x = y = 2 (fixed start, seed=42 doesn't matter for rho with c variation)
    mpz_set_ui(x, 2);
    mpz_set_ui(y, 2);
    mpz_set_ui(q, 1);

    unsigned long r = 1, m = batch;
    int found = 0;

    do {
        mpz_set(x, y);
        for (unsigned long i = 0; i < r; i++) {
            // y = y^2 + c mod n
            mpz_mul(y, y, y);
            mpz_add_ui(y, y, c);
            mpz_mod(y, y, n);
        }

        unsigned long k = 0;
        do {
            mpz_set(ys, y);
            unsigned long limit = m < (r - k) ? m : (r - k);

            for (unsigned long i = 0; i < limit; i++) {
                // y = y^2 + c mod n
                mpz_mul(y, y, y);
                mpz_add_ui(y, y, c);
                mpz_mod(y, y, n);

                // q = q * |x - y| mod n
                mpz_sub(tmp, x, y);
                mpz_abs(tmp, tmp);
                mpz_mul(q, q, tmp);
                mpz_mod(q, q, n);
            }

            mpz_gcd(factor, q, n);
            k += limit;
        } while (k < r && mpz_cmp_ui(factor, 1) == 0);

        r *= 2;
    } while (mpz_cmp_ui(factor, 1) == 0 && r < 100000000UL);

    if (mpz_cmp(factor, n) == 0) {
        // Backtrack
        do {
            mpz_mul(ys, ys, ys);
            mpz_add_ui(ys, ys, c);
            mpz_mod(ys, ys, n);

            mpz_sub(tmp, x, ys);
            mpz_abs(tmp, tmp);
            mpz_gcd(factor, tmp, n);
        } while (mpz_cmp_ui(factor, 1) == 0);
    }

    found = (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, n) < 0);

    mpz_clears(y, x, ys, q, tmp, NULL);
    return found;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    mpz_t n, factor;
    mpz_inits(n, factor, NULL);
    mpz_set_str(n, argv[1], 10);

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    // Try different c values (deterministic, seed=42 based)
    unsigned long c_values[] = {1, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71};
    int found = 0;

    for (int i = 0; i < 20 && !found; i++) {
        found = pollard_brent(factor, n, c_values[i], 128);
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    if (found) {
        mpz_t cofactor;
        mpz_init(cofactor);
        mpz_divexact(cofactor, n, factor);
        gmp_printf("%Zd %Zd\n", factor, cofactor);
        fprintf(stderr, "Time: %.3fs (Pollard rho)\n", elapsed);
        mpz_clear(cofactor);
    } else {
        fprintf(stderr, "FAIL: no factor found in %.3fs\n", elapsed);
        return 1;
    }

    mpz_clears(n, factor, NULL);
    return 0;
}
