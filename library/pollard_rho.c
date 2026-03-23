/*
 * Pollard's rho factoring with Brent's improvement.
 * Usage: ./pollard_rho <number>
 * Single-threaded, seed=42.
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <number>\n", argv[0]);
        return 1;
    }

    mpz_t n, x, y, d, q, ys, r, tmp;
    mpz_init(n); mpz_init(x); mpz_init(y); mpz_init(d);
    mpz_init(q); mpz_init(ys); mpz_init(r); mpz_init(tmp);

    if (mpz_set_str(n, argv[1], 10) != 0) {
        fprintf(stderr, "Invalid number\n");
        return 1;
    }

    struct timespec tstart;
    clock_gettime(CLOCK_MONOTONIC, &tstart);

    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, 42);

    int found = 0;
    for (int attempt = 0; attempt < 100 && !found; attempt++) {
        /* Random starting point and c */
        mpz_urandomm(x, rstate, n);
        mpz_t c;
        mpz_init(c);
        mpz_urandomm(c, rstate, n);
        if (mpz_cmp_ui(c, 0) == 0) mpz_set_ui(c, 1);

        mpz_set(y, x);
        mpz_set_ui(d, 1);

        /* Brent's cycle detection with batch GCD */
        unsigned long m = 128;
        unsigned long r_val = 1;

        while (mpz_cmp_ui(d, 1) == 0) {
            mpz_set(ys, y);
            for (unsigned long i = 0; i < r_val; i++) {
                mpz_mul(y, y, y);
                mpz_add(y, y, c);
                mpz_mod(y, y, n);
            }

            unsigned long k = 0;
            while (k < r_val && mpz_cmp_ui(d, 1) == 0) {
                mpz_set(ys, y);
                mpz_set_ui(q, 1);
                unsigned long batch = (m < r_val - k) ? m : r_val - k;

                for (unsigned long j = 0; j < batch; j++) {
                    mpz_mul(y, y, y);
                    mpz_add(y, y, c);
                    mpz_mod(y, y, n);
                    mpz_sub(tmp, x, y);
                    mpz_abs(tmp, tmp);
                    mpz_mul(q, q, tmp);
                    mpz_mod(q, q, n);
                }

                mpz_gcd(d, q, n);
                k += batch;
            }

            if (mpz_cmp(d, n) == 0) {
                /* Backtrack */
                mpz_set_ui(d, 1);
                mpz_set(y, ys);
                while (mpz_cmp_ui(d, 1) == 0) {
                    mpz_mul(y, y, y);
                    mpz_add(y, y, c);
                    mpz_mod(y, y, n);
                    mpz_sub(tmp, x, y);
                    mpz_abs(tmp, tmp);
                    mpz_gcd(d, tmp, n);
                }
                if (mpz_cmp(d, n) == 0) {
                    mpz_set_ui(d, 1);
                    break; /* try new c */
                }
            }

            mpz_set(x, y);
            r_val *= 2;

            /* Timeout check every so often */
            if (r_val > 1000000) {
                struct timespec tnow;
                clock_gettime(CLOCK_MONOTONIC, &tnow);
                double elapsed = (tnow.tv_sec - tstart.tv_sec) + (tnow.tv_nsec - tstart.tv_nsec) / 1e9;
                if (elapsed > 280) {
                    fprintf(stderr, "Pollard rho timeout after %.1fs\n", elapsed);
                    printf("FAIL\n");
                    return 1;
                }
            }
        }
        mpz_clear(c);

        if (mpz_cmp_ui(d, 1) > 0 && mpz_cmp(d, n) < 0) {
            struct timespec tnow;
            clock_gettime(CLOCK_MONOTONIC, &tnow);
            double elapsed = (tnow.tv_sec - tstart.tv_sec) + (tnow.tv_nsec - tstart.tv_nsec) / 1e9;
            gmp_printf("FACTOR: %Zd\n", d);
            mpz_t cofactor; mpz_init(cofactor);
            mpz_divexact(cofactor, n, d);
            gmp_printf("COFACTOR: %Zd\n", cofactor);
            fprintf(stderr, "Pollard rho found factor in %.3fs (attempt %d)\n", elapsed, attempt);
            mpz_clear(cofactor);
            found = 1;
        }
    }

    if (!found) {
        fprintf(stderr, "Pollard rho failed\n");
        printf("FAIL\n");
    }

    gmp_randclear(rstate);
    mpz_clear(n); mpz_clear(x); mpz_clear(y); mpz_clear(d);
    mpz_clear(q); mpz_clear(ys); mpz_clear(r); mpz_clear(tmp);
    return found ? 0 : 1;
}
