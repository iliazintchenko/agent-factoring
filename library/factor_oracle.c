/*
 * factor_oracle.c - Multi-strategy factoring oracle
 *
 * Tries cheap methods first before falling back to SIQS:
 * 1. Trial division (tiny factors)
 * 2. Fermat's method (close factors, |p-q| < N^1/4 or so)
 * 3. Pollard rho (factors up to ~10^13 quickly)
 * 4. SQUFOF (factors up to ~10^18 in O(N^1/4))
 * 5. ECM with escalating curves (medium factors)
 * 6. p-1 / p+1 (special structure)
 * 7. Fall through to YAFU SIQS
 *
 * For balanced semiprimes, most of these will fail quickly (O(seconds))
 * and we fall to SIQS. But occasionally one hits and saves huge time.
 *
 * Compile: gcc -O3 -march=native -o factor_oracle library/factor_oracle.c -lgmp -lecm -lm
 * Usage: ./factor_oracle <N>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gmp.h>
#include <ecm.h>

static struct timespec g_start;

static double elapsed_sec(void) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec - g_start.tv_sec) + (now.tv_nsec - g_start.tv_nsec) / 1e9;
}

/* Trial division up to limit */
static int trial_div(mpz_t N, mpz_t factor, unsigned long limit) {
    for (unsigned long p = 2; p <= limit; p++) {
        if (mpz_divisible_ui_p(N, p)) {
            mpz_set_ui(factor, p);
            return 1;
        }
    }
    return 0;
}

/* Fermat's method: search for s such that s^2 - N = t^2 */
static int fermat_method(mpz_t N, mpz_t factor, unsigned long max_iter) {
    mpz_t s, s2, diff, t;
    mpz_inits(s, s2, diff, t, NULL);

    mpz_sqrt(s, N);
    mpz_add_ui(s, s, 1); /* ceil(sqrt(N)) */

    for (unsigned long i = 0; i < max_iter; i++) {
        mpz_mul(s2, s, s);
        mpz_sub(diff, s2, N);

        mpz_sqrt(t, diff);
        mpz_mul(s2, t, t);
        if (mpz_cmp(s2, diff) == 0) {
            /* s^2 - N = t^2 => N = (s-t)(s+t) */
            mpz_sub(factor, s, t);
            if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0) {
                mpz_clears(s, s2, diff, t, NULL);
                return 1;
            }
        }
        mpz_add_ui(s, s, 1);
    }
    mpz_clears(s, s2, diff, t, NULL);
    return 0;
}

/* Pollard rho (Brent variant) */
static int pollard_rho(mpz_t N, mpz_t factor, unsigned long max_iter) {
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, 42);

    mpz_t x, y, d, diff, prod;
    mpz_inits(x, y, d, diff, prod, NULL);

    for (int attempt = 0; attempt < 5; attempt++) {
        mpz_urandomm(x, rstate, N);
        mpz_set(y, x);
        unsigned long c = 1 + attempt;

        unsigned long m = 128;
        unsigned long iter = 0;

        while (iter < max_iter) {
            mpz_set(d, y);
            mpz_set_ui(prod, 1);

            for (unsigned long i = 0; i < m && iter < max_iter; i++, iter++) {
                /* y = y^2 + c mod N */
                mpz_mul(y, y, y);
                mpz_add_ui(y, y, c);
                mpz_mod(y, y, N);

                /* Batch GCD: prod *= |x - y| mod N */
                mpz_sub(diff, x, y);
                mpz_abs(diff, diff);
                mpz_mul(prod, prod, diff);
                mpz_mod(prod, prod, N);
            }

            mpz_gcd(d, prod, N);
            if (mpz_cmp_ui(d, 1) > 0 && mpz_cmp(d, N) < 0) {
                mpz_set(factor, d);
                mpz_clears(x, y, d, diff, prod, NULL);
                gmp_randclear(rstate);
                return 1;
            }
            if (mpz_cmp(d, N) == 0) break; /* Backtrack needed */

            mpz_set(x, y);
        }
    }

    mpz_clears(x, y, d, diff, prod, NULL);
    gmp_randclear(rstate);
    return 0;
}

/* ECM with GMP-ECM library */
static int ecm_method(mpz_t N, mpz_t factor, double B1, int curves) {
    ecm_params params;
    ecm_init(params);

    for (int i = 0; i < curves; i++) {
        if (elapsed_sec() > 280.0) break;

        mpz_set_ui(params->sigma, 42 + i);
        params->method = ECM_ECM;

        int ret = ecm_factor(factor, N, B1, params);
        if (ret > 0 && mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0) {
            ecm_clear(params);
            return 1;
        }
    }

    ecm_clear(params);
    return 0;
}

/* P-1 method */
static int pm1_method(mpz_t N, mpz_t factor, double B1) {
    ecm_params params;
    ecm_init(params);
    params->method = ECM_PM1;
    mpz_set_ui(params->sigma, 42);

    int ret = ecm_factor(factor, N, B1, params);
    ecm_clear(params);

    return (ret > 0 && mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0);
}

/* P+1 method */
static int pp1_method(mpz_t N, mpz_t factor, double B1) {
    ecm_params params;
    ecm_init(params);
    params->method = ECM_PP1;
    mpz_set_ui(params->sigma, 42);

    int ret = ecm_factor(factor, N, B1, params);
    ecm_clear(params);

    return (ret > 0 && mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0);
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    clock_gettime(CLOCK_MONOTONIC, &g_start);

    mpz_t N, factor, cofactor;
    mpz_inits(N, factor, cofactor, NULL);
    mpz_set_str(N, argv[1], 10);

    int digits = mpz_sizeinbase(N, 10);

    /* Stage 1: Trial division (< 1ms) */
    if (trial_div(N, factor, 1000000)) {
        mpz_divexact(cofactor, N, factor);
        gmp_printf("%Zd = %Zd * %Zd\n", N, factor, cofactor);
        fprintf(stderr, "Trial division: %.3fs\n", elapsed_sec());
        goto done;
    }

    /* Stage 2: Fermat's method (close factors, ~1ms) */
    /* For balanced semiprimes, |p-q| < N^(1/4) means at most ~N^(1/4) iterations */
    {
        unsigned long fermat_limit = 100000; /* Quick check */
        if (fermat_method(N, factor, fermat_limit)) {
            mpz_divexact(cofactor, N, factor);
            gmp_printf("%Zd = %Zd * %Zd\n", N, factor, cofactor);
            fprintf(stderr, "Fermat: %.3fs\n", elapsed_sec());
            goto done;
        }
    }

    /* Stage 3: Pollard rho (~100ms for 10^13 factors) */
    if (digits <= 50) {
        if (pollard_rho(N, factor, 10000000)) {
            mpz_divexact(cofactor, N, factor);
            gmp_printf("%Zd = %Zd * %Zd\n", N, factor, cofactor);
            fprintf(stderr, "Pollard rho: %.3fs\n", elapsed_sec());
            goto done;
        }
    }

    /* Stage 4: P-1 with escalating bounds */
    {
        double pm1_bounds[] = {1e5, 1e6, 1e7, 5e7, 0};
        for (int i = 0; pm1_bounds[i] > 0; i++) {
            if (elapsed_sec() > 10.0) break; /* Don't spend more than 10s on p-1 */
            if (pm1_method(N, factor, pm1_bounds[i])) {
                mpz_divexact(cofactor, N, factor);
                gmp_printf("%Zd = %Zd * %Zd\n", N, factor, cofactor);
                fprintf(stderr, "P-1 (B1=%.0f): %.3fs\n", pm1_bounds[i], elapsed_sec());
                goto done;
            }
        }
    }

    /* Stage 5: P+1 */
    {
        double pp1_bounds[] = {1e5, 1e6, 1e7, 0};
        for (int i = 0; pp1_bounds[i] > 0; i++) {
            if (elapsed_sec() > 15.0) break;
            if (pp1_method(N, factor, pp1_bounds[i])) {
                mpz_divexact(cofactor, N, factor);
                gmp_printf("%Zd = %Zd * %Zd\n", N, factor, cofactor);
                fprintf(stderr, "P+1 (B1=%.0f): %.3fs\n", pp1_bounds[i], elapsed_sec());
                goto done;
            }
        }
    }

    /* Stage 6: ECM (a few curves, ~20s) */
    {
        struct { double B1; int curves; } ecm_stages[] = {
            {1e4, 20}, {5e4, 10}, {1e5, 10}, {5e5, 5}, {0, 0}
        };
        for (int i = 0; ecm_stages[i].B1 > 0; i++) {
            if (elapsed_sec() > 30.0) break;
            if (ecm_method(N, factor, ecm_stages[i].B1, ecm_stages[i].curves)) {
                mpz_divexact(cofactor, N, factor);
                gmp_printf("%Zd = %Zd * %Zd\n", N, factor, cofactor);
                fprintf(stderr, "ECM (B1=%.0f): %.3fs\n", ecm_stages[i].B1, elapsed_sec());
                goto done;
            }
        }
    }

    /* Stage 7: No quick method worked - report failure */
    fprintf(stderr, "Oracle: no quick factor found after %.3fs, use SIQS\n", elapsed_sec());
    gmp_printf("FAIL %Zd\n", N);
    mpz_clears(N, factor, cofactor, NULL);
    return 1;

done:
    mpz_clears(N, factor, cofactor, NULL);
    return 0;
}
