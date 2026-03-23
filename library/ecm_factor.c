/*
 * ECM factoring using GMP-ECM library.
 * Usage: ./ecm_factor <number> [B1] [curves]
 * Default: B1=auto-scaled by digit count, curves=100
 * Seed is always 42.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#include <ecm.h>

static double default_b1(int digits) {
    /* Optimal B1 bounds - tuned for balanced semiprimes where factor has 'digits' digits */
    if (digits <= 15) return 50000;
    if (digits <= 20) return 500000;
    if (digits <= 25) return 3000000;
    if (digits <= 30) return 11000000;
    if (digits <= 35) return 43000000;
    if (digits <= 40) return 110000000;
    if (digits <= 45) return 260000000;
    if (digits <= 50) return 850000000;
    return 2900000000.0;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <number> [B1] [curves]\n", argv[0]);
        return 1;
    }

    mpz_t n, f;
    mpz_init(n);
    mpz_init(f);
    mpz_set_str(n, argv[1], 10);

    int digits = strlen(argv[1]);
    /* For balanced semiprimes, each factor has ~digits/2 digits */
    int factor_digits = digits / 2;

    double b1 = (argc >= 3) ? atof(argv[2]) : default_b1(factor_digits);
    int curves = (argc >= 4) ? atoi(argv[3]) : 500;

    ecm_params params;
    ecm_init(params);
    params->method = ECM_ECM;

    /* Seed = 42 always - use param type 0 (original sigma) with values < 2^32 */
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, 42);

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    int found = 0;
    for (int i = 0; i < curves; i++) {
        /* Generate sigma in range [6, 2^32-1] from seeded RNG */
        mpz_t max_sigma;
        mpz_init(max_sigma);
        mpz_set_ui(max_sigma, 0xFFFFFFF0UL);
        mpz_urandomm(params->sigma, rstate, max_sigma);
        mpz_add_ui(params->sigma, params->sigma, 6);
        mpz_clear(max_sigma);
        params->param = ECM_PARAM_SUYAMA;

        mpz_set(f, n);
        int ret = ecm_factor(f, n, b1, params);
        if (ret > 0) {
            /* Found a factor */
            if (mpz_cmp_ui(f, 1) > 0 && mpz_cmp(f, n) < 0) {
                clock_gettime(CLOCK_MONOTONIC, &end);
                double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
                mpz_t q;
                mpz_init(q);
                mpz_divexact(q, n, f);
                gmp_printf("%Zd %Zd\n", f, q);
                fprintf(stderr, "ECM found factor after %d curves, B1=%.0f, time=%.3fs\n", i+1, b1, elapsed);
                mpz_clear(q);
                found = 1;
                break;
            }
        }
    }

    if (!found) {
        clock_gettime(CLOCK_MONOTONIC, &end);
        double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
        fprintf(stderr, "ECM failed after %d curves, B1=%.0f, time=%.3fs\n", curves, b1, elapsed);
        printf("FAIL\n");
    }

    ecm_clear(params);
    gmp_randclear(rstate);
    mpz_clear(n);
    mpz_clear(f);
    return found ? 0 : 1;
}
