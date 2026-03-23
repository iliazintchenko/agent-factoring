/*
 * ECM factoring wrapper using GMP-ECM library.
 * Usage: ./ecm_factor <number>
 * Single-threaded, seed = 42.
 */

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <gmp.h>
#include <ecm.h>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <number>\n", argv[0]);
        return 1;
    }

    struct timespec start_time, end_time;
    clock_gettime(CLOCK_MONOTONIC, &start_time);

    mpz_t n, f;
    mpz_init_set_str(n, argv[1], 10);
    mpz_init(f);

    int digits = (int)mpz_sizeinbase(n, 10);
    int factor_digits = digits / 2;

    // ECM parameters - progressive B1 strategy
    // Try increasing B1 values
    struct Stage { double B1; int curves; };
    Stage stages[] = {
        {2000, 50},
        {11000, 150},
        {50000, 400},
        {250000, 1000},
        {1000000, 2500},
        {3000000, 5000},
        {11000000, 10000},
        {43000000, 20000},
        {110000000, 50000},
    };

    // Determine max stage based on factor size
    int max_stage;
    if (factor_digits <= 15) max_stage = 1;
    else if (factor_digits <= 20) max_stage = 2;
    else if (factor_digits <= 25) max_stage = 3;
    else if (factor_digits <= 30) max_stage = 4;
    else if (factor_digits <= 35) max_stage = 5;
    else if (factor_digits <= 40) max_stage = 6;
    else if (factor_digits <= 45) max_stage = 7;
    else if (factor_digits <= 50) max_stage = 8;
    else max_stage = 9;

    fprintf(stderr, "ECM: %d-digit number (expecting ~%d-digit factor), max_stage=%d\n",
            digits, factor_digits, max_stage);

    // Use deterministic seed
    gmp_randstate_t rstate;
    gmp_randinit_mt(rstate);
    gmp_randseed_ui(rstate, 42);

    int found = 0;
    int total_curves = 0;

    for (int s = 0; s < max_stage && !found; s++) {
        double B1 = stages[s].B1;
        int num_curves = stages[s].curves;

        for (int i = 0; i < num_curves && !found; i++) {
            ecm_params params;
            ecm_init(params);
            params->method = ECM_ECM;
            params->B1done = 1.0;

            // Generate sigma in valid range [6, 2^32-1]
            unsigned long sigma_val = gmp_urandomm_ui(rstate, (1UL << 32) - 7) + 6;
            mpz_set_ui(params->sigma, sigma_val);

            int ret = ecm_factor(f, n, B1, params);
            total_curves++;

            ecm_clear(params);

            if (ret > 0 && mpz_cmp_ui(f, 1) > 0 && mpz_cmp(f, n) < 0) {
                found = 1;
                break;
            }

            if (total_curves % 500 == 0) {
                fprintf(stderr, "  ECM: %d curves done (stage %d, B1=%.0f)\n", total_curves, s+1, B1);
            }
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &end_time);
    double elapsed = (end_time.tv_sec - start_time.tv_sec) +
                     (end_time.tv_nsec - start_time.tv_nsec) / 1e9;

    if (found) {
        gmp_printf("%Zd\n", f);
        fprintf(stderr, "ECM: factored in %.3f seconds (%d curves)\n", elapsed, total_curves);
    } else {
        fprintf(stderr, "ECM: FAILED after %.3f seconds (%d curves)\n", elapsed, total_curves);
        mpz_clear(n);
        mpz_clear(f);
        gmp_randclear(rstate);
        return 1;
    }

    mpz_clear(n);
    mpz_clear(f);
    gmp_randclear(rstate);
    return 0;
}
