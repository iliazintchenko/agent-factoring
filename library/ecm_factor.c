/*
 * ECM factoring using GMP-ECM library.
 * Usage: ./ecm_factor <number> [B1_start] [B1_end]
 * Single-threaded, seed=42.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#include <ecm.h>

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <number> [B1_start] [B1_end]\n", argv[0]);
        return 1;
    }

    mpz_t n, f;
    mpz_init(n);
    mpz_init(f);

    if (mpz_set_str(n, argv[1], 10) != 0) {
        fprintf(stderr, "Invalid number: %s\n", argv[1]);
        return 1;
    }

    if (mpz_cmp_ui(n, 1) <= 0) { printf("TRIVIAL\n"); return 0; }
    if (mpz_even_p(n)) { printf("FACTOR: 2\n"); return 0; }

    double B1_start = 1e4;
    double B1_end = 1e11;
    if (argc >= 3) B1_start = atof(argv[2]);
    if (argc >= 4) B1_end = atof(argv[3]);

    struct timespec start, now;
    clock_gettime(CLOCK_MONOTONIC, &start);

    double B1 = B1_start;
    int curve = 0;
    int curves_per_B1 = 50;

    /* Use a deterministic sequence of A values for Montgomery curves */
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, 42);

    while (B1 <= B1_end) {
        for (int i = 0; i < curves_per_B1; i++) {
            ecm_params params;
            ecm_init(params);

            /* Use Suyama parametrization with deterministic sigma from seed=42 RNG */
            unsigned long sigma_val = gmp_urandomm_ui(rstate, 1000000000UL) + 6;
            params->B1done = 1.0;
            params->method = ECM_ECM;
            params->param = ECM_PARAM_SUYAMA;
            params->sigma_is_A = 0;
            mpz_set_ui(params->sigma, sigma_val);

            mpz_set_ui(f, 0);
            int ret = ecm_factor(f, n, B1, params);
            ecm_clear(params);

            if (ret > 0) {
                clock_gettime(CLOCK_MONOTONIC, &now);
                double elapsed = (now.tv_sec - start.tv_sec) + (now.tv_nsec - start.tv_nsec) / 1e9;
                gmp_printf("FACTOR: %Zd\n", f);
                fprintf(stderr, "ECM found factor with B1=%.0f, curve=%d, time=%.3fs\n", B1, curve, elapsed);

                mpz_t q, r;
                mpz_init(q); mpz_init(r);
                mpz_tdiv_qr(q, r, n, f);
                if (mpz_cmp_ui(r, 0) == 0 && mpz_cmp_ui(f, 1) > 0 && mpz_cmp(f, n) < 0) {
                    gmp_printf("COFACTOR: %Zd\n", q);
                    fprintf(stderr, "VERIFIED\n");
                } else {
                    fprintf(stderr, "VERIFICATION FAILED\n");
                }
                mpz_clear(q); mpz_clear(r);
                gmp_randclear(rstate);
                mpz_clear(n); mpz_clear(f);
                return 0;
            }
            curve++;
        }

        clock_gettime(CLOCK_MONOTONIC, &now);
        double elapsed = (now.tv_sec - start.tv_sec) + (now.tv_nsec - start.tv_nsec) / 1e9;
        fprintf(stderr, "B1=%.0f done (%d curves, %.1fs elapsed)\n", B1, curve, elapsed);
        B1 *= 3.0;
    }

    fprintf(stderr, "ECM failed to find factor with B1 up to %.0f (%d curves)\n", B1_end, curve);
    printf("FAIL\n");
    gmp_randclear(rstate);
    mpz_clear(n); mpz_clear(f);
    return 1;
}
