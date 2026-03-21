/*
 * ECM worker: runs ECM with a specific B1 and curve count.
 * Usage: ./ecm_worker <N> <B1> <curves>
 * On success: prints factor to stdout, exits 0
 * On failure: exits 1
 */
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <ecm.h>

int main(int argc, char *argv[]) {
    if (argc < 4) { fprintf(stderr, "Usage: %s <N> <B1> <curves>\n", argv[0]); return 1; }

    mpz_t n, f;
    mpz_inits(n, f, NULL);
    mpz_set_str(n, argv[1], 10);
    double b1 = atof(argv[2]);
    int curves = atoi(argv[3]);

    ecm_params params;
    ecm_init(params);
    params->verbose = 0;

    for (int i = 0; i < curves; i++) {
        mpz_set(f, n);
        int ret = ecm_factor(f, n, b1, params);
        if (ret > 0 && mpz_cmp_ui(f, 1) > 0 && mpz_cmp(f, n) < 0) {
            mpz_t other; mpz_init(other);
            mpz_divexact(other, n, f);
            if (mpz_cmp(f, other) > 0) mpz_swap(f, other);
            gmp_printf("%Zd\n", f);
            mpz_clear(other);
            ecm_clear(params);
            mpz_clears(n, f, NULL);
            return 0;
        }
        mpz_set_ui(params->sigma, 0);
    }
    ecm_clear(params);
    mpz_clears(n, f, NULL);
    return 1;
}
