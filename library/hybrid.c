/*
 * Hybrid Factoring Engine
 *
 * Combines trial division, Pollard's rho, P-1, P+1, and ECM.
 * Uses GMP-ECM library. All randomness seeded with 42.
 *
 * Output: "factor1 factor2" (smaller first)
 */

#include <gmp.h>
#include <ecm.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static struct timespec g_start;

static double elapsed_sec(void) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec - g_start.tv_sec) +
           (now.tv_nsec - g_start.tv_nsec) / 1e9;
}

static int trial_divide(mpz_t factor, const mpz_t N, unsigned long limit)
{
    if (mpz_divisible_ui_p(N, 2)) { mpz_set_ui(factor, 2); return 1; }
    if (mpz_divisible_ui_p(N, 3)) { mpz_set_ui(factor, 3); return 1; }
    for (unsigned long d = 5; d <= limit; d += (d % 6 == 5) ? 2 : 4) {
        if (mpz_divisible_ui_p(N, d)) {
            mpz_set_ui(factor, d);
            return 1;
        }
    }
    return 0;
}

static int pollard_rho(mpz_t factor, const mpz_t N, unsigned long max_iters,
                       gmp_randstate_t rng)
{
    mpz_t x, y, d, acc, temp, c;
    mpz_inits(x, y, d, acc, temp, c, NULL);
    int found = 0;

    for (int attempt = 0; attempt < 20 && !found; attempt++) {
        mpz_urandomm(c, rng, N);
        mpz_add_ui(c, c, 1);
        mpz_urandomm(x, rng, N);
        mpz_set(y, x);
        mpz_set_ui(acc, 1);
        unsigned long batch = 0;

        for (unsigned long i = 0; i < max_iters && !found; i++) {
            mpz_mul(x, x, x); mpz_add(x, x, c); mpz_mod(x, x, N);
            mpz_mul(y, y, y); mpz_add(y, y, c); mpz_mod(y, y, N);
            mpz_mul(y, y, y); mpz_add(y, y, c); mpz_mod(y, y, N);

            mpz_sub(temp, x, y);
            if (mpz_sgn(temp) == 0) break;
            mpz_mul(acc, acc, temp);
            mpz_mod(acc, acc, N);
            batch++;

            if (batch >= 256) {
                mpz_gcd(d, acc, N);
                if (mpz_cmp_ui(d, 1) > 0 && mpz_cmp(d, N) < 0) {
                    mpz_set(factor, d);
                    found = 1;
                } else if (mpz_cmp(d, N) == 0) {
                    break; /* accumulated to N, need finer steps */
                }
                mpz_set_ui(acc, 1);
                batch = 0;
            }
        }
        if (!found && batch > 0) {
            mpz_gcd(d, acc, N);
            if (mpz_cmp_ui(d, 1) > 0 && mpz_cmp(d, N) < 0) {
                mpz_set(factor, d);
                found = 1;
            }
        }
    }

    mpz_clears(x, y, d, acc, temp, c, NULL);
    return found;
}

int main(int argc, char *argv[])
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    clock_gettime(CLOCK_MONOTONIC, &g_start);

    mpz_t N, factor, cofactor;
    mpz_inits(N, factor, cofactor, NULL);
    mpz_set_str(N, argv[1], 10);

    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, 42);

    int digits = strlen(argv[1]);
    int found = 0;

    /* Stage 0: Trial division */
    found = trial_divide(factor, N, 1000000);

    /* Stage 1: Pollard's rho — 500K iters (~1 sec for 30-digit) */
    if (!found) {
        found = pollard_rho(factor, N, 500000, rng);
        if (found)
            fprintf(stderr, "rho: %.3fs\n", elapsed_sec());
    }

    /* Stage 2: ECM via GMP-ECM */
    if (!found) {
        /* For balanced semiprimes: factor ~digits/2 digits */
        int fd = (digits + 1) / 2;

        /* B1 levels and expected curves from GMP-ECM docs */
        struct { double b1; int curves; } levels[] = {
            {2e3,    25},    /* up to 15 digits */
            {1.1e4,  90},    /* up to 20 */
            {5e4,    300},   /* up to 25 */
            {2.5e5,  700},   /* up to 30 */
            {1e6,    1800},  /* up to 35 */
            {3e6,    5100},  /* up to 40 */
            {1.1e7,  10600}, /* up to 45 */
            {4.3e7,  19300}, /* up to 50 */
            {1.1e8,  49000}, /* up to 55 */
        };
        int nlevels = sizeof(levels)/sizeof(levels[0]);

        ecm_params ecm_par;
        ecm_init(ecm_par);
        ecm_par->verbose = 0;

        /* First try P-1 at each level */
        for (int lev = 0; lev < nlevels && !found; lev++) {
            double b1 = levels[lev].b1;

            /* P-1 */
            ecm_par->method = ECM_PM1;
            mpz_set_si(ecm_par->B2, -1);
            ecm_par->B1done = 0;
            int ret = ecm_factor(factor, N, b1, ecm_par);
            if (ret > 0 && mpz_cmp_ui(factor, 1) > 0 &&
                mpz_cmp(factor, N) < 0) {
                found = 1;
                fprintf(stderr, "P-1 B1=%.0f: %.3fs\n", b1, elapsed_sec());
                break;
            }

            /* P+1 */
            ecm_par->method = ECM_PP1;
            mpz_set_si(ecm_par->B2, -1);
            ecm_par->B1done = 0;
            ret = ecm_factor(factor, N, b1, ecm_par);
            if (ret > 0 && mpz_cmp_ui(factor, 1) > 0 &&
                mpz_cmp(factor, N) < 0) {
                found = 1;
                fprintf(stderr, "P+1 B1=%.0f: %.3fs\n", b1, elapsed_sec());
                break;
            }

            if (elapsed_sec() > 280.0) break;
        }

        /* ECM: run curves at each level progressively */
        if (!found) {
            ecm_par->method = ECM_ECM;
            unsigned long sigma_counter = 7; /* Starting sigma (>= 6 required) */

            for (int lev = 0; lev < nlevels && !found; lev++) {
                double b1 = levels[lev].b1;
                int nc = levels[lev].curves;

                for (int c = 0; c < nc && !found; c++) {
                    ecm_par->param = ECM_PARAM_SUYAMA;
                    mpz_set_ui(ecm_par->sigma, sigma_counter++);
                    mpz_set_si(ecm_par->B2, -1);
                    ecm_par->B1done = 0;

                    int ret = ecm_factor(factor, N, b1, ecm_par);
                    if (ret > 0 && mpz_cmp_ui(factor, 1) > 0 &&
                        mpz_cmp(factor, N) < 0) {
                        found = 1;
                        fprintf(stderr, "ECM B1=%.0f sigma=%lu curve#%d: %.3fs\n",
                                b1, sigma_counter-1, c, elapsed_sec());
                    }

                    if (elapsed_sec() > 280.0) break;
                }
                if (elapsed_sec() > 280.0) break;
            }
        }

        ecm_clear(ecm_par);
    }

    if (found) {
        mpz_divexact(cofactor, N, factor);
        if (mpz_cmp(factor, cofactor) > 0) mpz_swap(factor, cofactor);
        gmp_printf("%Zd %Zd\n", factor, cofactor);
        fprintf(stderr, "Total: %.3fs (%d digits)\n", elapsed_sec(), digits);
    } else {
        fprintf(stderr, "FAILED after %.3fs (%d digits)\n", elapsed_sec(), digits);
        return 1;
    }

    mpz_clears(N, factor, cofactor, NULL);
    gmp_randclear(rng);
    return 0;
}
