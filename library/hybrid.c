/*
 * Hybrid ECM+MPQS factoring.
 *
 * Strategy:
 * 1. Run ECM with small B1 bounds (a few curves) for quick wins
 * 2. If ECM doesn't find a factor, fall back to MPQS
 *
 * For balanced semiprimes, ECM has a probability of ~1/sqrt(B1) of finding
 * a factor per curve. With a few curves at B1=1000-50000, we get a fast
 * filter that catches "lucky" numbers where p-1 or |E(F_p)| is smooth.
 *
 * Usage: ./hybrid <number>
 * Single-threaded, seed=42.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gmp.h>
#include <ecm.h>

int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <number>\n", argv[0]); return 1; }

    mpz_t n, f;
    mpz_init(n); mpz_init(f);
    if (mpz_set_str(n, argv[1], 10) != 0) return 1;

    struct timespec tstart;
    clock_gettime(CLOCK_MONOTONIC, &tstart);

    /* Phase 1: Quick ECM with small B1 */
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, 42);

    size_t digits = mpz_sizeinbase(n, 10);

    /* ECM parameters: few curves at increasing B1 */
    double B1_values[] = {1000, 5000, 25000, 100000};
    int curves_per_B1[] = {10, 10, 5, 5};
    int num_levels = 4;

    /* For small numbers, ECM alone is fast enough */
    if (digits <= 35) {
        num_levels = 4;
        curves_per_B1[0] = 20;
        curves_per_B1[1] = 20;
        curves_per_B1[2] = 15;
        curves_per_B1[3] = 10;
    } else if (digits <= 50) {
        /* Spend less time on ECM for larger numbers */
        num_levels = 3;
        curves_per_B1[0] = 5;
        curves_per_B1[1] = 5;
        curves_per_B1[2] = 3;
    } else {
        num_levels = 2;
        curves_per_B1[0] = 3;
        curves_per_B1[1] = 3;
    }

    int ecm_found = 0;
    for (int level = 0; level < num_levels && !ecm_found; level++) {
        double B1 = B1_values[level];
        for (int c = 0; c < curves_per_B1[level] && !ecm_found; c++) {
            ecm_params params;
            ecm_init(params);
            unsigned long sigma = gmp_urandomm_ui(rstate, 1000000000UL) + 6;
            params->B1done = 1.0;
            params->method = ECM_ECM;
            params->param = ECM_PARAM_SUYAMA;
            params->sigma_is_A = 0;
            mpz_set_ui(params->sigma, sigma);

            mpz_set_ui(f, 0);
            int ret = ecm_factor(f, n, B1, params);
            ecm_clear(params);

            if (ret > 0 && mpz_cmp_ui(f, 1) > 0 && mpz_cmp(f, n) < 0) {
                struct timespec tnow;
                clock_gettime(CLOCK_MONOTONIC, &tnow);
                double elapsed = (tnow.tv_sec - tstart.tv_sec) + (tnow.tv_nsec - tstart.tv_nsec) / 1e9;
                gmp_printf("FACTOR: %Zd\n", f);
                mpz_t cof; mpz_init(cof);
                mpz_divexact(cof, n, f);
                gmp_printf("COFACTOR: %Zd\n", cof);
                fprintf(stderr, "Hybrid: ECM found factor at B1=%.0f, curve %d, time=%.3fs\n", B1, c, elapsed);
                mpz_clear(cof);
                ecm_found = 1;
            }
        }
    }

    if (ecm_found) {
        gmp_randclear(rstate);
        mpz_clear(n); mpz_clear(f);
        return 0;
    }

    struct timespec tnow;
    clock_gettime(CLOCK_MONOTONIC, &tnow);
    double ecm_time = (tnow.tv_sec - tstart.tv_sec) + (tnow.tv_nsec - tstart.tv_nsec) / 1e9;
    fprintf(stderr, "Hybrid: ECM phase done (%.1fs), falling back to MPQS\n", ecm_time);

    /* Phase 2: MPQS via exec */
    gmp_randclear(rstate);
    mpz_clear(f);

    char cmd[2048];
    snprintf(cmd, sizeof(cmd), "timeout %d %s/mpqs %s",
             (int)(295.0 - ecm_time),
             "library", /* relative path */
             argv[1]);

    /* Get the directory of our binary */
    char *dir = strdup(argv[0]);
    char *slash = strrchr(dir, '/');
    if (slash) {
        *slash = '\0';
        snprintf(cmd, sizeof(cmd), "timeout %d %s/mpqs %s",
                 (int)(295.0 - ecm_time), dir, argv[1]);
    } else {
        snprintf(cmd, sizeof(cmd), "timeout %d ./mpqs %s",
                 (int)(295.0 - ecm_time), argv[1]);
    }
    free(dir);

    FILE *fp = popen(cmd, "r");
    if (!fp) {
        fprintf(stderr, "Failed to run MPQS\n");
        printf("FAIL\n");
        mpz_clear(n);
        return 1;
    }

    char buf[4096];
    char last_line[4096] = "";
    while (fgets(buf, sizeof(buf), fp)) {
        strncpy(last_line, buf, sizeof(last_line) - 1);
        /* Pass through FACTOR: lines */
        if (strncmp(buf, "FACTOR:", 7) == 0) {
            printf("%s", buf);
        }
    }
    int status = pclose(fp);

    /* Check if MPQS output two numbers on last line */
    if (strlen(last_line) > 0) {
        char *space = strchr(last_line, ' ');
        if (space) {
            /* Two numbers on last line = factor + cofactor */
            *space = '\0';
            clock_gettime(CLOCK_MONOTONIC, &tnow);
            double total = (tnow.tv_sec - tstart.tv_sec) + (tnow.tv_nsec - tstart.tv_nsec) / 1e9;
            printf("FACTOR: %s\n", last_line);
            printf("COFACTOR: %s\n", space + 1);
            fprintf(stderr, "Hybrid: MPQS found factor, total time=%.3fs\n", total);
        }
    }

    if (status != 0 && strlen(last_line) == 0) {
        printf("FAIL\n");
    }

    mpz_clear(n);
    return 0;
}
