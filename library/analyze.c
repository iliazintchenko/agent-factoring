/*
 * Semiprime Factoring Difficulty Analyzer
 * Computes theoretical optimal SIQS parameters and estimated time.
 *
 * For balanced semiprime N = p*q:
 * - Uses Dickman's rho function to estimate smoothness probability
 * - Computes optimal factor base size and sieve interval
 * - Estimates expected number of polynomials and total time
 *
 * Compile: gcc -O2 -march=native library/analyze.c -lgmp -lm -o analyze
 * Usage: ./analyze <N>
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

/*
 * Dickman's rho function ρ(u):
 * For u >= 1, ρ(u) gives the probability that a random number N
 * has all prime factors <= N^(1/u).
 *
 * Uses the approximation ρ(u) ≈ u^(-u) for u > 1.
 * More precise: ρ(u) = 1 for u <= 1.
 * For u > 1, use the recurrence or a table of known values.
 */
static double dickman_rho(double u) {
    if (u <= 0) return 1.0;
    if (u <= 1.0) return 1.0;
    if (u <= 2.0) return 1.0 - log(u);

    /* For u > 2, use the asymptotic approximation:
     * ρ(u) ≈ exp(-u * (ln(u) + ln(ln(u)) - 1 + (ln(ln(u)) - 2) / ln(u)))
     * This is the de Bruijn approximation.
     */
    double lnu = log(u);
    double lnlnu = log(lnu);
    double exponent = -u * (lnu + lnlnu - 1.0 + (lnlnu - 2.0) / lnu);
    return exp(exponent);
}

/*
 * L(N) = exp(sqrt(ln(N) * ln(ln(N))))
 * The sub-exponential complexity function for QS/NFS.
 */
static double L_function(double ln_N) {
    return exp(sqrt(ln_N * log(ln_N)));
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    mpz_t N;
    mpz_init_set_str(N, argv[1], 10);
    int digits = (int)mpz_sizeinbase(N, 10);
    int bits = (int)mpz_sizeinbase(N, 2);
    double ln_N = bits * log(2.0);
    double log2_N = bits;

    printf("=== Semiprime Factoring Analysis ===\n");
    printf("N: %d digits, %d bits\n\n", digits, bits);

    /* L(N) value */
    double L_N = L_function(ln_N);
    double ln_L = sqrt(ln_N * log(ln_N));
    printf("L(N) = exp(%.1f) ≈ 10^%.1f\n", ln_L, ln_L / log(10));

    /* Optimal SIQS parameters (theoretical) */
    /* B (smoothness bound) ≈ L(N)^(1/sqrt(2)) */
    double opt_ln_B = ln_L / sqrt(2.0);
    double opt_B = exp(opt_ln_B);
    int opt_fb_size = (int)(opt_B / (log(opt_B) - 1)); /* π(B) ≈ B/ln(B) */

    /* M (sieve half-interval) ≈ L(N)^(sqrt(2)) / sqrt(2N) */
    /* But in practice, M is chosen so that max|Q(x)| ≈ sqrt(2N) * M is optimal */
    double opt_M_factor = exp(ln_L * sqrt(2.0));
    double opt_M = opt_M_factor / exp(ln_N / 2.0);
    if (opt_M < 32768) opt_M = 32768;

    printf("--- Theoretical Optimal SIQS Parameters ---\n");
    printf("Smoothness bound B ≈ %.0f (log2 = %.1f)\n", opt_B, log2(opt_B));
    printf("Factor base size π(B) ≈ %d primes\n", opt_fb_size);
    printf("Sieve half-interval M ≈ %.0f\n", opt_M);

    /* Smoothness probability */
    /* For Q(x) ≈ sqrt(N) * M, the u-value is: */
    double log_Q = log2_N / 2.0 + log2(opt_M);
    double u = log_Q * log(2.0) / opt_ln_B;
    double smooth_prob = dickman_rho(u);

    printf("\n--- Smoothness Analysis ---\n");
    printf("Max |Q(x)| ≈ 2^%.1f\n", log_Q);
    printf("u = ln(Q)/ln(B) = %.2f\n", u);
    printf("ρ(u) ≈ %.3e (probability of Q being B-smooth)\n", smooth_prob);

    /* With large primes (single LP), effective probability is ~4x higher */
    double lp_factor = 4.0; /* empirical for single LP with LP bound = 30*B */
    double effective_prob = smooth_prob * lp_factor;
    printf("With large primes: effective prob ≈ %.3e (×%.0f)\n", effective_prob, lp_factor);

    /* Relations needed ≈ FB_size + extra */
    int rels_needed = opt_fb_size + (int)(opt_fb_size * 0.1);
    printf("\n--- Time Estimation ---\n");
    printf("Relations needed ≈ %d\n", rels_needed);

    /* Sieve positions per polynomial = 2*M */
    double polys_per_rel = 1.0 / (effective_prob * 2 * opt_M);
    double total_polys = rels_needed * polys_per_rel;
    printf("Expected polys ≈ %.0f\n", total_polys);

    /* Each polynomial costs: FB_size root updates + 2*M sieve operations + scan */
    /* On modern hardware: sieve ≈ 2*M bytes touched at ~10GB/s */
    double sieve_time_per_poly = 2 * opt_M / 10e9; /* seconds */
    double root_time_per_poly = (double)opt_fb_size * 50e-9; /* ~50ns per root update */
    double poly_time = sieve_time_per_poly + root_time_per_poly;
    double est_total_time = total_polys * poly_time;

    printf("Time per poly ≈ %.3e s (sieve=%.3e, roots=%.3e)\n",
           poly_time, sieve_time_per_poly, root_time_per_poly);
    printf("Estimated total sieve time ≈ %.1f s\n", est_total_time);

    /* Block Lanczos time (empirical: ~2-8% of sieve time for 80-89d) */
    double la_time = est_total_time * 0.05;
    printf("Estimated LA time ≈ %.1f s\n", la_time);
    printf("Estimated TOTAL time ≈ %.1f s\n", est_total_time + la_time);

    /* YAFU empirical comparison */
    printf("\n--- YAFU Empirical Parameters ---\n");
    if (bits <= 100) printf("smallmpqs: ~0.01s\n");
    else if (bits <= 170) printf("siqs (default params): <1s\n");
    else if (bits <= 200) printf("siqs (default params): 1-5s\n");
    else if (bits <= 230) printf("siqs (-siqsNB 10): 5-20s\n");
    else if (bits <= 250) printf("siqs (-siqsNB 12): 20-55s\n");
    else if (bits <= 270) printf("siqs (-siqsNB 14): 55-100s\n");
    else if (bits <= 285) printf("siqs (-siqsNB 18 -siqsB 70-90K): 100-200s\n");
    else if (bits <= 295) printf("siqs (-siqsNB 18 -siqsB 80-100K): 200-295s\n");
    else printf("INFEASIBLE: >300s single-core for any SIQS configuration\n");

    /* For balanced semiprimes, check if factors could be close */
    printf("\n--- Balanced Semiprime Analysis ---\n");
    printf("Expected factor size: ~%d digits (%d bits)\n", digits / 2, bits / 2);
    printf("Factor range: [10^%d, 10^%d]\n", (digits - 1) / 2, (digits + 1) / 2);

    /* Fermat's method check: would need |p-q| < N^(1/4) for fast factoring */
    printf("Fermat boundary: |p-q| < N^(1/4) ≈ 2^%d ≈ 10^%d digits\n",
           bits / 4, (int)(bits * 0.075));
    printf("For random balanced semiprime: |p-q| ≈ √N ≈ 2^%d — Fermat useless\n", bits / 2);

    /* ECM analysis */
    printf("\n--- ECM Analysis (for balanced semiprimes) ---\n");
    double p_bits = bits / 2.0;
    double ecm_B1_needed = exp(sqrt(2.0 * p_bits * log(2.0) * log(p_bits * log(2.0))));
    printf("Factor p has ~%d bits\n", (int)p_bits);
    printf("ECM optimal B1 ≈ %.0e (too large for practical single-curve)\n", ecm_B1_needed);
    printf("ECM is NOT viable for balanced semiprimes (factors too large)\n");

    /* NFS crossover analysis */
    printf("\n--- NFS Crossover ---\n");
    double nfs_complexity = exp(1.923 * pow(ln_N, 1.0/3.0) * pow(log(ln_N), 2.0/3.0));
    double siqs_complexity = L_N;
    printf("SIQS complexity L(N) ≈ exp(%.1f)\n", ln_L);
    printf("GNFS complexity ≈ exp(%.1f)\n",
           1.923 * pow(ln_N, 1.0/3.0) * pow(log(ln_N), 2.0/3.0));
    if (nfs_complexity < siqs_complexity)
        printf("NFS is theoretically faster at %d digits\n", digits);
    else
        printf("SIQS is faster at %d digits (NFS crossover ~100-110d)\n", digits);

    mpz_clear(N);
    return 0;
}
