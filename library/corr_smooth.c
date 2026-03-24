/*
 * corr_smooth.c — Experiment: Smoothness correlations across polynomial pairs.
 *
 * HYPOTHESIS: For two QS-style polynomials P1(x) = (x+m)^2 - N and
 * P2(x) = (x+m)^2 - kN (for small k), the smoothness of P1(x) and P2(x)
 * might be CORRELATED because they share the factor (x+m)^2 and differ
 * only in the subtracted term.
 *
 * If Pr[P1 smooth AND P2 smooth] >> Pr[P1 smooth] * Pr[P2 smooth],
 * then simultaneously sieving both polynomials could yield more relations
 * per sieve step than independent sieving.
 *
 * More interesting pairing: P1(x) = (x+m)^2 - N and
 * P3(x) = (x+m+1)^2 - N = P1(x) + 2(x+m) + 1.
 * P1 and P3 differ by 2(x+m)+1, which is odd and ~2√N.
 * If P1(x) = a*b and P3(x) = c*d where a,c share factors,
 * we get additional relations.
 *
 * EXPERIMENT: For each x, compute P1(x) and P1(x+d) for several d.
 * Measure the conditional probability Pr[P1(x+d) smooth | P1(x) smooth].
 * Compare to the unconditional Pr[P1(x+d) smooth].
 * If conditional > unconditional, correlation exists.
 */

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

static int is_smooth(mpz_t n, unsigned long B) {
    mpz_t r;
    mpz_init(r);
    mpz_abs(r, n);
    if (mpz_cmp_ui(r, 0) == 0) { mpz_clear(r); return 1; }

    /* Trial divide by primes up to B */
    for (unsigned long p = 2; p <= B; p++) {
        /* Quick primality check not needed — composite factors will be caught */
        while (mpz_divisible_ui_p(r, p))
            mpz_divexact_ui(r, r, p);
    }

    int s = (mpz_cmp_ui(r, 1) == 0);
    mpz_clear(r);
    return s;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N> [B] [range]\n", argv[0]);
        return 1;
    }

    mpz_t N;
    mpz_init(N);
    mpz_set_str(N, argv[1], 10);

    unsigned long B = (argc >= 3) ? atol(argv[2]) : 5000;
    long range = (argc >= 4) ? atol(argv[3]) : 200000;

    size_t ndig = mpz_sizeinbase(N, 10);
    fprintf(stderr, "N=%zu digits, B=%lu, range=+-%ld\n", ndig, B, range);

    mpz_t sqN, val, res1, res2, tmp;
    mpz_init(sqN); mpz_init(val); mpz_init(res1); mpz_init(res2); mpz_init(tmp);
    mpz_sqrt(sqN, N);

    /* Test several shift values d */
    int shifts[] = {1, 2, 3, 5, 7, 10, 20, 50, 100};
    int n_shifts = sizeof(shifts) / sizeof(shifts[0]);

    /* For each x in range, compute P1(x) = (x+m)^2 - N
     * and check smoothness */
    int *smooth1 = calloc(2 * range + 1, sizeof(int));

    fprintf(stderr, "Computing smoothness for %ld values...\n", 2 * range + 1);

    int total_smooth1 = 0;
    for (long x = -range; x <= range; x++) {
        if (x >= 0) mpz_add_ui(val, sqN, x);
        else mpz_sub_ui(val, sqN, -x);
        mpz_mul(res1, val, val);
        mpz_sub(res1, res1, N);
        smooth1[x + range] = is_smooth(res1, B);
        if (smooth1[x + range]) total_smooth1++;
    }

    double p_smooth = (double)total_smooth1 / (2 * range + 1);
    fprintf(stderr, "Pr[smooth] = %d / %ld = %.6f\n", total_smooth1, 2*range+1, p_smooth);

    /* For each shift d, measure conditional probability */
    fprintf(stderr, "\nShift | Pr[both] | Pr[s1]*Pr[s2] | Ratio | Interpretation\n");
    fprintf(stderr, "------+----------+---------------+-------+---------------\n");

    for (int si = 0; si < n_shifts; si++) {
        int d = shifts[si];
        int both_smooth = 0;
        int count = 0;

        for (long x = -range; x <= range - d; x++) {
            int s1 = smooth1[x + range];
            int s2 = smooth1[x + range + d];
            if (s1 && s2) both_smooth++;
            count++;
        }

        double p_both = (double)both_smooth / count;
        double p_indep = p_smooth * p_smooth;
        double ratio = (p_indep > 0) ? p_both / p_indep : 0;

        fprintf(stderr, "  %3d | %.6f | %.6f       | %.3f | %s\n",
                d, p_both, p_indep, ratio,
                ratio > 1.5 ? "STRONG correlation!" :
                ratio > 1.1 ? "Weak correlation" :
                ratio > 0.9 ? "No correlation" :
                "Anti-correlation");
    }

    /* Test cross-multiplier correlation:
     * P1(x) = (x+m)^2 - N vs P2(x) = (x+m')^2 - 2N where m' = floor(sqrt(2N)) */
    fprintf(stderr, "\nCross-multiplier test (k=2):\n");
    {
        mpz_t sqkN, valk, resk;
        mpz_init(sqkN); mpz_init(valk); mpz_init(resk);
        mpz_mul_ui(tmp, N, 2);
        mpz_sqrt(sqkN, tmp);

        int both = 0, s1_count = 0, s2_count = 0;
        long cnt = 0;
        for (long x = -range; x <= range; x++) {
            if (x >= 0) { mpz_add_ui(val, sqN, x); mpz_add_ui(valk, sqkN, x); }
            else { mpz_sub_ui(val, sqN, -x); mpz_sub_ui(valk, sqkN, -x); }

            mpz_mul(res1, val, val); mpz_sub(res1, res1, N);
            mpz_mul(resk, valk, valk); mpz_sub(resk, resk, tmp);

            int s1 = is_smooth(res1, B);
            int s2 = is_smooth(resk, B);
            if (s1) s1_count++;
            if (s2) s2_count++;
            if (s1 && s2) both++;
            cnt++;
        }

        double p1 = (double)s1_count / cnt;
        double p2 = (double)s2_count / cnt;
        double pb = (double)both / cnt;
        double pi = p1 * p2;
        fprintf(stderr, "Pr[s1]=%.6f, Pr[s2]=%.6f, Pr[both]=%.6f, independent=%.6f, ratio=%.3f\n",
                p1, p2, pb, pi, pi > 0 ? pb / pi : 0);

        mpz_clear(sqkN); mpz_clear(valk); mpz_clear(resk);
    }

    /* Test GCD-based correlation:
     * For pairs (P1(x), P1(x+d)), compute gcd(P1(x), P1(x+d)).
     * If gcd > 1 frequently, shared factors exist. */
    fprintf(stderr, "\nGCD analysis for d=1:\n");
    {
        int gcd_nontrivial = 0;
        long cnt = 0;
        for (long x = -range; x <= range - 1; x++) {
            if (x >= 0) mpz_add_ui(val, sqN, x);
            else mpz_sub_ui(val, sqN, -x);
            mpz_mul(res1, val, val); mpz_sub(res1, res1, N);
            mpz_abs(res1, res1);

            mpz_add_ui(val, val, 1);
            mpz_mul(res2, val, val); mpz_sub(res2, res2, N);
            mpz_abs(res2, res2);

            mpz_gcd(tmp, res1, res2);
            if (mpz_cmp_ui(tmp, 1) > 0) gcd_nontrivial++;
            cnt++;
        }
        fprintf(stderr, "gcd(P(x), P(x+1)) > 1: %d / %ld = %.4f\n",
                gcd_nontrivial, cnt, (double)gcd_nontrivial / cnt);
    }

    free(smooth1);
    mpz_clear(N); mpz_clear(sqN); mpz_clear(val);
    mpz_clear(res1); mpz_clear(res2); mpz_clear(tmp);

    return 0;
}
