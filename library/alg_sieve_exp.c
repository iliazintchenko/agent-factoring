/*
 * alg_sieve_exp.c — Measure doubly-smooth yield in NFS vs QS search spaces
 *
 * For a given N, compare:
 * 1. QS: f(x) = (x+m_qs)^2 - N, x in [-A, A]
 * 2. NFS degree-3: (a-b*m, f_hom(a,b)), (a,b) in grid
 * 3. NFS degree-4: similar
 *
 * For each approach, measure:
 * - Value sizes (bits) of the evaluated norms
 * - Smoothness probability (trial divide + count)
 * - Relations per unit of search effort
 *
 * This determines whether NFS-style evaluation gives more relations
 * per candidate than QS, even without lattice sieving.
 *
 * Usage: ./alg_sieve_exp <N>
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

static int gcd_int(int a, int b) {
    a = a < 0 ? -a : a; b = b < 0 ? -b : b;
    while (b) { int t = a % b; a = b; b = t; } return a;
}

int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }

    mpz_t N;
    mpz_init(N);
    mpz_set_str(N, argv[1], 10);
    int digits = strlen(argv[1]);

    /* Smoothness bound */
    double lnN = digits * 2.302585;
    double lnlnN = log(lnN);
    double B_val = exp(0.5 * sqrt(lnN * lnlnN));
    unsigned long B = (unsigned long)B_val;
    if (B > 100000) B = 100000; /* cap for speed */
    printf("# %d-digit N, smoothness bound B = %lu\n", digits, B);

    mpz_t val, cof, tmp;
    mpz_init(val); mpz_init(cof); mpz_init(tmp);

    /* Sieve of small primes for smoothness testing */
    int n_primes = 0;
    unsigned long primes[10000];
    for (unsigned long p = 2; p <= B && n_primes < 10000; p++) {
        int ip = 1;
        for (unsigned long d = 2; d*d <= p; d++) if (p%d==0) { ip=0; break; }
        if (ip) primes[n_primes++] = p;
    }
    printf("# Factor base: %d primes up to %lu\n", n_primes, primes[n_primes-1]);

    /* QS: f(x) = (x + ceil(sqrt(N)))^2 - N */
    mpz_t sqrtN;
    mpz_init(sqrtN);
    mpz_sqrt(sqrtN, N);
    mpz_add_ui(sqrtN, sqrtN, 1);

    int qs_total = 0, qs_smooth = 0, qs_1lp = 0;
    long qs_val_bits_sum = 0;
    int range = 100000;

    for (int x = 1; x <= range; x++) {
        mpz_set(val, sqrtN);
        mpz_add_ui(val, val, x);
        mpz_mul(val, val, val);
        mpz_sub(val, val, N);
        mpz_abs(cof, val);

        qs_val_bits_sum += mpz_sizeinbase(cof, 2);

        for (int i = 0; i < n_primes; i++) {
            while (mpz_divisible_ui_p(cof, primes[i]))
                mpz_divexact_ui(cof, cof, primes[i]);
        }
        qs_total++;
        if (mpz_cmp_ui(cof, 1) == 0) qs_smooth++;
        else if (mpz_sizeinbase(cof, 2) <= 30) qs_1lp++;
    }

    printf("\nQS: %d candidates in [1, %d]\n", qs_total, range);
    printf("  avg value size: %.1f bits\n", (double)qs_val_bits_sum / qs_total);
    printf("  fully smooth: %d (%.4f%%)\n", qs_smooth, 100.0*qs_smooth/qs_total);
    printf("  1 large prime: %d (%.4f%%)\n", qs_1lp, 100.0*qs_1lp/qs_total);

    /* NFS degree 3: base-m decomposition */
    for (int degree = 3; degree <= 5; degree += 1) {
        mpz_t m, coeff[8];
        mpz_init(m);
        for (int i = 0; i <= degree; i++) mpz_init(coeff[i]);

        /* m = N^(1/(d+1)) */
        mpz_root(m, N, degree + 1);

        /* Base-m decomposition */
        mpz_t rem;
        mpz_init(rem);
        mpz_set(rem, N);
        for (int i = 0; i < degree; i++) {
            mpz_fdiv_qr(rem, coeff[i], rem, m);
        }
        mpz_set(coeff[degree], rem);
        mpz_clear(rem);

        /* Verify */
        mpz_set(tmp, coeff[degree]);
        for (int i = degree - 1; i >= 0; i--) {
            mpz_mul(tmp, tmp, m);
            mpz_add(tmp, tmp, coeff[i]);
        }
        int poly_ok = (mpz_cmp(tmp, N) == 0);

        double log2_m = log2(mpz_get_d(m));
        printf("\nNFS degree %d: m = %.0f (~2^%.1f), poly %s\n",
               degree, mpz_get_d(m), log2_m, poly_ok ? "OK" : "FAIL");
        if (!poly_ok) { printf("  SKIPPING due to polynomial failure\n"); continue; }

        /* Search over (a,b) pairs */
        int a_max = (int)(pow(mpz_get_d(m), 0.3));
        if (a_max < 100) a_max = 100;
        if (a_max > 5000) a_max = 5000;
        int b_max = a_max / 2;
        if (b_max < 10) b_max = 10;

        int nfs_total = 0, nfs_rat_smooth = 0, nfs_alg_smooth = 0, nfs_both = 0;
        long nfs_rat_bits_sum = 0, nfs_alg_bits_sum = 0;

        for (int b = 1; b <= b_max; b++) {
            for (int a = -a_max; a <= a_max; a++) {
                if (a == 0) continue;
                if (gcd_int(a, b) != 1) continue;

                /* Rational norm: |a - b*m| */
                mpz_set_si(val, a);
                mpz_submul_ui(val, m, b);
                mpz_abs(cof, val);

                int rat_bits = mpz_sizeinbase(cof, 2);
                nfs_rat_bits_sum += rat_bits;

                /* Trial divide rational side */
                for (int i = 0; i < n_primes; i++) {
                    while (mpz_divisible_ui_p(cof, primes[i]))
                        mpz_divexact_ui(cof, cof, primes[i]);
                }
                int rat_smooth = (mpz_cmp_ui(cof, 1) == 0);
                if (rat_smooth) nfs_rat_smooth++;

                /* Algebraic norm: f_hom(a, b) = sum c_i * a^i * b^(d-i) */
                mpz_set_ui(val, 0);
                for (int i = 0; i <= degree; i++) {
                    mpz_set(tmp, coeff[i]);
                    for (int j = 0; j < i; j++) mpz_mul_si(tmp, tmp, a);
                    for (int j = 0; j < degree - i; j++) mpz_mul_si(tmp, tmp, b);
                    mpz_add(val, val, tmp);
                }
                mpz_abs(cof, val);

                int alg_bits = mpz_sizeinbase(cof, 2);
                nfs_alg_bits_sum += alg_bits;

                /* Trial divide algebraic side */
                for (int i = 0; i < n_primes; i++) {
                    while (mpz_divisible_ui_p(cof, primes[i]))
                        mpz_divexact_ui(cof, cof, primes[i]);
                }
                int alg_smooth = (mpz_cmp_ui(cof, 1) == 0);
                if (alg_smooth) nfs_alg_smooth++;

                if (rat_smooth && alg_smooth) nfs_both++;
                nfs_total++;
            }
        }

        printf("  search range: a in [-%d,%d], b in [1,%d] -> %d candidates\n",
               a_max, a_max, b_max, nfs_total);
        printf("  avg rational norm: %.1f bits\n", (double)nfs_rat_bits_sum / nfs_total);
        printf("  avg algebraic norm: %.1f bits\n", (double)nfs_alg_bits_sum / nfs_total);
        printf("  rational smooth: %d (%.4f%%)\n", nfs_rat_smooth, 100.0*nfs_rat_smooth/nfs_total);
        printf("  algebraic smooth: %d (%.4f%%)\n", nfs_alg_smooth, 100.0*nfs_alg_smooth/nfs_total);
        printf("  BOTH smooth: %d (%.4f%%)\n", nfs_both, 100.0*nfs_both/nfs_total);

        /* Compare with QS */
        double qs_rate = (double)(qs_smooth + qs_1lp) / qs_total;
        double nfs_rate = (double)nfs_both / nfs_total;
        printf("  Comparison: QS smooth rate = %.4f%%, NFS both-smooth rate = %.4f%%\n",
               qs_rate * 100, nfs_rate * 100);
        printf("  NFS/QS ratio: %.2f\n", nfs_rate / qs_rate);

        for (int i = 0; i <= degree; i++) mpz_clear(coeff[i]);
        mpz_clear(m);
    }

    mpz_clear(N); mpz_clear(sqrtN);
    mpz_clear(val); mpz_clear(cof); mpz_clear(tmp);
    return 0;
}
