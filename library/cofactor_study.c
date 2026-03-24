/*
 * cofactor_study.c — Study cofactor distribution in QS-style sieving
 *
 * For a given N, generate many QS candidates (x² - N), extract the
 * smooth part (trial divide by primes up to B), and analyze:
 * 1. Distribution of cofactor sizes
 * 2. Collision rate among cofactors (shared prime factors)
 * 3. Effectiveness of different cofactor bounds
 *
 * This data guides development of novel factoring approaches.
 *
 * Usage: ./cofactor_study <N> [fb_size]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

#define MAX_FB 4096
#define MAX_CAND 200000
#define HIST_BINS 64

static unsigned long primes[MAX_FB];
static int nprimes = 0;

static int is_prime(unsigned long n) {
    if (n < 2) return 0;
    if (n < 4) return 1;
    if (n % 2 == 0 || n % 3 == 0) return 0;
    for (unsigned long d = 5; d * d <= n; d += 6)
        if (n % d == 0 || n % (d+2) == 0) return 0;
    return 1;
}

static unsigned long modsqrt(unsigned long a, unsigned long p) {
    if (a == 0 || p == 2) return a % p;
    unsigned long e = (p-1)/2, b = a, r = 1;
    while (e > 0) { if (e&1) r = (__uint128_t)r*b%p; b = (__uint128_t)b*b%p; e >>= 1; }
    if (r != 1) return 0;
    if (p % 4 == 3) {
        b = a; e = (p+1)/4; r = 1;
        while (e > 0) { if (e&1) r = (__uint128_t)r*b%p; b = (__uint128_t)b*b%p; e >>= 1; }
        return r;
    }
    /* Full Tonelli-Shanks */
    unsigned long Q = p-1, S = 0;
    while (!(Q&1)) { Q >>= 1; S++; }
    unsigned long z = 2;
    while (1) { b = z; e = (p-1)/2; r = 1; while (e > 0) { if (e&1) r = (__uint128_t)r*b%p; b = (__uint128_t)b*b%p; e >>= 1; } if (r == p-1) break; z++; }
    unsigned long M = S, c, t, R;
    b = z; e = Q; c = 1; while (e > 0) { if (e&1) c = (__uint128_t)c*b%p; b = (__uint128_t)b*b%p; e >>= 1; }
    b = a; e = Q; t = 1; while (e > 0) { if (e&1) t = (__uint128_t)t*b%p; b = (__uint128_t)b*b%p; e >>= 1; }
    b = a; e = (Q+1)/2; R = 1; while (e > 0) { if (e&1) R = (__uint128_t)R*b%p; b = (__uint128_t)b*b%p; e >>= 1; }
    while (t != 1) {
        unsigned long i = 0, tt = t;
        while (tt != 1) { tt = (__uint128_t)tt*tt%p; i++; }
        b = c; for (unsigned long j = 0; j < M-i-1; j++) b = (__uint128_t)b*b%p;
        M = i; c = (__uint128_t)b*b%p; t = (__uint128_t)t*c%p; R = (__uint128_t)R*b%p;
    }
    return R;
}

int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N> [fb_size]\n", argv[0]); return 1; }

    mpz_t N, sqrtN, x, val, cof;
    mpz_init(N); mpz_init(sqrtN); mpz_init(x); mpz_init(val); mpz_init(cof);
    mpz_set_str(N, argv[1], 10);
    int digits = strlen(argv[1]);

    int fb_target = 200;
    if (argc > 2) fb_target = atoi(argv[2]);
    if (fb_target > MAX_FB) fb_target = MAX_FB;

    /* Build factor base */
    nprimes = 0;
    for (unsigned long p = 2; nprimes < fb_target; p++) {
        if (!is_prime(p)) continue;
        unsigned long nmp = mpz_fdiv_ui(N, p);
        if (p > 2 && modsqrt(nmp, p) == 0 && nmp != 0) continue;
        primes[nprimes++] = p;
    }
    printf("# %d-digit N, factor base: %d primes up to %lu\n", digits, nprimes, primes[nprimes-1]);

    mpz_sqrt(sqrtN, N);
    mpz_t tmp; mpz_init(tmp);
    mpz_mul(tmp, sqrtN, sqrtN);
    if (mpz_cmp(tmp, N) < 0) mpz_add_ui(sqrtN, sqrtN, 1);
    mpz_clear(tmp);

    /* Generate candidates and study cofactors */
    int n_total = 0, n_full = 0, n_1lp = 0, n_2lp = 0, n_large = 0;
    int cofactor_bits_hist[HIST_BINS] = {0};

    /* Track cofactors for collision analysis */
    mpz_t *cofactors = malloc(MAX_CAND * sizeof(mpz_t));
    int n_cofactors = 0;

    int sieve_range = (fb_target < 500) ? 100000 : 500000;

    for (int t = 1; t <= sieve_range && n_cofactors < MAX_CAND; t++) {
        /* Positive side */
        mpz_set(x, sqrtN);
        mpz_add_ui(x, x, t);
        mpz_mul(val, x, x);
        mpz_sub(val, val, N);
        mpz_abs(cof, val);

        /* Trial divide */
        for (int i = 0; i < nprimes; i++) {
            while (mpz_divisible_ui_p(cof, primes[i]))
                mpz_divexact_ui(cof, cof, primes[i]);
        }

        n_total++;
        int bits = mpz_sizeinbase(cof, 2);
        if (bits <= 1) { n_full++; bits = 0; }
        else if (bits <= 30) n_1lp++;
        else if (bits <= 62) n_2lp++;
        else n_large++;

        int bin = bits * HIST_BINS / 200; /* 200 bits max */
        if (bin >= HIST_BINS) bin = HIST_BINS - 1;
        cofactor_bits_hist[bin]++;

        /* Save cofactors for collision analysis (only medium-sized ones) */
        if (bits > 1 && bits <= 62 && n_cofactors < MAX_CAND) {
            mpz_init_set(cofactors[n_cofactors], cof);
            n_cofactors++;
        }

        /* Negative side too */
        if (t <= sieve_range / 2) {
            mpz_set(x, sqrtN);
            mpz_sub_ui(x, x, t);
            mpz_mul(val, x, x);
            mpz_sub(val, val, N);
            mpz_abs(cof, val);
            for (int i = 0; i < nprimes; i++) {
                while (mpz_divisible_ui_p(cof, primes[i]))
                    mpz_divexact_ui(cof, cof, primes[i]);
            }
            n_total++;
            bits = mpz_sizeinbase(cof, 2);
            if (bits <= 1) { n_full++; bits = 0; }
            else if (bits <= 30) n_1lp++;
            else if (bits <= 62) n_2lp++;
            else n_large++;
            bin = bits * HIST_BINS / 200;
            if (bin >= HIST_BINS) bin = HIST_BINS - 1;
            cofactor_bits_hist[bin]++;
            if (bits > 1 && bits <= 62 && n_cofactors < MAX_CAND) {
                mpz_init_set(cofactors[n_cofactors], cof);
                n_cofactors++;
            }
        }
    }

    printf("# Candidates tested: %d\n", n_total);
    printf("# Fully smooth:     %d (%.2f%%)\n", n_full, 100.0*n_full/n_total);
    printf("# 1 large prime:    %d (%.2f%%)\n", n_1lp, 100.0*n_1lp/n_total);
    printf("# 2 large primes:   %d (%.2f%%)\n", n_2lp, 100.0*n_2lp/n_total);
    printf("# Too large:        %d (%.2f%%)\n", n_large, 100.0*n_large/n_total);

    printf("# Cofactor size distribution (bits -> count):\n");
    for (int i = 0; i < HIST_BINS; i++) {
        if (cofactor_bits_hist[i] > 0) {
            int lo = i * 200 / HIST_BINS;
            int hi = (i+1) * 200 / HIST_BINS - 1;
            printf("#   %d-%d bits: %d\n", lo, hi, cofactor_bits_hist[i]);
        }
    }

    /* Collision analysis: find cofactors that share a prime factor */
    printf("# Collision analysis (%d cofactors):\n", n_cofactors);
    if (n_cofactors > 1) {
        /* Count exact matches */
        int exact_matches = 0;
        /* Use batch GCD: compute product of all cofactors, then gcd each with product/self */
        /* For efficiency, use a simpler approach: product tree */

        /* Build product of all cofactors */
        mpz_t prod;
        mpz_init_set_ui(prod, 1);
        int batch_size = (n_cofactors > 50000) ? 50000 : n_cofactors;
        for (int i = 0; i < batch_size; i++)
            mpz_mul(prod, prod, cofactors[i]);

        /* For each cofactor, compute gcd(cof, prod/cof) */
        int shared_factor = 0;
        mpz_t g, rem;
        mpz_init(g); mpz_init(rem);
        for (int i = 0; i < batch_size; i++) {
            mpz_divexact(rem, prod, cofactors[i]);
            mpz_gcd(g, cofactors[i], rem);
            if (mpz_cmp_ui(g, 1) > 0)
                shared_factor++;
        }
        printf("#   Cofactors with shared factor: %d/%d (%.2f%%)\n",
               shared_factor, batch_size, 100.0*shared_factor/batch_size);

        /* Also count exact cofactor matches (same value) */
        /* Quick sort-based count */
        /* For simplicity, hash-based */
        int hash_size = batch_size * 3;
        unsigned long *hashes = calloc(hash_size, sizeof(unsigned long));
        exact_matches = 0;
        for (int i = 0; i < batch_size; i++) {
            unsigned long h = mpz_fdiv_ui(cofactors[i], hash_size);
            int found = 0;
            for (int p = 0; p < hash_size; p++) {
                int slot = (h + p) % hash_size;
                if (hashes[slot] == 0) { hashes[slot] = i + 1; break; }
                /* Very rough check: hash collision means possible exact match */
            }
        }
        free(hashes);

        mpz_clear(prod); mpz_clear(g); mpz_clear(rem);
    }

    /* Cleanup */
    for (int i = 0; i < n_cofactors; i++) mpz_clear(cofactors[i]);
    free(cofactors);
    mpz_clear(N); mpz_clear(sqrtN); mpz_clear(x); mpz_clear(val); mpz_clear(cof);
    return 0;
}
