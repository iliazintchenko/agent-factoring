/*
 * Smooth Number Amplification via Multiplicative Combinations (SNAMC)
 *
 * Novel approach: instead of sieving a polynomial, generate B-smooth
 * numbers near sqrt(N) by enumerating products of factor base primes.
 * For each smooth s ≈ sqrt(N), compute s^2 mod N and check smoothness.
 *
 * If s^2 - N is also smooth, we have a relation: s^2 ≡ (s^2 - N) (mod N).
 *
 * Build: gcc -O3 -o snamc snamc.c -lgmp -lm
 */

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef unsigned int u32;
typedef unsigned long long u64;

#define MAX_FB 4096

static mpz_t g_N, g_sqrtN;
static int g_digits;
static u32 g_primes[MAX_FB];
static int g_nprimes;
static u32 g_B;

/* Generate primes up to B */
static void gen_primes(void) {
    char *sieve = calloc(g_B + 1, 1);
    for (u32 i = 2; i <= g_B; i++) sieve[i] = 1;
    for (u32 i = 2; (u64)i * i <= g_B; i++)
        if (sieve[i])
            for (u32 j = i * i; j <= g_B; j += i) sieve[j] = 0;
    g_nprimes = 0;
    for (u32 i = 2; i <= g_B && g_nprimes < MAX_FB; i++)
        if (sieve[i]) g_primes[g_nprimes++] = i;
    free(sieve);
}

/* Check if val is B-smooth. Returns 1 if smooth, fills exponent vector. */
static int is_smooth(mpz_t val, u64 *vec, int vec_words) {
    mpz_t rem;
    mpz_init_set(rem, val);
    memset(vec, 0, vec_words * sizeof(u64));

    int neg = mpz_sgn(rem) < 0;
    if (neg) { mpz_neg(rem, rem); vec[0] |= 1ULL; }

    for (int i = 0; i < g_nprimes; i++) {
        int exp = 0;
        while (mpz_fdiv_ui(rem, g_primes[i]) == 0) {
            mpz_fdiv_q_ui(rem, rem, g_primes[i]);
            exp++;
        }
        if (exp & 1) {
            int bit = i + 1;
            vec[bit / 64] ^= 1ULL << (bit % 64);
        }
    }

    int ok = (mpz_cmp_ui(rem, 1) == 0);
    mpz_clear(rem);
    return ok;
}

/* Enumerate B-smooth numbers near sqrt(N) using DFS on prime products.
   For each smooth s with |s - sqrt(N)| < D, compute s^2 - N and test. */
static int found_smooth = 0, tested = 0;

static void enumerate(mpz_t current, int prime_idx, mpz_t target, mpz_t D,
                       int vec_words) {
    /* Check if current is within distance D of target (sqrt(N)) */
    mpz_t diff, sq, residual;
    mpz_init(diff); mpz_init(sq); mpz_init(residual);

    mpz_sub(diff, current, target);
    mpz_abs(diff, diff);

    if (mpz_cmp(diff, D) <= 0 && mpz_cmp_ui(current, 1) > 0) {
        tested++;
        /* Compute s^2 - N */
        mpz_mul(sq, current, current);
        mpz_sub(residual, sq, g_N);

        u64 vec[MAX_FB / 64 + 2];
        if (is_smooth(residual, vec, (g_nprimes + 1 + 63) / 64)) {
            found_smooth++;
            /* Found a relation! */
            gmp_fprintf(stderr, "  SMOOTH: s=%Zd, s^2-N=%Zd\n", current, residual);
        }

        if (tested % 100000 == 0)
            fprintf(stderr, "  tested=%d found=%d\n", tested, found_smooth);
    }

    /* Prune: if current > target + D, no point going higher */
    mpz_add(diff, target, D);
    if (mpz_cmp(current, diff) > 0) {
        mpz_clear(diff); mpz_clear(sq); mpz_clear(residual);
        return;
    }

    /* Try multiplying by each prime from prime_idx onward */
    for (int i = prime_idx; i < g_nprimes && i < 30; i++) {
        mpz_t next;
        mpz_init(next);
        mpz_mul_ui(next, current, g_primes[i]);

        /* If next > target + D, skip (pruning) */
        mpz_add(diff, target, D);
        if (mpz_cmp(next, diff) <= 0) {
            enumerate(next, i, target, D, vec_words);
        }
        mpz_clear(next);
    }

    mpz_clear(diff); mpz_clear(sq); mpz_clear(residual);
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    mpz_init_set_str(g_N, argv[1], 10);
    g_digits = (int)mpz_sizeinbase(g_N, 10);
    mpz_init(g_sqrtN);
    mpz_sqrt(g_sqrtN, g_N);

    /* Parameters */
    double logN = g_digits * log(10.0);
    g_B = (u32)exp(0.4 * sqrt(logN * log(logN)));
    if (g_B < 50) g_B = 50;
    if (g_B > 10000) g_B = 10000;

    gen_primes();

    /* D = distance around sqrt(N) to search */
    mpz_t D;
    mpz_init(D);
    /* D = B^2 so that |s^2 - N| ≈ 2*D*sqrt(N) ≈ 2*B^2*sqrt(N) */
    mpz_set_ui(D, (u64)g_B * g_B);

    fprintf(stderr, "SNAMC: %d digits, B=%u, primes=%d, D=",
            g_digits, g_B, g_nprimes);
    gmp_fprintf(stderr, "%Zd\n", D);

    /* Start enumeration from 1 */
    mpz_t start;
    mpz_init_set_ui(start, 1);

    int vec_words = (g_nprimes + 1 + 63) / 64;
    enumerate(start, 0, g_sqrtN, D, vec_words);

    fprintf(stderr, "Total: tested=%d, smooth=%d\n", tested, found_smooth);

    mpz_clear(start); mpz_clear(D);
    mpz_clear(g_N); mpz_clear(g_sqrtN);
    return (found_smooth > 0) ? 0 : 1;
}
