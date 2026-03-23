// Power Smooth Factoring - Novel approach
//
// Key idea: Instead of sieving x^2 mod N for smooth numbers (QS approach),
// use modular exponentiation to generate numbers that are GUARANTEED to have
// specific prime factors, then check the remaining cofactor for smoothness.
//
// Method:
// 1. Compute a = g^(B#) mod N where B# = product of primes up to B (primorial)
// 2. a-1 mod N is guaranteed to be divisible by p-1 (if p-1 | B#) for each prime p | N
// 3. But we need smooth numbers in the "QS sense": x^2 ≡ smooth (mod N)
//
// Modified approach: "Power-residue assisted QS"
// For each FB prime p, compute r_p = g^(k_p) mod N where k_p is chosen
// so that r_p^2 mod N is partially smooth. The partial smoothness from
// the power structure reduces the sieving needed.
//
// Even simpler: use the QS framework but with a modified polynomial:
// Instead of Q(x) = (x + sqrt(N))^2 - N, use
// Q(x) = (g^x mod N)^2 - something
// where the power structure helps with smoothness.
//
// Actually, the most promising variant: "Structured Quadratic Residues"
// For each x, compute y = x^2 mod N and check if y is smooth.
// But instead of random x, choose x from a carefully structured set
// where x^2 mod N tends to be smoother.
//
// Strategy: Let S be the set of B-smooth numbers up to some bound.
// For each s in S, s is already smooth by construction.
// Compute x = s mod N. Then x^2 mod N = s^2 mod N.
// s^2 mod N is a random-looking number mod N, so not useful.
//
// Better strategy: Use the factored form of the convergents.
// In CFRAC, the convergents p_k/q_k give p_k^2 ≡ small value (mod N).
// What if we use Lehmer's "nearby" trick: find numbers NEAR p_k that
// are also squares, and whose Q values share factors?
//
// This is actually the "smooth number chain" idea from the filename.
// Let me implement it properly:
//
// For consecutive CF convergents k and k+1:
// p_k^2 ≡ A_k (mod N)  and  p_{k+1}^2 ≡ A_{k+1} (mod N)
// Then (p_k * p_{k+1})^2 ≡ A_k * A_{k+1} (mod N)
// If A_k * A_{k+1} is smooth, we have a relation!
// The product A_k * A_{k+1} tends to share factors with A_k and A_{k+1} individually.
//
// More generally: for any subset S of convergent indices,
// (∏_{k in S} p_k)^2 ≡ ∏_{k in S} A_k (mod N)
// We need the product of A_k values to be smooth.
//
// If individual A_k values are ~2*sqrt(N) ~ N^{1/2}, then a product of t values
// is ~N^{t/2}. This is LARGER, so products are LESS smooth, not more.
//
// The trick: find SUBSETS where the product has many shared/cancelling factors.
// This is equivalent to the standard linear algebra step, just viewed differently.
//
// So this doesn't actually help beyond standard CFRAC. Let me try something else.
//
// NOVEL APPROACH: "Residue Amplification Sieve" (RAS)
//
// Observation: In QS, we sieve g(x) = a*x^2 + 2*b*x + c for smooth values.
// The size of g(x) is ~sqrt(2N)*M for |x| ~ M.
// The smoothness probability is ~u^(-u) where u = log(sqrt(2N)*M) / log(B).
//
// What if we could REDUCE the effective size of the numbers being tested?
//
// Key idea: Use GCD-based reduction. For each candidate x, compute
// g_reduced(x) = g(x) / gcd(g(x), P)
// where P is a precomputed product of "helper" values.
//
// If P contains many factors of g(x), then g_reduced(x) is much smaller
// and more likely to be smooth over the remaining factor base.
//
// How to construct P? Use product trees!
// P = product of g(y) for many y values.
// Then gcd(g(x), P) removes shared factors between g(x) and the g(y) values.
//
// This is actually Bernstein's batch smoothness idea, but used for PARTIAL smoothness.
//
// Implementation plan:
// 1. Use SIQS to generate many g(x) values
// 2. Compute P = product of g(y) for a batch of y values
// 3. For each g(x), compute gcd(g(x), P) to extract shared factors
// 4. Check if the remaining cofactor is smooth
//
// This should increase the yield of smooth numbers from each batch of candidates.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <gmp.h>

// This is a proof-of-concept. The actual implementation follows the SIQS
// framework but with batch GCD for partial factor extraction.

static struct timespec g_start;
static double elapsed_sec() {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (now.tv_sec - g_start.tv_sec) + (now.tv_nsec - g_start.tv_nsec) / 1e9;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <number>\n", argv[0]);
        return 1;
    }

    clock_gettime(CLOCK_MONOTONIC, &g_start);

    mpz_t N, sqrtN;
    mpz_init(N);
    mpz_init(sqrtN);

    if (mpz_set_str(N, argv[1], 10) != 0) {
        fprintf(stderr, "Invalid number\n");
        return 1;
    }

    mpz_sqrt(sqrtN, N);

    size_t digits = mpz_sizeinbase(N, 10);

    // Compute primorial P# (product of primes up to some bound)
    unsigned long prim_bound = 100000;
    if (digits <= 30) prim_bound = 10000;
    else if (digits <= 40) prim_bound = 50000;
    else if (digits <= 50) prim_bound = 200000;

    mpz_t primorial;
    mpz_init(primorial);
    mpz_set_ui(primorial, 1);

    // Sieve primes up to prim_bound
    std::vector<bool> is_prime(prim_bound + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (unsigned long i = 2; i * i <= prim_bound; i++)
        if (is_prime[i])
            for (unsigned long j = i*i; j <= prim_bound; j += i)
                is_prime[j] = false;

    for (unsigned long i = 2; i <= prim_bound; i++) {
        if (is_prime[i]) {
            // Include p^k where p^k <= prim_bound
            unsigned long pk = i;
            while (pk <= prim_bound) {
                mpz_mul_ui(primorial, primorial, i);
                pk *= i;
            }
        }
    }

    fprintf(stderr, "RAS: %zu digits, primorial %zu bits\n",
            digits, mpz_sizeinbase(primorial, 2));

    // Pollard p-1 style: compute g = 2^primorial mod N
    mpz_t g;
    mpz_init(g);
    mpz_set_ui(g, 2);
    mpz_powm(g, g, primorial, N);

    // gcd(g-1, N) might give a factor if p-1 is prim_bound-smooth
    mpz_t test;
    mpz_init(test);
    mpz_sub_ui(test, g, 1);
    mpz_t factor;
    mpz_init(factor);
    mpz_gcd(factor, test, N);

    if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0) {
        mpz_t cofactor;
        mpz_init(cofactor);
        mpz_divexact(cofactor, N, factor);
        if (mpz_cmp(factor, cofactor) > 0) mpz_swap(factor, cofactor);
        gmp_printf("%Zd %Zd\n", factor, cofactor);
        fprintf(stderr, "RAS: factored via p-1 in %.3fs\n", elapsed_sec());
        mpz_clear(cofactor);
        mpz_clear(N); mpz_clear(sqrtN); mpz_clear(primorial); mpz_clear(g);
        mpz_clear(test); mpz_clear(factor);
        return 0;
    }

    // p-1 didn't work (expected for random semiprimes).
    // Try Williams p+1: if p+1 is smooth, we can detect it.
    // Lucas sequences: V_0 = 2, V_1 = P, V_n = P*V_{n-1} - V_{n-2}
    // If p+1 | primorial, then V_{primorial} ≡ 2 (mod p)

    // Try a few P values
    for (int P = 3; P <= 20; P++) {
        mpz_t V_prev, V_curr, V_next;
        mpz_init_set_ui(V_prev, 2);
        mpz_init_set_ui(V_curr, P);
        mpz_init(V_next);

        // Compute V_{primorial} mod N using the doubling formulas:
        // V_{2k} = V_k^2 - 2
        // V_{2k+1} = V_k * V_{k+1} - P

        // Binary method for V_n:
        size_t prim_bits = mpz_sizeinbase(primorial, 2);
        mpz_set_ui(V_prev, 2); // V_0
        mpz_set_ui(V_curr, P); // V_1

        for (size_t bit = prim_bits - 1; bit > 0; bit--) {
            if (mpz_tstbit(primorial, bit - 1)) {
                // V_{2k+1}: V_prev = V_k * V_{k+1} - P, V_curr = V_{k+1}^2 - 2
                mpz_mul(V_next, V_prev, V_curr);
                mpz_sub_ui(V_next, V_next, P);
                mpz_mod(V_next, V_next, N);

                mpz_mul(V_curr, V_curr, V_curr);
                mpz_sub_ui(V_curr, V_curr, 2);
                mpz_mod(V_curr, V_curr, N);

                mpz_set(V_prev, V_next);
            } else {
                // V_{2k}: V_curr = V_k * V_{k+1} - P, V_prev = V_k^2 - 2
                mpz_mul(V_next, V_prev, V_curr);
                mpz_sub_ui(V_next, V_next, P);
                mpz_mod(V_next, V_next, N);

                mpz_mul(V_prev, V_prev, V_prev);
                mpz_sub_ui(V_prev, V_prev, 2);
                mpz_mod(V_prev, V_prev, N);

                mpz_set(V_curr, V_next);
            }
        }

        // Check gcd(V_primorial - 2, N)
        mpz_sub_ui(test, V_prev, 2);
        mpz_gcd(factor, test, N);

        if (mpz_cmp_ui(factor, 1) > 0 && mpz_cmp(factor, N) < 0) {
            mpz_t cofactor;
            mpz_init(cofactor);
            mpz_divexact(cofactor, N, factor);
            if (mpz_cmp(factor, cofactor) > 0) mpz_swap(factor, cofactor);
            gmp_printf("%Zd %Zd\n", factor, cofactor);
            fprintf(stderr, "RAS: factored via p+1 (P=%d) in %.3fs\n", P, elapsed_sec());
            mpz_clear(cofactor);
            mpz_clear(V_prev); mpz_clear(V_curr); mpz_clear(V_next);
            mpz_clear(N); mpz_clear(sqrtN); mpz_clear(primorial); mpz_clear(g);
            mpz_clear(test); mpz_clear(factor);
            return 0;
        }

        mpz_clear(V_prev); mpz_clear(V_curr); mpz_clear(V_next);
    }

    // Neither p-1 nor p+1 worked. Try ECM-style with multiple curves.
    // For balanced semiprimes, this is expected.
    fprintf(stderr, "RAS: p-1 and p+1 methods failed (expected for balanced semiprimes)\n");
    fprintf(stderr, "FAIL: falling back to SIQS needed\n");

    mpz_clear(N); mpz_clear(sqrtN); mpz_clear(primorial); mpz_clear(g);
    mpz_clear(test); mpz_clear(factor);
    return 1;
}
