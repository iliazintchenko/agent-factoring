/*
 * order_probe.c — Experimental: Multiplicative order probing for factoring
 *
 * Novel approach: Instead of smooth numbers, exploit the MULTIPLICATIVE
 * ORDER structure of (Z/NZ)*. For N = pq:
 *   ord_N(a) = lcm(ord_p(a), ord_q(a))
 *
 * If we can find a such that:
 *   a^k ≡ 1 (mod p) but a^k ≢ 1 (mod q)  [or vice versa]
 * then gcd(a^k - 1, N) gives a factor.
 *
 * Key insight: for random a, ord_p(a) | (p-1) and ord_q(a) | (q-1).
 * The 2-adic valuations v_2(ord_p(a)) and v_2(ord_q(a)) are usually
 * different. By computing a^((N-1)/2^k) for increasing k, we can
 * detect when the "p part" and "q part" diverge.
 *
 * This is essentially the Miller-Rabin witness idea applied to factoring.
 * It's O(log N) per attempt, but the success probability per attempt
 * depends on the structure of p-1 and q-1.
 *
 * NOVEL TWIST: Instead of random a, use a = product of small primes
 * (smooth bases). For smooth a, the multiplicative order ord(a) is
 * constrained by the smooth part of the group order. By combining
 * information from many smooth bases, we might reconstruct enough
 * of the group structure to factor N.
 *
 * This is related to Pohlig-Hellman but applied to FACTORING, not DLP.
 *
 * Usage: ./order_probe <N>
 * Output: FACTOR: <p> (if successful)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

static double walltime(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/*
 * Method 1: Power-of-2 descent
 * Compute a^((N-1)/2^k) for k = 0, 1, 2, ...
 * When the result transitions from ≢1 to ≡1, check gcd.
 *
 * For N = pq, (N-1) = pq - 1.
 * a^((N-1)/2) = a^((pq-1)/2).
 *
 * This doesn't directly relate to ord_p(a) because (N-1)/2 ≠ (p-1)/2.
 *
 * Better: use the "strong pseudoprime" idea.
 * Write N-1 = 2^s * d (d odd).
 * Compute a^d, then square repeatedly.
 * If a^(2^r * d) ≡ -1 (mod N) for some r, it's a "strong pseudoprime".
 * If a^d ≡ 1 (mod N), also strong pseudoprime.
 * Otherwise: gcd(a^(2^r * d) - 1, N) or gcd(a^(2^r * d) + 1, N) might factor.
 *
 * For N = pq: this works when v_2(ord_p(a)) ≠ v_2(ord_q(a)).
 */
static int method_strong_test(mpz_t N, mpz_t factor) {
    mpz_t d, a, x, g;
    mpz_init(d); mpz_init(a); mpz_init(x); mpz_init(g);

    /* N-1 = 2^s * d */
    mpz_sub_ui(d, N, 1);
    unsigned long s = 0;
    while (mpz_even_p(d)) { mpz_divexact_ui(d, d, 2); s++; }

    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, 42);

    int found = 0;
    for (int trial = 0; trial < 10000 && !found; trial++) {
        /* Random base a in [2, N-2] */
        mpz_urandomm(a, rng, N);
        if (mpz_cmp_ui(a, 2) < 0) mpz_set_ui(a, 2);

        /* x = a^d mod N */
        mpz_powm(x, a, d, N);

        if (mpz_cmp_ui(x, 1) == 0 || mpz_cmp(x, d) == 0) continue; /* trivial */

        /* Check if x ≡ -1 mod N */
        mpz_t nm1;
        mpz_init(nm1);
        mpz_sub_ui(nm1, N, 1);
        if (mpz_cmp(x, nm1) == 0) { mpz_clear(nm1); continue; }

        /* Square x repeatedly */
        for (unsigned long r = 1; r < s; r++) {
            mpz_t prev;
            mpz_init_set(prev, x);
            mpz_mul(x, x, x);
            mpz_mod(x, x, N);

            if (mpz_cmp_ui(x, 1) == 0) {
                /* prev^2 ≡ 1 but prev ≢ ±1 → non-trivial square root of 1 */
                mpz_sub_ui(prev, prev, 1);
                mpz_gcd(g, prev, N);
                if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                    mpz_set(factor, g);
                    found = 1;
                }
                mpz_clear(prev);
                break;
            }

            if (mpz_cmp(x, nm1) == 0) {
                mpz_clear(prev);
                break; /* -1 reached, no useful info */
            }
            mpz_clear(prev);
        }
        mpz_clear(nm1);
    }

    gmp_randclear(rng);
    mpz_clear(d); mpz_clear(a); mpz_clear(x); mpz_clear(g);
    return found;
}

/*
 * Method 2: Multi-prime order accumulation
 * For each small prime ℓ, compute a^((N^2-1)/ℓ) mod N.
 * If a has order divisible by ℓ in (Z/pZ)* but not (Z/qZ)*,
 * then a^((N^2-1)/ℓ) ≢ 1 (mod p) but ≡ 1 (mod q).
 * Then gcd(a^((N^2-1)/ℓ) - 1, N) = q.
 *
 * Why N^2-1? Because |GL(1, Z/NZ)| = φ(N) divides N^2-1...
 * Actually φ(N) = (p-1)(q-1) doesn't divide N^2-1 in general.
 *
 * Better: compute M = lcm(1, 2, ..., B) for increasing B.
 * For each B, compute a^M mod N and check gcd(a^M - 1, N).
 * This is Pollard's p-1 method! But it only works when p-1 is B-smooth.
 *
 * Novel twist: instead of a single base a, use MANY bases.
 * For each base a_i, compute a_i^M mod N. Then:
 * gcd(∏(a_i^M - 1), N) has a better chance of revealing p.
 *
 * Even more novel: use the PRODUCTS a_i^M to detect structure.
 * If a_1^M ≡ 1 (mod p) and a_2^M ≢ 1 (mod p), then
 * a_1^M - a_2^M ≡ -a_2^M ≢ 0 (mod p), so gcd(a_1^M - a_2^M, N)
 * doesn't help. But if we compute gcd(a_1^M * a_2^(-M) - 1, N)...
 * that's gcd((a_1/a_2)^M - 1, N), which is just the p-1 test with
 * base a_1/a_2.
 *
 * So multiple bases don't help for p-1 — each base is independent.
 * Unless we use a MORE SOPHISTICATED combination.
 */
static int method_multi_base_p1(mpz_t N, mpz_t factor) {
    /* Pollard p-1 with multiple bases and accumulation */
    mpz_t a, prod, g;
    mpz_init(a); mpz_init(prod); mpz_init(g);
    mpz_set_ui(prod, 1);

    int found = 0;

    /* Try with base 2, accumulating prime powers */
    mpz_set_ui(a, 2);
    for (unsigned long p = 2; p < 100000 && !found; p++) {
        /* Check primality (simple) */
        int is_p = 1;
        for (unsigned long d = 2; d * d <= p; d++)
            if (p % d == 0) { is_p = 0; break; }
        if (!is_p) continue;

        /* Compute a = a^p mod N (accumulate) */
        unsigned long pp = p;
        while (pp <= 100000) { /* include prime powers */
            mpz_powm_ui(a, a, p, N);
            pp *= p;
        }

        /* Periodically check gcd */
        if (p % 100 == 99 || p > 99000) {
            mpz_sub_ui(g, a, 1);
            mpz_gcd(g, g, N);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                mpz_set(factor, g);
                found = 1;
            } else if (mpz_cmp(g, N) == 0) {
                /* Overshot — need smaller step. p-1 is very smooth. */
                /* This shouldn't happen for balanced semiprimes */
            }
        }
    }

    mpz_clear(a); mpz_clear(prod); mpz_clear(g);
    return found;
}

/*
 * Method 3: Quadratic character divergence
 * For small primes ℓ, the Legendre symbol (ℓ/p) and (ℓ/q) can differ.
 * The Jacobi symbol (ℓ/N) = (ℓ/p)(ℓ/q).
 *
 * If (ℓ/N) = -1: then (ℓ/p) ≠ (ℓ/q) — we know they differ!
 *   But we can't extract (ℓ/p) or (ℓ/q) individually.
 *
 * If (ℓ/N) = 1: then either both are +1 or both are -1.
 *
 * Can we use the PATTERN of Jacobi symbols to factor?
 *
 * For each ℓ where (ℓ/N) = -1, we know (ℓ/p) ≠ (ℓ/q).
 * By quadratic reciprocity: (ℓ/p) relates to (p/ℓ).
 * (p/ℓ) is determined by p mod ℓ.
 *
 * If we could determine (ℓ/p) for enough primes ℓ, we could
 * reconstruct p using CRT + Jacobi constraints.
 *
 * This is the "genus theory" approach:
 * The genus of p mod 4d (where d = discriminant) determines the
 * Jacobi symbols. But extracting the genus requires factoring.
 *
 * Novel idea: use the CORRELATION between Jacobi symbols of
 * different ℓ values. For ℓ₁ and ℓ₂:
 * (ℓ₁/p)(ℓ₂/p) = (ℓ₁ℓ₂/p)
 *
 * If (ℓ₁/N) = (ℓ₂/N) = -1, then (ℓ₁/p) ≠ (ℓ₂/p) each.
 * But (ℓ₁ℓ₂/N) = (ℓ₁/N)(ℓ₂/N) = 1.
 * And (ℓ₁ℓ₂/p) = (ℓ₁/p)(ℓ₂/p) could be ±1.
 *
 * The product of two QNRs mod p is a QR. So (ℓ₁ℓ₂/p) = +1.
 * Similarly (ℓ₁ℓ₂/q) = +1.
 * So (ℓ₁ℓ₂/N) = 1, which we already knew. No new information.
 *
 * What about THREE primes? (ℓ₁ℓ₂ℓ₃/N) = (ℓ₁/N)(ℓ₂/N)(ℓ₃/N).
 * Still just a product of known values.
 *
 * Conclusion: Jacobi symbols alone can't factor.
 */

/*
 * Method 4: Order-based factoring via Pohlig-Hellman structure
 *
 * For each small prime ℓ, determine whether ℓ | ord_p(a) - ord_q(a).
 * If so, we can separate the p and q components.
 *
 * Compute a^(φ/ℓ) where φ = (p-1)(q-1). Problem: we don't know φ.
 *
 * But: for a known to have smooth order, we can try:
 * a^(M/ℓ) where M = a large smooth number likely to be a multiple of φ.
 *
 * If M = φ: a^(M/ℓ) ≡ a^(φ/ℓ) mod N. This is non-trivial mod p iff
 * ℓ | ord_p(a), and non-trivial mod q iff ℓ | ord_q(a).
 * If the two differ: gcd(a^(M/ℓ) - 1, N) gives a factor.
 *
 * This is exactly Pollard p-1 Stage 2 / Williams p+1 idea.
 */

int main(int argc, char *argv[]) {
    if (argc < 2) { fprintf(stderr, "Usage: %s <N>\n", argv[0]); return 1; }

    double start = walltime();
    mpz_t N, factor;
    mpz_init(N); mpz_init(factor);
    mpz_set_str(N, argv[1], 10);

    /* Trial division */
    for (unsigned long p = 2; p < 1000000; p++) {
        if (mpz_divisible_ui_p(N, p)) {
            printf("FACTOR: %lu\n", p);
            return 0;
        }
    }

    fprintf(stderr, "OrderProbe: trying strong pseudoprime method...\n");
    if (method_strong_test(N, factor)) {
        mpz_t cof; mpz_init(cof);
        mpz_divexact(cof, N, factor);
        if (mpz_cmp(factor, cof) > 0) mpz_set(factor, cof);
        gmp_printf("FACTOR: %Zd\n", factor);
        fprintf(stderr, "OrderProbe: success (strong test) in %.3fs\n", walltime() - start);
        return 0;
    }

    fprintf(stderr, "OrderProbe: trying multi-base p-1 method...\n");
    if (method_multi_base_p1(N, factor)) {
        mpz_t cof; mpz_init(cof);
        mpz_divexact(cof, N, factor);
        if (mpz_cmp(factor, cof) > 0) mpz_set(factor, cof);
        gmp_printf("FACTOR: %Zd\n", factor);
        fprintf(stderr, "OrderProbe: success (p-1) in %.3fs\n", walltime() - start);
        return 0;
    }

    fprintf(stderr, "OrderProbe: all methods failed (%.3fs)\n", walltime() - start);
    mpz_clear(N); mpz_clear(factor);
    return 1;
}
