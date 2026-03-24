/*
 * spectral_walk.c — Novel factoring approach: "Spectral Random Walk"
 *
 * HYPOTHESIS: Random walks on (Z/NZ)* produce different collision statistics
 * depending on the group structure Z_{p-1} × Z_{q-1}. By analyzing walk
 * statistics (collision rates, cycle lengths, autocorrelations) we might
 * extract factoring information without finding smooth numbers.
 *
 * APPROACH:
 * 1. Run many short random walks on (Z/NZ)* with different starting points
 * 2. Project walks onto small moduli (track x_i mod m for small m)
 * 3. Detect period structure in projections — periods divide ord(x) in Z_N*
 * 4. The period structure reveals information about lcm(p-1, q-1)
 * 5. From period information, reconstruct p and q
 *
 * KEY INNOVATION: Instead of finding ONE period (like Pollard rho), find
 * MANY short partial periods simultaneously and combine them. Each walk
 * gives O(1) bits of information about the group structure. With enough
 * walks, we accumulate enough information to factor.
 *
 * THEORETICAL QUESTION: Is the information rate per walk step better
 * than Pollard rho's O(1) bit per √N steps?
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gmp.h>

/* Walk function: x → x^2 * g^e mod N, where e depends on low bits of x */
static void walk_step(mpz_t x, mpz_t N, mpz_t *multipliers, int n_mult) {
    /* Partition step: use low bits of x to choose multiplier */
    unsigned long low = mpz_fdiv_ui(x, n_mult);
    mpz_mul(x, x, x);
    mpz_mul(x, x, multipliers[low]);
    mpz_mod(x, x, N);
}

/* Distinguished point detection */
static int is_distinguished(mpz_t x, unsigned long mask) {
    return (mpz_fdiv_ui(x, mask + 1) == 0);
}

/*
 * Attempt 1: Parallel random walks with cross-collision detection.
 *
 * Run K walks simultaneously. Track (walk_id, step_count, x_value).
 * When two walks from DIFFERENT starting points collide (x_i == x_j mod N),
 * we might have x_i ≡ x_j mod p but x_i ≢ x_j mod q (or vice versa).
 * Then gcd(x_i - x_j, N) = p.
 *
 * This is essentially multi-start Pollard rho with birthday collision.
 * Expected: O(N^{1/4} / sqrt(K)) steps total. Still O(N^{1/4}) overall.
 */

/*
 * Attempt 2: Walk on quotient groups.
 *
 * For a small prime ℓ, compute x^((N-1)/ℓ) mod N.
 * If ℓ | (p-1) but ℓ ∤ (q-1), then x^((N-1)/ℓ) mod p is an ℓ-th root of unity,
 * but x^((N-1)/ℓ) mod q = x^((N-1)/ℓ mod (q-1)) mod q, which is "random".
 *
 * The distribution of x^((N-1)/ℓ) mod N reveals whether ℓ | (p-1) and/or ℓ | (q-1).
 *
 * Specifically:
 * - If ℓ | (p-1) and ℓ | (q-1): x^((N-1)/ℓ) takes ℓ^2 distinct values
 * - If ℓ | (p-1) xor ℓ | (q-1): x^((N-1)/ℓ) takes ℓ distinct values
 * - If ℓ ∤ (p-1) and ℓ ∤ (q-1): x^((N-1)/ℓ) takes 1 value (always 1)
 *
 * Wait — (N-1)/ℓ might not be integer if ℓ ∤ (N-1). Also (N-1) = pq-1.
 * For ℓ | (p-1): ℓ | (N-1) iff ℓ | (pq-1) = p(q-1) + (p-1).
 * Since ℓ | (p-1), we need ℓ | p(q-1). If ℓ ∤ p (which is true for ℓ < p),
 * then ℓ | (q-1). So ℓ | (N-1) iff ℓ | (p-1) AND ℓ | (q-1), or neither.
 *
 * This means (N-1)/ℓ is integer only when BOTH p-1 and q-1 are divisible by ℓ.
 * That doesn't help distinguish them.
 *
 * BETTER: Use x^((N^2-1)/ℓ) mod N. Since N^2-1 = (N-1)(N+1) = (pq-1)(pq+1).
 * Or use Jacobi symbol computations.
 */

/*
 * Attempt 3: Character sum distinguisher.
 *
 * Compute S(a) = sum_{x=0}^{M-1} e^{2πi * a*x^2 / N} for various a.
 * This is a Gauss sum variant. For a that's a QR mod p but QNR mod q,
 * |S(a)| ≈ √(M*p) instead of √(M*N).
 *
 * BUT: Computing this sum classically requires M additions.
 * For |S| to be statistically distinguishable, M must be ≈ N.
 * Total cost: N operations. Same as trial division.
 *
 * Not useful.
 */

/*
 * Attempt 4: Subgroup order testing via GCD accumulation.
 *
 * For random g, compute g^k for k = 2, 3, 4, 5, ...
 * Accumulate product P = ∏(g^k - 1) mod N.
 * Periodically check gcd(P, N).
 *
 * If ord_p(g) = r, then g^r ≡ 1 mod p, so p | (g^r - 1), so p | P.
 * gcd(P, N) reveals p when we hit k = r.
 *
 * This is just Pollard p-1! Only works if p-1 is smooth.
 *
 * TWIST: Instead of g^k, use g^{prime_k} (only prime exponents).
 * Then p | P when prime_k | ord_p(g). The order divides p-1.
 * With many random g's, the orders cover different factors of p-1.
 * If p-1 = ∏ qi^ei, and for each qi we find a g with qi | ord(g),
 * then combining GCDs gives p.
 *
 * The probability that a random g has qi | ord(g) is 1 - 1/qi.
 * So with ~2 random g's per qi, we cover all qi.
 * But we don't know which qi to test.
 *
 * FURTHER TWIST: Use batch GCD. Compute g^k for many k simultaneously
 * and batch-test gcd(g^k - 1, N).
 */

/*
 * Let me implement something concrete and testable:
 *
 * "Multi-base order fragment detection" (MOFD)
 *
 * Algorithm:
 * 1. Pick many random bases g_1, ..., g_K
 * 2. For each g_i, compute g_i^m mod N for smooth m (m = product of small primes)
 * 3. If g_i^m ≡ 1 (mod p), then p | gcd(g_i^m - 1, N)
 * 4. Accumulate: P = ∏_i gcd(g_i^m - 1, N)
 * 5. If we're lucky, P = p for some accumulated value
 *
 * This is essentially multi-start Pollard p-1 with varying B1.
 * For balanced semiprimes, p-1 is unlikely smooth, so this fails.
 *
 * UNLESS we use ECM-style curves instead of the multiplicative group.
 * Each curve gives a different group order near p, and SOME might be smooth.
 *
 * This is literally ECM. Not novel.
 */

/*
 * OK, let me try something genuinely different:
 *
 * "Collision-based period extraction" (CPE)
 *
 * For a random base g, the sequence g, g^2, g^3, ... mod N has period r = ord(g).
 * In Z_{p-1} × Z_{q-1}, this is (g mod p)^1, (g mod p)^2, ... with period r_p = ord_p(g)
 * and similarly r_q = ord_q(g). So r = lcm(r_p, r_q).
 *
 * If r_p and r_q are coprime, we can extract r_p and r_q from r.
 * But we don't know r either (that's the whole problem).
 *
 * IDEA: Instead of finding r, find a MULTIPLE of r_p that's NOT a multiple of r_q.
 * Then g^{kr_p} ≡ 1 (mod p) but g^{kr_p} ≢ 1 (mod q), giving gcd(g^{kr_p} - 1, N) = p.
 *
 * How to find kr_p? If we compute gcd(g^m - 1, N) for various m, and get p
 * (not 1 or N), then m is a multiple of r_p but not r_q.
 *
 * Trying random m: the probability that m is a multiple of r_p is 1/r_p.
 * For r_p ≈ p ≈ √N, this is exponentially unlikely.
 *
 * Trying STRUCTURED m: if m = ∏(small primes), this is Pollard p-1.
 * If m = k! for growing k, same thing.
 *
 * No improvement over known methods.
 */

/*
 * Let me try implementing the most concrete testable idea:
 * "Birthday collision factoring with multiple hash functions"
 *
 * Use K different hash functions h_1, ..., h_K on Z/NZ.
 * For each function, compute h_i(x) for random x and look for collisions.
 * If h_i maps Z/NZ → Z/pZ × Z/qZ, collisions in the p-projection
 * reveal p via gcd.
 *
 * Hash functions: h_i(x) = x^i mod N for small i.
 * Collision: h_i(x) = h_i(y) mod p iff x^i ≡ y^i mod p iff (x/y)^i ≡ 1 mod p.
 * This means x/y is an i-th root of unity mod p.
 *
 * The number of i-th roots of unity mod p is gcd(i, p-1).
 * If gcd(i, p-1) > 1 but gcd(i, q-1) = 1, then collisions in the p-projection
 * are more likely than in the q-projection.
 *
 * For i = 2: gcd(2, p-1) = 2 always. So x^2 ≡ y^2 mod p iff (x/y)^2 ≡ 1 mod p,
 * i.e., x ≡ ±y mod p. Two classes, so 2x more collisions.
 * Similarly for q. No asymmetry.
 *
 * For i = 3: gcd(3, p-1) is 3 if 3|(p-1), else 1.
 * If 3|(p-1) but 3∤(q-1): 3x more collisions in p-projection.
 * Then gcd(x^3 - y^3, N) = p for a collision pair.
 *
 * But we need collisions: with M random x values, expected collisions ≈ M^2/(2p).
 * Need M ≈ √p ≈ N^{1/4} for ≈1 collision. Same as Pollard rho.
 *
 * The i=3 insight doesn't help because both projections have similar collision rates.
 * The factor of 3 only helps when we already have a collision.
 */

/* Let me just implement a basic experiment: measure walk statistics
 * for different iteration functions and see if any exhibit non-trivial
 * structure that could be exploited. */

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N> [num_walks] [walk_length]\n", argv[0]);
        return 1;
    }

    mpz_t N;
    mpz_init(N);
    mpz_set_str(N, argv[1], 10);

    int num_walks = (argc >= 3) ? atoi(argv[2]) : 100;
    int walk_len = (argc >= 4) ? atoi(argv[3]) : 10000;

    gmp_randstate_t rng;
    gmp_randinit_mt(rng);
    gmp_randseed_ui(rng, 42);

    mpz_t x, y, g, tmp;
    mpz_init(x); mpz_init(y); mpz_init(g); mpz_init(tmp);

    struct timespec t0;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    /* Experiment 1: Multi-base GCD accumulation
     * For each base g, compute g^(M!) mod N for growing M.
     * Track when gcd(g^(M!) - 1, N) becomes nontrivial. */
    fprintf(stderr, "=== Experiment: Multi-base GCD accumulation ===\n");

    int found = 0;
    for (int w = 0; w < num_walks && !found; w++) {
        /* Random base */
        mpz_urandomm(g, rng, N);
        if (mpz_cmp_ui(g, 2) < 0) mpz_set_ui(g, 2);

        mpz_set(x, g);

        /* g^(2*3*4*5*...*k) = ((g^2)^3)^4)^5... */
        for (int k = 2; k <= walk_len && !found; k++) {
            mpz_powm_ui(x, x, k, N);

            /* Periodically check GCD */
            if (k % 100 == 0 || k <= 20) {
                mpz_sub_ui(tmp, x, 1);
                mpz_gcd(tmp, tmp, N);
                if (mpz_cmp_ui(tmp, 1) > 0 && mpz_cmp(tmp, N) < 0) {
                    mpz_t q;
                    mpz_init(q);
                    mpz_divexact(q, N, tmp);
                    gmp_printf("%Zd %Zd\n", tmp, q);
                    struct timespec t1;
                    clock_gettime(CLOCK_MONOTONIC, &t1);
                    double el = (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)/1e9;
                    fprintf(stderr, "Factor found! base=%d, k=%d, time=%.3fs\n", w, k, el);
                    fprintf(stderr, "  This means %s has a factor with %s-1 being %d-smooth\n",
                            argv[1], mpz_get_str(NULL, 10, tmp), k);
                    mpz_clear(q);
                    found = 1;
                }
            }
        }
    }

    if (!found) {
        fprintf(stderr, "p-1 method: No factor found (p-1 not %d-smooth)\n", walk_len);
    }

    /* Experiment 2: Pollard rho with multiple polynomials
     * Track cycle lengths across polynomials to detect structure */
    fprintf(stderr, "\n=== Experiment: Multi-polynomial Pollard rho ===\n");
    found = 0;

    for (int c_val = 1; c_val <= num_walks && !found; c_val++) {
        mpz_set_ui(x, 2);
        mpz_set_ui(y, 2);
        mpz_t c;
        mpz_init_set_ui(c, c_val);

        for (int step = 0; step < walk_len && !found; step++) {
            /* x = x^2 + c mod N (tortoise) */
            mpz_mul(x, x, x); mpz_add(x, x, c); mpz_mod(x, x, N);
            /* y = (y^2+c)^2+c mod N (hare) */
            mpz_mul(y, y, y); mpz_add(y, y, c); mpz_mod(y, y, N);
            mpz_mul(y, y, y); mpz_add(y, y, c); mpz_mod(y, y, N);

            mpz_sub(tmp, x, y);
            mpz_gcd(tmp, tmp, N);
            if (mpz_cmp_ui(tmp, 1) > 0 && mpz_cmp(tmp, N) < 0) {
                mpz_t q; mpz_init(q); mpz_divexact(q, N, tmp);
                gmp_printf("%Zd %Zd\n", tmp, q);
                struct timespec t1; clock_gettime(CLOCK_MONOTONIC, &t1);
                double el = (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)/1e9;
                fprintf(stderr, "Rho factor! c=%d, step=%d, time=%.3fs\n", c_val, step, el);
                mpz_clear(q);
                found = 1;
            }
        }
        mpz_clear(c);
    }

    if (!found) {
        fprintf(stderr, "Rho: No factor in %d walks of %d steps\n", num_walks, walk_len);
    }

    /* Experiment 3: ECM-like approach with accumulation
     * For multiple elliptic curves, accumulate scalar multiplication results
     * and periodically check GCD */
    fprintf(stderr, "\n=== Experiment: Multi-curve ECM accumulation ===\n");
    found = 0;

    mpz_t acc; mpz_init_set_ui(acc, 1);

    for (int curve = 0; curve < num_walks && !found; curve++) {
        /* Use ECM library directly via ecm_factor */
        /* For now, just use multiplicative group with different bases */
        mpz_urandomm(g, rng, N);
        if (mpz_cmp_ui(g, 2) < 0) mpz_set_ui(g, 2);

        /* Compute g^B! for B = small bound */
        mpz_set(x, g);
        int B1 = 1000;
        for (int k = 2; k <= B1; k++) {
            mpz_powm_ui(x, x, k, N);
        }

        /* Accumulate (x - 1) into product */
        mpz_sub_ui(tmp, x, 1);
        mpz_mul(acc, acc, tmp);
        mpz_mod(acc, acc, N);

        /* Check GCD periodically */
        if ((curve + 1) % 10 == 0) {
            mpz_gcd(tmp, acc, N);
            if (mpz_cmp_ui(tmp, 1) > 0 && mpz_cmp(tmp, N) < 0) {
                mpz_t q; mpz_init(q); mpz_divexact(q, N, tmp);
                gmp_printf("%Zd %Zd\n", tmp, q);
                struct timespec t1; clock_gettime(CLOCK_MONOTONIC, &t1);
                double el = (t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)/1e9;
                fprintf(stderr, "Accumulated factor! curves=%d, B1=%d, time=%.3fs\n",
                        curve+1, B1, el);
                mpz_clear(q);
                found = 1;
            }
        }
    }

    if (!found) {
        fprintf(stderr, "ECM-accum: No factor from %d curves with B1=1000\n", num_walks);
    }

    mpz_clear(acc);
    mpz_clear(x); mpz_clear(y); mpz_clear(g); mpz_clear(tmp);
    gmp_randclear(rng);
    mpz_clear(N);

    return found ? 0 : 1;
}
