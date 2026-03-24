/*
 * structured_walk.c — Factoring via structured random walks in Z_N*
 *
 * NOVEL APPROACH: Instead of Pollard's rho (random walk with x → x²+c),
 * use a STRUCTURED walk that exploits the product structure Z_N* ≅ Z_{p-1} × Z_{q-1}.
 *
 * Key insight: a random walk on Z_N* is simultaneously a walk on Z_{p-1}
 * and Z_{q-1}. If the walks have DIFFERENT cycle lengths mod p vs mod q,
 * detecting this difference reveals a factor.
 *
 * Standard Pollard rho detects this via birthday collisions (O(N^{1/4})).
 * Can we do better with a STRUCTURED walk?
 *
 * IDEA 1: Multi-degree walk
 * Instead of x → x²+c (degree 2), use x → x^k for carefully chosen k.
 * The iteration x → x^k has period dividing (p-1)/gcd(k-1, p-1) mod p.
 * If we choose k such that gcd(k-1, p-1) is large but gcd(k-1, q-1)
 * is small, the periods differ maximally, enabling faster detection.
 * Problem: we don't know p-1 or q-1.
 *
 * IDEA 2: Smooth-exponent accumulation
 * Compute x → x^(product of small primes) iteratively.
 * After accumulating primes up to B, if p-1 is B-smooth, x^{lcm(1..B)} ≡ 1 (mod p).
 * This is Pollard p-1. For balanced semiprimes, p-1 is not B-smooth.
 *
 * IDEA 3: Multi-walk correlation (NOVEL)
 * Run K independent walks simultaneously:
 *   w_i(t+1) = f_i(w_i(t)) for i = 1..K
 * where f_i are different step functions.
 *
 * At each step, compute GCD of all pairwise PRODUCTS:
 *   gcd(w_i(t) * w_j(t) - w_k(t) * w_l(t), N)
 *
 * The product w_i * w_j mod p and mod q follows different distributions.
 * If any such product equals another mod p but not mod q, gcd reveals p.
 *
 * Expected collision time: with K walks, there are O(K^4) product pairs,
 * giving collision probability K^4/N per step. After T steps: K^4*T/N.
 * For K^4*T = N: T = N/K^4. Total work: K*T = N/K^3.
 * Minimizing over K: K should be as large as feasible (memory-limited).
 * With K = N^{1/8}: T = N^{1/2}, total = N^{3/8}. But K = N^{1/8} walks
 * requires N^{1/8} memory, which is exponential.
 *
 * For practical K (constant): total work = O(N/K^3) = O(N^{1/4}) × constant.
 * Same as Pollard rho. No asymptotic improvement.
 *
 * IDEA 4: Structured walk with prime-power steps (TESTING THIS)
 * Walk: x → x * g^{random_small_prime} where g is a fixed generator.
 * After T steps: x_T = x_0 * g^{sum_of_small_primes}.
 * The sum of T random small primes has a specific distribution.
 * If this distribution concentrates on values divisible by (p-1) or (q-1)
 * with higher probability than expected, we can detect period collisions faster.
 *
 * The "structured" part: instead of purely random steps, use steps
 * that are SMOOTH numbers. This biases the walk toward elements
 * whose order divides smooth numbers — which are exactly the elements
 * that would reveal factors.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

/* Multi-walk with product correlation */
int factor_multiwalk(mpz_t N, mpz_t factor, int K, long max_steps) {
    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, 42);

    /* Initialize K walks */
    mpz_t *w = malloc(K * sizeof(mpz_t));
    mpz_t *c = malloc(K * sizeof(mpz_t)); /* constants for each walk */

    for (int i = 0; i < K; i++) {
        mpz_init(w[i]);
        mpz_init(c[i]);
        mpz_urandomm(w[i], rng, N);
        if (mpz_sgn(w[i]) == 0) mpz_set_ui(w[i], 2);
        mpz_set_ui(c[i], i + 1);
    }

    mpz_t prod, diff, g;
    mpz_inits(prod, diff, g, NULL);

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    int found = 0;
    for (long step = 0; step < max_steps && !found; step++) {
        /* Advance each walk: w_i → w_i^2 + c_i mod N */
        for (int i = 0; i < K; i++) {
            mpz_mul(w[i], w[i], w[i]);
            mpz_add(w[i], w[i], c[i]);
            mpz_mod(w[i], w[i], N);
        }

        /* Check pairwise products for collisions */
        /* For efficiency, only check a subset of pairs */
        for (int i = 0; i < K && !found; i++) {
            for (int j = i + 1; j < K && !found; j++) {
                mpz_sub(diff, w[i], w[j]);
                mpz_gcd(g, diff, N);
                if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                    mpz_set(factor, g);
                    found = 1;
                }
            }
        }

        /* Also try product-based collisions every 100 steps */
        if (step % 100 == 0 && K >= 4 && !found) {
            for (int i = 0; i < K && !found; i++) {
                for (int j = i + 1; j < K && !found; j++) {
                    mpz_mul(prod, w[i], w[j]);
                    mpz_mod(prod, prod, N);
                    for (int k = j + 1; k < K && !found; k++) {
                        for (int l = k + 1; l < K && !found; l++) {
                            mpz_t prod2;
                            mpz_init(prod2);
                            mpz_mul(prod2, w[k], w[l]);
                            mpz_mod(prod2, prod2, N);
                            mpz_sub(diff, prod, prod2);
                            mpz_gcd(g, diff, N);
                            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                                mpz_set(factor, g);
                                found = 1;
                            }
                            mpz_clear(prod2);
                        }
                    }
                }
            }
        }

        if (step % 1000000 == 0 && step > 0) {
            clock_gettime(CLOCK_MONOTONIC, &t1);
            double elapsed = (t1.tv_sec - t0.tv_sec) +
                           (t1.tv_nsec - t0.tv_nsec) / 1e9;
            fprintf(stderr, "MultiWalk: %ld steps, %.1fs\n", step, elapsed);
        }
    }

    if (found) {
        clock_gettime(CLOCK_MONOTONIC, &t1);
        double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
        fprintf(stderr, "MultiWalk: Factor found at step %ld in %.3fs\n",
                max_steps, elapsed);
    }

    for (int i = 0; i < K; i++) { mpz_clear(w[i]); mpz_clear(c[i]); }
    free(w); free(c);
    mpz_clears(prod, diff, g, NULL);
    gmp_randclear(rng);

    return found;
}

/* Smooth-step walk: accumulate smooth exponents */
int factor_smooth_walk(mpz_t N, mpz_t factor, long max_steps) {
    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, 42);

    /* Small primes for step generation */
    unsigned long small_primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
                                     31, 37, 41, 43, 47, 53, 59, 61, 67, 71};
    int n_primes = sizeof(small_primes) / sizeof(small_primes[0]);

    /* Multiple walkers with different generators */
    int K = 8;
    mpz_t *x = malloc(K * sizeof(mpz_t));
    mpz_t *g_base = malloc(K * sizeof(mpz_t));

    for (int i = 0; i < K; i++) {
        mpz_init(x[i]);
        mpz_init(g_base[i]);
        mpz_set_ui(g_base[i], small_primes[i]);
        mpz_set_ui(x[i], 2 + i);
    }

    mpz_t g, diff, exp_val;
    mpz_inits(g, diff, exp_val, NULL);

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    int found = 0;
    for (long step = 0; step < max_steps && !found; step++) {
        /* Each walker takes a random smooth step */
        for (int i = 0; i < K; i++) {
            /* Pick a random small prime and multiply x by g^p */
            int pidx = step % n_primes;
            mpz_set_ui(exp_val, small_primes[(pidx + i) % n_primes]);
            mpz_powm(x[i], x[i], exp_val, N);
        }

        /* Check: is any x_i ≡ 1 (mod N)? (Would mean order divides accumulated exponent) */
        for (int i = 0; i < K; i++) {
            mpz_sub_ui(diff, x[i], 1);
            mpz_gcd(g, diff, N);
            if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                mpz_set(factor, g);
                found = 1;
                break;
            }
        }

        /* Also check pairwise differences */
        if (!found) {
            for (int i = 0; i < K && !found; i++) {
                for (int j = i + 1; j < K && !found; j++) {
                    mpz_sub(diff, x[i], x[j]);
                    mpz_gcd(g, diff, N);
                    if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, N) < 0) {
                        mpz_set(factor, g);
                        found = 1;
                    }
                }
            }
        }

        if (step % 100000 == 0 && step > 0) {
            clock_gettime(CLOCK_MONOTONIC, &t1);
            double elapsed = (t1.tv_sec - t0.tv_sec) +
                           (t1.tv_nsec - t0.tv_nsec) / 1e9;
            fprintf(stderr, "SmoothWalk: %ld steps, %.1fs\n", step, elapsed);
        }
    }

    if (found) {
        clock_gettime(CLOCK_MONOTONIC, &t1);
        double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;
        fprintf(stderr, "SmoothWalk: Factor found in %.3fs\n", elapsed);
    }

    for (int i = 0; i < K; i++) { mpz_clear(x[i]); mpz_clear(g_base[i]); }
    free(x); free(g_base);
    mpz_clears(g, diff, exp_val, NULL);
    gmp_randclear(rng);

    return found;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N> [mode: multi|smooth] [max_steps]\n", argv[0]);
        return 1;
    }

    mpz_t N, f;
    mpz_inits(N, f, NULL);
    mpz_set_str(N, argv[1], 10);

    const char *mode = argc > 2 ? argv[2] : "multi";
    long max_steps = argc > 3 ? atol(argv[3]) : 10000000;

    int ok = 0;
    if (strcmp(mode, "smooth") == 0) {
        ok = factor_smooth_walk(N, f, max_steps);
    } else {
        int K = 16;
        ok = factor_multiwalk(N, f, K, max_steps);
    }

    if (ok) {
        mpz_t c; mpz_init(c);
        mpz_divexact(c, N, f);
        if (mpz_cmp(f, c) > 0) mpz_swap(f, c);
        gmp_printf("%Zd %Zd\n", f, c);
        mpz_clear(c);
    } else {
        fprintf(stderr, "FAIL\n");
    }

    mpz_clears(N, f, NULL);
    return ok ? 0 : 1;
}
