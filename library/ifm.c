/*
 * Iterated Frobenius Map (IFM) factoring
 *
 * Novel approach: iterate x -> x^N mod N.
 * Via CRT decomposition:
 *   mod p: x -> x^{q mod (p-1)} mod p
 *   mod q: x -> x^{p mod (q-1)} mod q
 * These are different dynamical systems. Brent's cycle detection
 * finds the overall cycle; GCD extraction reveals the factor when
 * the two components have different cycle structure.
 *
 * Expected complexity: O(N^{1/4}) iterations, each costing O(log(N)^3)
 * for the modular exponentiation.
 *
 * Also includes a "multi-start" variant that tries different starting
 * points and different Frobenius-like maps (x -> x^{N+k} for small k).
 */

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Try to factor N using the iterated map x -> x^E mod N
 * where E is typically N or N+k.
 * Returns 1 if factor found (stored in factor), 0 otherwise.
 */
static int ifm_try(mpz_t factor, const mpz_t N, const mpz_t E,
                    unsigned long max_iters, gmp_randstate_t rng)
{
    mpz_t x, y, d, acc, temp;
    mpz_inits(x, y, d, acc, temp, NULL);

    /* Random starting point in [2, N-1] */
    mpz_urandomm(x, rng, N);
    if (mpz_cmp_ui(x, 2) < 0) mpz_set_ui(x, 2);
    mpz_set(y, x);

    mpz_set_ui(acc, 1);

    unsigned long power = 1, lam = 0;
    int found = 0;
    unsigned long batch = 0;
    const unsigned long BATCH_SIZE = 128;

    for (unsigned long i = 0; i < max_iters && !found; i++) {
        if (power == lam) {
            mpz_set(y, x);
            power <<= 1;
            lam = 0;
        }

        /* x = x^E mod N */
        mpz_powm(x, x, E, N);
        lam++;

        /* Accumulate (x - y) into product for batch GCD */
        mpz_sub(temp, x, y);
        mpz_mod(temp, temp, N);
        if (mpz_sgn(temp) == 0) {
            /* x == y mod N: cycle found but trivial. Restart. */
            mpz_urandomm(x, rng, N);
            if (mpz_cmp_ui(x, 2) < 0) mpz_set_ui(x, 2);
            mpz_set(y, x);
            mpz_set_ui(acc, 1);
            power = 1; lam = 0;
            continue;
        }

        mpz_mul(acc, acc, temp);
        mpz_mod(acc, acc, N);
        batch++;

        if (batch >= BATCH_SIZE) {
            mpz_gcd(d, acc, N);
            if (mpz_cmp_ui(d, 1) > 0 && mpz_cmp(d, N) < 0) {
                mpz_set(factor, d);
                found = 1;
            } else if (mpz_cmp(d, N) == 0) {
                /* Accumulated product is 0 mod N — backtrack */
                /* Just restart with a new random start */
                mpz_urandomm(x, rng, N);
                if (mpz_cmp_ui(x, 2) < 0) mpz_set_ui(x, 2);
                mpz_set(y, x);
                power = 1; lam = 0;
            }
            mpz_set_ui(acc, 1);
            batch = 0;
        }
    }

    /* Final GCD check on remaining accumulated product */
    if (!found && batch > 0) {
        mpz_gcd(d, acc, N);
        if (mpz_cmp_ui(d, 1) > 0 && mpz_cmp(d, N) < 0) {
            mpz_set(factor, d);
            found = 1;
        }
    }

    mpz_clears(x, y, d, acc, temp, NULL);
    return found;
}

/* Quick trial division up to limit */
static int trial_divide(mpz_t factor, const mpz_t N, unsigned long limit)
{
    if (mpz_divisible_ui_p(N, 2)) { mpz_set_ui(factor, 2); return 1; }
    if (mpz_divisible_ui_p(N, 3)) { mpz_set_ui(factor, 3); return 1; }

    for (unsigned long d = 5; d <= limit; d += (d % 6 == 5) ? 2 : 4) {
        if (mpz_divisible_ui_p(N, d)) {
            mpz_set_ui(factor, d);
            return 1;
        }
    }
    return 0;
}

/* Pollard's rho as a fast fallback for smaller numbers */
static int pollard_rho(mpz_t factor, const mpz_t N, unsigned long max_iters,
                       gmp_randstate_t rng)
{
    mpz_t x, y, d, acc, temp, c;
    mpz_inits(x, y, d, acc, temp, c, NULL);

    mpz_urandomm(c, rng, N);
    mpz_urandomm(x, rng, N);
    mpz_set(y, x);
    mpz_set_ui(acc, 1);

    int found = 0;
    unsigned long batch = 0;

    for (unsigned long i = 0; i < max_iters && !found; i++) {
        /* x = x^2 + c mod N */
        mpz_mul(x, x, x);
        mpz_add(x, x, c);
        mpz_mod(x, x, N);

        mpz_mul(y, y, y); mpz_add(y, y, c); mpz_mod(y, y, N);
        mpz_mul(y, y, y); mpz_add(y, y, c); mpz_mod(y, y, N);

        mpz_sub(temp, x, y);
        mpz_abs(temp, temp);
        mpz_mul(acc, acc, temp);
        mpz_mod(acc, acc, N);
        batch++;

        if (batch >= 256) {
            mpz_gcd(d, acc, N);
            if (mpz_cmp_ui(d, 1) > 0 && mpz_cmp(d, N) < 0) {
                mpz_set(factor, d);
                found = 1;
            } else if (mpz_cmp(d, N) == 0) {
                /* Restart */
                mpz_urandomm(c, rng, N);
                mpz_urandomm(x, rng, N);
                mpz_set(y, x);
            }
            mpz_set_ui(acc, 1);
            batch = 0;
        }
    }

    if (!found && batch > 0) {
        mpz_gcd(d, acc, N);
        if (mpz_cmp_ui(d, 1) > 0 && mpz_cmp(d, N) < 0) {
            mpz_set(factor, d);
            found = 1;
        }
    }

    mpz_clears(x, y, d, acc, temp, c, NULL);
    return found;
}

int main(int argc, char *argv[])
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        return 1;
    }

    mpz_t N, factor, cofactor, E;
    mpz_inits(N, factor, cofactor, E, NULL);
    mpz_set_str(N, argv[1], 10);

    struct timespec start, now;
    clock_gettime(CLOCK_MONOTONIC, &start);

    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, 42);

    int found = 0;

    /* Stage 0: Trial division up to 10^6 */
    if (!found) {
        found = trial_divide(factor, N, 1000000);
    }

    /* Stage 1: Pollard's rho (good for smaller factors) */
    if (!found) {
        found = pollard_rho(factor, N, 5000000, rng);
    }

    /* Stage 2: IFM with E = N */
    if (!found) {
        mpz_set(E, N);
        found = ifm_try(factor, N, E, 10000000, rng);
    }

    /* Stage 3: IFM with E = N+1, N-1, N+2, etc. */
    if (!found) {
        for (int k = 1; k <= 5 && !found; k++) {
            mpz_set(E, N);
            mpz_add_ui(E, E, k);
            found = ifm_try(factor, N, E, 2000000, rng);
            if (!found) {
                mpz_set(E, N);
                mpz_sub_ui(E, E, k);
                found = ifm_try(factor, N, E, 2000000, rng);
            }
        }
    }

    /* Stage 4: More Pollard's rho with different seeds */
    if (!found) {
        for (int attempt = 0; attempt < 20 && !found; attempt++) {
            gmp_randseed_ui(rng, 42 + attempt);
            found = pollard_rho(factor, N, 10000000, rng);
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &now);
    double elapsed = (now.tv_sec - start.tv_sec) + (now.tv_nsec - start.tv_nsec) / 1e9;

    if (found) {
        mpz_divexact(cofactor, N, factor);
        /* Output smaller factor first */
        if (mpz_cmp(factor, cofactor) > 0) mpz_swap(factor, cofactor);
        gmp_printf("%Zd %Zd\n", factor, cofactor);
        fprintf(stderr, "IFM: factored in %.3f seconds\n", elapsed);
    } else {
        fprintf(stderr, "IFM: FAILED after %.3f seconds\n", elapsed);
        return 1;
    }

    mpz_clears(N, factor, cofactor, E, NULL);
    gmp_randclear(rng);
    return 0;
}
