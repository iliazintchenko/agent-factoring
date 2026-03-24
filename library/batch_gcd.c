/*
 * batch_gcd.c — Bernstein product-tree batch GCD implementation.
 */

#include "batch_gcd.h"
#include <stdlib.h>
#include <string.h>

/* Round up to next power of 2 */
static size_t next_pow2(size_t n) {
    size_t p = 1;
    while (p < n) p <<= 1;
    return p;
}

void pt_build(product_tree_t *pt, const mpz_t *candidates, size_t n_orig) {
    size_t n = next_pow2(n_orig);
    pt->n = n;
    pt->orig_n = n_orig;
    pt->nodes = (mpz_t *)malloc(2 * n * sizeof(mpz_t));

    /* Initialize all nodes */
    for (size_t i = 0; i < 2 * n; i++)
        mpz_init(pt->nodes[i]);

    /* Set leaves */
    for (size_t i = 0; i < n_orig; i++)
        mpz_set(pt->nodes[n + i], candidates[i]);
    /* Pad with 1s */
    for (size_t i = n_orig; i < n; i++)
        mpz_set_ui(pt->nodes[n + i], 1);

    /* Build tree bottom-up: nodes[i] = nodes[2i] * nodes[2i+1] */
    for (size_t i = n - 1; i >= 1; i--)
        mpz_mul(pt->nodes[i], pt->nodes[2 * i], pt->nodes[2 * i + 1]);
}

void pt_remainders(const product_tree_t *pt, const mpz_t root, mpz_t *results) {
    size_t n = pt->n;

    /* Temporary array for remainder tree */
    mpz_t *rem = (mpz_t *)malloc(2 * n * sizeof(mpz_t));
    for (size_t i = 0; i < 2 * n; i++)
        mpz_init(rem[i]);

    /* Root remainder: root mod product_tree_root */
    mpz_mod(rem[1], root, pt->nodes[1]);

    /* Top-down: rem[i] = rem[parent] mod nodes[i] */
    for (size_t i = 2; i < 2 * n; i++) {
        if (mpz_sgn(pt->nodes[i]) > 0)
            mpz_mod(rem[i], rem[i / 2], pt->nodes[i]);
    }

    /* Extract leaf remainders */
    for (size_t i = 0; i < pt->orig_n; i++)
        mpz_set(results[i], rem[n + i]);

    for (size_t i = 0; i < 2 * n; i++)
        mpz_clear(rem[i]);
    free(rem);
}

void pt_free(product_tree_t *pt) {
    for (size_t i = 0; i < 2 * pt->n; i++)
        mpz_clear(pt->nodes[i]);
    free(pt->nodes);
    pt->nodes = NULL;
}

void compute_primorial(mpz_t result, unsigned long B) {
    /* Sieve primes up to B */
    unsigned char *is_prime = (unsigned char *)calloc(B + 1, 1);
    for (unsigned long i = 2; i <= B; i++) is_prime[i] = 1;
    for (unsigned long i = 2; i * i <= B; i++) {
        if (is_prime[i]) {
            for (unsigned long j = i * i; j <= B; j += i)
                is_prime[j] = 0;
        }
    }

    mpz_set_ui(result, 1);
    mpz_t pp;
    mpz_init(pp);

    /* For each prime p <= B, include p^k where p^k <= B^2 */
    unsigned long B2 = B * B;
    for (unsigned long p = 2; p <= B; p++) {
        if (!is_prime[p]) continue;
        unsigned long pk = p;
        while (pk <= B2 / p) pk *= p;  /* p^k <= B^2 */
        mpz_set_ui(pp, pk);
        mpz_mul(result, result, pp);
    }

    mpz_clear(pp);
    free(is_prime);
}

size_t batch_smooth_detect(const mpz_t *candidates, size_t n, unsigned long B,
                           mpz_t *smooth_part) {
    if (n == 0) return 0;

    /* Compute primorial */
    mpz_t P;
    mpz_init(P);
    compute_primorial(P, B);

    /* Build product tree */
    product_tree_t pt;
    pt_build(&pt, candidates, n);

    /* Compute remainders: rem[i] = P mod candidates[i] */
    mpz_t *rem = (mpz_t *)malloc(n * sizeof(mpz_t));
    for (size_t i = 0; i < n; i++) mpz_init(rem[i]);
    pt_remainders(&pt, P, rem);

    /* smooth_part[i] = gcd(candidates[i], rem[i]) = gcd(candidates[i], P) */
    size_t count = 0;
    for (size_t i = 0; i < n; i++) {
        mpz_gcd(smooth_part[i], candidates[i], rem[i]);
        /* Check if fully smooth: candidates[i] / smooth_part[i] == 1 */
        mpz_t cofactor;
        mpz_init(cofactor);
        mpz_divexact(cofactor, candidates[i], smooth_part[i]);
        /* Repeatedly extract smooth part */
        int changed = 1;
        while (changed) {
            changed = 0;
            mpz_gcd(smooth_part[i], cofactor, P);
            if (mpz_cmp_ui(smooth_part[i], 1) > 0) {
                mpz_divexact(cofactor, cofactor, smooth_part[i]);
                changed = 1;
            }
        }
        if (mpz_cmp_ui(cofactor, 1) == 0) count++;
        /* Store the actual cofactor for caller to inspect */
        mpz_set(smooth_part[i], cofactor);
        mpz_clear(cofactor);
    }

    for (size_t i = 0; i < n; i++) mpz_clear(rem[i]);
    free(rem);
    pt_free(&pt);
    mpz_clear(P);

    return count;
}
