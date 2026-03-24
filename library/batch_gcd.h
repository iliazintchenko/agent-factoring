/*
 * batch_gcd.h — Bernstein's product-tree batch GCD for smooth number detection.
 *
 * Given a set of candidates {c_1, ..., c_n} and a primorial P = prod(primes <= B),
 * computes gcd(c_i, P) for all i simultaneously in O(n log^2(nB)) time.
 *
 * This enables "sieve-free" smooth number detection on ARBITRARY candidate sets.
 */

#ifndef BATCH_GCD_H
#define BATCH_GCD_H

#include <gmp.h>
#include <stddef.h>

/* Product tree node */
typedef struct {
    mpz_t *nodes;   /* nodes[0..2*n-1], leaves at [n..2*n-1] */
    size_t n;        /* number of leaves (padded to power of 2) */
    size_t orig_n;   /* original number of candidates */
} product_tree_t;

/* Build product tree from leaves. n_orig candidates in candidates[]. */
void pt_build(product_tree_t *pt, const mpz_t *candidates, size_t n_orig);

/* Compute remainders: for each leaf i, result[i] = root mod candidates[i].
   root is typically the primorial P. */
void pt_remainders(const product_tree_t *pt, const mpz_t root, mpz_t *results);

/* Free product tree */
void pt_free(product_tree_t *pt);

/* Batch smooth detection:
   Given candidates[0..n-1] and smoothness bound B,
   sets smooth_part[i] = gcd(candidates[i], B#^k) where B# is the primorial
   and k is large enough to capture prime powers.
   Returns count of fully smooth candidates. */
size_t batch_smooth_detect(const mpz_t *candidates, size_t n, unsigned long B,
                           mpz_t *smooth_part);

/* Compute primorial: product of all primes <= B, with prime powers p^k <= B^2 */
void compute_primorial(mpz_t result, unsigned long B);

#endif /* BATCH_GCD_H */
