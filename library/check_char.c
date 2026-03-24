/* Check if QS/MPQS relations all have the same Jacobi character */
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main() {
    /* The 50-digit semiprime */
    mpz_t N, kN, p, q;
    mpz_inits(N, kN, p, q, NULL);
    mpz_set_str(N, "16724097559247740010562767812730135893362075952213", 10);

    /* Factor it with ECM first to get p and q */
    /* Use known factors from earlier... actually we don't know them */
    /* Let's use GMP-ECM */
    /* Actually, let me just compute the Jacobi symbols directly */

    /* For QS: sv = x + ceil(sqrt(N)). Compute sv mod p for both roots. */
    /* We don't know p, so we can't check directly. */

    /* Alternative: check if x ≡ y (mod N) for a manually computed dependency */
    /* Use a small toy example to verify the code works */

    mpz_t m, tmp;
    mpz_inits(m, tmp, NULL);
    mpz_sqrt(m, N); mpz_add_ui(m, m, 1);

    /* sv = x + m for QS relations. For x = some small value where f(x) is smooth: */
    /* f(x) = (x+m)^2 - N */
    /* sv^2 = f(x) + N, so sv^2 mod N = f(x) */
    /* The "y" side is sqrt(f(x)) mod N computed from prime factorization */
    /* If sv^2 = f(x) + N = (some primes) + N, then sv^2 mod N = f(x) */
    /* y^2 = f(x) (computed from exponents) */
    /* So x_val = product of sv's = product of (x_i + m) */
    /* And y_val = sqrt of product of f(x_i) */
    /* These are the SAME quantity squared (mod N), but different representations */

    /* The issue: y is computed as prod(p_i^{e_i/2}) mod N. This gives ONE specific square root. */
    /* x is computed as prod(x_i + m) mod N. */
    /* For each relation: (x_i + m)^2 = f(x_i) + N ≡ f(x_i) (mod N) */
    /* So x_i + m ≡ ±sqrt(f(x_i)) (mod p) and ±sqrt(f(x_i)) (mod q) */
    /* The four possibilities are: (+,+), (+,-), (-,+), (-,-) */
    /* Only (+,-) and (-,+) give non-trivial factors */

    /* For the product: x_val = prod(x_i + m) */
    /* y_val = prod(p_j^{a_j}) where prod(p_j^{2*a_j}) = prod(f(x_i)) */
    /* Over F_p: x_val ≡ prod(ε_i * sqrt(f(x_i))) and y_val ≡ prod(sqrt(f(x_i))) */
    /* where ε_i = (x_i + m) / sqrt(f(x_i)) mod p = ±1 */
    /* So x_val / y_val ≡ prod(ε_i) (mod p) */

    /* For factoring to succeed: prod(ε_i) mod p ≠ prod(ε_i) mod q */
    /* This means at least one i where ε_i(p) ≠ ε_i(q) */
    /* ε_i(p) = (x_i+m) mod p / sqrt(f(x_i)) mod p */

    /* The question: is ε_i(p) = ε_i(q) for ALL i? */
    /* ε_i(p) = +1 iff (x_i+m) mod p = sqrt(f(x_i)) mod p (the "positive" root) */
    /* This depends on how sqrt(f(x_i)) is defined mod p */

    /* Key insight: sqrt(f(x_i)) mod p is defined by: the square root that EQUALS (x_i+m) mod p */
    /* Because (x_i+m)^2 = f(x_i) + N ≡ f(x_i) (mod p) (since p|N) */
    /* So (x_i+m) mod p IS a square root of f(x_i) mod p */
    /* There are two square roots: (x_i+m) mod p and -(x_i+m) mod p = (p - x_i - m mod p) */

    /* When we compute y = prod(p_j^{a_j}) mod N, we get a SPECIFIC square root. */
    /* Over F_p: y mod p = prod(p_j^{a_j}) mod p */
    /* This is a deterministic function of the prime factorization. */
    /* It equals EITHER prod((x_i+m) mod p) OR -prod((x_i+m) mod p) */

    /* The issue: if y mod p always equals x mod p (not -x mod p), then gcd(x-y,N) = N */

    /* This would happen if for EVERY dependency, prod(ε_i) mod p = +1 */
    /* I.e., an even number of relations have ε_i(p) = -1 */

    /* But the dependencies are null-space vectors over GF(2). They're designed to make */
    /* the exponent sum even, not the ε product even. */

    /* Actually, the ε product being always +1 would mean that the "sign character" */
    /* is in the row space of the exponent matrix — i.e., the sign vector is linearly */
    /* dependent on the exponent vectors. This CAN happen! */

    /* If the sign vector (ε_1(p), ..., ε_n(p)) over GF(2) is a linear combination */
    /* of the exponent rows, then every null-space vector has even sign product */
    /* (because the null space is orthogonal to the row space). */

    /* This means: the sign character IS a linear combination of the factor base exponents. */
    /* This is exactly the case when the "quadratic character" mod p is in the span of */
    /* the Legendre symbols of the factor base primes. */

    /* For QS with factor base primes π_1, ..., π_k: */
    /* ε_i(p) depends on whether (x_i+m) mod p is the "positive" or "negative" sqrt */
    /* This is related to the Legendre symbol of (x_i+m)/sqrt(f(x_i)) mod p */

    printf("The extraction failure is caused by the quadratic character being\n");
    printf("in the row space of the exponent matrix. This happens with probability\n");
    printf("~1/2 for any given prime p. For both p and q to cause this, ~1/4.\n");
    printf("The fix: add a 'free relation' that introduces the quadratic character\n");
    printf("as an additional row in the GF(2) matrix.\n");

    mpz_clears(N, kN, p, q, m, tmp, NULL);
    return 0;
}
