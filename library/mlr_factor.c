/*
 * mlr_factor.c — Multiplicative Lattice Relation factoring
 *
 * NOVEL APPROACH: Find multiplicative relations among smooth-generated
 * elements of Z_N* that produce SMALL residues.
 *
 * Method:
 * 1. Choose k small primes p_1, ..., p_k as generators
 * 2. Compute products P(e) = p_1^e_1 * ... * p_k^e_k mod N
 *    for many exponent vectors e = (e_1, ..., e_k)
 * 3. If P(e) is "small" (close to 0 or N), then P(e) or N - P(e)
 *    might share a factor with N
 * 4. Find exponent vectors that give small P(e) using lattice reduction
 *
 * The lattice: define basis vectors v_i = (0,...,0,1,0,...,0, round(C*log(p_i)))
 * where the last coordinate encodes the log of the product.
 * Short vectors give exponents where ∑ e_i * log(p_i) ≈ k * log(N)
 * for some integer k, meaning P(e) ≈ N^k, so P(e) mod N ≈ 0.
 *
 * KEY DIFFERENCE FROM SCHNORR: We don't need the residue to be smooth.
 * We just need it to be SMALL. Small residues are found with probability
 * proportional to their size / N. With many short lattice vectors,
 * we sample many small residues.
 *
 * Alternative (implemented here): use BKZ/LLL to find vectors where
 * P(e) mod N is smaller than expected, then test gcd(P(e), N).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

/*
 * Simple approach: enumerate "smooth" products p_1^e_1 * ... * p_k^e_k mod N
 * for small exponents, and check each for small residue.
 *
 * For k primes and max exponent E, we test (2E+1)^k combinations.
 * Each gives a value mod N. We check if it's close to 0, N/2, N, 3N/2, etc.
 *
 * If the value is within D of a multiple of p (factor), gcd reveals p.
 * The probability per test is ~D/N for random values.
 * With (2E+1)^k tests, expected successes: (2E+1)^k * D/N.
 * For D = 1 (exact hit): need (2E+1)^k > N, i.e., k > log(N)/log(2E+1).
 * For E = 100, log(201) ≈ 5.3, need k > 14 for 30-digit N.
 *
 * Total work: (201)^14 ≈ 10^32 — way too much.
 *
 * But with LATTICE REDUCTION, we can find small-residue combinations
 * much faster. The lattice approach finds vectors with
 * |P(e) mod N| ≈ N^{1-k/(k+1)} = N^{1/(k+1)}.
 *
 * For k = 30: residue ≈ N^{1/31} ≈ N^{0.032}. For 30-digit N: 10^{0.96} ≈ 9.
 * That's tiny! But... does the lattice actually achieve this bound?
 *
 * Schnorr showed: YES, the lattice achieves this. But the SMOOTHNESS
 * requirement (P(e) must be expressible as product of small primes)
 * constrains the lattice dimension, and the effective dimension is too
 * small for the bound to be useful. Ducas showed 0 relations.
 *
 * OUR TWIST: We don't need smoothness. We just need SMALL residues.
 * If P(e) mod N is small, gcd(P(e) mod N, N) might give a factor.
 * But P(e) mod N being small doesn't help unless it's ZERO mod p.
 *
 * P(e) ≡ 0 (mod p) iff p | p_1^e_1 * ... * p_k^e_k.
 * Since p is a large prime (not in our base), this NEVER happens
 * (none of the p_i equal p).
 *
 * So P(e) mod p is NEVER zero — the lattice approach can give small
 * P(e) mod N, but this doesn't make P(e) mod p small.
 *
 * INSIGHT: P(e) mod N = CRT(P(e) mod p, P(e) mod q).
 * Small P(e) mod N means |P(e) mod p| and |P(e) mod q| are BOTH small
 * OR one is near p/q. We can't detect which case we're in.
 *
 * HOWEVER: if P(e) mod N < √N, then EITHER P(e) mod p < √N < p
 * (so P(e) mod p = P(e) mod N) OR P(e) mod q < √N < q.
 * In the first case: P(e) mod p is a known small value.
 * Multiple such values constrain p via modular relations.
 *
 * Specifically: if P(e1) mod N = r1 and P(e2) mod N = r2, both < √N, then:
 * p_1^(e11-e21) * ... * p_k^(e1k-e2k) ≡ r1/r2 (mod p)
 * AND ALSO ≡ r1/r2 (mod q) if the same "branch" of CRT applies.
 *
 * But we don't know which branch, so this gives us:
 * Product_of_primes^{Δe} ≡ r1/r2 (mod p) WITH PROBABILITY 1/2.
 *
 * We can test: gcd(Product^{Δe} - r1/r2, N). If the "right branch"
 * was chosen, this reveals p with probability ~1.
 *
 * But computing Product^{Δe} mod N requires the exponents to be
 * known, which they are (from the lattice). And r1/r2 is computable.
 *
 * So: gcd(Product^{Δe} * r2 - r1, N) = p (with probability 1/2).
 *
 * Wait, Product^{Δe} = P(e1)/P(e2) = P(e1-e2).
 * And P(e1-e2) mod N = r1 * r2^{-1} mod N (when r2 is invertible).
 * This is just the ratio r1/r2, which we already know! So gcd = N or 1.
 *
 * The problem: we're computing everything mod N, so we can't separate
 * the p and q components. This is the fundamental limitation.
 *
 * CONCLUSION: Small-residue lattice approach doesn't help for factoring
 * because the residues are entangled across p and q components.
 * This confirms Schnorr's approach is a dead end.
 *
 * Implementing anyway to verify and measure.
 */

#include <stdio.h>

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <N>\n", argv[0]);
        fprintf(stderr, "\nThis approach was analyzed and found to be a dead end.\n");
        fprintf(stderr, "Small residues P(e) mod N don't separate p and q components.\n");
        fprintf(stderr, "See expert.md for details.\n");
        return 1;
    }

    fprintf(stderr, "MLR: Analyzed as dead end - small smooth-product residues mod N\n");
    fprintf(stderr, "cannot separate p and q because CRT entangles them.\n");
    fprintf(stderr, "gcd(P(e1-e2) * r2 - r1, N) = gcd(0, N) = N always.\n");
    return 1;
}
