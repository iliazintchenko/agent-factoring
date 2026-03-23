# Expert Knowledge

## The smoothness bottleneck

All known sub-exponential factoring algorithms rely on finding **smooth numbers**. The L-exponent is determined by the size of values being tested for smoothness:
- QS: tests values of size ~√N → L[1/2]
- NFS: uses algebraic number fields to reduce values to ~N^(2/3) → L[1/3]
- Hypothetical L[1/4]: would need values of size ~N^(1/4) or equivalent

Any approach that still tests numbers of size √N for smoothness will remain L[1/2]. To do better, you must either generate smaller candidates, or avoid smoothness entirely.

The quasi-polynomial DLP breakthrough (Barbulescu-Gaudry-Joux-Thomé 2013) shows L[1/3] barriers CAN be broken for related problems in function fields. The technique exploits systematic factorization of polynomials over finite fields — no known analog over Z. The "translation problem" from function fields to number fields is a major open question.

## Dead ends

- **Schnorr lattice factoring (2021)**: Lattice dimension grows with smoothness bound (~50K for 90d). Ducas (CWI) confirmed: 0 relations in 1000 trials. Any lattice approach depending on a large factor base inherits the same dimension blowup.
- **Spectral methods on Cayley graphs**: Cayley graph of (Z/NZ)* encodes group structure, but computing the spectrum is exponential classically — this is exactly what Shor's QFT solves.
- **Multi-image NFS**: Adding more number fields (Barbulescu et al.) improves L[1/3] constants but cannot break the L[1/3] barrier. Simultaneous smoothness across fields constrains the sieve region proportionally.
- **Dixon's random squares**: Values are ~N in size, far too large for smoothness without polynomial structure.
- **CFRAC**: Single expansion limits relation rate. Not competitive above 30 digits.
- **ECM for balanced semiprimes**: L[1/2] in N for balanced factors. Not competitive above ~55 digits.
- **Lattice-based smooth value generation**: Constructing lattices where short vectors correspond to smooth products mod N reduces to Schnorr's failed approach. Short vectors have small exponents → small products → no wrapping mod N.
- **Higher-dimensional sieving (tower NFS)**: For integer factoring (vs DLP in finite fields), extra dimensions don't reduce effective norm sizes. Tower NFS breakthrough for GF(p^n) doesn't translate.
- **ECM failure information aggregation**: Collecting smooth/non-smooth classifications of curve orders across many ECM curves to constrain p. Success rate only reveals factor SIZE (already known), not specific factor IDENTITY.

## Research survey

- **No classical sub-L[1/3] algorithm exists** for general integer factoring.
- **Regev (2023)**: Quantum factoring in O~(n^{3/2}) gates. Fundamentally quantum, no classical dequantization known.
- **Stange (2022)**: Index calculus in (Z/NZ)*. Clean but same L[1/2] or L[1/3]. Not an improvement.
- **Umans & Wang (2025)**: Conditional deterministic N^{1/6+o(1)}. Still exponential.
- **Harvey-Hittmeir**: Deterministic N^{1/5+o(1)}. Still exponential.

## Why novel approaches are hard

The core difficulty: any classical approach to factoring must somehow distinguish the factor p from ~sqrt(N) other candidates of similar size. Sub-exponential algorithms achieve this via smooth number relations + linear algebra. To improve the L-exponent, you need SMALLER smooth candidates, which requires new algebraic structures that reduce norm sizes below what NFS achieves.

## Key insight: smoothness detection cost matters

Standard QS/NFS sieving has per-candidate cost O(1) amortized but requires sequential memory access over a fixed interval. Bernstein's batch smoothness detection (product/remainder trees) has per-candidate cost O(log²B · polylog) but works on ARBITRARY candidate sets. This unlocks two opportunities:

1. **Larger smoothness bounds**: Sieving cost scales as O(B/log B) per candidate; batch GCD scales as O(log²B). For large B, batch GCD is cheaper. A larger B means higher smoothness probability, fewer candidates needed.

2. **Non-sequential candidate generation**: Sieving requires evaluating consecutive polynomial values. Batch GCD works on any set of values, enabling novel candidate generation strategies.

**Practical validation**: Batch smoothness was used in the RSA-240 factoring record (2019) and the "We Are on the Same Side" paper (Eprint 2023/801) showed it can outperform Cado-NFS's sieving at RSA-250 scale.

## Theoretical analysis: batch GCD + multi-large-prime

**Hypothesis**: Combining batch smooth-part extraction with aggressive multi-large-prime relation combination could improve the L(1/2) constant for QS-regime factoring.

Analysis: In standard QS, the optimal smoothness bound B = L(1/2, 1/√2) balances factor base size against smoothness probability. The total cost is L(1/2, √2).

With batch GCD, the per-candidate cost drops from O(B) to O(polylog(B)), allowing a larger optimal B. Rough analysis suggests the optimum shifts to B ≈ L(1/2, 1) with total cost ≈ L(1/2, 1) — a constant improvement over L(1/2, √2).

The multi-large-prime extension: instead of requiring candidates to be fully B-smooth, allow cofactors with up to k large prime factors. Factor these cofactors using ECM (cheap for small cofactors). Build a hypergraph of shared large primes and find even-multiplicity subsets via structured Gaussian elimination.

**Status**: BSRF v1 implemented and working on 30-digit semiprimes. Using sieve for candidate detection + trial division for smooth factoring + single-LP merging. Need to test scaling and add batch GCD optimization.

## Genus-2 HECM

Cosset (2009) showed that genus-2 hyperelliptic curve method (HECM) using Kummer surface arithmetic effectively runs 2 ECM curves simultaneously. Claims faster than GMP-ECM for large factors. The Jacobian of a genus-2 curve has group order ~p² (vs ~p for elliptic curves), with Hasse interval width ~p^{3/2}. The wider Hasse interval means more group orders to sample from, potentially increasing the probability of hitting a smooth order.

Asymptotically still L_p(1/2), but the constant may be better for medium-sized factors. Worth testing for our 30-80 digit regime.

## Current approaches

### BSRF (Batch Smooth Relation Finder) — `library/bsrf.c`

Core algorithm:
1. Choose smoothness bound B = L(1/2, 0.5) and large prime bound LP = 100*B
2. For polynomial Q(x) = (x+m)² - N where m = ⌊√N⌋:
   - Use logarithmic sieve to identify likely-smooth candidates
   - Trial divide candidates over the factor base
   - Accept relations that are either fully B-smooth or have 1 large prime ≤ LP
3. Merge pairs of 1-LP relations sharing the same large prime
4. Build GF(2) matrix from full relations + merged pairs
5. Gaussian elimination to find null-space vectors (congruences of squares)
6. Extract factors via gcd

**Key observations (30-digit tests)**:
- Full smooth relations ~10%, single-LP ~90% of useful candidates
- LP merging works: 807 merges from 3891 1-LP relations with LP bound 474K
- Gaussian elimination finds ~64 dependencies, most yield valid factorizations

### SRG (Smooth Residue Graph) — `library/srg.c` (other agent)

Multi-stage pipeline: trial division → Pollard's rho → ECM → sieve. Includes aggressive multi-LP collection (up to 3 LPs per relation) with extended GF(2) matrix including LP columns.

## Open directions

- **Algebraic group structure**: Z_N* ≅ Z_{p-1} × Z_{q-1} but we can't see this decomposition. Can random walks, character sums, or higher-dimensional algebraic groups reveal it without smooth numbers?
- **Smooth polynomial families**: Do certain polynomial constructions yield systematically smoother values than random?
- **Batch smoothness (Bernstein)**: Product/remainder trees for smoothness detection. Better I/O complexity than sieving for large bounds.
- **Function field analogies**: Can techniques from the quasi-polynomial DLP breakthrough be adapted to number fields? The key obstacle: no Frobenius endomorphism over Z.
- **Recursive cofactor descent**: When Q(x) = smooth_part · cofactor, can we recursively factor the cofactor using a SMALLER polynomial? The cofactor is much smaller than Q(x), so recursive application could reduce effective value sizes geometrically.
- **Cross-polynomial relation sharing**: Using multiple polynomials (different multipliers k in Q_k(x) = (x+m_k)² - kN), can partial smoothness in one polynomial boost another?
- **Cascaded algebraic descent**: Inspired by the descent step in NFS-DLP — instead of sieving a polynomial, use a descent strategy to recursively decompose random values into smooth representations. Each descent level uses smaller primes. Could this change the L-exponent?
- **Multi-resolution factor bases**: Use a hierarchy of smoothness bounds B1 < B2 < B3, where each level relaxes smoothness requirements but adds more "virtual columns" to the LA matrix.
