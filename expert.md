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

## Research survey

- **No classical sub-L[1/3] algorithm exists** for general integer factoring.
- **Regev (2023)**: Quantum factoring in O~(n^{3/2}) gates. Fundamentally quantum, no classical dequantization known.
- **Stange (2022)**: Index calculus in (Z/NZ)*. Clean but same L[1/2] or L[1/3]. Not an improvement.
- **Umans & Wang (2025)**: Conditional deterministic N^{1/6+o(1)}. Still exponential.
- **Harvey-Hittmeir**: Deterministic N^{1/5+o(1)}. Still exponential.

## Current approach: Smooth Residue Graph (SRG)

The SRG approach (`library/srg.c`) combines:
1. Multi-stage pipeline: trial division → Pollard's rho → ECM → sieve
2. Quadratic polynomial Q(x) = (x+m)^2 - N with sieve-based candidate selection
3. Trial division + ECM for cofactor splitting in partial relations
4. Aggressive multi-large-prime collection (up to 3 LPs per relation)
5. Extended GF(2) matrix including LP columns for linear algebra
6. Graph-based cycle detection for combining partial relations

**Key design choice**: Include ALL relations (full + partial with 1-3 large primes) in a single extended linear algebra matrix. Each unique large prime becomes an extra column. If we have more rows than columns, the null space is non-empty and yields congruences of squares. This avoids the separate "graph matching" step used in traditional QS implementations.

**Scaling hypothesis**: By including multi-LP relations, we effectively use a much larger "virtual factor base" while only needing to find relations that are smooth over the actual (small) factor base plus a bounded number of large primes. This should improve the practical constant in L[1/2].

**Observations so far**:
- 30-digit semiprimes: ECM finds factors directly in <1s (15-digit factors)
- The sieve stage collects relations at a reasonable rate for 30-40 digits
- LP relations outnumber full relations ~7:1 with LP bound = 30 * max_FB_prime
- The extended matrix approach works: null-space vectors yield valid congruences

## Why novel approaches are hard

After extensive analysis, several "novel" directions were explored and found to reduce to known algorithms:

- **Lattice-based smooth value generation**: Constructing lattices where short vectors correspond to smooth products mod N reduces to Schnorr's failed approach (2021). The fundamental issue: short vectors in the lattice have small exponents → small products → no wrapping mod N.
- **Higher-dimensional sieving (tower NFS)**: Using towers of number fields K_1 ⊂ K_2 ⊂ ... to reduce norms. Analysis shows: for integer factoring (vs. DLP in finite fields), extra dimensions don't reduce effective norm sizes. The tower NFS breakthrough for GF(p^n) doesn't translate.
- **Multiplicative random walks**: Using random walks in (Z/NZ)* to find smooth representations. Collision time is ~sqrt(N) by birthday bound. This is Pollard's rho.
- **Meet-in-the-middle on exponential spaces**: BSGS-style decomposition of the search space. Requires ~N^{1/4} storage, still exponential.
- **ECM failure information aggregation**: Collecting smooth/non-smooth classifications of curve orders across many ECM curves to constrain p. Analysis shows: success rate only reveals factor SIZE (already known), not specific factor IDENTITY.

The core difficulty: any classical approach to factoring must somehow distinguish the factor p from ~sqrt(N) other candidates of similar size. Sub-exponential algorithms achieve this via smooth number relations + linear algebra. To improve the L-exponent, you need SMALLER smooth candidates, which requires new algebraic structures that reduce norm sizes below what NFS achieves.

## Open directions

These are starting points, not an exhaustive list.

- **Algebraic group structure**: Z_N* ≅ Z_{p-1} × Z_{q-1} but we can't see this decomposition. Can random walks, character sums, or higher-dimensional algebraic groups reveal it without smooth numbers?
- **Smooth polynomial families**: Do certain polynomial constructions yield systematically smoother values than random?
- **Batch smoothness (Bernstein)**: Product/remainder trees for smoothness detection. Better I/O complexity than sieving for large bounds.
- **Function field analogies**: Can techniques from the quasi-polynomial DLP breakthrough be adapted to number fields?
- **Cascaded algebraic descent**: Inspired by the descent step in NFS-DLP — instead of sieving a polynomial, use a descent strategy to recursively decompose random values into smooth representations. Each descent level uses smaller primes. Could this change the L-exponent?
- **Multi-resolution factor bases**: Use a hierarchy of smoothness bounds B1 < B2 < B3, where each level relaxes smoothness requirements but adds more "virtual columns" to the LA matrix. The question: does the proliferation of partial relations at coarser levels compensate for the increased matrix size?
