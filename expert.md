# Expert Knowledge

## The smoothness bottleneck

All known sub-exponential factoring algorithms rely on finding **smooth numbers**. The L-exponent is determined by the size of values being tested for smoothness:
- QS: tests values of size ~√N → L[1/2]
- NFS: uses algebraic number fields to reduce values to ~N^(2/3) → L[1/3]
- Hypothetical L[1/4]: would need values of size ~N^(1/4) or equivalent

Any approach that still tests numbers of size √N for smoothness will remain L[1/2]. To do better, you must either generate smaller candidates, or avoid smoothness entirely.

The quasi-polynomial DLP breakthrough (Barbulescu-Gaudry-Joux-Thomé 2013) shows L[1/3] barriers CAN be broken for related problems in function fields. The technique exploits systematic factorization of polynomials over finite fields — no known analog over Z. The "translation problem" from function fields to number fields is a major open question.

## Why the Pell equation is special (degree-2 uniqueness)

The continued fraction of √N produces convergents p_k/q_k with **bounded residues**: |p_k² - N·q_k²| < 2√N, regardless of k. This infinite family of small residues is what makes CFRAC and QS possible.

For higher-degree analogs (e.g., N^{1/3}), the situation is fundamentally different:
- CF of N^{1/3}: convergent residues |p_k³ - N·q_k³| ≈ N^{2/3}·q_k, which **grows** with q_k
- Thue's theorem: the equation |x^d - N·y^d| < C has **finitely many** solutions for d ≥ 3
- Consequence: there is NO infinite family of small "cubic residues" analogous to the Pell equation

This means degree ≥ 3 polynomial evaluation over Z cannot produce an unbounded family of small values. The ONLY way to use higher-degree polynomials effectively is via **number fields** (NFS), where the "norm" provides a different notion of size that isn't subject to Thue's finiteness.

## The L-exponent formula for sieve algorithms

For a sieve testing values of size N^β for B-smoothness:
- Optimal smoothness bound: B = L[1/2, √(β/2)]
- Smooth probability: ρ(u) where u = β·ln(N)/ln(B)
- Total complexity: **L[1/2, √(2β)]**

| Approach | Value size (β) | L-constant √(2β) |
|----------|---------------|-------------------|
| Dixon's  | 1             | √2 ≈ 1.414       |
| QS       | 1/2           | 1.000             |
| Hypothetical β=1/3 | 1/3 | √(2/3) ≈ 0.816  |
| NFS      | ~2/3 but split | L[1/3] (different formula) |

NFS achieves L[1/3] by splitting the norm across TWO polynomials (algebraic + rational), exploiting the number field structure to make each part smaller than what a single polynomial could achieve.

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
- **Cubic/higher-degree continued fractions**: Residues grow (Thue's theorem), giving values LARGER than QS, not smaller. A degree-d polynomial near N^{1/d} gives values ~d·N^{(d-1)/d}·M, which is worse than QS for d > 2.
- **Group-order methods (p-1, p+1, Gaussian/Eisenstein integers)**: All give fixed group orders (p±1 variants). For balanced semiprimes where p±1 are not smooth, these are ineffective. ECM randomizes the group order but is still L[1/2].
- **Character sum / autocorrelation approaches**: Detecting the period p in χ_N(n) requires Ω(√N) samples (information-theoretic lower bound). Same complexity as trial division.
- **Birthday-paradox collision methods (Pollard rho variants)**: O(N^{1/4}) for balanced semiprimes regardless of the map used. Not sub-exponential.

## Research survey

- **No classical sub-L[1/3] algorithm exists** for general integer factoring.
- **Regev (2023)**: Quantum factoring in O~(n^{3/2}) gates. Fundamentally quantum, no classical dequantization known.
- **Stange (2022)**: Index calculus in (Z/NZ)*. Clean but same L[1/2] or L[1/3]. Not an improvement.
- **Umans & Wang (2025)**: Conditional deterministic N^{1/6+o(1)}. Still exponential.
- **Harvey-Hittmeir**: Deterministic N^{1/5+o(1)}. Still exponential.
- **Tower NFS (Barbulescu-Kim 2016)**: Achieves sub-L[1/3] for DLP in specific finite field extensions via tower of number fields. Key ingredient: Frobenius action allows systematic polynomial splitting at each tower level. No analog over Z for integer factoring — the "descent" works because polynomials over finite fields factor efficiently, but integers don't.

## Why novel approaches are hard

The core difficulty: any classical approach to factoring must somehow distinguish the factor p from ~sqrt(N) other candidates of similar size. Sub-exponential algorithms achieve this via smooth number relations + linear algebra. To improve the L-exponent, you need SMALLER smooth candidates, which requires new algebraic structures that reduce norm sizes below what NFS achieves.

## Key insight: smoothness detection cost matters

Standard QS/NFS sieving has per-candidate cost O(1) amortized but requires sequential memory access over a fixed interval. Bernstein's batch smoothness detection (product/remainder trees) has per-candidate cost O(log²B · polylog) but works on ARBITRARY candidate sets. This unlocks two opportunities:

1. **Larger smoothness bounds**: Sieving cost scales as O(B/log B) per candidate; batch GCD scales as O(log²B). For large B, batch GCD is cheaper. A larger B means higher smoothness probability, fewer candidates needed.

2. **Non-sequential candidate generation**: Sieving requires evaluating consecutive polynomial values. Batch GCD works on any set of values, enabling novel candidate generation strategies.

## Genus-2 HECM

Cosset (2009) showed that genus-2 hyperelliptic curve method (HECM) using Kummer surface arithmetic effectively runs 2 ECM curves simultaneously. Claims faster than GMP-ECM for large factors. The Jacobian of a genus-2 curve has group order ~p² (vs ~p for elliptic curves), with Hasse interval width ~p^{3/2}. The wider Hasse interval means more group orders to sample from, potentially increasing the probability of hitting a smooth order.

Asymptotically still L_p(1/2), but the constant may be better for medium-sized factors. Worth testing for our 30-80 digit regime.

## Current approaches

### DBRM (Descent-Based Relation Mining) — `library/dbrm.c`

Aggressive multi-large-prime sieve with graph-based relation combining. Uses single polynomial Q(x) = (x+m)²-N with large smoothness bound B and LP matching.

**Key observations (30-digit tests)**:
- With B=50000 (fb=2539): 5165 smooth + 4352 LP-merged = 9517 rels for 2540 cols
- Factors 30-digit semiprimes in ~7s
- LP merging provides ~45% of useful relations

### BSRF (Batch Smooth Relation Finder) — `library/bsrf.c`

Core algorithm: sieve + trial division + single-LP merging + GF(2) elimination.

**Key observations (30-digit tests)**:
- Full smooth relations ~10%, single-LP ~90% of useful candidates
- LP merging works: 807 merges from 3891 1-LP relations with LP bound 474K

### SRG (Smooth Residue Graph) — `library/srg.c`

Multi-stage pipeline: trial division → Pollard's rho → progressive ECM → sieve with multi-LP relations + extended GF(2) LA.

**Key design**: Progressive ECM with increasing B1 bounds (2K→43M) dramatically reduces worst-case variance. Also includes a sieve stage with aggressive multi-LP collection (up to 3 LPs per relation) — all relations (full + partial) go into an extended GF(2) matrix where each unique large prime is an extra column. Null-space vectors yield congruences of squares.

**Scaling data (worst-case across 5 semiprimes per size)**:
- 30-36d: <1s (Pollard's rho or quick ECM)
- 37-43d: 1-3s (ECM level 1-2)
- 44-50d: 3-7s (ECM level 3)
- 51-55d: 7-75s (ECM level 4-5, high variance)
- 56-60d: 12-88s (ECM level 5-6)
- 61-65d: 44-215s (ECM level 6-7, approaching limits)
- 66d+: ECM timeout on some inputs, sieve kicks in but too slow (needs multi-polynomial)

**Bottleneck at 66+ digits**: The sieve collects ~165 full relations vs ~11K needed. Partial-1 relations (12K) outnumber full relations 77:1, but each unique LP adds a column to the matrix, so rows < columns. Need either (a) much longer sieve time, (b) SIQS-style multi-polynomial to generate smaller Q(x) values, or (c) tighter LP bound to increase LP collisions.

### Hybrid (P-1/P+1/ECM engine) — `library/hybrid.c`

Uses GMP-ECM library with progressive parameter scheduling. Pipeline:
1. Trial division to 10^6
2. Pollard's rho (500K iterations)
3. P-1 at multiple B1 levels (2K to 110M)
4. P+1 at same levels
5. ECM with Suyama parameterization, progressive B1 schedule

**Key observations**:
- Surprisingly effective on balanced semiprimes: P-1/P+1 alone factor most 30-50 digit semiprimes
- This suggests the test semiprimes have somewhat smooth p-1 or p+1 values
- Scaling 30-55 digits: {30:1.7s, 35:2.6s, 40:3.1s, 45:3.5s, 50:30s, 55:45s}
- The jump at 50 digits suggests P-1/P+1 reaches its limits and ECM kicks in
- For truly hard semiprimes (both p-1 and p+1 have large prime factors), would need ECM with many more curves

### IFM (Iterated Frobenius Map) — `library/ifm.c`

Novel iteration function: x → x^N mod N (instead of rho's x → x^2+c). Via CRT, decomposes into independent power maps mod p and mod q with different dynamics. **Dead end**: each iteration costs O(log N) multiplications (100x more than rho per step), but cycle length is similar (O(N^{1/4})). Net result is ~100x slower than Pollard's rho for no compensating advantage.

## Active approaches

### Multi-Multiplier Sieve (MMS)

**Core idea**: For each sieve point x, evaluate *multiple* polynomials simultaneously:
  f_k(x) = (x + ⌊√(kN)⌋)² − kN,  for k = 1, 2, 3, …, K

Since (x + m_k)² ≡ f_k(x) (mod N) for each k, *products* across multipliers give valid congruences of squares:
  ∏(x + m_k)² ≡ ∏ f_k(x) (mod N)

**Why this differs from known methods**:
- Not QS/MPQS: those use different polynomials for different x ranges; MMS uses multiple polynomials at the *same* x and combines across them.
- Not NFS: no number fields involved. Values are still ~√(kN)·M, i.e., L[1/2]-class.
- The novelty is *cross-multiplier combination*: partial factorizations from different k values for the same x can be combined. If f_1(x) has large cofactor L and f_3(x) also has L as a factor, their product has L² (even power) — a full relation emerges from two partial ones.

**Expected benefits**:
1. Each prime p with Legendre symbol (kN|p) = 1 for *any* k enters the factor base. More k values → more primes → denser sieving.
2. Cross-k large-prime collisions add "free" relations via birthday paradox. With K multipliers, O(K²) collision opportunities per x.
3. The Knuth-Schroeppel multiplier effect is captured automatically — the best k for each prime is used.

**Theoretical expectation**: Still L[1/2] asymptotically (value sizes unchanged), but potentially much better constants. The practical scaling curve may bend differently than single-polynomial methods for 30-70 digit numbers.

**Key implementation detail**: Multipliers MUST be square-free. If k₂ = m² · k₁ for integer m, then f_{k₂}(x) = m² · f_{k₁}(x') for appropriate x'. The factor m² has even exponent, so both relations have identical GF(2) parity vectors. Combined relations from such pairs are trivially zero — useless for linear algebra. Using square-free k values (1, 2, 3, 5, 6, 7, 10, ...) avoids this.

**Implementation**: `library/mms.c`. Sieve + trial division + single-LP matching + GF(2) elimination.

**Scaling results (30-35 digits)**:
| Digits | Worst-of-5 time |
|--------|----------------|
| 30     | 0.074s         |
| 31     | 0.095s         |
| 32     | 0.103s         |
| 33     | 0.205s         |
| 34     | 0.236s         |
| 35     | 0.278s         |

Larger sizes in progress.

## Open directions

- **Algebraic group structure**: Z_N* ≅ Z_{p-1} × Z_{q-1} but we can't see this decomposition. Can random walks, character sums, or higher-dimensional algebraic groups reveal it without smooth numbers?
- **Smooth polynomial families**: Do certain polynomial constructions yield systematically smoother values than random?
- **Batch smoothness (Bernstein)**: Product/remainder trees for smoothness detection. Better I/O complexity than sieving for large bounds.
- **Function field analogies**: Can techniques from the quasi-polynomial DLP breakthrough be adapted to number fields? The key obstacle: no Frobenius endomorphism over Z.
- **Recursive cofactor descent**: When Q(x) = smooth_part · cofactor, can we recursively factor the cofactor using a SMALLER polynomial? The cofactor is much smaller than Q(x), so recursive application could reduce effective value sizes geometrically.
- **Cross-polynomial relation sharing**: Using multiple polynomials (different multipliers k in Q_k(x) = (x+m_k)² - kN), can partial smoothness in one polynomial boost another?
- **Cascaded algebraic descent**: Inspired by the descent step in NFS-DLP — instead of sieving a polynomial, use a descent strategy to recursively decompose random values into smooth representations. Each descent level uses smaller primes. Could this change the L-exponent?
- **Multi-resolution factor bases**: Use a hierarchy of smoothness bounds B1 < B2 < B3, where each level relaxes smoothness requirements but adds more "virtual columns" to the LA matrix.
- **Aggressive descent for factoring**: Explore whether multi-level descent (as used in DLP) can improve the factoring sieve beyond the standard large-prime variation.
