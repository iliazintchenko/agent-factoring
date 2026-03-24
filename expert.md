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
- **Bilinear smoothness decomposition**: Tried to split smoothness testing into two independent N^{1/3}-sized tests: find a,b ~ N^{1/3} with a*b mod N smooth. But a*b ~ N^{2/3} (no modular reduction since a*b < N), so the "third part" is N^{2/3}, not N^{1/3}. Can't achieve NFS-like norm splitting without actual number field structure.
- **LLL polynomial selection for degree 2**: For f(x) = ax²+bx+c with f(m) ≡ 0 (mod N), LLL finds (a,b) with small coefficients. But the OPTIMAL degree-2 polynomial is a=1, b=-2m (the standard QS polynomial). LLL can't improve on this — it already IS the shortest vector. For degree 3+, coefficients ~N^{1/(d+1)} give values ~N^{d/(d+1)} which is WORSE than QS for integer smoothness. Only NFS's algebraic norm makes higher degree useful.
- **Smooth number enumeration near targets (SNAMC)**: Generating B-smooth numbers near √N by multiplying factor base primes. This is a subset-sum problem — exponentially hard. The Dickman-based density of smooth numbers is high, but FINDING specific-valued ones is as hard as enumerating ALL smooth numbers up to the target.
- **Extended matrix with ECM-descended cofactors (CDS)**: Using ECM to split sieve cofactors into medium primes, then adding each medium prime as a column in an extended GF(2) matrix. **Fails** because the number of unique medium primes grows faster than the number of relations — each ECM descent produces new, previously-unseen primes, rapidly expanding the column count. With B=2201 on 30-digit N: 1764 relations but 2487 unique medium primes = 2798 columns. The matrix is always underdetermined. The correct approach for cofactor descent is LP-matching (combine partials sharing a medium prime), not extended matrix.

## Research survey

- **No classical sub-L[1/3] algorithm exists** for general integer factoring.
- **Regev (2023)**: Quantum factoring in O~(n^{3/2}) gates. Fundamentally quantum, no classical dequantization known.
- **Stange (2022)**: Index calculus in (Z/NZ)*. Clean but same L[1/2] or L[1/3]. Not an improvement.
- **Umans & Wang (2025)**: Conditional deterministic N^{1/6+o(1)}. Still exponential.
- **Harvey-Hittmeir**: Deterministic N^{1/5+o(1)}. Still exponential.
- **Tower NFS (Barbulescu-Kim 2016)**: Achieves sub-L[1/3] for DLP in specific finite field extensions via tower of number fields. Key ingredient: Frobenius action allows systematic polynomial splitting at each tower level. No analog over Z for integer factoring — the "descent" works because polynomials over finite fields factor efficiently, but integers don't.

## Why novel approaches are hard — comprehensive analysis

The core difficulty: any classical approach to factoring must somehow distinguish the factor p from ~√N other candidates. Sub-exponential algorithms do this via smooth number relations + linear algebra.

**The L-exponent hierarchy and what determines it:**

1. **L[1/2]** (QS-class): Sieve a degree-2 polynomial f(x) = x² - N. Values ≈ M√N at sieve boundary M. LLL analysis shows THIS IS OPTIMAL for degree-2 polynomials over Z — the standard QS polynomial IS the shortest vector in the coefficient lattice. No lattice trick can improve it.

2. **L[1/3]** (NFS-class): Use algebraic number fields to decompose the norm into TWO parts (algebraic + rational), each smaller than the single polynomial value. This "norm splitting" is the KEY INGREDIENT that improves the exponent. It requires:
   - A polynomial f(x) of degree d ≥ 3 with a root m mod N
   - An algebraic number field Q(α) where f(α) = 0
   - The ability to define "smoothness" in this number field (via prime ideals)

3. **Sub-L[1/3]** (Tower NFS for DLP, no analog for factoring): Tower of field extensions enables recursive descent where polynomials factor efficiently. Works for DLP in GF(p^n) because the Frobenius endomorphism provides a systematic splitting. No analog over Z.

**What WOULD be needed to beat L[1/3] for factoring:**
- A new algebraic structure over Z (not number fields, not function fields) where "norms" are smaller than N^{2/3} per component
- OR a way to generate smooth numbers of size < N^{1/3} that still relate to N's factorization
- OR a completely non-smoothness-based approach (no known sub-exponential method of this type exists)

**Approaches that CAN'T help (proven by our experiments):**
- Multi-multiplier combination (MMS): only improves constants, not the exponent
- Bilinear decomposition: can't split norms without number fields
- LLL polynomial selection: degree-2 is already optimal; degree-3+ needs NFS
- Smooth number enumeration: subset-sum hard
- ECM cofactor descent (CDS): too many unique medium primes for extended matrix; LP-matching doesn't change the exponent
- Birthday/collision methods: O(N^{1/4}), not sub-exponential

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
- Scaling data (worst-case across 5 semiprimes per size):
  - 30-40d: 1.7-28.9s (mostly P-1/P+1, some ECM)
  - 41-50d: 28-35s (ECM with progressive bounds)
  - 51-55d: 8-45s (ECM level 3-5)
  - 56-61d: 49-221s (ECM level 5-6, high variance)
- The approach covers 30-65d range but is outperformed by sieve methods (MMS) for <55d
- Best niche: 55-65d where ECM's progressive schedule handles balanced semiprimes well

### CAD (Cascaded Algebraic Descent) — `library/cad.c`

QS-variant with aggressive cofactor splitting via Pollard's rho and ECM. Uses log sieve for candidate identification + trial division + single-LP merging + GF(2) elimination.

**Key observations**:
- 30-digit: ~10s (much slower than MMS at 0.074s)
- The bottleneck is trial division, not the sieve or LA
- Cofactor splitting via rho/ECM adds ~10-15% more usable relations but the per-candidate cost is high
- **Dead end for small digits**: standard sieve (like MMS) is much faster because the per-candidate cost of trial division + cofactor splitting outweighs the benefit of additional relations
- **Potential for larger digits**: at 50+ digits, the ratio of LP-resolvable cofactors increases, and each resolved cofactor saves significant sieve time

### BatchQS (Batch GCD QS variant) — `library/batchqs.c`

Uses Bernstein's product tree approach for batch smoothness detection instead of a sieve. Computes gcd(primorial^2, Q(x)) for batches of candidates using remainder trees.

**Status**: Prototype working but low yield — the batch GCD correctly identifies smooth candidates but misses ~80% compared to a standard sieve. Likely cause: the primorial power isn't high enough to extract all smooth prime powers. Needs further investigation.

### MMCFRAC (Multi-Multiplier CF) — `library/mmcfrac.c`

Novel: combines CFRAC with MMS-style multi-multiplier. For K square-free multipliers, runs CF expansion of sqrt(kN) simultaneously. CF convergents naturally minimize |p^2 - kNq^2| < 2√(kN). Cross-multiplier LP merging adds "free" relations.

**Scaling**: 30d: 2.2s, 35d: 20s, 40d: 191s, 45d: FAIL. Bottleneck: trial division O(fb) per CF step with no pre-filtering.

### IFM (Iterated Frobenius Map) — `library/ifm.c`

**Dead end**: x → x^N mod N iteration. O(N^{1/4}) cycle length but each step costs O(log N) multiplications (100x more than rho). Net ~100x slower for no compensating advantage.

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

**Key finding on multiplier selection**: Smallest square-free multipliers (1,2,3,5) outperform KS-scored larger multipliers (22,43,30,42) by 10-20%. Reason: KS scoring optimizes for factor base density but ignores that larger k gives larger polynomial values f_k(x) ≈ 2x√(kN). The value size increase outweighs the FB density gain.

**Key implementation detail**: Multipliers MUST be square-free. If k₂ = m² · k₁ for integer m, then f_{k₂}(x) = m² · f_{k₁}(x') for appropriate x'. The factor m² has even exponent, so both relations have identical GF(2) parity vectors. Combined relations from such pairs are trivially zero — useless for linear algebra. Using square-free k values (1, 2, 3, 5, 6, 7, 10, ...) avoids this.

**Implementation**: `library/mms.c`. Sieve + trial division + single-LP matching + GF(2) elimination.

**Scaling results (worst-of-5 per size)**:
| Digits | Time   | Digits | Time    |
|--------|--------|--------|---------|
| 30     | 0.07s  | 41     | 1.97s   |
| 31     | 0.10s  | 42     | 2.86s   |
| 32     | 0.11s  | 43     | 3.66s   |
| 33     | 0.21s  | 44     | 5.06s   |
| 34     | 0.24s  | 45     | 6.79s   |
| 35     | 0.28s  | 46     | 10.0s   |
| 36     | 0.40s  | 47     | 14.3s   |
| 37     | 0.69s  | 48     | 16.9s   |
| 38     | 0.76s  | 49     | 19.1s   |
| 39     | 1.15s  | 50     | 20.0s   |
| 40     | 1.41s  | 51     | 30.6s   |

Comparison with SRG (ECM-based): MMS is faster at 30-42 digits (e.g., 30d: 0.07s vs 0.83s), but slower at 50+ digits (50d: 20s vs 2.9s). The SRG approach benefits from ECM's subexponential scaling, while MMS is a pure sieve method with L[1/2] scaling.

55-digit test: ~75s (truncated FB) to ~132s (full FB). Beyond 55 digits impractical with current parameters.

**Scaling analysis**: L[1/2] and L[1/3] fits comparably (SSE 0.477 vs 0.440 over 30-53d). Growth ~1.3x per digit.

**Key finding**: MMS is competitive only at 30-35 digits. BSRF-v3 (uses SIQS-style smaller polynomials) is 3-20x faster across all sizes. The multi-multiplier LP matching provides ~60% more relations (constant factor), but doesn't shrink the polynomial values — the fundamental bottleneck. Lesson: improving the L-constant requires smaller polynomial values, not just more relations from the same-sized values.

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
