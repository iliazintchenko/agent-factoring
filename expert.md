# Expert Knowledge

## The smoothness bottleneck

All known sub-exponential factoring algorithms rely on finding **smooth numbers**. The L-exponent is determined by the size of values being tested for smoothness:
- QS: tests values of size ~√N → L[1/2]
- NFS: uses algebraic number fields to reduce values to ~N^(2/3) → L[1/3]
- Hypothetical L[1/4]: would need values of size ~N^(1/4) or equivalent

Any approach that still tests numbers of size √N for smoothness will remain L[1/2]. To do better, you must either generate smaller candidates, or avoid smoothness entirely.

The quasi-polynomial DLP breakthrough (Barbulescu-Gaudry-Joux-Thomé 2013) shows L[1/3] barriers CAN be broken for related problems in function fields. The technique exploits systematic factorization of polynomials over finite fields — no known analog over Z.

## Why degree-2 is special (Thue's theorem)

The CF of √N produces convergents with **bounded residues**: |p_k² - N·q_k²| < 2√N regardless of k. This infinite family of small residues is what makes QS possible.

For degree ≥ 3, Thue's theorem proves |x^d - N·y^d| < C has **finitely many** solutions. No infinite family of small higher-degree residues exists over Z. The ONLY way to use higher-degree polynomials is via **number fields** (NFS), where the norm provides a different notion of size.

## The L-exponent formula

For a sieve testing values of size N^β: total complexity is **L[1/2, √(2β)]**.

| Approach | Value size (β) | Complexity |
|----------|---------------|------------|
| Dixon's  | 1             | L[1/2, √2] |
| QS       | 1/2           | L[1/2, 1]  |
| NFS      | ~2/3 but split across two number fields | L[1/3] |

LLL analysis confirms QS's degree-2 polynomial IS the shortest vector in the coefficient lattice — no lattice trick can improve it. For degree 3+, you need number fields (NFS).

## Why higher-dimensional algebras don't help

Quadratic extensions Z[√d]: norm is a² - d·b², degree 2 in two variables. Same value sizes as QS regardless of d. Can't improve beyond L[1/2].

Quaternion algebras: norm is a² + b² + c² + d², degree 4 overall but still **degree 2 in each variable**. Values are ~max(a,b,c,d)² — four variables don't reduce the norm size. More variables don't help because the norm is quadratic per variable.

**Key constraint**: NFS achieves L[1/3] not by having a better norm, but by splitting the problem across two different norm maps (rational and algebraic) and requiring simultaneous smoothness. Adding more number fields (multi-image NFS) only improves L[1/3] constants, not the exponent. To beat L[1/3], you'd need either a fundamentally different splitting mechanism or an approach that avoids smoothness-based norms entirely.

## What would be needed to beat L[1/3]

- A new algebraic structure over Z where "norms" are smaller than N^{2/3} per component
- OR smooth numbers of size < N^{1/3} that relate to N's factorization
- OR a completely non-smoothness-based sub-exponential approach (none known)

Tower NFS achieves sub-L[1/3] for DLP in GF(p^n) via Frobenius-enabled polynomial splitting. No analog over Z for integer factoring — the descent works because polynomials over finite fields factor efficiently, but integers don't.

## Dead ends (tested or analyzed)

- **Schnorr lattice factoring**: Lattice dimension grows with smoothness bound. Ducas (CWI): 0 relations in 1000 trials.
- **Spectral methods on Cayley graphs**: Computing the spectrum of (Z/NZ)* is exponential classically — exactly what Shor's QFT solves.
- **Multi-image NFS**: More number fields improve L[1/3] constants but can't break the barrier. Simultaneous smoothness across fields constrains the sieve region proportionally.
- **Iterated Frobenius Map** (x → x^N mod N): O(N^{1/4}) cycle length but 100x cost per step vs Pollard rho. Net slower.
- **ECM cofactor descent (CDS)**: Splitting cofactors with ECM produces too many unique medium primes — matrix is always underdetermined.
- **Smooth number enumeration near √N**: Equivalent to subset-sum — exponentially hard.
- **LLL polynomial selection**: Degree-2 already optimal (QS polynomial IS the shortest vector). Degree-3+ needs NFS.
- **Bilinear smoothness decomposition**: Can't split norms without number fields — a·b mod N doesn't reduce since a·b < N.
- **Character sum / autocorrelation**: Detecting period p requires Ω(√N) samples.
- **Cubic/higher-degree CF**: Residues grow (Thue), worse than QS.
- **Group-order methods for balanced semiprimes**: p±1 methods fail when p±1 aren't smooth. ECM randomizes but stays L[1/2].
- **Multi-multiplier sieve (MMS)**: Cross-multiplier LP matching gives ~60% more relations (constant factor) but doesn't shrink polynomial values.
- **Quaternion norms for smoothness testing**: Norm is still quadratic per variable — more variables don't reduce value sizes. However, quaternion algebras have other relevant structure (Deuring correspondence with supersingular elliptic curves, Brandt matrices encoding Hecke operators, Pizer's Ramanujan graphs) that could theoretically leak factoring information through non-smoothness-based channels. Unexplored.
- **BatchQS (naive implementation)**: Bernstein product-tree batch GCD misses ~80% of smooth candidates vs standard sieve. The technique itself is sound but the implementation used insufficient primorial powers.
- **Berlekamp polynomial splitting over Z/NZ**: Compute x^N mod (f(x), N) for random degree-2 polynomials f. The gcd(x^N - x, f) over Z/NZ should split f differently mod p vs q, revealing a factor. FAILS for balanced semiprimes because: over F_p with root r of f, x^N maps r → r^q, and gcd(x^N-x, f) is non-trivial only when r^{q-1} ≡ 1 (mod p), which requires ord(r) | (q-1). For random r ∈ F_p*, this has probability gcd(q-1, p-1)/(p-1) ≈ log(p)/p ≈ 0. This is essentially the Pollard p-1 condition — the approach degenerates to p-1 factoring for balanced semiprimes.
- **Multiplicative Lattice Descent (MLD)**: Extending the LP bound from B to B^α (α>2) increases partial relation rate but the cofactor factor base grows as B^(α-1)/log(B). Tested for 30-digit N with α=3: MLD needs 2.16x MORE sieve work than standard QS because the matrix size explosion (100k × 100k) outweighs the partial relation gain (143x). The cofactor-to-matrix tradeoff is inherently bounded: any LP expansion that multiplies the relation rate by R also multiplies the matrix dimension by ≥ R, so the L-exponent cannot improve. This is the fundamental reason LP variations (single, double, triple) improve constants but not the exponent.
- **Class group 2-Sylow exploration**: For disc -4N, the class group has ambiguous forms including (p,0,q). Finding this via repeated squaring g^(2^k) gives order-2 elements, but most are "trivial" (a=2 type, from the discriminant factor 4). Needs BSGS or Shanks-style search to find the factoring form — O(N^{1/4}) time, same as Pollard rho. Not better than QS.
- **Multi-base GCD accumulation**: Running Pollard p-1 with many random bases and B1 up to 50000 doesn't factor balanced semiprimes (p-1 not smooth). Multi-curve ECM accumulation with B1=1000 also fails for 30-digit semiprimes. Confirmed: group-order methods fundamentally fail for balanced random semiprimes.
- **Smoothness correlations across polynomial pairs**: Tested whether P1(x)=(x+m)²-N and P1(x+d) have correlated smoothness (would allow joint sieving). For 15-digit N with B=2000, measured conditional vs. unconditional smooth probability across shifts d=1..100 and cross-multiplier k=2. Result: no statistically significant correlation (ratios 0.85-1.63, all within 2σ of 1.0). Theoretical analysis confirms: for each prime p, the sieve conditions at x and x+d are independent (different residue classes), so joint smoothness equals the product of marginal probabilities. Cross-multiplier smoothness is also independent for the same reason. Constant-factor improvements from LP variations are possible but L-exponent is unaffected.
- **Polynomial splitting mod composite N** (poly_split.c): For monic f(x) of degree d, compute gcd(x^{(N-1)/2} - 1, f(x)) mod N. If f splits differently mod p vs mod q, intermediate GCD steps encounter non-invertible leading coefficients, revealing a factor. HOWEVER: for monic polynomials, the probability of hitting a non-invertible coefficient is O(1/√N) per trial — essentially zero for large N. This is because all intermediate coefficients are sums/products of N-bit numbers mod N, which are essentially random. 500 trials on 30-digit N found zero factors. Dead end.
- **CRT-structured candidate generation**: Use CRT to force Q(x) = (x+m)² - N to be divisible by several chosen primes simultaneously. Testing shows higher raw smoothness rate (0.59% vs 0.24% for sequential at B=50000, 30 digits), but this is illusory: the CRT forces x to be of size ∏primes, making Q(x) much larger. The cofactor Q(x)/∏primes ≈ 2√N regardless — same as standard QS. No asymptotic improvement.
- **Galois action on NFS relations**: For Galois polynomials, the d automorphisms permute roots but the algebraic norm (= product of all conjugates) is Galois-invariant by definition. Applying automorphisms gives the SAME norm, not new relations. The Galois structure provides d representations of one relation, not d independent relations. At best a constant factor of d improvement in candidate coverage, which doesn't change the L-exponent.
- **Multiplicative lattice relations (MLR)**: Find exponents e_i such that p_1^{e_1} * ... * p_k^{e_k} mod N is small, using LLL on the log-embedding lattice. Even when the residue P(e) mod N is small (achievable via lattice reduction), this doesn't help: P(e) mod p is never zero (since p doesn't divide any p_i), and small P(e) mod N doesn't imply small P(e) mod p due to CRT entanglement. Computing gcd(P(e1)/P(e2) - r1/r2, N) always gives N because the ratio is computed mod N, erasing the p-vs-q distinction. This confirms Schnorr's lattice approach is fundamentally limited — not just practically but theoretically.
- **Best-of-K polynomial selection** (surface_factor.c): For each sieve position, evaluate K different polynomials Q_k(x) and pick the one with smallest |Q_k(x)|. Tested with K=10 on 30-digit N, B=50000, 500K candidates: smoothness rate improved 1.92x but at 10x computational cost → net 5.3x WORSE per smooth value. Theory: the smoothness improvement from minimum of K values scales as K^{u/ln(B)} ≈ K^{0.5} (where u = ln(Q)/ln(B)), which is always SUBLINEAR in K. Best-of-K cherry-picking never beats standard sequential evaluation for any K, B, or N.

## Key insight: smoothness detection cost

Sieving: O(1) amortized per candidate but requires sequential memory access. Bernstein's batch GCD: O(log²B) per candidate but works on ARBITRARY candidate sets. For large B, batch GCD is cheaper and enables non-sequential candidate generation. This is unexploited — no implementation has combined batch GCD with a non-sequential candidate generation strategy effectively.

## Research survey

- **No classical sub-L[1/3] algorithm exists** for general integer factoring.
- **Regev (2023)**: Quantum O~(n^{3/2}) gates. No classical dequantization.
- **Stange (2022)**: Index calculus in (Z/NZ)*. Same L[1/2] or L[1/3].
- **Tower NFS (Barbulescu-Kim 2016)**: Sub-L[1/3] for DLP in GF(p^n) only. No factoring analog.
- **Umans & Wang (2025)**: Conditional deterministic N^{1/6+o(1)}. Exponential.
- **Cosset (2009)**: Genus-2 HECM on Kummer surfaces — 2 ECM curves simultaneously. Still L_p(1/2); only useful if constant improvement hints at a different scaling regime.

## Implementation status

### MPQS baseline (lgsh.c)
Working MPQS implementation with:
- Standard QS + MPQS polynomial generation (single-prime MPQS: A=q², Hensel-lifted B)
- LP matching with correct cofactor handling (divide sv product by lp⁻¹ mod N)
- Knuth multiplier fallback (tries k=1,3,5,7,11,... until extraction succeeds)
- GF(2) Gaussian elimination with up to 128 dependencies

**Measured scaling (LGSH-v5, worst-of-5, MPQS + multiplier + ECM fallback):**
| Digits | Worst (s) | Avg (s) |
|--------|-----------|---------|
| 30 | 0.17 | 0.13 |
| 34 | 0.25 | 0.23 |
| 38 | 3.6 | 1.5 |
| 42 | 10.1 | 8.1 |
| 46 | 24.7 | 15.5 |
| 50 | 42.1 | 33.0 |
| 54 | ~80+ | (running) |

Scaling: roughly consistent with L[1/2]. Key optimizations: smaller FB (exp(0.50*Lexp)) reduces matrix solver time significantly. ECM fallback handles 50+ digit extraction failures. MPQS extraction degeneracy (all dependencies give trivial x≡±y mod N) affects ~25% of 50-digit semiprimes, requiring ECM fallback.

### ECM baseline (ecm_baseline.c)
GMP-ECM with auto-scaled B1 (Suyama parametrization, sigma = 6 + curve_index).

**Worst-of-5 measurements:**

| Digits | Worst-of-5 (s) |
|--------|----------------|
| 30 | 0.07 |
| 32 | 0.05 |
| 34 | 0.10 |
| 36 | 0.12 |
| 38 | 0.47 |
| 40 | 0.41 |
| 42 | 2.21 |
| 44 | 1.69 |
| 46 | 2.33 |
| 48 | 6.67 |
| 50 | 13.97 |
| 52 | 14.20 |
| 54 | 12.68 |
| 56 | 32.75 |
| 60 | 8.20 |

ECM is highly variable (factor of 10x between lucky/unlucky runs). Non-monotonic timing (60-digit faster than 56-digit) is due to favorable p-1 structure in specific test numbers. ECM is L_p[1/2] = L_N[1/2] for balanced semiprimes (since p ~ √N, ln(p) ~ ln(N)/2).

### Known bugs / lessons
- **LP matching must divide by cofactor**: When combining two single-LP relations sharing cofactor c, the combined sqrt_val must be sv1·sv2·c⁻¹ (mod N), not just sv1·sv2. Otherwise x = y for all dependencies.
- **Sign handling in extraction**: When total sign count / 2 is odd, negate y. Without this, some dependencies give incorrect x²≡y² congruences.
- **Integer overflow in parameter selection**: `(int)pow(L, 1.1)` overflows when L > 2³¹. Use `exp(1.1 * log(L))` with cap checks.
- **Degenerate extraction**: Some N values cause ALL dependencies to give trivial factors (x≡±y mod N). Knuth multiplier k·N fixes this by changing the algebraic structure.

### Batch GCD cofactor matching (siqs.c, experimental)
Implemented Bernstein product-tree batch GCD for finding shared factors among cofactors. The idea: instead of exact LP matching (same cofactor), find ANY common prime between cofactors of different partial relations. In O(n log²n) vs O(n²) for pairwise GCD. Allows much larger cofactors (up to B² instead of B) since we can decompose them after the batch GCD.

Status: infrastructure built, needs testing at scale to measure whether the extra relations from shared-factor matching actually improve throughput.

### Cofactor distribution study (cofactor_study.c)
Measured cofactor sizes and batch GCD collision rates across digit ranges (with optimal factor base sizes):

| Digits | FB size | Full smooth | 1LP    | 2LP     | Batch GCD collision |
|--------|---------|-------------|--------|---------|---------------------|
| 40     | 1854    | 0.004%      | 0.25%  | 21.5%   | 11.3%               |
| 50     | 4096*   | 0.001%      | 0.04%  | 3.8%    | 4.0%                |
| 60     | 4096*   | 0.0001%     | 0.003% | 0.57%   | 0.93%               |

*FB capped at 4096, suboptimal for 50-60 digit.

**Key finding**: Batch GCD collision rate drops from 11.3% → 4.0% → 0.93% as N grows. This confirms batch GCD matching gives constant-factor improvement in the 2LP regime but cannot improve the L-exponent. The 2LP pool is 5000x larger than the fully smooth pool, confirming that aggressive large prime handling is essential for practical performance, but this is already well-known.

### CCD factoring (ccd_factor.c)
QS-style sieve with batch cofactor matching. Critical lesson: partial relation matching must mark partners as "used" to prevent duplicate combined relations. Duplicates cause ALL GF(2) null vectors to be trivially degenerate (X ≡ ±Y mod N always).

## Open directions

These are starting points, not an exhaustive list.

- **Algebraic group structure**: Z_N* ≅ Z_{p-1} × Z_{q-1} but we can't see this decomposition. Can random walks, character sums, or higher-dimensional algebraic groups reveal it?
- **Novel norm structures**: NFS beats L[1/2] by splitting the norm across two maps, not by having a better single norm. What algebraic structure would enable a third splitting or a fundamentally different norm reduction to beat L[1/3]?
- **Function field analogies**: Can the quasi-polynomial DLP breakthrough technique be adapted? Key obstacle: no Frobenius over Z.
- **Batch GCD + non-sequential candidates**: Bernstein's batch smoothness works on arbitrary candidate sets. What candidate generation strategy (not sieving) would best exploit this?
- **Recursive cofactor descent**: CDS (extended matrix) failed because unique medium primes grow faster than relations. But LP-matching based descent is unexplored — could the collision rate be improved with structured cofactor generation?
- **Correlated norms across number fields**: In multi-image NFS, the norms from different fields are essentially independent. If we could design fields where smoothness of one norm CORRELATES with smoothness of another (e.g., via Galois structure or isogenies), the simultaneous smoothness probability would increase, potentially improving the exponent. No known construction achieves this.
- **High-dimensional lattice descent on cofactors**: For k cofactors c₁...cₖ, the lattice {(a₁,...,aₖ) : Σ aᵢcᵢ ≡ 0 (mod N)} has determinant N and shortest vector ~N^{1/(k+1)}. For k=3, vectors are ~N^{1/4}. If these short vectors could be used to build multiplicative (not just additive) cofactor relations, the effective smoothness bound would decrease. Obstacle: additive lattice relations don't directly give multiplicative relations needed for congruence-of-squares.
