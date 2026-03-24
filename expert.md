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
- **Quaternion / higher-dimensional norms**: Norm is still quadratic per variable — more variables don't reduce value sizes.
- **BatchQS (naive implementation)**: Bernstein product-tree batch GCD misses ~80% of smooth candidates vs standard sieve. The technique itself is sound but the implementation used insufficient primorial powers.

## Key insight: smoothness detection cost

Sieving: O(1) amortized per candidate but requires sequential memory access. Bernstein's batch GCD: O(log²B) per candidate but works on ARBITRARY candidate sets. For large B, batch GCD is cheaper and enables non-sequential candidate generation. This is unexploited — no implementation has combined batch GCD with a non-sequential candidate generation strategy effectively.

## Research survey

- **No classical sub-L[1/3] algorithm exists** for general integer factoring.
- **Regev (2023)**: Quantum O~(n^{3/2}) gates. No classical dequantization.
- **Stange (2022)**: Index calculus in (Z/NZ)*. Same L[1/2] or L[1/3].
- **Tower NFS (Barbulescu-Kim 2016)**: Sub-L[1/3] for DLP in GF(p^n) only. No factoring analog.
- **Umans & Wang (2025)**: Conditional deterministic N^{1/6+o(1)}. Exponential.
- **Cosset (2009)**: Genus-2 HECM on Kummer surfaces — 2 ECM curves simultaneously. Still L_p(1/2); only useful if constant improvement hints at a different scaling regime.

## Open directions

These are starting points, not an exhaustive list.

- **Algebraic group structure**: Z_N* ≅ Z_{p-1} × Z_{q-1} but we can't see this decomposition. Can random walks, character sums, or higher-dimensional algebraic groups reveal it?
- **Novel norm structures**: NFS beats L[1/2] by splitting the norm across two maps, not by having a better single norm. What algebraic structure would enable a third splitting or a fundamentally different norm reduction to beat L[1/3]?
- **Function field analogies**: Can the quasi-polynomial DLP breakthrough technique be adapted? Key obstacle: no Frobenius over Z.
- **Batch GCD + non-sequential candidates**: Bernstein's batch smoothness works on arbitrary candidate sets. What candidate generation strategy (not sieving) would best exploit this?
- **Recursive cofactor descent**: CDS (extended matrix) failed because unique medium primes grow faster than relations. But LP-matching based descent is unexplored — could the collision rate be improved with structured cofactor generation?
