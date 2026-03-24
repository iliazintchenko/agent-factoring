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

### Why L[1/3] appears to be a fundamental classical barrier (analysis summary)

Every known approach to integer factoring falls into one of three categories, each with a proven or conjectured barrier:

1. **Smoothness-based** (QS, NFS, CFRAC, Dixon): The L-exponent is determined by the ratio log(norm_value)/log(factor_base_bound). NFS achieves L[1/3] by splitting norms across two number fields. Tested improvements — LP expansion, multi-polynomial, batch GCD, multipliers, correlated norms, Galois action, CRT candidates, Best-of-K — all give at most constant-factor improvements. The exponent is locked by the polynomial degree and Thue's theorem (no small high-degree residues over Z).

2. **Group-order-based** (ECM, Pollard p-1/p+1, rho, class group methods): Complexity depends on factor size p, not N. For balanced semiprimes (p ≈ √N), these match QS at L[1/2]. ECM is the best in class. Division polynomial, Schoof-like, and Deuring correspondence approaches all reduce to this category.

3. **Algebraic structure** (Hecke operators, quaternion algebras, isogeny graphs, lattice methods): Either require knowing the factorization to compute (circular), or reduce to birthday bounds O(N^{1/4}), or reduce to smoothness methods. No non-smoothness channel has been found despite extensive search.

The function field sieve breaks L[1/3] for DLP because it exploits TWO properties absent over Z: (a) Frobenius endomorphism provides systematic degree reduction, and (b) polynomials over finite fields factor efficiently (Berlekamp). Without analogs of either property, integer factoring appears stuck at L[1/3].

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
- **Quaternion norms for smoothness testing**: Norm is still quadratic per variable — more variables don't reduce value sizes.
- **Deuring correspondence / isogeny graph factoring** (analyzed): Supersingular curves over F_p correspond to quaternion orders via Deuring. The ℓ-isogeny graph mod N is a product of graphs mod p and mod q. THREE approaches analyzed: (1) Detecting supersingularity mod p would reveal #E = p+1, enabling Pollard p+1 — but detecting supersingularity requires knowing p. (2) Spectral decomposition of the product graph would separate p and q components — but the graph has N/12 vertices, requiring exponential work. (3) Random walks to measure mixing time — gives O(N^{1/4}) birthday bound, same as Pollard rho. All three channels are blocked. The Deuring correspondence is NOT a viable non-smoothness-based approach to factoring.
- **BatchQS (naive implementation)**: Bernstein product-tree batch GCD misses ~80% of smooth candidates vs standard sieve. The technique itself is sound but the implementation used insufficient primorial powers.
- **Multiplicative Lattice Descent (MLD)**: Extending the LP bound from B to B^α (α>2) increases partial relation rate but the cofactor factor base grows as B^(α-1)/log(B). Tested for 30-digit N with α=3: MLD needs 2.16x MORE sieve work than standard QS because the matrix size explosion (100k × 100k) outweighs the partial relation gain (143x). The cofactor-to-matrix tradeoff is inherently bounded: any LP expansion that multiplies the relation rate by R also multiplies the matrix dimension by ≥ R, so the L-exponent cannot improve. This is the fundamental reason LP variations (single, double, triple) improve constants but not the exponent.
- **Class group 2-Sylow exploration**: For disc -4N, the class group has ambiguous forms including (p,0,q). Finding this via repeated squaring g^(2^k) gives order-2 elements, but most are "trivial" (a=2 type, from the discriminant factor 4). Needs BSGS or Shanks-style search to find the factoring form — O(N^{1/4}) time, same as Pollard rho. Not better than QS.
- **Quadratic character columns for extraction degeneracy**: QC columns in the GF(2) matrix DON'T help with degenerate extraction. Reason: QC bits are REDUNDANT — for any fully smooth relation, the product being a perfect square automatically satisfies all QC conditions. Tested with QC primes where kN is QNR mod q: all 128 dependencies still give trivial GCD. Also tried XOR combinations of dependency pairs — still degenerate. The problem is that the ENTIRE null space maps to x ≡ ±y (mod N), indicating the FB doesn't separate p and q for these specific N values. Known fix: Block Lanczos (different null-space sampling) or simply fall back to ECM.
- **Multi-base GCD accumulation**: Running Pollard p-1 with many random bases and B1 up to 50000 doesn't factor balanced semiprimes (p-1 not smooth). Multi-curve ECM accumulation with B1=1000 also fails for 30-digit semiprimes. Confirmed: group-order methods fundamentally fail for balanced random semiprimes.
- **Smoothness correlations across polynomial pairs**: Tested whether P1(x)=(x+m)²-N and P1(x+d) have correlated smoothness (would allow joint sieving). For 15-digit N with B=2000, measured conditional vs. unconditional smooth probability across shifts d=1..100 and cross-multiplier k=2. Result: no statistically significant correlation (ratios 0.85-1.63, all within 2σ of 1.0). Theoretical analysis confirms: for each prime p, the sieve conditions at x and x+d are independent (different residue classes), so joint smoothness equals the product of marginal probabilities. Cross-multiplier smoothness is also independent for the same reason. Constant-factor improvements from LP variations are possible but L-exponent is unaffected.
- **Schoof-like torsion factoring**: Compute x^N mod (x³+ax+b, N) for random curves E. Tested: works on small N (3, 19 digits) but FAILS on all 30-40 digit semiprimes. Since x³+ax+b is MONIC, no inversions during polynomial reduction — factors found only from final coefficient GCDs with probability ~1/p per curve. Total work O(√N) — same as trial division. NOT polynomial time despite initial appearance.
- **Division polynomial factoring**: Evaluate ψ_ℓ(P) for random points P on random curves, checking gcd(ψ_ℓ(P), N) for small primes ℓ. This is ECM with single-prime stage 1: tests if ℓ | #E(F_p) one prime at a time, vs ECM's simultaneous B1! test. Strictly worse than ECM because: ECM tests all primes up to B1 in one O(B1 log N) exponentiation, while divpoly tests each ℓ separately in O(ℓ²) per prime. Works on 3-digit N, fails on 30-digit.
- **Hecke operators on modular forms**: For Eisenstein series E_k, T_N(E_k) = σ_{k-1}(N) · E_k where σ_{k-1}(N) = (1+p^{k-1})(1+q^{k-1}). Knowing σ_1(N) = 1+p+q+N gives p+q and factors N. But computing T_N REQUIRES knowing divisors of N (for the Fourier coefficient recurrence), making it circular. Li (2025, ePrint 2025/1681) explicitly notes this: "Eisenstein series trivialize the Hecke problem into a factoring problem." Hecke operators REDUCE TO factoring, they don't solve it.
- **CM-ECM (near-perfect-square group orders)**: For CM disc D with |D|=O(1), #E ≈ (√p±1)², so smoothness depends on √p-sized factors (u halves). BUT finding curves with t ≈ ±2√p requires O(√p) trials (Sato-Tate). The √p penalty exceeds the smoothness advantage: √p × ρ(u/2) >> ρ(u). Still L_p[1/2].
- **QC columns for extraction degeneracy**: QC bits are REDUNDANT — perfect squares automatically satisfy QC conditions. Tested on 49-digit semiprimes: all 128 deps + XOR pairs degenerate. Entire null space maps to x ≡ ±y (mod N). Fix: Block Lanczos or ECM fallback.
- **Polynomial splitting mod composite N** (poly_split.c): For monic f(x) of degree d, compute gcd(x^{(N-1)/2} - 1, f(x)) mod N. If f splits differently mod p vs mod q, intermediate GCD steps encounter non-invertible leading coefficients. HOWEVER: for monic polynomials, probability of non-invertible coefficient is O(1/√N) per trial. 500 trials on 30-digit N found zero factors. Dead end.
- **CRT-structured candidate generation**: Force Q(x) to be divisible by chosen primes via CRT. Higher raw smoothness rate (0.59% vs 0.24% at B=50000) is illusory: CRT forces x to be large, so cofactor Q(x)/∏primes ≈ 2√N regardless. No asymptotic improvement.
- **Best-of-K polynomial selection** (surface_factor.c): Pick the smallest |Q_k(x)| from K polynomials per position. K=10 gives 1.92x rate improvement but 10x cost → 5.3x worse net. Improvement scales as K^{u/ln(B)} ≈ K^{0.5}, always sublinear. Never beats sequential.
- **CM curve ECM** (cm_ecm.c): For CM disc D, #E ≈ (√p ± 1)² — near-perfect square. Appears to improve smoothness since √p is smaller than p. BUT: ECM needs B1-**powersmooth** orders, and the squared exponents require primes ≤ √B1 instead of ≤ B1. This gives u_CM = ln(√p)/ln(√B1) = ln(p)/ln(B1) = u_standard. **Exactly zero advantage.** The near-square structure is perfectly canceled by the doubled prime power requirement.
- **Approximate DL relations** (approx_dl.c): Find e where 2^e mod N is close to a small prime. P(gcd reveals factor) ≈ 2/p ≈ 10^{-15} per hit. Exponential.
- **Galois action on NFS relations**: Algebraic norm is Galois-invariant (product of all conjugates). Automorphisms give same norm, not new relations. At best constant factor d improvement.

## Key insight: smoothness detection cost

Sieving: O(1) amortized per candidate but requires sequential memory access. Bernstein's batch GCD: O(log²B) per candidate but works on ARBITRARY candidate sets. For large B, batch GCD is cheaper and enables non-sequential candidate generation. This is unexploited — no implementation has combined batch GCD with a non-sequential candidate generation strategy effectively.

## Research survey

- **No classical sub-L[1/3] algorithm exists** for general integer factoring (confirmed by comprehensive 2020-2025 literature survey).
- **van Leeuwen (2020)**: Rigorous proof that randomized NFS runs in L[1/3] (previously only heuristic).
- **Harvey & Hittmeir (2020-2022)**: Best deterministic factoring improved from N^{1/4} to N^{1/5+o(1)}. Still exponential.
- **Boudot et al. (2022)**: Comprehensive survey confirming NFS at L[1/3, (64/9)^{1/3}] still SOTA. RSA-250 factored.
- **Bouillaguet et al. (2023)**: Alternative NFS sieving strategies give ~5% speedup. Constant factor only.
- **Henry Cohn (MIT)**: Argues no known barrier prevents progress below L[1/3], noting historical 1→1/2→1/3 progression. Conjectures L[ε] for all ε > 0 may be achievable. No concrete algorithm proposed.
- **Factoring via multiplicative relations mod n (2022)**: L[1/2] without NFS, L[1/3] with it. No improvement.
- **Integer factoring via CF and quadratic forms (2024)**: L[1/2] class. Worse than NFS.
- **Tower NFS**: Applied ONLY to DLP in GF(p^n) with composite n. NOT applicable to integer factoring. No analog exists.
- **Regev (2023)**: Quantum O~(n^{3/2}) gates. No classical dequantization.
- **Cosset (2009)**: Genus-2 HECM on Kummer surfaces. Still L_p(1/2).

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
| 54 | 80.5 | 73.1 |
| 58 | 244.1 | 133.1 |

Fitting log(time) vs digits suggests L[1/2] scaling as expected for QS variants.

### Known bugs / lessons
- **LP matching must divide by cofactor**: When combining two single-LP relations sharing cofactor c, the combined sqrt_val must be sv1·sv2·c⁻¹ (mod N), not just sv1·sv2. Otherwise x = y for all dependencies.
- **Sign handling in extraction**: When total sign count / 2 is odd, negate y. Without this, some dependencies give incorrect x²≡y² congruences.
- **Integer overflow in parameter selection**: `(int)pow(L, 1.1)` overflows when L > 2³¹. Use `exp(1.1 * log(L))` with cap checks.
- **Degenerate extraction**: Some N values cause ALL dependencies to give trivial factors (x≡±y mod N). Knuth multiplier k·N fixes this by changing the algebraic structure.

### MPQS extraction degeneracy: root cause analysis
Single-prime MPQS (A = q² for prime q) has a fundamental character diversity problem. For each polynomial, sv = (Ax+B)·q⁻¹ mod N. The Jacobi symbol (q⁻¹ / p) is fixed for all relations from the same q, and consecutive q values (which are close to √(√(2N)/M)) tend to share the same (q/p) character. This means the quadratic character vector lies in the row space of the exponent matrix with high probability, making ALL GF(2) null-space vectors produce trivial x ≡ ±y (mod N).

**Fix**: Multi-prime A (proper SIQS): A = q₁·q₂·...·qₖ where qᵢ are factor base primes. Different subsets give different A values with guaranteed diverse characters. Block Lanczos (instead of Gaussian elimination) also helps by sampling the null space more uniformly.

**Current workaround**: ECM fallback with optimal B1 for balanced semiprimes.

### ECM with optimal B1 for balanced semiprimes (ecm_recycle.c)
**KEY FINDING**: Pure ECM with B1 = exp(√(0.5·ln(N)·ln(ln(N))))/5 dramatically outperforms MPQS for balanced semiprimes:

| Digits | ECM-optimal worst | LGSH-v5 worst | Speedup |
|--------|-------------------|---------------|---------|
| 30 | 0.1s | 0.2s | 2x |
| 38 | 1.0s | 3.6s | 3.6x |
| 42 | 3.1s | 10.1s | 3.3x |
| 46 | 6.3s | 24.7s | 3.9x |
| 50 | 17.0s | 42.1s | 2.5x |
| 54 | 43.8s | 80.5s | 1.8x |
| 58 | 103.0s | 244.1s | 2.4x |

The key insight: for balanced semiprimes (p ≈ q ≈ √N), the optimal ECM B1 is proportional to L_p[1/2] ≈ L_N[1/4]. Starting B1 near optimal/5 and ramping up by 1.3x every 10 curves gives excellent performance. Most 50-digit semiprimes are factored on curve 0-2 (B1 already sufficient).

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

### NFS brute-force evaluation experiment (alg_sieve_exp.c)
Compared doubly-smooth yield of QS vs NFS (degree 3, 4, 5) for 40-digit N with B=27031:

| Method | Norms (bits) | Both-smooth rate | Ratio vs QS |
|--------|-------------|------------------|-------------|
| QS     | 81.8        | 0.273%           | 1.00        |
| NFS d=3| 40+59       | 0.022%           | 0.08        |
| NFS d=4| 32+53       | 0.289%           | 1.06        |
| NFS d=5| 27+50       | 0.717%           | **2.63**    |

**KEY RESULT**: NFS degree 5 produces 2.63x more doubly-smooth relations per candidate than QS at 40 digits, even with brute-force (a,b) enumeration (no lattice sieving). This is because the norms are substantially smaller: 27+50 = 77 total bits vs 82 bits for QS. The double-smoothness penalty is more than offset by smaller norms at higher degree. This advantage should GROW with N (since NFS is L[1/3] vs QS's L[1/2]). The practical challenge is implementing efficient line sieving or special-q sieving for the algebraic side. The 50+ digit experiments timed out with brute-force evaluation, confirming that some form of sieving is necessary.

### Hybrid QS-NFS (hybrid_qs_nfs.c)
QS-style 1D sieving on NFS polynomials with b=1,2,3,... For 30-digit N with degree 4: 54 doubly-smooth relations found in 120s (need 710). Rate ~10/s. The approach WORKS but the search space (b=1..5000 × a=-50000..50000) is exhausted before collecting enough relations.

**Analysis**: The NFS per-candidate advantage (2.63x at d=5) is offset by the SMALLER search space: algebraic norms grow as a^d, limiting A_MAX. For 30-digit, the two effects roughly cancel. For larger N, the NFS advantage grows (L[1/3] vs L[1/2] asymptotically).

**QS-NFS crossover**: Theoretical crossover at ~9 digits (from L-function analysis) but practical crossover at ~70-100 digits due to NFS overhead (double-smoothness, algebraic square root, polynomial selection).

### Candidate generation strategies (candidate_strategies.c)
Compared 5 strategies for 30-digit N, B=5000, 50K candidates:

| Strategy    | Smooth rate | Notes |
|-------------|------------|-------|
| Near-zero   | 0.022%    | Best: smallest |x²-N| |
| Sequential  | 0.018%    | Standard QS baseline |
| Multi-k     | 0.016%    | Knuth multipliers, comparable |
| Random      | 0.008%    | Worst: large random offsets |
| CRT-skip    | 0%/34     | Too few candidates (CRT forces large x) |

Confirms: candidates closest to √N (smallest values) are best. This is the foundation of QS. No candidate generation strategy can beat sequential-near-√N for QS-type polynomials.

### Combined approach performance (best_factor.sh)
Tested combined ECM→LGSH→HSD pipeline:

| Digits | Time (worst-of-5) | Method |
|--------|-------------------|--------|
| 30     | 0.02s            | ECM    |
| 40     | 0.68s            | ECM    |
| 50     | 19.0s            | ECM    |
| 55     | 36.7s            | ECM    |
| 60     | FAIL             | ECM timeout |
| 62     | 36s (optimized)  | ECM with B1 skip |
| 64     | FAIL             | ECM needs ~3300 curves at B1=43M, timeout at 295s |

**Gap at 60-70 digits**: ECM runs out of curves in time limit. MPQS degeneracy affects 25% of cases. NFS not competitive without lattice sieve. This 60-70 digit range is the hardest for current implementations.

## Open directions

These are starting points, not an exhaustive list.

- **Algebraic group structure**: Z_N* ≅ Z_{p-1} × Z_{q-1} but we can't see this decomposition. Can random walks, character sums, or higher-dimensional algebraic groups reveal it?
- **Novel norm structures**: NFS beats L[1/2] by splitting the norm across two maps, not by having a better single norm. What algebraic structure would enable a third splitting or a fundamentally different norm reduction to beat L[1/3]?
- **Function field analogies**: Can the quasi-polynomial DLP breakthrough technique be adapted? Key obstacle: no Frobenius over Z.
- **Batch GCD + non-sequential candidates**: Bernstein's batch smoothness works on arbitrary candidate sets. What candidate generation strategy (not sieving) would best exploit this?
- **Recursive cofactor descent**: CDS (extended matrix) failed because unique medium primes grow faster than relations. But LP-matching based descent is unexplored — could the collision rate be improved with structured cofactor generation?
- **Correlated norms across number fields**: In multi-image NFS, the norms from different fields are essentially independent. If we could design fields where smoothness of one norm CORRELATES with smoothness of another (e.g., via Galois structure or isogenies), the simultaneous smoothness probability would increase, potentially improving the exponent. No known construction achieves this.
- **Class group index calculus**: Hafner-McCurley (1989) showed class group computation for disc D takes L_D[1/2]. For D = -4N, this is L_N[1/2] — same exponent as QS. The class number formula h(-4N) = (2/π)√N · L(1, χ_{-4N}) connects to the factorization but computing h(-4N) requires O(√N) work. Real quadratic fields don't help either. Not an improvement path.
- **High-dimensional lattice descent on cofactors**: For k cofactors c₁...cₖ, the lattice {(a₁,...,aₖ) : Σ aᵢcᵢ ≡ 0 (mod N)} has determinant N and shortest vector ~N^{1/(k+1)}. For k=3, vectors are ~N^{1/4}. If these short vectors could be used to build multiplicative (not just additive) cofactor relations, the effective smoothness bound would decrease. Obstacle: additive lattice relations don't directly give multiplicative relations needed for congruence-of-squares.
