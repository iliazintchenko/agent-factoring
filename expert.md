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

## Why L[1/3] appears to be a fundamental classical barrier

Every known approach falls into one of three categories:

1. **Smoothness-based** (QS, NFS): L-exponent locked by polynomial degree and Thue's theorem. LP expansion, multi-polynomial, batch GCD, multipliers, Galois action, CRT candidates, Best-of-K — all give at most constant-factor improvements.

2. **Group-order-based** (ECM, p-1/p+1, rho, class group): Complexity depends on factor size p. For balanced semiprimes (p ≈ √N), matches QS at L[1/2]. Division polynomial, Schoof-like, and Deuring correspondence approaches all reduce to this.

3. **Algebraic structure** (Hecke operators, quaternion algebras, isogeny graphs, lattice methods): Either require knowing the factorization to compute (circular), reduce to birthday bounds O(N^{1/4}), or reduce to smoothness methods.

The function field sieve breaks L[1/3] for DLP because of TWO properties absent over Z: (a) Frobenius endomorphism provides systematic degree reduction, (b) polynomials over finite fields factor efficiently (Berlekamp).

## Information-theoretic structure of factoring

Factoring an n-bit semiprime requires ~n/2 bits. Each smooth relation provides exactly **1 useful bit** (one GF(2) vector) despite carrying ~7.7 bits of "surprise" at u=4. Most information in the smooth factorization is wasted — only parity of exponents matters. The 1-bit-per-relation bottleneck is fundamental to the GF(2) linear algebra framework. Breaking it would require a new algebraic framework beyond exponent-parity relations.

Small-subgroup DL projection (g^{(N-1)/l} mod N for l | N-1): recovers O(1) bits total regardless of N's size, because gcd(N-1, φ(N)) is typically tiny. Strictly weaker than Pollard p-1.

Factorizations of N±k for small k: yield O(k·ln B) bits about p, but O(log N) needed. Information gap is exponential.

## Why L[1/3] IS the fundamental classical barrier

**The L-exponent formula**: α = 1/(k+1) where k = number of independent reduction stages in the algorithm.
- QS: k=1 (smoothness only) → α = 1/2
- NFS: k=2 (polynomial degree + smoothness) → α = 1/3
- FFS: k=3+ (+ tower descent via Frobenius) → α = 1/4 → 0

**Why NFS is stuck at k=2**: Adding a 3rd stage requires an iterable "compression map" Φ that reduces element size while preserving smoothness structure. Integer size is ARCHIMEDEAN — it cannot be iteratively reduced by algebraic maps. Polynomial degree is NON-ARCHIMEDEAN — Frobenius + quotient reduction decreases it freely. This asymmetry is the deep reason FFS reaches quasi-polynomial while NFS cannot break L[1/3].

**Constants**: GNFS c = (64/9)^{1/3} ≈ 1.923 (tight). MNFS floor: c = (32/9)^{1/3} ≈ 1.526 (tight, uniquely determined, no slack). No known optimization changes the leading constant for GNFS. All improvements (polynomial selection, lattice sieving, large primes, batch methods) affect only sub-leading terms.

All known sub-exponential factoring approaches are index-calculus methods hitting the same Dickman-de Bruijn tradeoff. Optimizing smoothness bound B, sieving range M, and relation count inherently yields α=1/3 through the three-way balance. This is robust across:

- Number field sieve and all known variants
- Any reduction to lattice problems (SVP requires dimension ≥ O(log N/log log N), encoding multiplicative structure)
- Group algebra / Hopf algebra decompositions (reduce to multi-polynomial NFS)
- Precomputation/non-uniform models (S·T ≥ L[2/3, c] conjectured, N-dependent sieve unavoidable)
- No known classical reduction to a problem with better algorithms

The quasi-polynomial DLP breakthrough in small characteristic does NOT transfer — it requires Frobenius endomorphisms and tower field structure absent over Z.

## What would be needed to beat L[1/3]

- A new algebraic structure over Z where "norms" are smaller than N^{2/3} per component
- OR smooth numbers of size < N^{1/3} that relate to N's factorization
- OR a completely non-smoothness-based sub-exponential approach (none known)
- OR a new algebraic framework beyond GF(2) exponent-parity that extracts >1 bit per relation
- OR a "descent" mechanism for Z: an endomorphism-like map that reduces "complexity" of elements recursively. Z has no non-trivial ring endomorphism (it's the initial object in Ring), which is WHY Frobenius descent fails over Z.

**Concrete target**: L[1/4] requires norm sizes N^{1/d²} instead of N^{1/d}. Verified by optimization: plugging N^{1/d²} into the NFS balance shifts cubic→quartic. Nested number field norms don't achieve this (composition = full norm, an invariant). This is the PRECISE mathematical barrier.

## Explored directions and conclusions

**Smoothness-based approaches explored:**
- **Schnorr lattice factoring**: Lattice dimension grows with smoothness bound. Ducas (CWI): 0 relations in 1000 trials.
- **LLL polynomial selection**: Degree-2 already optimal (shortest vector). Degree-3+ needs NFS.
- **Multi-image NFS**: Improves L[1/3] constants only, not the exponent.
- **Multi-multiplier sieve (MMS)**: ~60% more relations (constant factor) but doesn't shrink polynomial values.
- **Multiplicative Lattice Descent (MLD)**: LP expansion from B to B^α multiplies matrix dimension by ≥ R for R× more relations. L-exponent cannot improve — this is why LP variations improve constants but not the exponent.
- **Best-of-K polynomial selection**: Smoothness improvement scales as K^{u/ln(B)} ≈ K^{0.5}, always sublinear in K. Never beats sequential.
- **Smoothness correlations across polynomials**: No significant correlation found — sieve conditions at different x or different multipliers are independent.
- **CRT-structured candidates**: Higher raw smoothness rate is illusory — CRT forces large x, cofactor ≈ 2√N regardless.
- **Batch GCD collision rate**: Drops from 11.3% (40d) to 0.93% (60d). Constant-factor improvement only.
- **Multiplicative lattice relations (MLR)**: LLL finds small P(e) mod N, but CRT entanglement means small mod N ≠ small mod p. Confirms Schnorr is theoretically limited.

**Group-order approaches explored:**
- **ECM for balanced semiprimes**: L[1/2] in N. Not competitive above ~55 digits.
- **p-1/p+1 for balanced semiprimes**: Fail when p±1 aren't smooth.
- **CM-ECM**: Near-square group orders (#E ≈ (√p±1)²) but finding such curves costs O(√p). The penalty exactly cancels the smoothness advantage.
- **Iterated Frobenius Map** (x → x^N mod N): 100x per-step cost vs Pollard rho. Net slower.
- **Class group 2-Sylow**: O(N^{1/4}), same as Pollard rho.
- **Division polynomial factoring**: ECM with single-prime stage 1. Strictly worse than ECM.
- **Schoof-like torsion factoring**: Monic polynomial → no inversions → O(√N) per curve.

**Algebraic structure approaches explored:**
- **Spectral methods on Cayley graphs**: Exponential classically — exactly what Shor's QFT solves.
- **Deuring correspondence / isogeny graphs**: Three channels analyzed (supersingularity detection, spectral decomposition, random walk mixing). All blocked — require knowing p, exponential graph size, or O(N^{1/4}) birthday bound.
- **Hecke operators on modular forms**: Computing T_N requires knowing divisors of N. Circular.
- **Quaternion norms for smoothness**: Quadratic per variable.
- **Brandt matrices mod N**: For primes m < N, the quaternion norm form a²+b²+Nc²+Nd² forces c=d=0, making representation counts N-independent. The Brandt matrix is identical for all N — no factoring information leaks through spectral structure.
- **Character sum / autocorrelation**: Detecting period p requires Ω(√N) samples. Jacobi symbol autocorrelation tested extensively: semiprime vs prime statistics indistinguishable for M << p. The only "loophole" (small divisors of p-1) is exactly Pollard's p-1 method.
- **Polynomial splitting mod N**: Probability O(1/√N) per trial for monic polynomials.
- **Berlekamp over Z/NZ**: Degenerates to p-1 condition for balanced semiprimes — requires ord(r) | (q-1) in F_p, probability ≈ log(p)/p ≈ 0.
- **Approximate DL relations**: P(gcd reveals factor) ≈ 2/p. Exponential.
- **Bilinear smoothness decomposition**: Can't split norms without number fields.
- **Galois action on NFS relations**: Norm is Galois-invariant. Same relation, not new ones.
- **Sum-product phenomena over Z/NZ**: Zero divisors are too sparse (only ~2√N out of N). Detection requires k ~ √N samples — exponential in factor bit-length. E*(A)/E+(A) ratio indistinguishable from prime for N > 10^6.
- **Tensor decomposition of Z/NZ**: CRT gives hidden rank-2 structure, but mod-N wrapping destroys it. Signal-to-noise ratio ~1/√N. Catch-22: samples ≤ √N avoid wrapping but give rank-1 (trivial); samples > √N wrap and give full rank.
- **p-adic/Newton lifting over Z/NZ**: Three approaches tested (inversion failure, convergence divergence, resultants). All fail: (1) inversion failure probability ~1/p per step (exponential), (2) Newton iteration over finite fields is chaotic, not convergent — orbits don't reach roots for p > ~1000, (3) resultant approach reduces to Pollard p-1.
- **Groebner bases over Z/NZ**: Per-inversion failure rate is ~2/√N regardless of system size k. Buchberger coefficients behave as pseudorandom elements. Need ~√N operations to expect one non-invertible coefficient — no better than random GCD.
- **Batched Pollard p-1**: Batch product P=∏(a_i^k-1) mod N gives zero additional factoring power. By Fermat's little theorem, if (p-1)|k then ALL bases satisfy a^k≡1 mod p simultaneously; if not, none do. No intermediate case. Correlated bases (consecutive) also provide no advantage.
- **CF-base recursive descent**: Expressing integers in base r (from CF convergents of √N where r²≡s mod N) gives digit-products of size ~√N — identical to CFRAC. Rediscovers Morrison-Brillhart from a different angle. No Frobenius analog exists over Z to enable true descent.
- **Isogeny walks over Z/NZ**: Computing isogenies requires polynomial root-finding mod N, which IS the factoring leak (reduces to Berlekamp over Z/NZ). The Ramanujan graph structure is irrelevant — walk divergence is undetectable because isogeny computation itself reveals factors or doesn't. Confirmed: 0/140 balanced semiprimes factored.
- **Non-abelian HSP embeddings**: (Z/NZ)*/(Z/NZ)*² has order 4, structurally trivial. Hardness is in computing the hiding function (distinguishing cosets), not group theory. Faithful embeddings require exponential dimension. All re-embeddings encode factoring circularly.
- **Subfield lifting / tower extensions**: R = Z/NZ[x]/(f(x)) decomposes asymmetrically mod p vs q, but every computation operates on both components. Tower NFS works for DLP because the field tower is explicitly known; for factoring, the CRT decomposition IS the secret. Zero balanced semiprimes factored across 7 strategies.
- **Schnorr lattice tower / BKZ**: Minkowski bound gives α→0 as k→∞, but finding short vectors in dimension k takes 2^Ω(k) time. BKZ with full block size gives L[1, c] — worse than NFS. Polynomial-time reduction (LLL) cannot achieve α<1/2. The lattice geometric possibility is blocked by SVP hardness.
- **Three number field NFS**: All-3-smooth is 42% worse constant; any-2-of-3-smooth: 3× probability gain exactly cancelled by 3× larger matrix. Coppersmith MNFS improves constant toward 1.526 but α stays 1/3. The L[1/3] exponent is fundamental to the three-way d/B/M balance.
- **Batch GCD + non-sequential candidates**: Sequential sieving wins because it naturally produces the smallest candidates (a²-N ≈ 2a·δ near √N). Random/lattice candidates have larger values → exponentially lower smoothness. Batch GCD advantageous only for external candidate sources, not single-N factoring.
- **Reed-Solomon decoding over Z/NZ**: Each arithmetic operation has probability ~2/√N of hitting zero-divisor. RS algebraic structure does NOT increase this above random baseline. Need n≈N^{1/4} code length — still exponential.
- **Virtual Frobenius (x→x^N mod N)**: Acts as double Frobenius on CRT components. 100% success but O(√N) worst-case iterations, each costing O(log²N). Worse than Pollard rho. AKS-style polynomial ring variants degenerate for large N. No nontrivial field extensions over Z/NZ to exploit.
- **L-functions / class numbers**: Computing h(-N) requires O(√N) terms of L(1,χ_N). h(-N) ≠ h(-p)·h(-q) — no useful factorization. h(-N) mod l not determined by (p mod l, q mod l). Class group computation and factoring are computationally equivalent (McCurley).

## Key insight: smoothness detection cost

Sieving: O(1) amortized per candidate but requires sequential memory access. Bernstein's batch GCD: O(log²B) per candidate but works on ARBITRARY candidate sets. For large B, batch GCD is cheaper and enables non-sequential candidate generation. This is unexploited — no implementation has combined batch GCD with a non-sequential candidate generation strategy effectively.

## Research survey

- **No classical sub-L[1/3] algorithm exists** for general integer factoring (confirmed by 2020-2025 literature survey).
- **Henry Cohn (MIT)**: Argues no known barrier prevents progress below L[1/3], noting historical 1→1/2→1/3 progression. No concrete algorithm proposed.
- **Regev (2023)**: Quantum O~(n^{3/2}) gates. No classical dequantization.
- **Tower NFS (Barbulescu-Kim 2016)**: Sub-L[1/3] for DLP in GF(p^n) only. No factoring analog.
- **Harvey & Hittmeir (2020-2022)**: Best deterministic factoring N^{1/5+o(1)}. Exponential.
- **Regev (2023)**: Improved Shor's circuit to O(n^{3/2}) gates (from O(n²)). Follow-ups by Ragavan-Vaikuntanathan, Pilatte, Ekerå-Gärtner. Hhan (EUROCRYPT 2025) proved matching lower bound in generic ring model.
- **Schnorr lattice (2021 claim)**: Refuted experimentally by Ducas (CWI), 0/1000 relations. Never published at peer-reviewed venue.
- **No classical dequantization of Shor**: Unlike quantum ML algorithms, factoring's HSP structure resists classical simulation. The gap O(√N) vs O(polylog N) is robust.
- **GNFS L[1/3, 1.923] unchanged since 1990s**: RSA-250 (Feb 2020) remains factoring record. All progress is in practical constants (better sieving, batch smoothness).
- **NFS L[1/3] is HEURISTIC, not proven**: Best rigorous bound is L[1/2] (Lenstra-Pomerance 1992). Lee-Venkatesan (2017) proved L[1/3] only for randomized variant finding square congruences. Key unproven assumptions: smoothness heuristic ("norms behave like random integers"), independence of smoothness events, monogenic field, square root step.

## Open directions (genuinely unexplored)

After ~250 investigations, the remaining open space is extremely narrow. These are the directions that have NOT been conclusively eliminated:

- **Correlated norms across number fields**: Measured ~15% positive correlation at 40 digits, but deep analysis shows it's entirely SIZE BIAS (both norms are functions of the same (a,b) coordinates), NOT algebraic correlation. Effect weakens at larger N (20%→5-7%). Already fully exploited by lattice sieving. Factor bases are independent conditional on norm size. No residual algebraic correlation found. Cannot be amplified beyond current NFS techniques.
- **Non-smoothness-based sub-exponential relations**: Every known index-calculus method uses smooth number relations (1 bit per relation via GF(2) linear algebra). A fundamentally different relation type — extracting >1 bit per relation, or using a different algebraic structure than exponent-parity — could change the game. No candidate is known.
- **Computational number theory surprises**: The factoring landscape has been stable since 1993. But mathematical breakthroughs are inherently unpredictable. Possible sources: new connections to algebraic geometry (motives, étale cohomology), advances in additive combinatorics, or unexpected consequences of the Langlands program.
- **Approximate/relaxed algebraic structures**: Most approaches seek EXACT algebraic properties (exact smoothness, exact homomorphisms). Relaxed versions with controlled error might open new avenues. No concrete proposal exists.

## Explored directions (comprehensive dead-end catalog)

See above sections for details. Summary of ~280 investigated approaches across categories:

**Smoothness-based (all L[1/2] or L[1/3], no improvement):** Schnorr lattice, LLL polynomial selection, multi-image NFS, 3-field NFS, MMS, MLD, Best-of-K, CRT candidates, batch GCD, MLR, K-large-prime, CF-base descent, multi-d CF, sparse polynomial NFS, group algebra Z[C_k], lattice tower/BKZ, lattice on O_K, batch GCD + non-sequential, higher-order residues GF(k), smooth bit-patterns.

**Group-order-based (all L[1/2] in p or worse):** ECM, p-1/p+1, CM-ECM, iterated Frobenius, class group 2-Sylow, division polynomials, Schoof-like, algebraic tori T_n, genus-2 HECM, T_6/XTR compression, Weil/Tate pairings, dynamical systems (Chebyshev/Lattès), curve point-counting.

**Algebraic/structural (all blocked):** Spectral/Cayley, Deuring/isogeny, Hecke operators, Brandt matrices, character sums, polynomial splitting, Berlekamp, DL relations, bilinear decomposition, Galois action, sum-product, tensor decomposition, p-adic lifting, Groebner, batched p-1, Newton/Hensel, isogeny walks, non-abelian HSP, subfield lifting, virtual Frobenius, higher residue symbols, exterior algebra/SNF, counting arguments, modular forms/theta, carry propagation, algebraic varieties, RS codes over Z/NZ, synthetic Frobenius, tropical geometry, higher reciprocity laws, Cayley graph word metric, Kloosterman/exponential sums, Newton identities/power sums (s_k=p^k+q^k needs s_1=p+q = circular), zeta function of Z/NZ (evaluation at any s encodes factorization).

**Geometric/cohomological (all reduce to classical):** Arakelov geometry (product formula = same norms as NFS, arithmetic RR gives same counts), condensed mathematics (Z/NZ is discrete, all frameworks trivialize), prismatic cohomology (requires choosing p first = circular), derived algebraic geometry, Shimura varieties, motivic integration, weight 3/2 modular forms/Shimura correspondence (coefficients encode factoring info but extraction has interpretation/computation/decomposition barriers).

**Analytic/approximation (all vacuous or blocked):** abc conjecture (N=pq squarefree makes all bounds trivial), Baker theory (|log(p/q)| bound gives 0 search-space bits), stochastic resonance (smoothness has no threshold dynamics; noise hurts via Jensen), simulated annealing (cofactor landscape has no spatial correlation), Gauss sum spectra (|G(χ)|² reveals factors but requires O(N) character enumeration), wavelet/multiresolution analysis (detects structure but signal weak, equivalent to trial division), Dedekind zeta of Q(√N) (requires O(√N) terms).

**Information-theoretic/complexity (barriers confirmed):** 1-bit/relation bottleneck, small-subgroup DL, nearby factorizations, oracle classification, interactive proofs, ML prediction, Coppersmith + partial info, precomputation S·T bound, analytic function impossibility, classical QFT analog, period-finding dequantization, non-standard computation, hybrid quantum, Z-descent barrier, L[1/3] robustness, Kolmogorov complexity (K(p|N)=O(1) but K^poly(p|N) large iff factoring hard — restates problem), Pell/class groups (genus theory gives only O(1) bits, CF period indistinguishable SP vs prime), communication complexity (CC bounds O(log log N) only, Coppersmith n/4 certificate tight for lattice methods), entropy/MI analysis (each trial yields ~10^{-3} bits, MI/cost optimum matches NFS parameter selection), SDP/LP relaxations (weak due to global carry structure — Lasserre level O(1) insufficient), descriptive complexity (factoring definable in ESO∩USO but FO+LFP definition = P algorithm = open).

**Dynamical/ergodic (all equivalent to known methods):** Ergodic theory of multiplication maps (mixing/orbit properties indistinguishable SP vs prime for generic observables, symbolic dynamics also indistinguishable), proof/circuit complexity (factoring likely not in AC0/TC0 but natural proofs barrier blocks proving it, OBDD exponential for mult bits, circuit inversion hits treewidth barrier).

**Parameterized/structural (no improvement):** FPT with k=log(p) via ECM at exp(√(k log k)), Coppersmith FPT with gap parameter, NFS NOT FPT in factor size. Matroid theory of relations (standard GF(2) linear matroid, no useful second matroid). Additive combinatorics / sumsets (QR doubling constant encodes CRT but information flows factorization→structure only). AG codes (need base field, Z/NZ not a field). Waring / r_4(N) (computing divisor sums = factoring). Compressed sensing (2-sparse divisor indicator, O(log N) measurements suffice info-theoretically but computable measurements circular). Fixed-point / PPAD (factoring likely in PPAD/PPP, finding fixed points still hard). Topos theory (functoriality barrier, effective topos respects Turing complexity). FI-modules / representation stability (wrong algebraic side). Derandomization (NFS randomness not the barrier, smoothness heuristic is). Model theory / ultraproducts (zero divisors not FO-definable without naming). Diophantine geometry / rational points (genus 0, Faltings inapplicable, Chabauty needs genus>rank). Tensor networks (Shor has O(n) entanglement, MPS/PEPS insufficient). FHE analogy (CRT is natural encryption, bootstrapping=factoring). Spin glass (glassy landscape, BP diverges, SA fails above ~15 bits). Incidence geometry / ST over Z/NZ (no statistical difference SP vs prime). PIT (universal vs existential statement mismatch). Bio-inspired (flat fitness landscape, no gradient).
