# Expert Knowledge

## The smoothness bottleneck

All known sub-exponential factoring algorithms rely on finding **smooth numbers**. The L-exponent is determined by the size of values tested for smoothness:
- QS: values ~√N → L[1/2]
- NFS: values ~N^(2/3) via number fields → L[1/3]
- Hypothetical L[1/4]: would need values ~N^(1/4)

To beat L[1/2], generate smaller candidates or avoid smoothness entirely.

## Thue's theorem and degree-2

CF of √N produces convergents with bounded residues: |p_k² - N·q_k²| < 2√N for all k. This infinite family of small residues makes QS possible.

For degree ≥ 3, Thue's theorem: |x^d - N·y^d| < C has **finitely many** solutions. No infinite family of small higher-degree residues exists over Z. Higher-degree polynomials require number fields (NFS).

## The L-exponent formula

For a sieve testing values of size N^β: complexity is **L[1/2, √(2β)]**.

| Approach | Value size (β) | Complexity |
|----------|---------------|------------|
| Dixon's  | 1             | L[1/2, √2] |
| QS       | 1/2           | L[1/2, 1]  |
| NFS      | ~2/3, split across two number fields | L[1/3] |

LLL confirms QS's degree-2 polynomial IS the shortest vector in the coefficient lattice. Degree 3+ requires number fields.

## Higher-dimensional algebras don't help

Z[√d]: norm a² - d·b², degree 2 per variable. Same value sizes as QS regardless of d. Can't beat L[1/2].

Quaternion algebras: norm a² + b² + c² + d², degree 2 per variable despite 4 variables. More variables don't reduce norm size.

NFS achieves L[1/3] by splitting across two norm maps (rational + algebraic) requiring simultaneous smoothness — not by having a better norm. Adding more number fields (MNFS) improves constants only.

## Three families of factoring algorithms

Every known approach falls into one of three categories:

1. **Smoothness-based** (QS, NFS): L-exponent locked by polynomial degree and Thue's theorem. All optimizations (LP expansion, multi-polynomial, batch GCD, multipliers, Galois action, CRT candidates, Best-of-K) give at most constant-factor improvements.

2. **Group-order-based** (ECM, p-1/p+1, rho, class group): Complexity depends on factor size p. For balanced semiprimes (p ≈ √N), matches QS at L[1/2].

3. **Algebraic structure** (Hecke operators, quaternion algebras, isogeny graphs, lattice methods): Either circular (require knowing factorization), reduce to birthday bounds O(N^{1/4}), or reduce to smoothness methods.

The function field sieve breaks L[1/3] for DLP via two properties absent over Z: (a) Frobenius endomorphism for systematic degree reduction, (b) polynomials over finite fields factor efficiently (Berlekamp). This does NOT transfer to integer factoring.

## Information-theoretic structure

Factoring an n-bit semiprime requires ~n/2 bits. Each smooth relation provides exactly **1 useful bit** (one GF(2) vector). The 1-bit-per-relation bottleneck is fundamental to the congruence-of-squares framework. Breaking it requires a new algebraic framework beyond exponent-parity.

Small-subgroup DL projection (g^{(N-1)/l} mod N): O(1) bits total regardless of N, because gcd(N-1, φ(N)) is typically tiny. Strictly weaker than Pollard p-1.

Factorizations of N±k for small k: O(k·ln B) bits about p, but O(log N) needed. Gap is exponential.

## Why NFS is stuck at k=2

**α = 1/(k+1)** where k = number of independent reduction stages:
- QS: k=1 → α = 1/2
- NFS: k=2 → α = 1/3
- FFS: k=3+ (tower descent via Frobenius) → α = 1/4 → 0

Adding a 3rd stage requires an iterable compression map that reduces element size while preserving smoothness structure. Integer size is ARCHIMEDEAN — cannot be iteratively reduced by algebraic maps. Polynomial degree is NON-ARCHIMEDEAN — Frobenius + quotient reduction decreases it freely. This asymmetry is the deep reason FFS reaches quasi-polynomial while NFS cannot break L[1/3].

Tower descent axiomatics confirms: FFS-style descent requires (A) non-trivial ring endomorphism (impossible: End_Ring(Z)={id,0}), (B) independent "degree" function additive on products (no such function), (C) Frobenius orbits on primes (primes algebraically independent), (D) smoothness testing without factoring (equivalent to factoring). In function fields, "degree" and "complexity" are INDEPENDENT parameters; in number fields they are the SAME (bit-length).

**Constants**: GNFS c = (64/9)^{1/3} ≈ 1.923 (tight). MNFS floor: c = (32/9)^{1/3} ≈ 1.526 (tight, uniquely determined). No known optimization changes the GNFS leading constant. All improvements affect sub-leading terms only. However, at cryptographic sizes (1024-4096 bits), constant improvements matter enormously — a factor-of-2 in c is worth many doublings of hardware.

All sub-exponential factoring approaches are index-calculus methods hitting the same Dickman-de Bruijn tradeoff. The three-way B/M/d balance inherently yields α=1/3. Robust across:
- NFS and all known variants
- Lattice reductions (SVP requires dimension ≥ O(log N/log log N))
- Group algebra / Hopf algebra decompositions (reduce to multi-polynomial NFS)
- Precomputation/non-uniform models (S·T ≥ L[2/3, c] conjectured)

**L[1/3] universality**: Any sieve-based algorithm with k parameters achieves L[1/(k+1)]. NFS has k=2 → 1/3. A third degree of freedom (k=3 → L[1/4]) requires a third independent smoothness source — all candidates systematically eliminated: MNFS improves c only, class/unit groups already in NFS, multiplicative orders wrong domain, ECM sub-dominant (L[1/6]), isogenies no smoothness connection. Only function fields have a third source (Frobenius). GNFS is conjectured asymptotically optimal within the smooth number paradigm.

## What would be needed to beat L[1/3]

- A new algebraic structure over Z where "norms" are smaller than N^{2/3} per component
- OR smooth numbers of size < N^{1/3} that relate to N's factorization
- OR a completely non-smoothness-based sub-exponential approach (none known)
- OR a new algebraic framework beyond GF(2) exponent-parity that extracts >1 bit per relation
- OR a "descent" mechanism for Z: an endomorphism-like map that reduces "complexity" recursively. Z has no non-trivial ring endomorphism (initial object in Ring), which is WHY Frobenius descent fails over Z.

**Concrete target**: L[1/4] requires norm sizes N^{1/d²} instead of N^{1/d}. Nested number field norms don't achieve this (composition = full norm, an invariant).

**Barrier structure**: The five barriers decompose into approach-dependent (L[1/3] wall, archimedean gap, GF(2) bottleneck — could be evaded by novel framework) and structural (CRT opacity, Z-rigidity — any classical algorithm must confront). Minimal barrier set: {CRT opacity, Z-rigidity} — the others are consequences within sieve methods.

**No one has proved L[1/3] is a hard lower bound.** The barrier is empirical — every known approach hits it, but every *possible* approach need not. The progression L[1] → L[1/2] → L[1/3] took decades per step; each breakthrough came from an unexpected direction.

## Probability estimate

**P(classical algorithm ever beats L[1/3]) ≈ 2-3%** (80% CI [0.5%, 8%]). Updated after ~530 systematic investigations.

Why not 0%:
- No proof that L[1/3] is a hard lower bound — the barrier is empirical, not information-theoretic
- The function field DLP breakthrough shows analogous barriers CAN be broken
- L[1/3] is a THEOREM within the sieve framework (Dickman function + 2-parameter optimization), NOT a proven lower bound for factoring in general
- Mathematics has surprised us before (AKS, Babai GI, fast matrix multiplication)

Why 2-3% (lowered from initial estimate):
- ~530 systematic investigations across every plausible mathematical avenue — NONE moved the exponent
- L[1/3] emerges independently in NFS, FFS analogues, and every index-calculus variant — it is a universal property
- ~35 years of effort by hundreds of mathematicians have improved only c, never the exponent 1/3
- Every approach that looked promising hit the Dickman-function wall when pushed to the general case
- The five structural barriers are deeply interrelated and mutually reinforcing
- Complete absence of "near miss" algorithms that would signal a breakthrough is close

## Research survey

- **No classical sub-L[1/3] algorithm exists** for general integer factoring (2020-2025 literature survey).
- **Henry Cohn (MIT)**: No known barrier prevents progress below L[1/3]. No concrete algorithm proposed.
- **Regev (2023)**: Quantum O~(n^{3/2}) gates (from Shor's O(n²)). Hhan (EUROCRYPT 2025) proved matching lower bound in generic ring model. No classical dequantization — the gap O(√N) vs O(polylog N) is robust.
- **Tower NFS (Barbulescu-Kim 2016)**: Sub-L[1/3] for DLP in GF(p^n) only. No factoring analog.
- **Harvey & Hittmeir (2020-2022)**: Best deterministic factoring N^{1/5+o(1)}. Exponential.
- **Schnorr lattice (2021 claim)**: Refuted by Ducas (CWI), 0/1000 relations. Never peer-reviewed.
- **GNFS L[1/3, 1.923] unchanged since 1990s**: RSA-250 (Feb 2020) remains record. All progress in practical constants.
- **NFS L[1/3] is HEURISTIC, not proven**: Best rigorous bound is L[1/2] (Lenstra-Pomerance 1992). Lee-Venkatesan (2017) proved L[1/3] only for finding square congruences. Key unproven assumptions: smoothness heuristic, independence of smoothness events, monogenic field, square root step.
- **2024-2026 literature survey**: No algorithm has improved the classical L-exponent since GNFS (early 1990s). Every sub-L[1/3] claim debunked or retracted (Schnorr 2021, Yan et al., Chen quantum LWE, JVG quantum). Notable new work: Stange 2022 (index calculus bridge), Gao et al. 2025 (rank-3 lattice second vector — elegant but no path to sub-L[1/3]), Gorodetsky 2023 (Dickman phase transition at y=(log x)^{3/2}), Pascadi 2025 (smooth number equidistribution to moduli x^{5/8}), Tao 2025 (max-entropy framework for smooth number anatomy), Friedlander-Granville smooth number improvements, Pinnacle quantum architecture. None change the L-exponent.
- **Rigorous L[1/3] gap**: Lee-Venkatesan proved L[1/3] for sieving phase. Remaining gaps: (1) algebraic square root extraction in L[1/3] time, (2) independence of smoothness events, (3) class group/unit group computation. Most promising path: complete Lee-Venkatesan by proving congruence-to-factor extraction in L[1/3] time. Medium-term: extend smooth number equidistribution to polynomial values. Deep obstruction: parity problem in sieve theory.
- **Regev lattice classically intractable**: Dimension O(√n), indistinguishable from random q-ary lattices by BKZ. Classical cost 2^{Θ(√n)} strictly worse than GNFS.
- **GF(2) reduction is optimal**: Z-lattice dependencies are a SUBSET of GF(2) dependencies. SNF invariant factors all 1 or 2. Mod-2 reduction RELAXES constraints to increase the solution space. 1-bit/relation is a feature within congruence-of-squares.

## Explored directions

~530 approaches investigated. None improved the L-exponent.

### Smoothness-based (all L[1/2] or L[1/3])

- **Schnorr lattice**: Lattice dimension grows with smoothness bound. Ducas: 0 relations in 1000 trials.
- **Lattice crypto techniques (BDD/dual/hybrid/primal)**: Kannan embedding already in NFS. No useful BDD formulation (Coppersmith needs partial info). Dual attack degenerates to Gaussian elimination (no noise unlike LWE). Hybrid BKZ+MITM adds cost at GNFS saddle point. Fundamental: p is not a short vector in any known poly-dim lattice; factoring needs π(B) relations not one vector.
- **LLL polynomial selection**: Degree-2 already optimal (shortest vector). Degree-3+ needs NFS.
- **Multi-image NFS**: Improves L[1/3] constants only.
- **Degree-4+ norms**: LLL gives α=0.45 (d=4), 0.31 (d=6), below N^{2/3}. But N^{2/3} is specific to degree-3 base-m; L[1/3, c] preserved because NFS optimization converges regardless of d. LLL improves c by ~5-10%. Cyclotomic/thin fields give smaller norms but need special-form N (SNFS).
- **Multi-multiplier sieve (MMS)**: ~60% more relations (constant factor), doesn't shrink polynomial values.
- **Multiplicative Lattice Descent (MLD)**: LP expansion multiplies matrix dimension proportionally. L-exponent cannot improve.
- **Best-of-K polynomial selection**: Smoothness improvement scales as K^{0.5}, always sublinear. Never beats sequential.
- **Smoothness correlations across polynomials**: None found — sieve conditions at different x or multipliers are independent.
- **CRT-structured candidates**: Higher raw smoothness rate illusory — CRT forces large x, cofactor ≈ 2√N.
- **Algebraic curve families** (Pell CF, Cornacchia, Lehman): Cornacchia shows 50-500x raw smoothness BUT rank only 20-190 out of 1229 needed — relations cluster on few primes. Pell values small but sparse (O(log N) convergents). Lehman = MPQS. No family achieves both high rate AND full-rank relations.
- **Batch GCD collision rate**: 11.3% at 40d → 0.93% at 60d. Constant-factor only. Sequential sieving wins because it produces the smallest candidates (near √N). Batch GCD advantageous only for external candidate sources.
- **Multiplicative lattice relations (MLR)**: LLL finds small P(e) mod N, but CRT entanglement: small mod N ≠ small mod p. Confirms Schnorr limited.
- **Correlated norms across number fields**: ~15% positive correlation at 40 digits — entirely SIZE BIAS (shared (a,b) coordinates), not algebraic. Weakens at larger N. Already exploited by lattice sieving.
- **Three number field NFS**: All-3-smooth 42% worse constant; any-2-of-3: 3× probability cancelled by 3× matrix. MNFS improves toward c=1.526, α stays 1/3.
- **K-large-prime**: Degenerates to L-smooth, same L[1/2]. Practical 5-10x only.
- **CF-base recursive descent**: Digit-products ~√N, identical to CFRAC. No Frobenius analog over Z.
- **Multi-d CF / cube root CF**: Degree-3 Thue barrier confirmed. Multi-d = MPQS in disguise.
- **Sparse polynomial NFS**: Fewer roots hurts sieve. Sparsity saves only a constant. Kleinjung already optimal.
- **Group algebra Z[C_k]**: Equivalent to MPNFS with cyclotomic fields.
- **Lattice on O_K**: Norms ~N^{1/2} (exponential), only (ln N)^{1/3} vectors. Rank-2 wins.
- **Schnorr lattice tower / BKZ**: Minkowski bound gives α→0 as k→∞, but SVP in dimension k takes 2^Ω(k). BKZ → L[1, c] (worse than NFS). LLL cannot achieve α<1/2.
- **Higher-order residues GF(k)**: Factor base shrinks by 1/(k-1), info/relation peaks at k=2. GF(2) optimal.
- **k-th power congruences (GF(k) LA)**: GF(k) LA still needs π(B)+1 smooth relations regardless of k. GF(k) kernel has same dimension as GF(2); extra bits/relation don't reduce relation count. No improvement over k=2.
- **Non-maximal order norms Z[cα]**: Conductor c inflates norms. Z[cα] ⊂ O_K means fewer elements, not more. Strictly worse than standard NFS.
- **Quaternion algebra reduced norms for sieving**: Degree-2 in 4 variables. Encoding factoring info (norm ≡ 0 mod N) forces norms ≥ N. Reduces to QS-class L[1/2].
- **Resultant smoothness**: Smoothness probability is the PRODUCT (not sum) of individual smoothness probabilities, making resultants strictly WORSE than norms directly.
- **Smooth number cascade/avalanche structure**: Full exponent vectors carry more info than GF(2) but extra info doesn't help factoring. The 1-bit-per-relation bound is AIRTIGHT within congruence-of-squares.
- **Smooth bit-patterns**: Boolean classifier finds shallow artifacts (2nd MSB bias) that diminish at crypto sizes.
- **LLL on cofactor log-lattice**: Smoothness is discrete; continuous log-lattice optimization targets wrong objective. 0 smooth combinations found.
- **Cofactor correlation in QS**: Cofactors at nearby positions completely independent. GCD rate = 0. Log-size correlation r < 0.02.
- **Cross-polynomial smoothness correlation (MPQS)**: Independent conditioned on x mod p. Apparent correlation is purely size effect.
- **Heavy-tail analysis**: Smoothness matches Dickman (ratio ≈ 1.08). Cofactor tail exponential not power-law.
- **Smooth value structure mining**: Gaps follow Poisson. APs enriched 1.5-2.5x but IS the sieve. All structure identical for SP vs prime.
- **QS cofactor product smoothness**: For 1-LP QS, cofactors are prime by necessity; products not smooth.
- **QS cofactor distribution**: Non-uniform (98.4% chi²) BUT entirely explained by public polynomial structure. 0/400 peaks match p mod l.
- **Additive structure of smooth numbers**: Tautology barrier (a + b ≡ 0 mod p vacuous when p|N). σ(N) = (1+p)(1+q) computing IS factoring. Additive energy indistinguishable SP vs prime.
- **Super-Dickman sieving regions**: Super-smooth windows exist (C ≈ 2-12x Dickman) and are findable, but: QS sieve already exploits this, enhancement bounded constant, large-prime divisibility structureless. Cannot improve L-exponent.
- **Product lattice / geometry of numbers**: L_rel dimension π(B) — too large for LLL. Bootstrapping paradox: constructing L_rel requires the smooth relations. NFS already optimally deploys lattice methods.
- **Sieve in alternative number systems**: Gaussian Z[i], Eisenstein Z[ω], Hurwitz quaternions. All produce norms ≥ NFS norms. NFS polynomial selection IS the optimal sieve across all known number systems.
- **ML for sieve acceleration**: AUC 0.80-0.87 with full features but features cost as much as the sieve itself. The sieve IS the optimal ML model for smoothness detection.
- **Lattice reductions of factoring (survey)**: Coppersmith needs n/4 bits of p. Schnorr debunked. Regev quantum-only. Classical lattice factoring needs dim ≈ n^{1/3} but no known construction without quantum help.
- **Coppersmith without partial info**: 0 bits recovered across all tested dimensions. n/4 threshold is sharp.
- **Multidimensional CF (Jacobi-Perron)**: Produces LARGER norms than standard CFRAC. Complex conjugates ruin the degree-2 estimate.
- **Selberg sieve for smoothness**: Optimized for EXCLUDING primes. Running it "in reverse" doesn't improve on Dickman.
- **Stange index calculus bridge (arXiv:2211.06821)**: Clean framework unifying factoring and DLP. Complexity matchable to L[1/3] with worse constants. No new algebraic structure.
- **Special Galois groups for NFS**: Galois structure does NOT change L-exponent — only affects constant c.
- **DGS for NFS square root**: Doesn't improve asymptotic because bottleneck is class group computation, not lattice sampling.
- **Simultaneous diophantine approximation (LLL)**: Norm sizes match or exceed standard NFS.
- **Recursive lattice descent for NFS**: Three failure modes: log-to-multiplicative gap, bootstrap impossibility, smoothness bottleneck. Quantitative: 2^{157} WORSE than NFS at 1024-bit.
- **PIT / derandomized smooth detection**: Smoothness polynomials have MAXIMAL circuit complexity. B-smooth integers don't form algebraic subgroups. Sieving IS essentially optimal.
- **Smooth number construction via algebraic dependencies**: Products of known relations are DEPENDENT. Independence barrier: new evaluations revert to Dickman rates.
- **Short interval smooth numbers / analytic NT**: Smooth values follow Poisson process. The sieve ALREADY fully exploits CRT structure. CADO-NFS, msieve already embody all relevant insights.
- **Number field tower descent**: Tower descent produces NO new independent relations. Size reduction bounded by extension degree (constant), unlike FFS.
- **Third homomorphism / third NFS channel**: Systematic analysis of all 5 candidates: (1) second number field = MNFS, constant only; (2) p-adic = circular; (3) function field = circular; (4) adelic = product formula constrains; (5) modular forms = circular. Z provably cannot support a third independent smoothness channel due to product formula and characteristic barrier.
- **Sparse polynomials / lacunary series**: Carry propagation destroys factor sparsity. Cyclotomic factoring reduces to Pollard p-1.
- **Algebraic statistics of sieve matrices**: BP fails (short cycles). ML/toric warm starts fail. O(n) bits of summary cannot determine O(d·n) bits of null space.
- **Gorodetsky phase transition at y=(log x)^{3/2}**: NFS operates FAR above threshold. Moving toward it catastrophically increases sieve cost. Tao max-entropy validates NFS optimizations.
- **Cohomological obstructions for smooth numbers**: Class number h correlates with sieve yield but with CONSTANTS not EXPONENTS. Smoothness is not a local-to-global property.
- **Height theory / Arakelov for NFS elements**: Height spectrum determined by |N|, not factorization. LLL/BKZ already approaches Arakelov bound. CONFIRMS NFS optimality.

### Group-order-based (all L[1/2] in p or worse)

- **ECM for balanced semiprimes**: L[1/2] in N. Wall at ~55-56 digits (28-digit factors with B1=10^7).
- **p-1/p+1**: Fail when p±1 aren't smooth.
- **CM-ECM**: Near-square group orders but finding such curves costs O(√p). Penalty cancels advantage.
- **Iterated Frobenius** (x → x^N mod N): 100x per-step cost vs Pollard rho. Net slower.
- **Class group 2-Sylow**: O(N^{1/4}), same as Pollard rho.
- **Division polynomial factoring**: ECM with single-prime stage 1. Strictly worse than ECM.
- **Schoof-like torsion**: Monic polynomial → no inversions → O(√N) per curve.
- **Algebraic tori T_n** (n∈{3,4,5,6,8,10,12}): Independent smoothness channels (Φ_n(p)), 2/30 unique wins. Real but limited.
- **Genus-2 HECM**: Smoothness advantage for u<5 but 4-5x cost neutralizes. No asymptotic gain.
- **T_6/XTR compression**: ρ(2u)/ρ(u) collapses, 34-89x slower ops, 0% smoothness at B=1000 for 30-bit primes.
- **Weil/Tate pairings**: 0/60 pairing-specific factorizations, all ECM-equivalent.
- **Dynamical systems (Chebyshev/Lattès)**: No improvement over Pollard rho. Combined p-1/p+1 real but limited.
- **Curve point-counting**: All genera degenerate to ECM. O(N^{3/4}) inversions = ECM exactly.
- **GL(2, Z/NZ)**: Element orders combine p-1 and p+1 channels. Same L[1/2] barrier.
- **Hilbert class polynomials** (h(D)=1,2,3,4,5): 0% success across 125+ pairs. Reduces to Pollard p+1.
- **CM endomorphisms**: Decomposes into known methods — ζ test = Solovay-Strassen, supersingular = Williams p+1, GLV already in GMP-ECM.
- **Formal group laws for EC factoring**: [l]-endomorphism reduces to ECM-like operations. Formal logarithm convergence radius depends on p but detecting it requires knowing p.
- **Special semiprimes bootstrap**: Pr[p-1 is L_p[1/3,1]-smooth] ≈ 10^{-12} to 10^{-34} at crypto sizes. ECM already optimally exploits all group-order special cases.
- **Random multi-curve point counts**: CRT multiplicativity holds perfectly. Reduces to running ECM on each curve independently.
- **EC L-function analytic rank**: L(E,1) computation requires Euler product factorization at p — circular.
- **Analytic rank jumps of EC families**: No statistically significant SP vs prime differences on any non-circular metric.
- **Genus-g Jacobian ECM scaling**: Efficiency PEAKS AT g=1. L-exponent scales as √g. Lenstra's 1987 choice was optimal.
- **Quadratic form / Bhargava composition**: Ambiguous forms encode factoring form (p,0,q) but finding it requires ~√N forms. Class group infrastructure = SQUFOF at O(N^{1/4}).
- **Analytic rank of quadratic twists E^{(N)}**: Root numbers and analytic ranks are functions of N ALONE. Decomposition chi_N = chi_p·chi_q requires knowing p,q.
- **Sato-Tate distribution over Z/NZ**: Point counts decompose via CRT. Normalized traces violate Hasse bound for semiprimes — trivially detectable but this is primality testing, not factoring.

### Algebraic/structural (all blocked)

- **Spectral methods on Cayley graphs**: Exponential classically — what Shor's QFT solves.
- **Deuring correspondence / isogeny graphs**: Three channels (supersingularity, spectral, mixing). All blocked — require knowing p, exponential graph, or O(N^{1/4}) birthday.
- **Hecke operators on modular forms**: Computing T_N requires knowing divisors. Circular.
- **Quaternion norms**: Degree 2 per variable. Smoothness only 6% above random for matched magnitude.
- **Non-commutative ring norms (M_2(Z), quaternions, Z[S_3])**: Apparent smoothness advantage (25-37%) is size bias. Any multiplicative norm N:R→Z projects factoring in R to factoring in Z. Every polynomial-time computable invariant of matrices over Z/NZ is a polynomial function of entries, which respects CRT. Non-decomposing invariants (orders, word length) are computationally intractable. No "moderately hard" middle ground.
- **Brandt matrices mod N**: Norm form a²+b²+Nc²+Nd² forces c=d=0 for primes m<N. Matrix is N-independent.
- **Character sum / autocorrelation**: Period p requires Ω(√N) samples. Only loophole = Pollard p-1.
- **Jacobi symbol deconvolution / ICA**: Binary ICA on ±1 signals info-theoretically impossible. All paths require Ω(√N).
- **Polynomial splitting mod N**: O(1/√N) per trial for monic polynomials.
- **Berlekamp over Z/NZ**: Probability ≈ log(p)/p ≈ 0. 0/3600 factorizations.
- **Approximate DL relations**: P(gcd reveals factor) ≈ 2/p. Exponential.
- **Galois action on NFS relations**: Norm is Galois-invariant. Same relation, not new ones.
- **Sum-product over Z/NZ**: Zero divisors too sparse (~2√N out of N). CRT isomorphism preserves sum-product ratios (ring-isomorphism invariant). Extended to all orbit types — zero separation.
- **Tensor decomposition of Z/NZ**: CRT rank-2 structure IS detectable but mod-N wrapping destroys it. SNR ~1/√N.
- **p-adic/Newton lifting**: (1) inversion failure ~1/p, (2) Newton chaotic for p>~1000, (3) resultant = Pollard p-1.
- **Groebner bases over Z/NZ**: Per-inversion failure ~2/√N. Need ~√N operations — no better than random GCD.
- **Batched Pollard p-1**: By FLT, all bases succeed or fail simultaneously. Zero batch advantage.
- **Isogeny walks over Z/NZ**: Computing isogenies requires polynomial root-finding mod N = the factoring leak. 0/140 factored.
- **Non-abelian HSP embeddings**: Faithful embeddings require exponential dimension. Computing hiding function IS factoring.
- **Subfield lifting / tower extensions**: CRT decomposition IS the secret. 0 balanced semiprimes factored across 7 strategies.
- **Reed-Solomon over Z/NZ**: Per-operation failure ~2/√N (random baseline). Need n≈N^{1/4} code length.
- **Virtual Frobenius (x→x^N)**: O(√N) worst-case iterations. Worse than Pollard rho.
- **L-functions / class numbers**: h(-N) requires O(√N) terms. Class group computation ≡ factoring (McCurley).
- **Pizer Ramanujan graphs via supersingular j-invariants**: 6 approaches tested. Approaches 1-5: O(1/√N) success. Approach 6 IS Pollard rho. Isogeny differences real but computationally inaccessible.
- **Additive structure of power residues via CRT**: Reconstructing CRT decomposition requires Jacobi sums mod N = factoring.
- **Eigenvalue structure of matrices over Z/NZ**: Eigenvalues leak via CRT but computing them → p-1 condition.
- **GF(3)/ternary relations**: QS produces quadratic relations; GF(3) kernel gives mixed a²≡b³ (unusable). GF(k) for k>2 strictly worse.
- **Higher residue symbols / reciprocity (cubic/quartic)**: 0/500+ factored. Quadratic reciprocity computable without factoring (Jacobi symbol); higher reciprocity CANNOT be computed mod N without factors — norm map involves factorization. Unique computational gift of degree 2. Multi-curve EC: underdetermination GROWTH theorem (more curves make it harder). Sato-Tate product moments depend only on N.
- **Bilinear smoothness decomposition**: Can't split norms without number fields.
- **Counting arguments** (r₂, r₄, class number, partitions): Factor-independent counts reveal nothing; factor-dependent costs ≥ factoring.
- **Modular forms / theta / tau**: r₄(N)=8(1+p)(1+q) elegant but requires O(N) work.
- **Newton identities / power sums**: s_k=p^k+q^k needs s_1=p+q = circular.
- **Zeta function of Z/NZ**: Z_N(1)=(p+1)(q+1)/N. Evaluation at any s encodes factorization.
- **Carry propagation**: Branching factor 2.0, ~1.25×2^{n/2} paths. MQ-hard mod 2.
- **MR witness fingerprinting**: Extended MR chain IS Pollard p-1 restricted to 2-part. General factoring requires probing ALL prime components.
- **Randomized ring embeddings** (GL(2), polynomial quotient, tensor): GL(2) = p±1/Williams (1982). Polynomial splitting = simplified NFS.
- **CF structure of √N**: Partial quotients follow Gauss-Kuzmin for both SP and prime (indistinguishable). No shortcut beyond CFRAC.
- **Multiplicative structure of small multiples kN**: gcd(kN+m, N)=gcd(m,N), multiplier irrelevant.
- **NFS iterated descent obstruction**: FFS descent uses 3 absent properties: (1) Frobenius, (2) norm linear in degree, (3) PID.
- **Batch GCD + non-sequential candidates**: Sequential sieving wins — smallest candidates, locality advantage.
- **Approximate Frobenius x→x^l mod N**: Orbit structure decomposes via CRT. Witt-Frobenius reduces to independent power maps per component. Frobenius error only useful when l IS a factor (circular).
- **Adaptive structured sampling for CRT probing**: 11 strategies tested. NO strategy accumulates CRT info faster than random. MI per sample ~0.04 bits regardless.
- **p-adic descent for factoring**: All info about p lives at archimedean place and place l=p, neither accessible without knowing p. No "partial Frobenius" exists.
- **CRT tree search**: Branching factor multiplicative across primes. Reduces to exhaustive search of O(√N).
- **Matrix Frobenius normal form gap**: Factor leakage during Gaussian elimination = random GCD at O(1/p). No structural advantage.
- **Representation theory / partial Gauss sums**: Partial Gauss sums DO encode factor info but extracting requires knowing which characters are p-only — circular. ANOVA shows NO useful signal.
- **Number wall / Padé approximants**: Zero patterns encode factoring info but only in rows ≥ p+1 (O(p) barrier).
- **Floor-division sequence discrepancy**: Perfectly structured when m|N, quasi-random otherwise — no intermediate regime.
- **Power residue codes / QR codes mod N**: d_min directly encodes factoring info but computing it requires knowing p,q.
- **Adelic methods**: Product formula ensures gains at one place offset by losses at others. Hardness persists in full adelic setting.
- **Expander graph combinatorics**: GLOBAL properties encode factoring but need O(N) computation. LOCAL properties reduce to Pollard rho/p-1/ECM. No middle ground.
- **Additive-multiplicative relations (a·b+c≡0 mod N)**: Smooth c is size ~N (worse than QS's ~√N). Ring structure fully mined by existing methods.
- **CRT descent with inequality/primality pruning**: Cumulative pruning only 1/ln(B) — logarithmic, not geometric. NO asymptotic improvement.
- **Randomized PIT for factoring**: O(N/d) evaluations to detect — exponential.
- **Non-abelian CFT**: Representation-theoretic structure faces same multiplicative mixing barrier. Langlands doesn't compute more efficiently.
- **Multilinear maps**: Even ideal MLMs would give Shor's structure, confirming quantum advantage.
- **Categorical Galois theory (Janelidze)**: Descent data = idempotents = factoring.
- **CFSG / group structure of (Z/NZ)***: All decompose via CRT. The Jacobi symbol's "shadow has exactly one dimension fewer, and that missing dimension IS the factorization."
- **Fourier analysis on (Z/NZ)***: Characters require knowing φ(N) = factoring. Partial Fourier carries O(1/√N) signal.
- **Hasse-Minkowski / local-to-global**: xy=N is genus 0 — Hasse principle holds trivially.
- **Smooth number graph spectral theory**: Good expansion (consistent with random), no SP/prime anomaly.
- **F_1 geometry / Frobenius homomorphism test**: Birthday collision on Frobenius failure = Pollard rho. F₁-geometry provides clean framework but no computational advantage.
- **Thin groups / super-approximation**: Apollonian curvatures have essentially random smoothness. Super-approximation is adversarial to factoring.
- **Random polynomial systems over Z/NZ**: Signal is O(1/N) per sample — undetectable with poly evaluations.
- **Multi-variable Diophantine**: All approaches either exponential or cost dominated by the factoring step.
- **Noetherian / ideal structure of smooth number monoid**: QS, Shor, p-1, ECM are all strategies for computing aspects of the same primary decomposition.
- **Modular representation theory of (Z/NZ)***: Detects l|(p-1) asymmetry — IS Pollard p-1 in representation-theoretic language.
- **Brumer-Stark units / explicit CFT**: Dasgupta-Kakde dominated by class group computation: L[1, 1/5], worse than L[1/3].
- **Valuation theory of Q(√N) / Pell units**: Computing fundamental unit requires O(√N). h(4N) IS as hard as factoring.
- **Multi-dimensional diophantine approximation**: LLL gives norms WORSE than NFS (the "folding barrier"). Three degrees of freedom collapse.

### Geometric/cohomological (all reduce to classical)

- **Arakelov geometry**: Product formula = same norms as NFS. Optimal sieving region = Minkowski ellipsoid (already used).
- **Condensed mathematics**: Z/NZ is discrete — all frameworks trivialize.
- **Prismatic cohomology**: Requires choosing p first = circular.
- **Derived algebraic geometry**: Derived invariants = classical invariants.
- **Shimura varieties**: Same CRT-decomposed info. Cost grows with dimension.
- **Motivic integration**: Motivic volumes decompose via CRT.
- **Weight 3/2 modular forms / Shimura correspondence**: Extraction has interpretation/computation/decomposition barriers.
- **Algebraic K-theory / cyclotomic trace / THH/TC**: K₁(Z/NZ) = (Z/NZ)*, K₂ decomposes via CRT. Dennis trace, cyclotomic trace, TR tower, Witt vectors, Nil K-theory, lambda/Adams operations, chromatic localizations ALL decompose. NK_i terms trivial for finite rings. Only non-decomposable info from RELATIVE K_*(Z, NZ) which reduces to φ(N). K-theory unifies known algorithms (p-1 = K₁, ECM = K₀ of curves, QS/NFS = K₁ via smooth relations) but provides no new shortcut.
- **Sheaf cohomology / étale fundamental groups**: H^i decomposes via CRT.
- **Hochschild / cyclic homology**: HH₀=Z/NZ, HH₁=Z/NZ (no factoring info). Higher HH decomposes.
- **Nerve complex of factor base**: Non-trivial topology BUT indistinguishable SP vs prime.
- **Persistent homology** (1D and 2D): CRT decomposition invisible to homology.
- **Modular forms at composite level N=pq**: Old/new decomposition encodes factorization but dimension formula uses divisor sums of N. Computing the space requires O(N) minimum.
- **Algebraic geometry of xy=N**: Genus 0, no abelian structure. ECM works for number-theoretic reasons, not algebraic geometry per se.
- **Deformation theory of Z/NZ**: HH^2 SPLITS via CRT. All deformation invariants decompose. No computational shortcut.
- **Arakelov intersection theory for NFS poly selection**: Gives same bounds as classical lattice geometry without Arakelov machinery.
- **Maass forms**: Less structured than holomorphic forms. Same barrier: spectrum requires knowing geometry.
- **Weil conjectures for N-dependent varieties**: Point counts give N mod l (trivial). Weil bounds are N-independent.
- **Vector bundles / Chern classes over Spec(Z/NZ)**: All bundles are free (Quillen-Suslin). Chern classes trivial. Idempotent-finding IS factoring.
- **Motivic homotopy theory / A¹-invariants**: Computing for Spec(Z/NZ) reduces to K-groups which decompose via CRT.
- **Persistent cohomology of (Z/NZ)***: Betti numbers indistinguishable SP vs prime at feasible sample sizes.
- **δ-rings and prismatic cohomology**: Z is the INITIAL δ-ring. δ-map causes EXPONENTIAL size increase — opposite of FFS. Prismatic framework formalizes Frobenius lifts but does not create new endomorphisms.
- **Perfectoid tilting**: Tilting is LOCAL; factoring is GLOBAL. The tilt requires computing the FULL inverse limit of Frobenius (classically O(N)). Shor's algorithm is precisely "efficient computation of the perfectoid tilt."
- **Approximate algebraic geometry over Z/NZ**: Epsilon-variety density is UNIVERSAL regardless of factorization type. Aggregate statistics cannot distinguish factorization structures.
- **Fargues-Scholze geometrization**: Geometric Langlands does NOT provide algorithmic shortcuts. All objects decompose via CRT or require factoring. A structural theorem, not a computational tool.
- **Galois cohomology of Z/NZ**: Fundamental trichotomy: every invariant is either (1) computable but CRT-decomposed, (2) encodes factoring but costs ≥ factoring, or (3) circular. No sweet spot.
- **Étale K-theory of Z/NZ**: K^ét₀ = Z^{ω(N)} (computing IS factoring). All decompose via CRT.

### Analytic/approximation (all vacuous or blocked)

- **abc conjecture / Baker theory / effective abc**: N=pq squarefree makes abc trivially satisfied (quality ~1/3). Baker bound astronomically weaker than useful. S-unit equation approach IS QS. Direction mismatch: abc gives LOWER bounds on radical while factoring needs UPPER bounds on smoothness.
- **Stochastic resonance**: Smoothness has no threshold dynamics. Noise reduces expected smoothness (Jensen).
- **Simulated annealing**: Cofactor landscape has no spatial correlation. SA fails above ~15 bits.
- **Gauss sum spectra**: |G(χ)|² reveals factors but requires O(N) character enumeration.
- **Wavelet/multiresolution**: Detects structure but signal weak, equivalent to trial division.
- **Dedekind zeta of Q(√N)**: Requires O(√N) terms.
- **Analytic class number at rational s**: Precision needed = O(√N) terms.
- **Stochastic resonance for LLL**: Only works for polynomial selection (already standard).
- **Analytic continuation of divisor Dirichlet series**: F(s) = (1+p^{-s})(1+q^{-s}) — evaluating requires enumerating divisors. All analytic approaches achieve O(√N) at best.
- **Bombieri-Vinogradov / primes in APs**: Distributional results validate NFS but don't improve it.
- **Zeta zeros near 2π/ln(p)**: Computing ζ(1/2+it) requires O(t^{1/2}) = O(N^{1/4}) operations = Pollard rho complexity.
- **Circle method**: Singular series encodes factors but computing it requires the factors. Additive-to-multiplicative bridge runs wrong direction.
- **Iwasawa theory / cyclotomic towers / l-adic**: All l-adic computations yield symmetric aggregates (AND/MIN/SUM). Recovering individual values IS factoring. Legendre symbol sequence indistinguishable SP vs prime.
- **RMT / arithmetic statistics / L-function zeros**: Zeros are union of two independent sets with detectably different correlation. BUT evaluation requires Ω(√N) coefficients. RMT describes the OUTPUT of an expensive computation, not a shortcut to it.
- **Selberg/Arthur trace formula**: Both sides encode factoring info. Both require factors to evaluate. Zero computational shortcut.

### Information-theoretic/complexity (barriers confirmed)

- **1-bit/relation bottleneck**: Fundamental to GF(2) framework. Proved: 1-bit-per-relation is a theorem within congruence-of-squares via universality of order-2 elements.
- **Oracle classification**: No known sub-exp function gives super-constant bits.
- **Interactive proofs / sum-check / algebraic proof systems**: IP amplifies verification not search. Divisibility function has no low-degree polynomial representation. Sumcheck prover work O(N) or O(N²).
- **ML prediction**: Accuracy → chance by 32 bits. Carry chain destroys locality.
- **Coppersmith + partial info**: n/4 threshold. Sources self-defeating.
- **Precomputation S·T bound**: No amortization. T≥L[1/3, 1.923-o(1)] for any subexp S.
- **Analytic function impossibility**: Evaluability + convergence + sensitivity impossible simultaneously.
- **Classical QFT analog**: Ω(r^{1/3}) classical lower bound.
- **Period-finding dequantization**: O(√r) tight. Tang-style inapplicable. Exponential gap robust.
- **Non-standard computation**: All classical models hit 2^n info barrier. Only quantum escapes.
- **Hybrid quantum**: O(log log N) qubits classically simulable. Shor's speedup is qualitative. Gap provably tight in generic group model.
- **Kolmogorov complexity**: K(p|N)=O(1) but K^poly(p|N) large iff factoring hard (restates problem).
- **Pell/class groups**: Genus theory gives O(1) bits only.
- **Communication complexity**: Upper O(n), lower Ω(log n). Coppersmith n/4 certificate tight.
- **Entropy/MI analysis**: MI/cost optimum matches NFS parameter selection.
- **SDP/LP relaxations**: Lasserre level O(1) insufficient. Level n/2 exact but exponential.
- **Descriptive complexity**: Factoring in ESO∩USO. FO+LFP definition = P algorithm = open.
- **Proof complexity / feasible interpolation**: Under crypto assumptions, proofs must be superpolynomial. No algorithmic leverage.
- **Bounded arithmetic**: NFS provable in S^1_2 iff factoring in P/poly. Non-constructive.
- **Combining weak signals**: All signals alpha-decay. Only algebraic (constant-SNR) signals work = QS/NFS.
- **First principles beyond GF(2)**: 4 axioms for alternative algebraic structures. NFS already exploits optimally. Beating it requires super-logarithmic info density AND polynomial combining.
- **Randomized LA for F₂ kernel**: Sketching, Kaczmarz, sparse recovery all fail over F₂. Block Lanczos appears essentially optimal.
- **Algebraic complexity theory (VP/VNP, tau conjecture)**: Factor polynomial falls outside VP/VNP. Algebraic natural proofs barrier obstructs.
- **MPC / oblivious transfer complexity**: MPC cannot resolve factoring's complexity status.
- **Compressed exponent vector search**: Exponent vectors without algebraic context carry zero information about N.
- **Succinct data structures**: S·T ≥ L[2/3]. No amortization across N values.
- **ML on exponent matrix structure**: Kernel is a global property. No local pattern predictive. Reduces to Block Lanczos.
- **Computability-theoretic / algebraic computation tree**: Steele-Yao Ω(n log n). Algebrization barrier blocks proving factoring ∉ BPP.
- **Random self-reducibility**: N → N·r² preserves factors; iterating doesn't make factoring easier.
- **Property testing for ring structure**: Full factorization requires L[1/3] complexity.
- **Classical QEC analog**: Computational advantage lives in superposition, not error correction.
- **AIT / sophistication**: Sophistication O(log n). Measures regularity structure, not accessibility.
- **Oracle separation / query complexity**: NFS uses 6 operations beyond generic ring model explaining the gap. DT model captures input complexity, not computational.
- **Hypercontractivity / global functions**: CRT opacity means no efficiently computable operation independently affects the two coordinates. Deep principle: algebraic structure ≠ computational access.
- **Non-malleable extractors / two-source extraction**: Extracted bits carry zero factoring information. GF(2) LA already extracts the MAXIMUM.
- **Beyond GF(2) exponent exploitation**: Same π(B)+1 relations needed for all k. GCD success for k>2 strictly worse (non-universal gcd(k,p-1)).
- **Meta-analysis of barrier structure**: Barriers decompose into approach-dependent (B1-B3) and structural (B4-B5). B4 (CRT opacity) is the MASTER BARRIER. B5 (Z-rigidity) is deepest. Most plausible breakthrough direction: high-dimensional abelian varieties with rich endomorphism rings over Z/NZ.
- **BBS pseudorandomness distinguishers**: BBS security EQUIVALENT to factoring. Every distinguisher reduces to detecting period.

### Dynamical/ergodic (all equivalent to known methods)

- **Ergodic theory of multiplication maps**: Mixing/orbit properties indistinguishable SP vs prime.
- **Higher-dim Pollard rho (2D, 3D)**: No improvement.
- **PCF maps over Z/NZ**: Reduce to Pollard rho variants.
- **Symbolic dynamics of multiplication-by-N map**: All invariants decompose via CRT or require O(N) computation.
- **Nonlinear non-algebraic maps (digits, Collatz, bit reversal)**: Noise-floor mutual information with x mod p. Non-algebraic operations cannot break CRT opacity.
- **Arithmetic dynamics / arboreal Galois**: Pollard rho O(N^{1/4}) OPTIMAL (Shoup lower bound). All polynomial schemes within generic group model. Preimage trees require solving x²≡a mod N = factoring.
- **Skolem-Mahler-Lech linear recurrences mod N**: GCD hit rate O(1/p). Recurrences have NO birthday analog. Period extraction = Pollard p-1.

### Parameterized/structural (no improvement)

- **FPT**: k=log(p) via ECM at exp(√(k log k)). NFS NOT FPT in factor size.
- **Matroid theory**: Standard GF(2) linear matroid, no useful second matroid.
- **AG codes**: Need base field; Z/NZ not a field.
- **Waring / r₄(N)**: Computing divisor sums = factoring.
- **Compressed sensing**: Computable measurements circular.
- **Fixed-point / PPAD**: Factoring likely in PPAD/PPP. Newton needs s=p+q.
- **Topos theory**: Functoriality barrier; effective topos respects Turing complexity.
- **FI-modules**: Wrong algebraic side.
- **Derandomization**: Smoothness heuristic is the barrier, not randomness.
- **Model theory / ultraproducts**: Zero divisors not FO-definable without naming.
- **Diophantine geometry**: xy=N is genus 0 (Faltings inapplicable).
- **Tensor networks**: Shor has O(n) entanglement. MPS/PEPS insufficient.
- **FHE analogy**: CRT is natural encryption. Bootstrapping = factoring.
- **Spin glass**: Glassy landscape, BP diverges, SA fails above ~15 bits.
- **Incidence geometry / ST over Z/NZ**: No SP vs prime difference.
- **PIT**: Universal vs existential mismatch.
- **Bio-inspired**: Flat fitness landscape, no gradient.
- **HoTT / cubical type theory**: All homotopy trivial. Univalence: computing the CRT isomorphism IS factoring.
- **Reverse math**: FTA in RCA₀, no large cardinals needed.
- **Sheaf theory on divisibility poset**: Descriptive not prescriptive.
- **Game theory**: Query complexity n/2 optimal.
- **Analytic interpolation / Carlson**: Power sums need s₁=p+q.
- **Abstract interpretation**: All domains give known bounds.
- **Additive combinatorics / Freiman/BSG**: Smooth x-positions have Freiman dim ~7 (IS the sieve). No SP distinction.
- **Szemeredi regularity**: Partition construction requires factor structure.
- **Category theory of CRT**: Idempotent splitting = factoring.
- **Differential algebra over Z/NZ**: All invariants decompose via CRT.
- **Tropical geometry**: Tropicalization erases arithmetic content.
- **Knot theory / braid groups**: Jones polynomial circular.
- **Soliton / inverse scattering**: Requires continuous measurements = trial division.
- **Quantum annealing / QUBO**: Gap closes exponentially. D-Wave best: 23-bit.
- **SAT/SMT**: Works but exponential. Wall at 40-50 bits.
- **Amortized factoring of K semiprimes**: 0/41 relations transferable. Each N independent.
- **NCG/Connes**: Spectral action = character theory. All computable quantities N-independent or circular.
- **Constraint propagation / belief propagation**: BP below condensation threshold. CSP exponential.
- **Error-correcting code / syndrome decoding**: Factor graph violates all LDPC conditions. Fundamentally broken.
- **Integer programming for xy=N**: Cannot avoid exhaustive search over √N candidates.
- **Digit distribution of √N**: Pseudorandom. Carry propagation dissolves factoring info.
- **AD / gradient methods**: Smoothness is CONJUNCTIVE, gradient optimizes ADDITIVE — fundamental mismatch.
- **MCTS / heuristic search**: Landscape information-theoretically flat until exact answer.
- **Matroid intersection**: Sieving phase lies outside matroid theory.
- **GI transfer**: Automorphism approach doesn't avoid CRT barrier.
- **Extremal graph theory on factor base**: Well-modeled by random graphs.
- **DP / noisy GCD**: Quadratic gap insurmountable.
- **PPAD / Nash equilibria**: Constant-rank bimatrix games still PPAD-hard.
- **Lattice LA for null space**: Lattice worse at all scales.
- **Sparse congruences via ISD**: 2^{O(n/log n)} worse than Gauss O(n³).
- **PFR / additive combinatorics for smooth numbers**: PFR bounds too large. Kelley-Meka 3-AP IS the sieve. Multiplication-addition barrier.

### Witt-Frobenius near-miss

p-th power is an approximate ring endomorphism with error O(p). Witt formalism makes this precise. But for composite N, no single polynomial works — CRT decomposition required. **Closest known analog to Z-Frobenius**, but precision loss is fundamental. Computational experiments confirmed: reduces to independent power maps per component. All information about primes lives at the archimedean place, making descent impossible.

## Open directions

Five structural barriers identified so far (see docs/barrier_synthesis.txt): (1) L[1/3] structural wall, (2) archimedean/non-archimedean gap, (3) GF(2) bottleneck, (4) CRT opacity, (5) Z-rigidity. Every investigated approach ran into at least one. A breakthrough could violate one, or could come from an entirely unexpected direction that sidesteps them all. These are the open directions we see:

**Super-smoothness (barrier #3 — GF(2) bottleneck):**
Can a single smooth relation yield >1 independent constraint? The algebraic object must be richer than "exponent vector mod k" for ANY k, since changing the target field doesn't help. Would need a fundamentally different algebraic coincidence.

**Near-endomorphisms of Z (barrier #5 — Z-rigidity):**
Z has no non-trivial ring endomorphism. What about endomorphisms of structures that faithfully encode factoring but are not rings? The Deuring correspondence connects to modular polynomials, class polynomials, and j-invariant arithmetic — explored through specific channels but not exhaustively. The Witt-Frobenius near-miss suggests there may be other approximate endomorphisms worth investigating.
