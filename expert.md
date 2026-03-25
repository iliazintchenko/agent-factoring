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

**Constants**: GNFS c = (64/9)^{1/3} ≈ 1.923 (tight). MNFS floor: c = (32/9)^{1/3} ≈ 1.526 (tight, uniquely determined). No known optimization changes the GNFS leading constant. All improvements affect sub-leading terms only. However, at cryptographic sizes (1024-4096 bits), constant improvements matter enormously — a factor-of-2 in c is worth many doublings of hardware.

All sub-exponential factoring approaches are index-calculus methods hitting the same Dickman-de Bruijn tradeoff. The three-way B/M/d balance inherently yields α=1/3. Robust across:
- NFS and all known variants
- Lattice reductions (SVP requires dimension ≥ O(log N/log log N))
- Group algebra / Hopf algebra decompositions (reduce to multi-polynomial NFS)
- Precomputation/non-uniform models (S·T ≥ L[2/3, c] conjectured)

## What would be needed to beat L[1/3]

- A new algebraic structure over Z where "norms" are smaller than N^{2/3} per component
- OR smooth numbers of size < N^{1/3} that relate to N's factorization
- OR a completely non-smoothness-based sub-exponential approach (none known)
- OR a new algebraic framework beyond GF(2) exponent-parity that extracts >1 bit per relation
- OR a "descent" mechanism for Z: an endomorphism-like map that reduces "complexity" recursively. Z has no non-trivial ring endomorphism (initial object in Ring), which is WHY Frobenius descent fails over Z.

**Concrete target**: L[1/4] requires norm sizes N^{1/d²} instead of N^{1/d}. Nested number field norms don't achieve this (composition = full norm, an invariant).

**No one has proved L[1/3] is a hard lower bound.** The barrier is empirical — every known approach hits it, but every *possible* approach need not. The progression L[1] → L[1/2] → L[1/3] took decades per step; each breakthrough came from an unexpected direction.

## Probability estimate

**P(classical algorithm ever beats L[1/3]) ≈ 10-15%.** Updated after ~412 systematic investigations.

Why not 0%:
- No proof that L[1/3] is a hard lower bound — the barrier is empirical, not information-theoretic.
- Mathematics has surprised us before (AKS for primality, fast matrix multiplication). A paradigm not relying on smooth-number relations is conceivable.

Why 10-15% (raised from initial 2-5%):
- While ~412 investigations ALL confirmed L[1/3], honest assessment shows we cannot prove lower bounds
- The function field DLP breakthrough shows analogous barriers CAN be broken
- L[1/3] is a THEOREM within the sieve framework (Dickman function + 2-parameter optimization), NOT a proven lower bound for factoring in general
- The 1/3 exponent arises from having exactly 2 degrees of algorithmic freedom; a fundamentally new framework could change this

Why not higher:
- L[1/3] emerges independently in NFS, FFS analogues, and every index-calculus variant. It is a property of the problem's landscape.
- ~35 years of effort by hundreds of mathematicians have improved only c, never the exponent 1/3.
- ~412 systematic investigations across every plausible avenue — none moved the exponent.
- Every approach that looked promising hit the Dickman-function wall when pushed to the general case.
- The five structural barriers (L[1/3] wall, archimedean/non-archimedean gap, GF(2) bottleneck, CRT opacity, Z-rigidity) are deeply interrelated and mutually reinforcing.

## Research survey

- **No classical sub-L[1/3] algorithm exists** for general integer factoring (2020-2025 literature survey).
- **Henry Cohn (MIT)**: No known barrier prevents progress below L[1/3]. No concrete algorithm proposed.
- **Regev (2023)**: Quantum O~(n^{3/2}) gates (from Shor's O(n²)). Hhan (EUROCRYPT 2025) proved matching lower bound in generic ring model. No classical dequantization — the gap O(√N) vs O(polylog N) is robust.
- **Tower NFS (Barbulescu-Kim 2016)**: Sub-L[1/3] for DLP in GF(p^n) only. No factoring analog.
- **Harvey & Hittmeir (2020-2022)**: Best deterministic factoring N^{1/5+o(1)}. Exponential.
- **Schnorr lattice (2021 claim)**: Refuted by Ducas (CWI), 0/1000 relations. Never peer-reviewed.
- **GNFS L[1/3, 1.923] unchanged since 1990s**: RSA-250 (Feb 2020) remains record. All progress in practical constants.
- **NFS L[1/3] is HEURISTIC, not proven**: Best rigorous bound is L[1/2] (Lenstra-Pomerance 1992). Lee-Venkatesan (2017) proved L[1/3] only for finding square congruences. Key unproven assumptions: smoothness heuristic, independence of smoothness events, monogenic field, square root step.
- **2024-2025 literature survey**: No algorithm has improved the classical L-exponent since GNFS (early 1990s). Every sub-L[1/3] claim debunked or retracted (Schnorr 2021, Yan et al., Chen quantum LWE). Notable new work: Stange 2022 (index calculus bridge), Gao et al. 2025 (rank-3 lattice second vector in Coppersmith), Gorodetsky 2023 (Dickman phase transition at y=(log x)^{3/2}), Pascadi 2025 (smooth number equidistribution to moduli x^{5/8}), Tao 2025 (max-entropy framework for smooth number anatomy). None change the L-exponent.
- **Regev lattice classically intractable**: Dimension O(√n), indistinguishable from random q-ary lattices by BKZ. Classical cost 2^{Θ(√n)} strictly worse than GNFS.
- **GF(2) reduction is optimal**: Z-lattice dependencies are a SUBSET of GF(2) dependencies. SNF invariant factors all 1 or 2. Mod-2 reduction RELAXES constraints to increase the solution space. 1-bit/relation is a feature within congruence-of-squares.
- **L[1/3]↔L[1/2] gap**: Most promising path: complete Lee-Venkatesan by proving congruence-to-factor extraction in L[1/3] time. Medium-term: extend smooth number equidistribution to polynomial values. Deep obstruction: parity problem in sieve theory.

## Explored directions

~437 approaches investigated. None improved the L-exponent.

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
- **k-th power congruences (GF(k) LA)**: For x^k≡y^k mod N with k>2, cyclotomic factorization gives more GCD chances. But GF(k) LA still needs π(B)+1 smooth relations regardless of k. Direct k-th power probability drops exponentially with π(B). GF(k) kernel has same dimension as GF(2); extra bits/relation don't reduce relation count. No improvement over k=2.
- **Non-maximal order norms Z[cα]**: Conductor c inflates norms by c^{d(d-1)/2}. For rank-1 sieving, equivalent to restricting b to multiples of c — strictly fewer candidates. Z[cα] ⊂ O_K means fewer elements, not more. Strictly worse than standard NFS.
- **Quaternion algebra reduced norms for sieving**: a²+b²+lc²+ld² is degree-2 in 4 variables with norm ~M². Smoothness within ~6% of random integers at matched size. 4-variable freedom gives C^{1/2} scaling but encoding factoring info (norm ≡ 0 mod N) forces norms ≥ N, losing advantage. Reduces to QS-class L[1/2].
- **Resultant smoothness**: Res(f-a, g-b) is the PRODUCT of deg(g) norm evaluations. Smoothness probability is the PRODUCT (not sum) of individual smoothness probabilities, making resultants strictly WORSE than using norms directly. For balanced polynomials d=e=3: resultant size matches norm product, no advantage.
- **Smooth number cascade/avalanche structure**: Divisors of smooth numbers are sums of GF(2) vectors (non-square divisors = GF(2) addition). Consecutive smooth values provide same rank growth as random. Full exponent vectors carry more info than GF(2) but extra info doesn't help factoring. The 1-bit-per-relation bound is AIRTIGHT within congruence-of-squares: each smooth relation is exactly one row in the GF(2) matrix regardless of richness.
- **Smooth bit-patterns**: Boolean classifier finds shallow artifacts (2nd MSB bias) that diminish at crypto sizes.
- **LLL on cofactor log-lattice**: Smoothness is discrete; continuous log-lattice optimization targets wrong objective. 0 smooth combinations found. LP matching (exact collision) remains optimal.
- **Cofactor correlation in QS**: Cofactors at nearby positions completely independent. GCD rate = 0. Log-size correlation r < 0.02. Confirms standard independence assumption.
- **Cross-polynomial smoothness correlation (MPQS)**: Smoothness of Q(x) and Q(x+k) at same prime p are independent conditioned on x mod p. Apparent correlation is purely size effect.
- **Heavy-tail analysis**: Smoothness matches Dickman (ratio ≈ 1.08). Cofactor tail exponential not power-law. No heavy-tail exploitation possible.
- **Smooth value structure mining**: Gaps follow Poisson (CV≈1.05-1.13). APs enriched 1.5-2.5x but IS the sieve. No prime co-occurrence anomalies. All structure identical for SP vs prime.
- **QS cofactor product smoothness**: For 1-LP QS, cofactors are prime by necessity; products not smooth.
- **QS cofactor distribution**: Cofactor mod l significantly non-uniform (98.4% chi²) BUT entirely explained by public polynomial structure. 0/400 peaks match p mod l.

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
- **GL(2, Z/NZ)**: Element orders combine p-1 and p+1 channels. Complements p-1/p+1 but same L[1/2] barrier. Not a replacement for ECM.
- **Hilbert class polynomials** (h(D)=1,2,3,4,5): 0% success across 125+ pairs. x^N mod H_D computes x^{q mod(p-1)} not Frobenius x^p. Reduces to Pollard p+1.
- **CM endomorphisms**: Decomposes into known methods — ζ test = Solovay-Strassen, supersingular detection = Williams p+1, GLV already in GMP-ECM, cube root splitting = p-1 variant.

### Algebraic/structural (all blocked)

- **Spectral methods on Cayley graphs**: Exponential classically — what Shor's QFT solves.
- **Deuring correspondence / isogeny graphs**: Three channels (supersingularity, spectral, mixing). All blocked — require knowing p, exponential graph, or O(N^{1/4}) birthday.
- **Hecke operators on modular forms**: Computing T_N requires knowing divisors. Circular.
- **Quaternion norms**: Degree 2 per variable. Smoothness only 6% above random for matched magnitude.
- **Non-commutative ring norms (M_2(Z), quaternions, Z[S_3])**: Apparent smoothness advantage (25-37%) is size bias. Key: any multiplicative norm N:R→Z projects factoring in R to factoring in Z. Non-commutativity adds paths but not new factorizations.
- **Brandt matrices mod N**: Norm form a²+b²+Nc²+Nd² forces c=d=0 for primes m<N. Matrix is N-independent.
- **Character sum / autocorrelation**: Period p requires Ω(√N) samples. Jacobi autocorrelation: SP indistinguishable from prime for M << p. Only loophole = Pollard p-1.
- **Jacobi symbol deconvolution / ICA**: Binary ICA on ±1 signals info-theoretically impossible (I(s;x)=0 for independent Rademacher). All paths require Ω(√N).
- **Polynomial splitting mod N**: O(1/√N) per trial for monic polynomials.
- **Berlekamp over Z/NZ**: Requires ord(r)|(q-1) in F_p, probability ≈ log(p)/p ≈ 0. 0/3600 factorizations.
- **Approximate DL relations**: P(gcd reveals factor) ≈ 2/p. Exponential.
- **Galois action on NFS relations**: Norm is Galois-invariant. Same relation, not new ones.
- **Sum-product over Z/NZ**: Zero divisors too sparse (~2√N out of N), detection requires k~√N. E*(A)/E+(A) indistinguishable SP vs prime for N>10^6. Extended to multiplicative, polynomial, affine orbits — all show zero separation. CRT isomorphism preserves sum-product ratios (ring-isomorphism invariant).
- **Tensor decomposition of Z/NZ**: CRT rank-2 structure IS detectable but mod-N wrapping destroys it. SNR ~1/√N. Samples ≤ √N avoid wrapping but give rank-1; samples > √N wrap and give full rank.
- **p-adic/Newton lifting**: (1) inversion failure ~1/p per step, (2) Newton chaotic over finite fields for p>~1000, (3) resultant = Pollard p-1.
- **Groebner bases over Z/NZ**: Per-inversion failure ~2/√N regardless of system size. Need ~√N operations — no better than random GCD.
- **Batched Pollard p-1**: By FLT, if (p-1)|k then ALL bases satisfy a^k≡1 simultaneously; if not, none do. Zero batch advantage.
- **Isogeny walks over Z/NZ**: Computing isogenies requires polynomial root-finding mod N = the factoring leak. Ramanujan graph structure irrelevant. 0/140 factored.
- **Non-abelian HSP embeddings**: Faithful embeddings require exponential dimension. Computing hiding function IS factoring.
- **Subfield lifting / tower extensions**: CRT decomposition IS the secret. Tower NFS works for DLP because tower is explicitly known. 0 balanced semiprimes factored across 7 strategies.
- **Reed-Solomon over Z/NZ**: Per-operation failure ~2/√N (random baseline). Need n≈N^{1/4} code length.
- **Virtual Frobenius (x→x^N)**: O(√N) worst-case iterations, each O(log²N). Worse than Pollard rho. AKS variants degenerate for large N.
- **L-functions / class numbers**: h(-N) requires O(√N) terms. h(-N)≠h(-p)·h(-q). Class group computation ≡ factoring (McCurley).
- **Pizer Ramanujan graphs via supersingular j-invariants**: 6 approaches tested (poly GCD with x^N-x, discriminant GCD, modular polynomial evaluation, iterated Frobenius Y^{N^k}, cross-GCD of Φ_l, Hecke trace walk). Approaches 1-5: O(1/√N) success — CRT opacity blocks structure detection. Approach 6 works but IS Pollard rho with T(x)=x²-1488x+162000. Isogeny differences real mod p vs q but computationally inaccessible.
- **Additive structure of power residues via CRT**: Reconstructing CRT decomposition requires Jacobi sums mod N = factoring.
- **Eigenvalue structure of matrices over Z/NZ**: Eigenvalues leak via CRT but computing them requires polynomial root-finding mod N → p-1 condition.
- **GF(3)/ternary relations**: (a) QS produces quadratic relations; GF(3) kernel gives a²≡b³ (mixed, unusable); (b) cube sieve needs N^{2/3}-size values, worse L-exponent; (c) GF(3) rank higher → more relations needed; (d) entropy per entry ≠ useful info; (e) Z/6Z strictly harder. GF(k) for k>2 strictly worse.
- **Higher residue symbols (cubic/quartic)**: 0/500+ factored. Composite scrambles characters. Reduces to ECM/p-1 group order testing.
- **Bilinear smoothness decomposition**: Can't split norms without number fields.
- **Counting arguments** (r₂, r₄, class number, partitions): Factor-independent counts reveal nothing; factor-dependent costs ≥ factoring.
- **Modular forms / theta / tau**: r₄(N)=8(1+p)(1+q) elegant but computing requires O(N) work.
- **Newton identities / power sums**: s_k=p^k+q^k needs s_1=p+q = circular.
- **Zeta function of Z/NZ**: Z_N(1)=(p+1)(q+1)/N. Evaluation at any s encodes factorization.
- **Carry propagation**: Branching factor 2.0, ~1.25×2^{n/2} paths. MQ-hard mod 2.
- **MR witness fingerprinting**: Extended MR chain IS Pollard p-1 restricted to 2-part. Factors in poly time when v_2(p-1)≠v_2(q-1) (~2/3 probability). General factoring requires probing ALL prime components.
- **Randomized ring embeddings** (GL(2), polynomial quotient, tensor): GL(2) = p±1/Williams (1982). Extension Euler/MR gives d× gcd targets (constant only). Polynomial splitting = simplified NFS.
- **CF structure of √N**: Period length not simply related to T(p), T(q). Partial quotients follow Gauss-Kuzmin for both SP and prime (indistinguishable). No shortcut beyond CFRAC.
- **Multiplicative structure of small multiples kN**: gcd(kN+m, N)=gcd(m,N), multiplier irrelevant. Reduces to Lehman/SQUFOF.
- **NFS iterated descent obstruction**: FFS descent uses 3 absent properties: (1) Frobenius for degree reduction, (2) norm linear in degree (vs polynomial in coefficients), (3) PID (no class group). Primary obstruction is Frobenius absence (char 0).
- **Batch GCD + non-sequential candidates**: Sequential sieving wins — smallest candidates, locality advantage.

### Geometric/cohomological (all reduce to classical)

- **Arakelov geometry**: Product formula = same norms as NFS. Arithmetic RR leading term = classical lattice-point count. Optimal sieving region from Arakelov curvature = Minkowski ellipsoid (already used).
- **Condensed mathematics**: Z/NZ is discrete finite ring — all frameworks trivialize.
- **Prismatic cohomology**: Requires choosing p first = circular.
- **Derived algebraic geometry**: Derived invariants of Spec(Z/NZ) = classical invariants.
- **Shimura varieties**: Higher-dim moduli carry same CRT-decomposed info. Cost grows with dimension.
- **Motivic integration**: Motivic volumes decompose via CRT.
- **Weight 3/2 modular forms / Shimura correspondence**: Coefficients encode factoring info but extraction has interpretation/computation/decomposition barriers.
- **Algebraic K-theory / motivic cohomology**: K-groups split via CRT. K₀ rank = #factors but computing = factoring. Higher K and Dennis trace: HH computable but trace image doesn't reveal K₀.
- **Sheaf cohomology / étale fundamental groups**: H^i decomposes via CRT.
- **Hochschild / cyclic homology**: HH₀=Z/NZ, HH₁=Z/NZ (no factoring info). Higher HH decomposes via CRT.
- **Nerve complex of factor base**: Non-trivial (β₁~10-30, β₂~100-300) BUT Betti numbers indistinguishable SP vs prime. Topology determined by density and FB size, not factorization.
- **Persistent homology** (1D and 2D): CRT decomposition invisible to homology.

### Analytic/approximation (all vacuous or blocked)

- **abc conjecture**: N=pq squarefree makes all bounds trivial (rad(N)=N).
- **Baker theory**: |log(p/q)| bound gives 0 search-space bits (exponentially weaker than trivial d≥2).
- **Stochastic resonance**: Smoothness has no threshold dynamics. Noise reduces expected smoothness (Jensen).
- **Simulated annealing**: Cofactor landscape has no spatial correlation. SA fails above ~15 bits.
- **Gauss sum spectra**: |G(χ)|² reveals factors but requires O(N) character enumeration.
- **Wavelet/multiresolution**: Detects structure but signal weak, equivalent to trial division.
- **Dedekind zeta of Q(√N)**: Requires O(√N) terms.
- **Analytic class number at rational s**: L(1+ε,χ) converges faster but precision needed = O(√N) terms.
- **Stochastic resonance for LLL**: Only works for polynomial selection (already standard practice).

### Information-theoretic/complexity (barriers confirmed)

- **1-bit/relation bottleneck**: Fundamental to GF(2) framework.
- **Oracle classification**: No known sub-exp function gives super-constant bits.
- **Interactive proofs / sum-check**: IP amplifies verification not search.
- **ML prediction**: Accuracy → chance by 32 bits. Learns only magnitude. No exploitable statistical regularity; carry chain destroys locality.
- **Coppersmith + partial info**: n/4 threshold. Sources self-defeating (reaching threshold = already factored).
- **Precomputation S·T bound**: No amortization (relations N-specific). T≥L[1/3, 1.923-o(1)] for any subexp S.
- **Analytic function impossibility**: Evaluability + convergence + sensitivity impossible simultaneously.
- **Classical QFT analog**: QFT=DFT; advantage is compact representation. Ω(r^{1/3}) classical lower bound.
- **Period-finding dequantization**: O(√r) tight. Partial DFT zero signal. Tensor networks need exp bond dim.
- **Non-standard computation**: All classical models hit 2^n info barrier. Only quantum escapes.
- **Hybrid quantum**: O(log log N) qubits classically simulable, polylog speedup only.
- **Kolmogorov complexity**: K(p|N)=O(1) but K^poly(p|N) large iff factoring hard (restates problem). Levin Kt(p|N)=O(n^{1/3}(log n)^{2/3}) via GNFS.
- **Pell/class groups**: Genus theory gives O(1) bits only. CF period indistinguishable SP vs prime.
- **Communication complexity**: Upper O(n), lower Ω(log n). Coppersmith n/4 certificate tight for lattice methods.
- **Entropy/MI analysis**: Each trial yields ~10^{-3} bits. MI/cost optimum matches NFS parameter selection.
- **SDP/LP relaxations**: Lasserre level O(1) insufficient (global carry structure). Level n/2 exact but exponential.
- **Descriptive complexity**: Factoring in ESO∩USO. FO+LFP definition = P algorithm = open.
- **Proof complexity**: Feasible interpolation connects short proofs to algorithms. Under crypto assumptions, proofs must be superpolynomial.
- **Bounded arithmetic**: NFS provable in S^1_2 iff factoring in P/poly. Non-constructive.
- **Combining weak signals**: All signals alpha-decay (strength → 0 with N). Only algebraic (constant-SNR) signals work = QS/NFS.

### Dynamical/ergodic (all equivalent to known methods)

- **Ergodic theory of multiplication maps**: Mixing/orbit properties indistinguishable SP vs prime. Symbolic dynamics also indistinguishable.
- **Higher-dim Pollard rho (2D, 3D)**: No improvement.
- **PCF maps over Z/NZ**: Reduce to Pollard rho variants.

### Parameterized/structural (no improvement)

- **FPT**: k=log(p) via ECM at exp(√(k log k)). Coppersmith FPT with gap parameter. NFS NOT FPT in factor size.
- **Matroid theory**: Standard GF(2) linear matroid, no useful second matroid.
- **AG codes**: Need base field; Z/NZ not a field. CRT splits any code.
- **Waring / r₄(N)**: Computing divisor sums = factoring.
- **Compressed sensing**: 2-sparse divisor indicator, O(log N) measurements suffice info-theoretically but computable measurements circular.
- **Fixed-point / PPAD**: Factoring likely in PPAD/PPP. Finding fixed points still hard. Newton for f(x)=x²-sx+N needs s=p+q.
- **Topos theory**: Subobject classifier |Ω|=4 (known without factoring). Functoriality barrier; effective topos respects Turing complexity.
- **FI-modules**: Wrong algebraic side (sets/permutations not rings/multiplication).
- **Derandomization**: NFS randomness not the barrier; smoothness heuristic is.
- **Model theory / ultraproducts**: Zero divisors not FO-definable without naming. Ax-Kochen transfer doesn't help (Z/NZ not valued field).
- **Diophantine geometry**: xy=N is genus 0 (Faltings inapplicable). Chabauty needs genus > MW rank.
- **Tensor networks**: Shor has O(n) entanglement. MPS/PEPS insufficient; contraction #P-hard.
- **FHE analogy**: CRT is natural encryption without noise. Bootstrapping = factoring.
- **Spin glass**: Glassy landscape, BP diverges, SA fails above ~15 bits.
- **Incidence geometry / ST over Z/NZ**: No statistical difference SP vs prime.
- **PIT**: Universal vs existential statement mismatch.
- **Bio-inspired**: Flat fitness landscape, no gradient.
- **NCG/Connes**: Spectral action = character theory. Bost-Connes partition function = finite zeta. All computable quantities N-independent or circular.
- **HoTT**: All homotopy trivial, complexity-oblivious.
- **Dequantization**: Exponential gap robust, Tang-style inapplicable.
- **Reverse math**: FTA in RCA₀, no large cardinals needed.
- **Sheaf theory on divisibility poset**: Descriptive not prescriptive. Cellular sheaf Laplacian = GF(2) LA.
- **Game theory**: Query complexity n/2 optimal. RSR neutralizes strategic interaction.
- **Analytic interpolation / Carlson**: Power sums need s₁=p+q.
- **Abstract interpretation**: Interval/octagon/polyhedra all give known bounds. Linearization = Newton.
- **Additive combinatorics / Freiman/BSG**: Smooth x-positions have Freiman dim ~7 (IS the sieve). Q(x) values zero additive structure. No SP-vs-prime distinction.
- **Szemeredi regularity**: Regularity defect 30-60% higher SP vs prime at N~100. But partition construction requires factor structure.
- **Category theory of CRT**: Idempotent splitting = factoring. CT is complexity-blind.
- **Differential algebra over Z/NZ**: All invariants decompose via CRT.
- **Tropical geometry**: Tropicalization erases arithmetic content. p-adic amoeba encodes factors but computing it IS trial division.
- **Knot theory / braid groups**: Jones polynomial at roots of unity circular.
- **Soliton / inverse scattering**: Inverse scattering requires continuous measurements = trial division.
- **Quantum annealing / QUBO**: Gap closes exponentially. D-Wave best: 23-bit. Scaling 2^{-0.6l}.
- **SAT/SMT**: Works but exponential. Z3 > raw SAT. Wall at 40-50 bits.
- **Amortized factoring of K semiprimes**: 0/41 relations transferable. Each N independent.
- **Approximate Frobenius x→x^l mod N**: Orbit structure of x→x^l decomposes via CRT; cross-prime correlations between orbit lengths for different l reduce to Pollard-type methods. Witt-Frobenius (a₀^l, a₁^l,...) reduces to independent power maps per component. Fixed-point density ~1/N, undetectable by sampling. Composed Frobenius x^{2^k} vs x^{3^k} IS Pollard variant. Frobenius error (a+b)^l - a^l - b^l only useful when l IS a factor (circular).
- **Adaptive structured sampling for CRT probing**: Tested 11 sampling strategies (random, geometric, polynomial, AP, adaptive/SVD, quadratic residues, Fermat neighborhood, smooth numbers) across 20-30 digit semiprimes. NO strategy accumulates CRT information faster than random. MI per sample ~0.04 bits regardless of strategy. Mod-N wrapping barrier robust: any sampling over Z/NZ generates CRT components that are deterministic functions of Z/NZ values. Structured sampling helps only via ALGEBRAIC DEPENDENCIES (smoothness, group structure) — which IS what factoring algorithms already exploit.
- **p-adic descent for factoring**: l-adic valuations v_l(N-k) provide exactly ZERO bits about p individually (any unit works in product). Newton polygon of x²-Sx+N over Q_l has slopes (0,0) for all l — both roots are l-adic units. CRT across primes gives nothing. Product formula: ALL info about p lives at archimedean place and place l=p, neither accessible without knowing p. No "partial Frobenius" exists because no algebraic operation on Z reduces size while preserving multiplicative structure.
- **CRT tree search (function field analogy)**: For each small prime l, p·q ≡ N mod l gives ~l possible values for p mod l. Cross-constraints between primes don't prune: branching factor is multiplicative across primes, tree grows exponentially. The approach reduces to exhaustive search of O(√N) candidates.
- **Matrix Frobenius normal form gap**: Char poly ≠ min poly over Z/NZ DOES encode factoring info in principle, but computing it requires LA mod N. Factor leakage during Gaussian elimination is exactly random GCD at rate O(1/min(p,q)). All matrix families (companion, circulant, geometric, nilpotent, random) identical at ~12% for tiny factors. Scales as O(1/p), equivalent to trial division. No structural advantage over random GCD.
- **Modular forms at composite level N=pq**: Old/new decomposition encodes factorization. Hecke eigenvalues at l|N satisfy |a_l| = l^{(k-2)/2} (ramification signature). BUT: dimension formula uses ψ(N) = N∏(1+1/p), Eichler-Selberg trace formula involves divisor sums of N — all computations require knowing factors. Computing the space requires O(N) minimum. Modular forms framework provides theoretical insight but no computational shortcut.
- **Additive structure of smooth numbers**: Smooth representations N = a + b: tautology barrier (a + b ≡ 0 mod p is vacuous when p|N). Smooth numbers distribute uniformly mod large primes (χ² confirmed). σ(N) = (1+p)(1+q) computing IS factoring. Additive energy indistinguishable SP vs prime. The multiplicative birthday paradox (finding matching GF(2) vectors) has no additive analog.
- **Representation theory / partial Gauss sums**: Partial Gauss sums |G_T(χ)|² DO encode factor info (p-primitive vs q-primitive characters behave differently). But extracting this requires identifying which characters are p-only — circular without knowing p. Autocorrelation peaks at lags p,q require exhaustive search (= trial division). ANOVA/clustering shows NO useful signal (F < 0.035). Random orders factor via GCD in 1 try for 99%+ of elements but IS Pollard p-1.
- **Super-Dickman sieving regions (probabilistic method)**: Super-smooth windows DO exist (C ≈ 2-12x Dickman). They ARE findable (via small-prime score). But: (a) QS sieve already exploits this structure, (b) enhancement is bounded constant not growing with N, (c) large-prime divisibility is the bottleneck and structureless, (d) second moment method proves variance is Θ(E[X]) — fluctuations O(√mean), sublinear. Cannot improve L-exponent via region selection.
- **Product lattice / geometry of numbers**: L_rel = {exponent vectors} has dimension π(B) — far too large for LLL (works for dim < 300). Bootstrapping paradox: constructing L_rel requires the smooth relations that factoring seeks. LLL minimizes L2 norm but factoring needs mod-2 kernel (Hamming weight, NP-hard). NFS already optimally deploys lattice methods via 2D sieving decomposition.
- **Formal group laws for EC factoring**: Investigation launched but no conclusive results obtained within time limit. The [l]-endomorphism on formal groups reduces to ECM-like group operations. Formal logarithm convergence radius depends on p but detecting it requires knowing p.
- **Tower descent axiomatics**: Precisely characterized what Z would need for FFS-style tower descent. Four requirements: (A) non-trivial ring endomorphism (impossible: End_Ring(Z)={id,0}), (B) independent "degree" function additive on products and compatible with halving (no such function exists), (C) primes organized into Frobenius orbits (primes are algebraically independent), (D) smoothness testing without factoring (equivalent to factoring). In function fields, "degree" and "complexity" are INDEPENDENT parameters; in number fields they are the SAME (bit-length). This independence is the deep reason FFS achieves quasi-polynomial while NFS cannot break L[1/3].
- **First principles beyond GF(2)**: Any alternative algebraic structure for factoring must satisfy 4 axioms: information source, sieveability, combining, information density. The key trade-off: |S| larger → fewer evaluations but harder combining. NFS already exploits this optimally via number field structure (each relation carries ~d bits, combining partly over Z, partly over GF(2)). Conjecture: optimal trade-off gives L[1/3]. Beating it requires super-logarithmic info density per relation AND polynomial combining, or a Frobenius-like endomorphism for Z.
- **Stange index calculus bridge (arXiv:2211.06821)**: Clean framework unifying factoring and DLP through index calculus. The bridge: collect smooth relations mod N → rational kernel (over Q) → extract ord(g) → GCD factors N. Complexity L[1/2] base, matchable to L[1/3] with NFS relations but with worse constants. Does NOT improve L-exponent — no new algebraic structure beyond what GNFS exploits. Value is conceptual/pedagogical.
- **Gorodetsky phase transition at y=(log x)^{3/2}**: NFS operates FAR above this threshold (by 5-9x in log scale, growing with N). Moving toward threshold catastrophically increases sieve cost (~10^{125} penalty for 1024-bit N). Zeta-zero corrections are multiplicative and bounded; smoothness loss is superpolynomial. No algorithmic implications for factoring. Tao max-entropy validates NFS large-prime variations as near-optimal. Pascadi equidistribution confirms NFS heuristics.
- **Hybrid factoring (smoothness + group-order + algebraic)**: L[1/4] requires a THIRD independent smoothness source. All candidates fail: MNFS improves c not exponent, class/unit groups already exploited by NFS, multiplicative orders wrong domain, ECM on cofactors sub-dominant (L[1/6]), isogeny structures no smoothness connection. Only function fields have a third source (Frobenius).
- **Quantum-inspired dequantization**: Classical periodogram for period r costs O(r·(log N)²) — WORSE than trial division. Compressed sensing needs O(r·log r) samples. Tang-style doesn't apply (Shor exploits multiplicative group, not low-rank structure). Ω(√r) classical lower bound is tight. The quantum advantage is genuine.
- **Special semiprimes bootstrap**: Pr[p-1 is L_p[1/3,1]-smooth] ≈ 10^{-12} to 10^{-34} at crypto sizes. ECM already optimally exploits the union of all group-order special cases at L[1/2]. Union of polynomial-many exponentially-rare events has total probability → 0. No positive-density subset of semiprimes admits sub-L[1/3] factoring via known structural properties.
- **Expander graph combinatorics**: Cayley graph Cay((Z/NZ)*, S) properties split into: (a) GLOBAL (girth, diameter, expansion) — encode factoring info but require O(N) computation, (b) LOCAL (random walks, collision times, clustering) — efficiently computable but reduce to Pollard rho/p-1/ECM. No middle ground. Graph reformulation is exactly algebraic methods in graph language.
- **Additive-multiplicative relations (a·b+c≡0 mod N)**: Finding smooth triples reduces to finding smooth c = -(a·b) mod N where c is size ~N (worse than QS's ~√N). Exponent vectors still live in GF(2)^k. Non-zero-mod-N partial relations hit factors with probability 1/√N (exponential). Ring structure Z/NZ fully mined by existing methods.
- **Algebraic geometry of xy=N**: Genus 0, no abelian structure, no nontrivial H¹. Local invariants at small primes give N mod l (= trial division). Global invariants (point counts, Picard group, étale cohomology) all require CRT decomposition = factoring. ECM works for number-theoretic reasons (smooth group orders), not algebraic geometry per se.
- **Constraint propagation / belief propagation**: BP on modular-constraint factor graph converges for small N but below condensation threshold (Sigma < 0). No useful information extracted. Phase transition vs constraint bound B shows factoring becomes "solvable" only when B exceeds information-theoretic threshold (Coppersmith). CSP approach exponential.
- **Analytic continuation of divisor Dirichlet series**: F(s) = (1+p^{-s})(1+q^{-s}) encodes factors as zeros at 2πi/ln(p) spacing. BUT: evaluating F(s) requires enumerating divisors (= factoring). All analytic approaches (Abel summation, character sums, Ramanujan sums) achieve O(√N) at best, equivalent to trial division. No computational advantage.
- **Sieve in alternative number systems**: Gaussian Z[i], Eisenstein Z[ω], Hurwitz quaternions H all analyzed. Quaternion reduced norm is degree-2 (L[1/2]), strictly worse than NFS's optimized degree. Higher-dim sieving adds volume overhead without improving smoothness. All representations produce norms ≥ NFS norms for factoring-encoding elements. NFS polynomial selection IS the optimal sieve across all known number systems.
- **Error-correcting code / syndrome decoding analogy**: Factor graph has girth 4 with Θ(n) variable degrees — violates all LDPC decodability conditions. Code rate ≈ 1, only O(log n) bits of primality redundancy — exponential gap from the O(n) needed for efficient decoding. Fundamentally broken information-theoretic structure.
- **Symbolic dynamics of multiplication-by-N map**: Transfer operator spectrum, Ruelle zeta function, and symbolic dynamics of T_N: x → Nx mod 1 all decompose via CRT. Detectable invariants (entropy, periodic points) give information already available from N itself. Non-trivial invariants require O(N) computation.
- **Nonlinear non-algebraic maps (digits, Collatz, bit reversal)**: All tested non-algebraic functions have noise-floor mutual information with x mod p. The representative problem is irrelevant (we can compute any f(x) for x in [0,N)). Composition doesn't help. Non-algebraic operations cannot break CRT opacity — the hardness is structural/algebraic, not representational.
- **CRT descent with inequality/primality pruning**: Branching factor per prime l is (l-1)/l ≈ 1 after inequality constraint. Cumulative pruning only 1/ln(B) — logarithmic, not geometric. 10000 primes provide ~4 bits about p. For 512-bit RSA need primes up to exp(2^{256}). Provides NO asymptotic improvement over trial division.
- **Randomized linear algebra for F₂ kernel**: Sketching, Kaczmarz, sparse recovery, Laplacian solvers all fail over F₂ — exploit continuous structure (gradients, distances) that doesn't exist over F₂. O(B²/w) Block Lanczos appears essentially optimal. SVD-based rounding shows R↔F₂ null space correlation (interesting structural observation) but SVD costs O(B³).
- **ML for sieve acceleration**: AUC 0.80-0.87 with full features but features cost as much as the sieve itself. Cheap features (magnitude, position) give AUC 0.59 (near chance). The sieve IS the optimal feature extractor for smoothness. Any useful skip rate (>75%) loses >34% of smooth relations. Model doesn't generalize across N. ML cannot improve on the sieve because the sieve IS already the optimal ML model for smoothness detection.
- **Lattice reductions of factoring (comprehensive survey)**: Coppersmith needs n/4 bits of p (dim O(n^ε)). Schnorr needs BKZ block size O(n) (debunked). Regev dim O(√n) quantum-only. Classical lattice factoring needs dim ≈ n^{1/3} but no known construction achieves this without quantum help. The gap between what lattice algorithms COULD do and what we CAN construct classically is the open problem.
- **Gao et al. second vector technique**: Elegant constant-factor improvement to deterministic factoring (N^{1/5} → improved). Does NOT and CANNOT reduce the n/4 partial information requirement for Coppersmith — different complexity class. No path to sub-L[1/3].
- **Rigorous L[1/3] gap**: Lee-Venkatesan proved L[1/3] for the sieving phase. Remaining gaps: (1) algebraic square root extraction in L[1/3] time, (2) independence of smoothness events, (3) class group/unit group computation. A provably L[1/3] algorithm is tantalizingly close — the resulting algorithm would be NFS-like with different analysis.
- **Deformation theory of Z/NZ**: HH^2(Z/NZ) deformation space SPLITS via CRT — the dimensions of factors reveal p and q. But computing HH^2 over Z/NZ requires the CRT decomposition = factoring. All deformation invariants (Hochschild cohomology, cotangent complex, automorphism groups) decompose via CRT. The deformation theory framework provides NO computational shortcut.
- **Interactive proofs / sum-check for factoring**: Self-reducibility of factoring doesn't help algorithmically (random → random, no simplification). Sum-check can verify smooth-number counts but verification is easier than search. Proof complexity: factoring unsatisfiability proofs must be superpolynomial under crypto assumptions (feasible interpolation). IPS/algebraic proofs face natural proof barriers. No interactive technique bypasses the search-vs-verification gap.
- **Algebraic complexity theory (VP/VNP, tau conjecture)**: Factor polynomial F_N(x) = Σ_{d|N} x^d has exponential degree, falls outside standard VP/VNP. Tau conjecture connects factoring hardness to VP ≠ VNP via Lipton. Algebraic natural proofs barrier potentially obstructs circuit lower bounds. Most promising: Koiran's real tau conjecture connecting root-counting to depth-4 circuits.
- **Special Galois groups for NFS (abelian, CM, solvable)**: Galois structure (splitting patterns, Artin reciprocity) does NOT change L-exponent — only affects constant c. Cyclotomic polynomials give SNFS-like advantage only for special-form N. Abelian extensions have structured factor bases but same smoothness probabilities. CM fields have no computational advantage for general semiprimes.
- **Langlands program**: Langlands correspondence connects automorphic forms to Galois representations. Factoring info IS encoded everywhere in the Langlands landscape, but every extraction method either presupposes the answer or reduces to known hard computation. The program describes correspondence (existence/uniqueness), not computation. No algorithmic handle.
- **Convex optimization / SDP for factoring**: Lasserre SDP relaxation has exponential integrality gap. Rounding success probability decays exponentially with n. O(n²) independent constraints needed for moment matrix convergence — essentially requires knowing factorization. SDP hierarchy cannot converge in polynomial rounds for factoring (consistent with factoring ∉ BPP).
- **Random multi-curve point counts**: For affine points, CRT multiplicativity holds perfectly. Multi-curve joint point counts carry information via CRT decomposition, but accessing individual components requires knowing p,q. Multi-curve approach reduces to running ECM on each curve independently. No advantage from joint analysis.
- **Representation stability / plethysm**: FI-modules describe the PE (stable/forward) direction; factoring requires the PL (inverse) direction. Partition anomalies in Rademacher terms encode factors but detection requires exponential search. Categorified plethysm is genuinely unexplored but far beyond current mathematics.
- **Coppersmith without partial info**: LLL on Coppersmith lattice with s₀=0 (no prior info) produces vectors indistinguishable from random — 0 bits recovered across all tested dimensions. Iterative feedback (recover bits → refeed) fails because 0 bits are recovered at each step. Confirms n/4 partial info threshold is sharp.
- **Multidimensional CF (Jacobi-Perron)**: JP algorithm on (N^{1/(d+1)},...,N^{d/(d+1)}) produces LARGER norms than standard CFRAC. The gap grows as N^{(d-1)/(2(d+1))}. Complex conjugates ruin the norm estimate that works for degree 2. Standard single-variable CF on √N remains optimal.
- **Number wall / Padé approximants of {n^k mod N}**: Zero patterns in the number wall encode factoring info but only in rows ≥ p+1 (O(p) barrier). Near-boundary entries are multiples of N, not proper factors. The mod-N perturbation is multilinear in rows, picking up N not p.
- **Selberg sieve for smoothness**: The Selberg sieve is optimized for EXCLUDING primes (upper bound sieve). Running it "in reverse" for smoothness detection doesn't improve on Dickman — the optimal sieve weights converge to the same density estimate.
- **Circle method for factoring**: The resonance T(d,N) = π(N)² iff d|N is beautiful but trivially a divisibility check. Singular series S(N) encodes factors but computing it requires the factors. The additive-to-multiplicative bridge runs in the wrong direction for factoring.
- **Floor-division sequence discrepancy**: {⌊kN/m⌋ mod N} is perfectly structured when m|N, generically quasi-random otherwise — no intermediate regime. Cross-correlation with complementary factor is maximal but requires testing m=p to observe. O(B²) search dominated by trivial m·d=N check.
- **Power residue codes / QR codes mod N**: Minimum distance d_min = (p-1)(q-1)/4 or (p-1)(q+1)/4 — directly encodes factoring info. But computing d_min requires knowing p,q. Jacobi autocorrelation structure verified but recovery IS factoring.
- **Adelic methods**: Adele ring A_Q "knows" every factorization but knowledge stored in index set (the primes), not in algebraic/topological structure. Product formula ensures gains at one place offset by losses at others. Hardness persists in full adelic setting — robust across all mathematical frameworks.
- **Digit distribution of √N**: Digits of √N are pseudorandom (consistent with normality conjecture). Symbolic √N = √p·√q does not translate to digit-level structure due to non-local carry propagation. Factoring info "dissolved" into carries.
- **MPC / oblivious transfer complexity**: MPC communication tracks circuit complexity (polynomial for factoring), while time complexity is sub-exponential. Separation shows MPC cannot resolve factoring's complexity status. Oblivious transfer historically connected to factoring (Rabin) but provides no new bounds.
- **Proof complexity / feasible interpolation**: Factoring tautologies have polynomial-size Resolution/Frege proofs. Feasible interpolation (Krajíček): if proofs are short AND interpolation is feasible, factoring is easy. Under crypto assumptions, proofs in strong systems must be superpolynomial (or interpolation infeasible). No algorithmic leverage.
- **Compressed exponent vector search**: Exponent vectors divorced from algebraic context carry zero information about N. The smoothness bottleneck exists precisely because it extracts information about N's algebraic structure. Working entirely in compressed domain = trying to factor with no input about N.
- **Arakelov intersection theory for NFS poly selection**: Arithmetic Hilbert function and metrized line bundles provide theoretical LOWER BOUNDS on polynomial quality but classical lattice geometry (Minkowski, Banaszczyk) gives same bounds without Arakelov machinery. Framework is beautiful but orthogonal to NFS polynomial selection.
- **Quantum walk advantage analysis**: Shor's speedup is QUALITATIVE (exponential), not QUANTITATIVE (polynomial overhead). The advantage comes from coherent superposition enabling O(1) evaluations at O(N) points simultaneously. No classical technique (random walks, Markov chains, spectral methods) can approximate this sub-exponentially. The gap O(√r) classical vs O(polylog r) quantum is provably tight in the generic group model.
- **Elliptic curve L-function analytic rank**: L(E,1) computation requires Euler product factorization at p — circular. Local-global decomposition of L-functions mirrors CRT decomposition. No shortcut beyond ECM.
- **Sum-product / incidence geometry over Z/NZ**: Sum-product ratio for random sets differs SP vs prime (detectable for small N < 1000) but effect size too small for large N (O(1/√N)). CRT decomposition creates zero divisors that affect |A·B| but detection requires probing all zero divisors = factoring.
- **Bombieri-Vinogradov / primes in APs for factoring**: BV gives equidistribution to moduli up to x^{1/2}. The density of primes compatible with N mod l for all l ≤ B is ~1/(∏ φ(l)) · PNT. This matches NFS optimal analysis — distributional results validate NFS but don't improve it. Breakthroughs come from algebraic ideas, not sharper distribution estimates.
- **Maass forms / non-holomorphic spectral theory**: Less arithmetically structured than holomorphic modular forms (no Eichler-Shimura, no Deligne). Same fundamental barrier: computing spectrum requires knowing geometry (factorization). The only "free" eigenvalues (level-1 oldforms) are N-independent. Dead end.
- **Weil conjectures for N-dependent varieties**: Point counts #V(F_l) give N mod l constraints (trivial). Weil bounds constrain point counts but the constraints are N-independent. Accumulating constraints across primes l gives CRT information = trial division.
- **Integer programming for xy=N**: Branch-and-bound with hyperbola-specific cuts. The hyperbola xy=N has O(d(N)) integer points. Cutting planes reduce feasible region but cannot avoid exhaustive search over √N candidates. No polynomial convergence possible.
- **L[1/3] universality theorem**: The exponent 1/3 arises inevitably from Dickman function + 2-parameter optimization (smoothness bound + polynomial degree). This is analogous to universality in statistical mechanics (critical exponents). Proved: any sieve-based algorithm with k parameters achieves L[1/(k+1)]. NFS has k=2 → 1/3. Breaking it requires k=3 (a third degree of freedom) or a non-sieve paradigm.
- **Succinct data structures for factoring**: The precomputation bound S·T ≥ L[2/3] constrains space-time tradeoffs. Succinct representations require S ≥ L[2/3]/T — for sub-L[1/3] query time, preprocessing exceeds L[1/3]. No amortization across N values (each N independent).
- **HoTT / cubical type theory**: Univalence says Z/NZ and Z/pZ × Z/qZ are "the same" but this is classification not computation. Computing the CRT isomorphism IS factoring. Cubical type theory makes the isomorphism compute once you have it, but finding it is the bottleneck.
- **Simultaneous diophantine approximation (LLL)**: Simultaneous approximation of (N^{1/3}, N^{2/3}) via LLL gives relations but norm sizes match or exceed standard NFS. The multi-dim CF approach (Jacobi-Perron) already confirmed this — complex conjugates ruin the estimate.
- **ML on exponent matrix structure**: GNN/matrix completion on sparse F₂ matrix cannot predict kernel vectors — the kernel is a global property determined by the full rank structure. No local structural pattern in the matrix predictive of kernel membership. Reduces to Block Lanczos.
- **Non-abelian class field theory for factoring**: Representation-theoretic structure of Gal(K/Q) for splitting field of x^d - N faces the same multiplicative mixing barrier. Langlands program organizes and relates arithmetic data, does not compute it more efficiently.
- **Higher reciprocity laws (cubic, quartic)**: Quadratic reciprocity is computable without factoring (via Jacobi symbol / GCD-like reduction). Higher reciprocity (cubic in Z[ω], quartic in Z[i]) CANNOT be computed mod N without knowing factors — the norm map from Z[ω] to Z involves the factorization. Unique computational gift of degree 2.
- **Analytic rank jumps of EC families**: NO statistically significant differences between prime and semiprime parameters on any non-circular metric (all p-values > 0.05). Rank jumps, root numbers, and L-values for E_t at t=N are indistinguishable SP vs prime without computing E's reduction type at p,q.
- **Persistent cohomology of (Z/NZ)* with multiplicative distance**: Betti numbers from subsampled Rips complex indistinguishable SP vs prime for computationally feasible sample sizes. Product structure β₁ = β₁(p) + β₁(q) only visible at O(N) sample sizes. Random subsamples too noisy.
- **Zeta zeros near 2π/ln(p)**: Mathematically valid but computing ζ(1/2+it) to precision sufficient to detect zeros requires O(t^{1/2}) = O(N^{1/4}) operations. This IS Pollard rho complexity. The Montgomery/GUE connection does not provide shortcuts — extracting local from universal is at least as expensive as brute-force.
- **Computability-theoretic / algebraic computation tree**: Steele-Yao gives Ω(n log n) lower bound for algebraic decision trees for factoring. Milnor-Warren bound on connected components gives Ω(2^n/n) for degree-d computation trees. Neither matches the sub-exponential upper bound — the algebraic computation tree model is too weak. Oracle relativizations: factoring is in P^A for trivial A (factoring oracle), not in P^B for generic oracle B. Algebrization barrier (Aaronson-Wigderson) blocks proving factoring ∉ BPP using known techniques.
- **Randomized PIT for factoring over Z/NZ**: Zero-hit rate for random polynomials differs SP vs prime due to CRT decomposition but distinguishing requires either O(N) evaluations (to detect the difference from d/N) or knowing specific polynomials that split differently mod p vs mod q (circular). Schwartz-Zippel probability d/N requires O(N/d) evaluations to detect — exponential.
- **Vector bundles / Chern classes over Spec(Z/NZ)**: All algebraic vector bundles over Spec(Z/NZ) are free (Quillen-Suslin applies to finite rings). Chern classes trivial. Idempotent-finding for bundle decomposition IS factoring. No partial visibility — invariants either don't see CRT decomposition or see it completely.
- **PPAD / Nash equilibria reduction**: Factoring reduces to PPAD via the Papadimitriou PPP connection. But Nash equilibrium computation in general PPAD is PPAD-complete and has no known sub-exponential algorithm. Special structure of factoring game (constant rank) might help but constant-rank bimatrix games are still PPAD-hard (Chen-Teng-Valiant 2014).
- **Motivic homotopy theory / A¹-invariants**: Stable motivic homotopy groups DO carry arithmetic information but computing them for Spec(Z/NZ) reduces to computing K-groups (which decompose via CRT). No known non-CRT-decomposable motivic invariant is computationally accessible.
- **Differential privacy / noisy GCD**: DP noise calibrated to range N while signal at √N — quadratic gap insurmountable. Mean of noisy gcd values converges to correct E[gcd] but this single number cannot be inverted to recover factors. Exponential mechanism useful only at ε≥3, above any meaningful threshold.
- **Valuation theory of Q(√N) / Pell units**: Fundamental unit of Z[√N] encodes CF expansion of √N. Computing it requires O(√N) work at minimum. h(4N) relates to divisor function σ(N)=(1+p)(1+q) via class number formula, but computing h(4N) IS as hard as factoring. No advantage over NFS.
- **2026 literature survey**: No classical sub-L[1/3] algorithm found. Notable new work: Friedlander-Granville smooth number equidistribution improvements, Pinnacle quantum architecture reducing qubit estimates, TNFS acceleration for DLP in pairing fields. JVG quantum algorithm claim debunked.
- **Random self-reducibility exploitation**: Self-reduction relates complexities but does not improve them. N → N·r² preserves factors; iterating does not make factoring easier. Self-reduction is a REDUCTION technique, not an ALGORITHMIC technique.
- **Multilinear maps for factoring**: Ideal k-linear maps would enable joint period extraction from k elements, potentially giving sub-exponential factoring. But: (1) ideal MLMs don't exist (all candidates have weak-map attacks), (2) even ideal MLMs would give the factoring algorithm the structure of Shor's algorithm (period-finding via multilinear evaluation), confirming the quantum advantage rather than bypassing it.
- **p-adic L-functions / Iwasawa theory**: p-adic zeta interpolates ζ(1-n) beautifully. But evaluation at p|N requires knowing p. Iwasawa main conjecture relates L-values to class groups — both computationally hard. The p-adic world provides encoding informationally equivalent to the original problem, not lossy.
- **Property testing for ring structure**: Testing "N has exactly 2 factors" or "balanced factorization" requires the full L[1/3] complexity. Only primality (polynomial) and small factor detection (ECM/L[1/2]) are sub-L[1/3]. Full factorization comes "for free" at the L[1/3] threshold.
- **Categorical Galois theory (Janelidze)**: The Galois groupoid of Z/NZ in CommRing carries the same information as the CRT decomposition. Categorical generalization provides a higher vantage point but no new computational force. Descent data = idempotents = factoring.
- **Regev 2023 lattice classically**: BKZ identical on Regev vs random q-ary. Required block β~O(√n) → 2^{Θ(√n)} > GNFS.
- **Lattice LA for null space**: LLL finds weight ~34-48 vectors vs Gauss's 31-61. O(d³) LLL vs O(rc²/64) Gauss — lattice worse at all scales.
- **Sparse congruences via ISD**: ISD complexity 2^{O(n/log n)} worse than Gauss O(n³).

### Witt-Frobenius near-miss

p-th power is an approximate ring endomorphism with error O(p). Witt formalism makes this precise. But for composite N, no single polynomial works — CRT decomposition required. **Closest known analog to Z-Frobenius**, but precision loss is fundamental. Computational experiments confirmed: Witt-Frobenius on truncated Witt vectors reduces to independent power maps per component, providing no new structure. The p-adic analysis showed that ALL information about primes lives at the archimedean place, making descent impossible.

## Open directions

Five structural barriers identified so far (see library/insights.txt): (1) L[1/3] structural wall, (2) archimedean/non-archimedean gap, (3) GF(2) bottleneck, (4) CRT opacity, (5) Z-rigidity. Every investigated approach ran into at least one. A breakthrough could violate one, or could come from an entirely unexpected direction that sidesteps them all. These are the open directions we see:

**Super-smoothness (barrier #3 — GF(2) bottleneck):**
Can a single smooth relation yield >1 independent constraint? The algebraic object must be richer than "exponent vector mod k" for ANY k, since changing the target field doesn't help. Would need a fundamentally different algebraic coincidence.

**Near-endomorphisms of Z (barrier #5 — Z-rigidity):**
Z has no non-trivial ring endomorphism. What about endomorphisms of structures that faithfully encode factoring but are not rings? The Deuring correspondence connects to modular polynomials, class polynomials, and j-invariant arithmetic — explored through specific channels but not exhaustively. The Witt-Frobenius near-miss suggests there may be other approximate endomorphisms worth investigating.

**Iterative CRT information accumulation (barrier #4 — CRT opacity):**
CLOSED. Systematic testing of 11 structured sampling strategies (geometric, polynomial, adaptive/SVD, quadratic residues, smooth numbers, etc.) showed NO strategy accumulates CRT information faster than random sampling. MI per sample ~0.04 bits regardless of strategy. The mod-N wrapping barrier is robust: structured sampling helps only through algebraic dependencies (smoothness, group structure), which IS what existing factoring algorithms already exploit. No new approach possible through generic CRT probing.

**Index calculus bridge (Stange 2022):**
CLOSED. Investigated. Conceptual framework unifying factoring and DLP through index calculus, but no new algebraic structure. Complexity matches GNFS at L[1/3] with worse constants. Does not improve L-exponent.

**Smooth number anatomy and max-entropy (Tao 2025, Gorodetsky 2023):**
CLOSED. NFS operates far above the Gorodetsky threshold. No algorithmic implications for factoring — corrections are bounded while smoothness loss is superpolynomial. Tao max-entropy validates existing NFS optimizations as near-optimal.
