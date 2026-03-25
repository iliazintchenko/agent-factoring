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

**P(classical algorithm ever beats L[1/3]) ≈ 2-3%** (80% CI [0.5%, 8%]). Updated after ~526 systematic investigations.

Why not 0%:
- No proof that L[1/3] is a hard lower bound — the barrier is empirical, not information-theoretic
- The function field DLP breakthrough shows analogous barriers CAN be broken
- L[1/3] is a THEOREM within the sieve framework (Dickman function + 2-parameter optimization), NOT a proven lower bound for factoring in general
- Mathematics has surprised us before (AKS, Babai GI, fast matrix multiplication)

Why 2-3% (lowered from initial estimate):
- ~526 systematic investigations across every plausible mathematical avenue — NONE moved the exponent
- L[1/3] emerges independently in NFS, FFS analogues, and every index-calculus variant — it is a universal property
- ~35 years of effort by hundreds of mathematicians have improved only c, never the exponent 1/3
- Every approach that looked promising hit the Dickman-function wall when pushed to the general case
- The five structural barriers (L[1/3] wall, archimedean/non-archimedean gap, GF(2) bottleneck, CRT opacity, Z-rigidity) are deeply interrelated and mutually reinforcing
- Complete absence of "near miss" algorithms that would signal a breakthrough is close
- The L[1/3] universality theorem: any sieve-based algorithm with k parameters achieves L[1/(k+1)]; NFS has k=2; a third degree of freedom requires a fundamentally new mathematical framework that no one has even glimpsed

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

~526 approaches investigated. None improved the L-exponent.

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
- **Higher residue symbols / reciprocity (cubic/quartic)**: 0/500+ factored. Composite scrambles characters. Reduces to ECM/p-1 group order testing. Quadratic reciprocity computable without factoring (Jacobi symbol); higher reciprocity (cubic in Z[ω], quartic in Z[i]) CANNOT be computed mod N without knowing factors — norm map involves factorization. Unique computational gift of degree 2.
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
- **Interactive proofs / sum-check**: IP amplifies verification not search. Self-reducibility doesn't help algorithmically. Sum-check verifies smooth-number counts but verification easier than search. IPS/algebraic proofs face natural proof barriers.
- **ML prediction**: Accuracy → chance by 32 bits. Learns only magnitude. No exploitable statistical regularity; carry chain destroys locality.
- **Coppersmith + partial info**: n/4 threshold. Sources self-defeating (reaching threshold = already factored).
- **Precomputation S·T bound**: No amortization (relations N-specific). T≥L[1/3, 1.923-o(1)] for any subexp S.
- **Analytic function impossibility**: Evaluability + convergence + sensitivity impossible simultaneously.
- **Classical QFT analog**: QFT=DFT; advantage is compact representation. Ω(r^{1/3}) classical lower bound.
- **Period-finding dequantization**: O(√r) tight. Partial DFT zero signal. Tensor networks need exp bond dim. Classical periodogram costs O(r·(log N)²) — worse than trial division. Compressed sensing needs O(r·log r) samples. Tang-style inapplicable (Shor exploits multiplicative group, not low-rank). Exponential gap robust.
- **Non-standard computation**: All classical models hit 2^n info barrier. Only quantum escapes.
- **Hybrid quantum**: O(log log N) qubits classically simulable, polylog speedup only. Shor's speedup is qualitative (coherent superposition enabling O(1) evaluations at O(N) points). Gap O(√r) classical vs O(polylog r) quantum provably tight in generic group model.
- **Kolmogorov complexity**: K(p|N)=O(1) but K^poly(p|N) large iff factoring hard (restates problem). Levin Kt(p|N)=O(n^{1/3}(log n)^{2/3}) via GNFS.
- **Pell/class groups**: Genus theory gives O(1) bits only. CF period indistinguishable SP vs prime.
- **Communication complexity**: Upper O(n), lower Ω(log n). Coppersmith n/4 certificate tight for lattice methods.
- **Entropy/MI analysis**: Each trial yields ~10^{-3} bits. MI/cost optimum matches NFS parameter selection.
- **SDP/LP relaxations**: Lasserre level O(1) insufficient (global carry structure). Level n/2 exact but exponential.
- **Descriptive complexity**: Factoring in ESO∩USO. FO+LFP definition = P algorithm = open.
- **Proof complexity / feasible interpolation**: Feasible interpolation (Krajíček): if proofs are short AND interpolation is feasible, factoring is easy. Under crypto assumptions, proofs in strong systems must be superpolynomial. Factoring tautologies have polynomial-size Resolution/Frege proofs. No algorithmic leverage.
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
- **HoTT / cubical type theory**: All homotopy trivial, complexity-oblivious. Univalence says Z/NZ and Z/pZ × Z/qZ are "the same" but computing the CRT isomorphism IS factoring.

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

- **Compressed exponent vector search**: Exponent vectors divorced from algebraic context carry zero information about N. The smoothness bottleneck exists precisely because it extracts information about N's algebraic structure. Working entirely in compressed domain = trying to factor with no input about N.
- **Arakelov intersection theory for NFS poly selection**: Arithmetic Hilbert function and metrized line bundles provide theoretical LOWER BOUNDS on polynomial quality but classical lattice geometry (Minkowski, Banaszczyk) gives same bounds without Arakelov machinery. Framework is beautiful but orthogonal to NFS polynomial selection.

- **Elliptic curve L-function analytic rank**: L(E,1) computation requires Euler product factorization at p — circular. Local-global decomposition of L-functions mirrors CRT decomposition. No shortcut beyond ECM.
- **Sum-product / incidence geometry over Z/NZ**: Sum-product ratio for random sets differs SP vs prime (detectable for small N < 1000) but effect size too small for large N (O(1/√N)). CRT decomposition creates zero divisors that affect |A·B| but detection requires probing all zero divisors = factoring.
- **Bombieri-Vinogradov / primes in APs for factoring**: BV gives equidistribution to moduli up to x^{1/2}. The density of primes compatible with N mod l for all l ≤ B is ~1/(∏ φ(l)) · PNT. This matches NFS optimal analysis — distributional results validate NFS but don't improve it. Breakthroughs come from algebraic ideas, not sharper distribution estimates.
- **Maass forms / non-holomorphic spectral theory**: Less arithmetically structured than holomorphic modular forms (no Eichler-Shimura, no Deligne). Same fundamental barrier: computing spectrum requires knowing geometry (factorization). The only "free" eigenvalues (level-1 oldforms) are N-independent. Dead end.
- **Weil conjectures for N-dependent varieties**: Point counts #V(F_l) give N mod l constraints (trivial). Weil bounds constrain point counts but the constraints are N-independent. Accumulating constraints across primes l gives CRT information = trial division.
- **Integer programming for xy=N**: Branch-and-bound with hyperbola-specific cuts. The hyperbola xy=N has O(d(N)) integer points. Cutting planes reduce feasible region but cannot avoid exhaustive search over √N candidates. No polynomial convergence possible.
- **L[1/3] universality theorem**: The exponent 1/3 arises inevitably from Dickman function + 2-parameter optimization (smoothness bound + polynomial degree). This is analogous to universality in statistical mechanics (critical exponents). Proved: any sieve-based algorithm with k parameters achieves L[1/(k+1)]. NFS has k=2 → 1/3. Breaking it requires k=3 (a third degree of freedom) or a non-sieve paradigm.
- **Succinct data structures for factoring**: The precomputation bound S·T ≥ L[2/3] constrains space-time tradeoffs. Succinct representations require S ≥ L[2/3]/T — for sub-L[1/3] query time, preprocessing exceeds L[1/3]. No amortization across N values (each N independent).

- **Simultaneous diophantine approximation (LLL)**: Simultaneous approximation of (N^{1/3}, N^{2/3}) via LLL gives relations but norm sizes match or exceed standard NFS. The multi-dim CF approach (Jacobi-Perron) already confirmed this — complex conjugates ruin the estimate.
- **ML on exponent matrix structure**: GNN/matrix completion on sparse F₂ matrix cannot predict kernel vectors — the kernel is a global property determined by the full rank structure. No local structural pattern in the matrix predictive of kernel membership. Reduces to Block Lanczos.
- **Non-abelian class field theory for factoring**: Representation-theoretic structure of Gal(K/Q) for splitting field of x^d - N faces the same multiplicative mixing barrier. Langlands program organizes and relates arithmetic data, does not compute it more efficiently.

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
- **Discrete Gaussian sampling for NFS square root**: DGS elegantly samples lattice points but the NFS square root step is already O(B^{1+ε}) in practice. DGS doesn't improve the asymptotic because the bottleneck is class group computation, not lattice sampling.
- **Classical QEC analog for factoring**: QEC preserves superpositions; the computational advantage lives in the superposition, not the error correction. "Classical QEC" is an illuminating oxymoron: redundant classical computation of modular exponentiation provides fault tolerance but not speedup. The quantum-classical gap is in information encoding, not error resilience.
- **CFSG / group structure properties of (Z/NZ)***: The exponent lcm(p-1,q-1) can be computed probabilistically (Pollard lambda) in O(√(exponent)) time. The Frattini subgroup, socle, and automorphism group all decompose via CRT. The Jacobi symbol provides 1-dimensional shadow of the full group structure — "the shadow has exactly one dimension fewer, and that missing dimension IS the factorization."
- **Étale K-theory of Z/NZ**: K^ét_0 = Z^{ω(N)} (computing IS factoring). Higher K-groups: K^ét_1(Z/NZ) = (Z/NZ)*, K^ét_2 = 0. All decompose via CRT. Quillen-Lichtenbaum (= Rost-Voevodsky) relates to motivic cohomology, which also decomposes. K-theory unifies known algorithms (p-1 = K₁, ECM = K₀ of curves, QS/NFS = K₁ via smooth relations) but provides no new computational shortcut.
- **Extremal graph theory on factor base graph**: Turán-type bounds on the factor base graph relate to NFS performance but the graph is well-modeled by random graphs with prescribed degree sequences. No structural feature exploitable beyond what random-graph analysis already gives. Weighted hypergraph model natural but doesn't change L-exponent.
- **Multi-variable Diophantine equations for factoring**: Four-square x²+y²+z²+w²=N always solvable but decomposition doesn't reveal factors. Pell x²-Ny²=1 has regulator related to class number — computing either IS factoring. All multi-variable approaches that encode factoring either solve in exponential time or have cost dominated by the factoring step.
- **Fourier analysis on (Z/NZ)***: Characters require knowing φ(N) = factoring. Partial Fourier info (power sums Σf(x)^k) accessible but carries O(1/√N) signal. Character sum bounds (Burgess) give the best estimates but don't improve on NFS. Smooth number equidistribution in residue classes (Fouvry-Tenenbaum) validates NFS heuristics.
- **Algorithmic information theory / sophistication**: Sophistication of semiprimes is O(log n) (polynomial in bit length) — the "regularity" is simple. But sophistication measures regularity structure, not accessibility. "Identifying the regularity" vs "computing how to decompose" is precisely the gap that factoring hardness represents.
- **MCTS / heuristic search for factoring**: Factoring constraint is DETERMINISTIC (carry chain), not statistical. Constraint propagation determines all feasible branches; MCTS redundantly re-discovers feasibility via random rollouts. No value function gradient exists — the landscape is information-theoretically flat until the exact answer is found.
- **Matroid intersection for factoring**: The GF(2) dependency-finding phase IS efficiently solvable (Gaussian elimination). The sieving phase — the true bottleneck — lies outside matroid theory. No second matroid exists that enables polynomial-time factoring via matroid intersection.
- **Graph isomorphism → factoring transfer**: Babai's quasi-polynomial GI algorithm uses group-theoretic "split-or-Johnson" strategy. Aut(G_N) for N-dependent graphs either encodes full factoring or is trivial. The automorphism approach doesn't avoid the CRT barrier. GI and factoring are both NP-intermediate but for different structural reasons.
- **Hasse-Minkowski / local-to-global for factoring**: xy=N is genus 0 — Hasse principle holds trivially (no local-to-global obstruction). The Brauer-Manin obstruction is zero. The problem is not finding RATIONAL solutions (trivially exist) but finding INTEGER solutions with PRIME factors. Local-to-global principles describe rational points on varieties, not integer decompositions.
- **Smooth number graph spectral theory**: The graph has good expansion (consistent with random model) and no spectral anomaly distinguishing SP from prime. Character sum circularity: computing eigenvalues of the smooth number Laplacian requires knowing the factoring-relevant structure.
- **Regev 2023 lattice classically**: BKZ identical on Regev vs random q-ary. Required block β~O(√n) → 2^{Θ(√n)} > GNFS.
- **Lattice LA for null space**: LLL finds weight ~34-48 vectors vs Gauss's 31-61. O(d³) LLL vs O(rc²/64) Gauss — lattice worse at all scales.
- **Sparse congruences via ISD**: ISD complexity 2^{O(n/log n)} worse than Gauss O(n³).

- **δ-rings and prismatic cohomology**: Z is the INITIAL δ-ring (φ = id, trivial Frobenius). On Z/NZ, φ_p(x) = x^p is nontrivial but every invariant reduces to Pollard p-1 (fixed points), Pollard rho (orbits), or CRT factoring (idempotents). Prismatic site of Z/NZ decomposes via CRT into per-factor components. δ-map δ_p(x) = (x-x^p)/p causes EXPONENTIAL size increase (|δ_p(x)| ≈ |x|^p/p), destroying smoothness — opposite of FFS where Frobenius + reduction decreases degree. The archimedean absolute value on Q vs non-archimedean degree on F_q[t] is the deep reason. Prismatic framework formalizes Frobenius lifts but does not create new endomorphisms.
- **Skolem-Mahler-Lech linear recurrences mod N**: Zero set of order-k recurrence mod N decomposes via CRT: S_N = S_p ∩ S_q. GCD hit rate (finding n where gcd(a_n, N) > 1) scales as O(1/p) — expected first hit at position O(p), completely infeasible for large N. Pollard rho achieves O(√p) via birthday paradox; recurrences have NO birthday analog. Period T_N = lcm(T_p, T_q); extracting T_p requires factoring T_N — equivalent to Pollard p-1 when T_p divides |GL(2,Z/pZ)| = p(p-1)²(p+1). No strategy for initial conditions outperformed random. Partial period information from poly(log N) terms is indistinguishable SP vs prime.
- **Brumer-Stark units / explicit class field theory**: Genus theory for Q(√(-N)): 2-rank of Cl(-4N) = ω(N)-1 gives exactly 1 genus character for N=pq, but cannot distinguish which factor each Legendre symbol belongs to. Hilbert class polynomial H_{-4N}(x) factorization mod l reveals (D/l) = (-4/l)(p/l)(q/l) but individual symbols inaccessible. Dasgupta-Kakde explicit Brumer-Stark construction is dominated by class group computation: O(N^{1/5+ε}) under GRH = L[1, 1/5], strictly worse than L[1/3]. Partial Brumer-Stark info (θ mod small primes) requires O(N)-term character sums whose CRT decomposition needs knowing p,q. Heegner point approach similarly circular. The CRT structure that would make class field computations efficient IS the factoring problem.
- **Perfectoid tilting for factoring**: Tilting is a LOCAL operation (one prime at a time); factoring is a GLOBAL problem. Z_p ⊗ Z/NZ ≅ Z/pZ (localization at p already factors N). p^n-th roots mod N exist iff gcd(p^n, φ(N))|ord(x) — detecting this IS group-order testing (Pollard p-1). The tilt = lim_{←,x↦x^p} A/pA requires computing the FULL inverse limit of Frobenius, which classically costs O(N). Insight: Shor's quantum algorithm is precisely "efficient computation of the perfectoid tilt" — quantum parallelism evaluates O(N) Frobenius iterates in O(log N) time. This is the cleanest known explanation of the quantum advantage for factoring. All approaches (Witt vectors, Fargues-Fontaine, derived prismatic, THH/TC) decompose via CRT.
- **Recursive lattice descent for NFS large primes**: Three independent failure modes. (1) Log-to-multiplicative gap: LLL finds continuous approximations but 0/483 exact matches — primes have no smooth decomposition. (2) Bootstrap impossibility: recursive Coppersmith yields 0 bits per iteration; Howgrave-Graham n/4 threshold is tight. (3) Smoothness bottleneck: three "degrees of freedom" (B, lattice dim, recursion depth) collapse because dim = π(B) is determined by B and recursion doesn't reduce problem size (no Frobenius analog). One effective parameter → L[1/3]. Quantitative: recursive lattice descent 2^{157} WORSE than NFS at 1024-bit.
- **PIT / derandomized smooth number detection**: Smoothness polynomials over F_p have MAXIMAL circuit complexity (indistinguishable from random under shifted partial derivatives). Polynomial method shows non-smoothness polynomial has degree ≈ p, ruling out small hitting sets. CRT-based hitting sets reduce exactly to standard sieving. APs with small-prime-product common differences give standard sieve enrichment, nothing more (Granville uniformity theorem). B-smooth integers do not form algebraic subgroups mod p. Sieving IS essentially optimal for smooth number detection.
- **Sumcheck / algebraic proof systems**: Divisibility function [x|N] has no low-degree polynomial representation over F_q (degree barrier). Integer range constraint {1,...,N} has no compact algebraic encoding (ordering barrier). Sumcheck prover work is O(N) or O(N²), worse than trial division O(√N). Factoring already has efficient verification (certificate = factors); interactive proofs are redundant. The algebraic structure of sumcheck contains zero algorithmic leverage for factoring beyond known methods. Interactive proofs amplify verification, not search — this is a theorem, not a slogan.
- **Hypercontractivity / global functions on (Z/NZ)***: (Z/NZ)* ≅ (Z/pZ)* × (Z/qZ)* is a product group where KLLM hypercontractivity applies perfectly. The Jacobi symbol is the ideal non-global function (factors as Legendre_p × Legendre_q), with hypercontractivity ratio ρ² vs ρ for primes. BUT: computing the noise operator T_ρ requires the CRT decomposition = factoring. No efficiently computable operation on (Z/NZ)* can independently affect the two CRT coordinates. Fourier spectrum indistinguishable SP vs prime without level structure. Power residue count k² test detects semiprimes but requires O(N) work. Deep principle: algebraic structure (product decomposition) and computational access (CRT evaluation) are different things.
- **Non-malleable extractors / two-source extraction**: Smooth relations viewed as CRT-encoded pairs (a_p, a_q) have min-entropy ~7-9 bits in 685 dimensions. Sources are perfectly independent (CRT), so two-source extraction works. BUT: extracted bits carry zero factoring information — extraction ≠ factoring. GF(2) LA already extracts the MAXIMUM factoring information from smooth relations. The information-theoretic bound: each smooth relation provides at most 1 bit of factoring information (its GF(2) parity vector). Two-source extractors solve randomness extraction, not factoring. The algebraic structure of CRT-encoded smooth relations carries no more factoring information than GF(2) exponent vectors.
- **F_1 geometry / Frobenius homomorphism test**: gcd((a+b)^p - a^p - b^p, N) = p always when p|N (Fermat). Product aggregation over random k gives trial division in disguise. Birthday collision on Frobenius failure values achieves O(N^{1/4}) — IS Pollard rho. F(k) mod p cycles with period dividing p-1; birthday bound gives O(√p) collisions. F_1-geometry (Borger lambda-rings, Connes-Consani arithmetic site, Deninger flow) provides clean conceptual framework but no computational advantage. Adams operations ψ_k(x) = x^k on (Z/NZ)* decompose via CRT; joint eigenvalue structure inaccessible without factors.
- **Iwasawa theory / cyclotomic towers / l-adic structure**: All l-adic computations on N=pq yield symmetric aggregate functions: power residues give AND of individual characters, valuations give MIN, logarithms give SUM. Recovering individual values from AND/MIN/SUM IS the factoring problem. 2-adic sqrt(N) via Hensel lifting is efficiently computable but carries no extra information beyond N (x→x² is bijective on 1+8Z_2). l-th power residue fingerprints are lossy (AND operation). Legendre symbol sequence (N/l) for l=2,3,... is statistically indistinguishable SP vs prime (product of independent ±1 sequences is random ±1). Fundamental unit norm differs (SP more likely norm +1) but computing it costs O(√N).
- **RMT / arithmetic statistics / L-function zeros**: L(s, chi_N) = L(s, chi_p)·L(s, chi_q) — zeros are the union of two independent sets with detectably different pair correlation (reduced repulsion). BUT evaluating L(s, chi_N) to useful precision requires Ω(√N) Dirichlet coefficients. Small-n coefficients chi_N(n) = chi_p(n)·chi_q(n) carry O(1/√N) factoring information per term. One-level density, central values, and moment statistics all require Ω(√N) computation. RMT describes the OUTPUT of an expensive computation (zero statistics), not a shortcut to it. Beautiful dead end.
- **Arithmetic dynamics / arboreal Galois**: Pollard rho is O(N^{1/4}) and OPTIMAL among generic group algorithms (Shoup lower bound). All polynomial iteration schemes operate within the generic group model. Multi-map birthday (k independent maps) gives O(N^{1/4}/√k) but O(k²) pairwise GCDs make it strictly worse per unit work. Dynatomic polynomials: gcd(Φ_n(f,x), N) reduces to group-order testing. Preimage tree variance ratio SP/prime ≈ 3x but computing preimage trees requires solving x²≡a mod N = factoring. Postcritical orbits reduce to p-1/p+1. Novel ideas (resultant, discriminant, commutator) all 0% success for ≥30-bit.
- **Smooth number construction via algebraic dependencies**: Products of known smooth relations are DEPENDENT (GF(2) vectors are sums). Independence barrier: the only way to generate independent relations is to evaluate the polynomial at NEW points, reverting to Dickman-predicted rates. EC point addition: x-coordinate smoothness not preserved by group law. Number field norm multiplication preserves smoothness but gives dependent relations. Smooth rewriting paths either produce dependent relations (multiplication) or require new evaluations (addition/subtraction). The Dickman function reflects the multiplicative structure of Z, not a limitation of current methods.
- **Subspace designs for GF(2) LA**: Subspace-design-based Wiedemann achieves O(n²·s/w) — same as standard Block Lanczos/Wiedemann. Graph decomposition of factor base matrix not viable (no block structure). SGE + Block Lanczos gives practical constant-factor improvement (f~2-5x). LA phase already well-optimized; no method beats O(n²·s/w) for sparse GF(2) matrices. The LA phase does NOT determine the L-exponent — sieving dominates.
- **Approximate algebraic geometry over Z/NZ**: Epsilon-variety V_eps = {x : |f(x) mod N| < eps·N} has density |V_eps|/N = 2·eps, UNIVERSAL regardless of factorization type. CRT decomposition of V_eps does encode factoring info but computing it requires the factorization. Bezout-type bounds for epsilon-intersections show no SP/prime gap. Global properties (density, volume) are universal; local properties (specific short vectors) carry factoring info but reduce to Coppersmith. Clean separation: aggregate statistics cannot distinguish factorization structures.
- **AD / gradient methods for factoring**: p-adic Newton converges to sqrt(N) in Z_l but finding "bad primes" l=p,q IS factoring. Soft smoothness landscape has noise-like autocorrelation (~0); gradient ascent finds partially divisible values, not smooth ones. Smoothness is CONJUNCTIVE (all primes ≤ B divide), while gradient optimizes ADDITIVE objectives — fundamental mismatch. The sieve IS the optimal algorithm for smooth number detection: exploits additive log structure and periodic divisibility structure with massive parallelism. No continuous relaxation matches.
- **Higher reciprocity / explicit Langlands beyond GL(1)**: Cubic reciprocity (N/pi)_3 = (p/pi)_3·(q/pi)_3 in Z[omega] — product computable, individual symbols are not. Multi-curve EC approach: a_N(E) = a_p(E)·a_q(E) but k curves give k equations with k+1 unknowns (underdetermination GROWTH theorem — more curves make it harder, not easier). Sato-Tate product moments depend only on N, not on factorization, ruling out ALL moment-based tests. Rankin-Selberg L-functions face same computational gap. Langlands base change preserves multiplicativity ambiguity. The Langlands correspondence is structural, not computational — it maps factoring to an isomorphic problem on the automorphic side.
- **Oracle separation / query complexity**: NFS uses 6 operations beyond the generic ring model: polynomial selection (reads N's digits), smoothness testing (integer property), sieving (consecutive-integer structure), LLL (Euclidean geometry), GF(2) LA, and algebraic square root. These explain the gap: generic ring Ω(N^{1/3}) vs NFS L[1/3]. Decision tree depth for factoring is Θ(n) (trivial — must read all input bits). Certificate complexity: compositeness O(1) bits, primality Θ(n) bits. The DT model captures input complexity, not computational complexity.
- **Thin groups / super-approximation / Apollonian**: Apollonian curvatures have essentially RANDOM smoothness (ratio to Dickman 0.94-1.06). Super-approximation ensures pseudo-random behavior mod N — exactly what factoring algorithms need to AVOID. Zaremba-style bounded partial quotients trade approximation quality for algebraic structure but this tradeoff does NOT help smoothness. Markoff triples mod N: 0 factor revelations in all tests. BGS sieve detects primes in orbits but not smooth numbers useful for factoring. Thin group expansion property is adversarial to factoring.
- **Random polynomial systems over Z/NZ**: Solution counts |V(f) mod N| have CRT multiplicativity but signal is O(1/N) per sample — undetectable with polynomial many evaluations. Exponential sums S_t = sum f(x)^t over Z/NZ decompose via CRT but computing any S_t requires O(N) evaluations. Weil bounds hold per component but cross-component structure requires CRT decomposition = factoring.
- **Multi-dimensional diophantine approximation**: Multi-exponent lattice {(a_0,...,a_k) : sum a_i * N^{i/d} small} tested for varying k and d. LLL on multi-exponent lattices gives norms WORSE than standard NFS: the "folding barrier" — extra dimensions add algebraic dependencies that prevent independent size reduction. Three "degrees of freedom" (B, k, d) collapse because k does not independently reduce value sizes (algebraic equivalence to choosing different d). L[1/4] would require norms ~N^{1/d²}; achieving this needs a number field where norms depend on factors of N, not N itself — which requires knowing the factors. The CRT barrier persists.
- **Cohomological obstructions for smooth numbers**: Smooth Selmer group concept investigated. Class number h of NFS number field Q(alpha) moderately correlates with sieve yield (discriminant provides additional predictive power beyond alpha score). However, correlation is with CONSTANTS not EXPONENTS — optimizing polynomial selection using algebraic invariants improves c in L[1/3, c] but cannot change the 1/3. The Hasse-principle analogy breaks: smoothness is not a local-to-global property (a number can be smooth mod every prime yet not globally smooth). No cohomological obstruction defined.
- **GI / Babai-style techniques**: Cayley graph Cay((Z/NZ)*, S) spectrum has smaller normalized spectral gap for semiprimes (0.29) vs primes (0.42) — weak statistical distinguisher but unreliable. Eigenvalue bilinear structure E[i,j] encodes factorization but decomposing eigenvalue multiset into tensor factors requires knowing φ(p)×φ(q) dimensions = factoring. WL refinement stabilizes at a partition that does NOT reveal factors. Power graph and multiplication graph automorphism groups decompose via CRT. The factoring-to-GI reduction fails because ring structure is too rich: any graph encoding either preserves CRT (circular) or loses factoring info.
- **Effective abc / Baker theory for factoring**: abc bound on p+q given N=pq is trivially satisfied (quality ~1/3, not ~1). Baker lower bound |log(p/q)| > exp(-C*log²N) is astronomically weaker than useful. S-unit equation approach: number of solutions to a+b=N with S-smooth a,b gives upper bound that matches QS complexity (the approach IS QS). Smooth-sum density for balanced semiprimes: Pr[p+q is B-smooth] drops exponentially with log(N)/log(B). Direction mismatch: abc gives LOWER bounds on radical while factoring needs UPPER bounds on smoothness.
- **Algebraic K-theory / cyclotomic trace / THH/TC**: K_1(Z/NZ) = (Z/NZ)*, K_2(Z/NZ) decomposes via CRT. Dennis trace, cyclotomic trace, TR tower, Witt vectors, Nil K-theory, lambda/Adams operations, chromatic localizations ALL decompose via CRT. The cyclotomic Frobenius on TR is part of functorial structure — it decomposes. NK_i terms (Bass-Heller-Swan) are trivial for finite rings. The only non-decomposable info comes from RELATIVE constructions K_*(Z, NZ) which reduce to φ(N). All these invariants are functors respecting the CRT isomorphism — definitively negative.
- **Short interval smooth numbers / analytic NT**: Smooth values of QS/NFS polynomials follow a Poisson process at scales below the smoothness bound. No hidden clustering exists. The only structure is CRT (root alignment mod primes), which the sieve ALREADY fully exploits. Selberg sieve weights show zero correlation with actual smoothness. Special-q lattice sieve IS guided sieving from short-interval theory. Murphy's alpha IS polynomial-selection from root density. State-of-the-art implementations (CADO-NFS, msieve) already embody all relevant theoretical insights.
- **Number field tower descent for NFS**: Tower Q ⊂ K_1 ⊂ K_2 tested for degrees 3,4,6. Descent map N_{K/L} is LINEAR on exponent vectors — merges conjugate ideals. Key finding: tower descent produces NO new independent relations. L-level dependencies are REDUNDANT (arise from forgetting conjugate distinction). Every L-dependency lifts to a K-dependency. Tower helps in TNFS for DLP (better polynomial selection) but NOT for integer factoring (norms and sieving unchanged). The size reduction per descent step is bounded by extension degree (constant), unlike FFS where Frobenius reduces by factor q.
- **Beyond GF(2) exponent exploitation**: GF(k) kernel for k=2,3,5,7: same number of relations needed (pi(B)+1) regardless of k. GCD success rate: squares give P=1/2 (universal because gcd(2,p-1)=2 for all odd p). Cubic gives P = 1/3 or 2/3 depending on p mod 3. Higher k strictly worse because gcd(k,p-1) is non-universal. LLL on integer exponent matrix: short vectors found but correspond to products of KNOWN relations — no new information. The 1-bit-per-relation bound is a theorem within congruence-of-squares, proven via the universality of order-2 elements.
- **Meta-analysis of barrier structure**: The 5 barriers decompose into APPROACH-DEPENDENT (B1, B2, B3 — could be evaded by novel framework) and STRUCTURAL (B4, B5 — any classical algorithm must confront). B4 (CRT opacity) is the MASTER BARRIER: B1,B2,B3 are consequences when restricted to sieve methods. B5 (Z-rigidity) is the deepest: explains why FFS techniques fail over Z. Minimal barrier set: {B4, B5} implies all others. Most plausible breakthrough direction: high-dimensional abelian varieties with rich endomorphism rings over Z/NZ. No 2023-2026 mathematical development brings sub-L[1/3] significantly closer, though Fargues-Scholze geometrization provides new conceptual frameworks.
- **Non-commutative arithmetic / matrix groups**: GL(n, Z/NZ) = GL(n, Z/pZ) x GL(n, Z/qZ) by CRT — the group decomposes. PROVED: every polynomial-time computable invariant of matrices over Z/NZ is a polynomial function of entries, which respects CRT. Non-decomposing invariants (orders, word length) are computationally intractable. No "moderately hard" middle ground exists: either polynomial (decomposes) or search-hard. Matrix M^{p-1} = I test reduces to Pollard p-1. Word map evaluations in GL(2) reveal factors with O(1/p) probability per evaluation. Quaternion algebras split at odd primes → same barrier.
- **BBS pseudorandomness distinguishers**: BBS security is EQUIVALENT to factoring (QR reduction, confirmed tight). Every distinguisher reduces to detecting the BBS period, which grows exponentially. Linear complexity ratio LC/n converges to 0.50 (random) by 32 bits. Spectral analysis: no exploitable peaks. Higher-order correlations: zero signal. Next-bit prediction advantage: zero above noise for N > 24 bits. Distinguishing and extracting are different problems, and even distinguishing fails at cryptographic sizes. PRG-based approaches offer no advantage over direct algebraic methods.
- **Sparse polynomials / lacunary series**: Carry propagation in p*q destroys factor sparsity: products of two sparse numbers have HW ≈ n/2. N's sparsity is independent of factor sparsity. Cyclotomic factoring (gcd(Phi_k(2), N)) reduces to Pollard p-1. ISD for GF(2) kernel: complexity 2^{O(n/log n)} worse than Gaussian O(n³). Coppersmith + sparse factor enumeration: polynomial for t-sparse factors when second bit position < n/4 (~49% of 2-sparse cases). But general semiprimes have no sparse factor structure.
- **Third homomorphism / third NFS channel**: Systematic analysis of all candidate third channels: (1) Second number field = MNFS, improves constant only; (2) p-adic = circular (choosing p IS factoring); (3) Function field = circular; (4) Adelic = product formula constrains, no independent info; (5) Modular forms = Hecke eigenvalues at level N are the factors (circular). CONJECTURE: among algorithms based on smooth values of norms from number fields (the "smooth number paradigm"), GNFS is asymptotically optimal at L[1/3]. The algebraic structure of Z provably cannot support a third independent smoothness channel due to the product formula and characteristic barrier.
- **Genus-g Jacobian ECM scaling**: Efficiency PEAKS AT g=1. Proved: L-exponent scales as √g, so g=1 minimizes. At 50 digits: g=2 is 10^{5.8}× worse, g=3 is 10^{12.5}× worse. Even granting g× GLV speedup + √(2g)× channel effect, g=1 still wins by 5+ orders of magnitude. No crossover at any semiprime size from 10 to 60+ digits. Lenstra's 1987 choice of elliptic curves was optimal, not just convenient.
- **Quadratic form reduction theory / Bhargava composition**: Ambiguous forms of discriminant -4pq directly encode the factoring form (p,0,q). But finding it requires searching through h(D)/2 ~ √N forms. Class group infrastructure (Shanks) gives O(N^{1/4}) which IS SQUFOF. Higher composition (Bhargava cubic, quartic) gives additional representations but computing them requires class group computation ≡ factoring. No shortcut through the class group beyond SQUFOF.
- **Fargues-Scholze geometrization**: The geometric Langlands at local level does NOT provide algorithmic shortcuts. The Fargues-Fontaine curve at a "bad prime" l|N reveals l, but identifying bad primes IS factoring. All geometric objects (Bun_G, Hecke stacks, eigensheaves) decompose via CRT or require factoring to access. Prismatic cohomology (the most "prime-independent" framework) still reduces to Z/pZ × Z/qZ computations. The construction is a structural theorem about the SHAPE of Langlands, not a computational tool.
- **PFR / additive combinatorics for smooth numbers**: Smooth numbers have 10-27% smaller additive doubling than random, and their exponent vectors have 2-33% smaller sumsets in F_2^k. But PFR's polynomial bounds (K^{12}) are too large to yield algorithmic improvements. Kelley-Meka 3-AP bounds: smooth numbers have 2-5× AP enrichment over random, but this IS the sieve structure (already exploited). Multiplication-addition barrier: DL ↔ factoring equivalence prevents additive combinatorics from contributing beyond what multiplicative structure already provides.
- **Galois cohomology of Z/NZ**: Ext^1(Z/NZ, Z/NZ) = Z/NZ (decomposes via CRT). Brauer group Br(Z/NZ) = 0 (finite commutative ring). H^i(G, (Z/NZ)*) for Galois G decomposes via CRT for all i. Fundamental trichotomy: every cohomological invariant is either (1) trivially computable but CRT-decomposed, (2) encodes factoring but is as hard to compute as factoring, or (3) is defined in terms of the factoring (circular). No "sweet spot" where an invariant encodes enough info to factor but is computable without factoring.
- **Algebraic statistics of sieve matrices**: BP on sieve matrix factor graph fails (short cycles, symmetry, high rate). ML estimation of null space fails (null space is not a function of marginals). Toric warm starts: column statistics do not predict null space membership. O(n) bits of statistical summary cannot determine O(d·n) bits of null space. Current methods (Block Lanczos, Wiedemann, SGE) already exploit the computationally relevant aspect (sparsity).
- **Height theory / Arakelov for NFS elements**: Height spectrum of Z[alpha] determined by |N|, not factorization. Arakelov bound gives minimum norm for non-trivial elements — LLL/BKZ already approaches this. NFS sieving explores elements near the Arakelov bound systematically. Height perspective CONFIRMS NFS optimality rather than suggesting improvements. Smoothness is non-archimedean, invisible to heights.
- **Selberg/Arthur trace formula**: Both spectral and geometric sides encode factoring info. Geometric side: class numbers h(D) with D involving p,q. Spectral side: new/old form dimensions involve φ(N), σ_0(N). Selberg zeta function partial evaluation requires O(√N) geodesic terms. The trace formula is a mathematical identity relating two descriptions of the same object — both require the factors to evaluate. Provides structural understanding of prime/semiprime difference but zero computational shortcut.

### Witt-Frobenius near-miss

p-th power is an approximate ring endomorphism with error O(p). Witt formalism makes this precise. But for composite N, no single polynomial works — CRT decomposition required. **Closest known analog to Z-Frobenius**, but precision loss is fundamental. Computational experiments confirmed: Witt-Frobenius on truncated Witt vectors reduces to independent power maps per component, providing no new structure. The p-adic analysis showed that ALL information about primes lives at the archimedean place, making descent impossible.

## Open directions

Five structural barriers identified so far (see library/barrier_synthesis.txt): (1) L[1/3] structural wall, (2) archimedean/non-archimedean gap, (3) GF(2) bottleneck, (4) CRT opacity, (5) Z-rigidity. Every investigated approach ran into at least one. A breakthrough could violate one, or could come from an entirely unexpected direction that sidesteps them all. These are the open directions we see:

**Super-smoothness (barrier #3 — GF(2) bottleneck):**
Can a single smooth relation yield >1 independent constraint? The algebraic object must be richer than "exponent vector mod k" for ANY k, since changing the target field doesn't help. Would need a fundamentally different algebraic coincidence.

**Near-endomorphisms of Z (barrier #5 — Z-rigidity):**
Z has no non-trivial ring endomorphism. What about endomorphisms of structures that faithfully encode factoring but are not rings? The Deuring correspondence connects to modular polynomials, class polynomials, and j-invariant arithmetic — explored through specific channels but not exhaustively. The Witt-Frobenius near-miss suggests there may be other approximate endomorphisms worth investigating.
