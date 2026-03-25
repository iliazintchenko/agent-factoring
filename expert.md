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

## Why L[1/3] is the current classical barrier (but likely not the final one)

Every known approach falls into one of three categories:

1. **Smoothness-based** (QS, NFS): L-exponent locked by polynomial degree and Thue's theorem. LP expansion, multi-polynomial, batch GCD, multipliers, Galois action, CRT candidates, Best-of-K — all give at most constant-factor improvements.

2. **Group-order-based** (ECM, p-1/p+1, rho, class group): Complexity depends on factor size p. For balanced semiprimes (p ≈ √N), matches QS at L[1/2]. Division polynomial, Schoof-like, and Deuring correspondence approaches all reduce to this.

3. **Algebraic structure** (Hecke operators, quaternion algebras, isogeny graphs, lattice methods): Either require knowing the factorization to compute (circular), reduce to birthday bounds O(N^{1/4}), or reduce to smoothness methods.

The function field sieve breaks L[1/3] for DLP because of TWO properties absent over Z: (a) Frobenius endomorphism provides systematic degree reduction, (b) polynomials over finite fields factor efficiently (Berlekamp).

## Information-theoretic structure of factoring

Factoring an n-bit semiprime requires ~n/2 bits. Each smooth relation provides exactly **1 useful bit** (one GF(2) vector) despite carrying ~7.7 bits of "surprise" at u=4. Most information in the smooth factorization is wasted — only parity of exponents matters. The 1-bit-per-relation bottleneck is fundamental to the GF(2) linear algebra framework. Breaking it would require a new algebraic framework beyond exponent-parity relations.

Small-subgroup DL projection (g^{(N-1)/l} mod N for l | N-1): recovers O(1) bits total regardless of N's size, because gcd(N-1, φ(N)) is typically tiny. Strictly weaker than Pollard p-1.

Factorizations of N±k for small k: yield O(k·ln B) bits about p, but O(log N) needed. Information gap is exponential.

## Why NFS is stuck at k=2

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

**Concrete target**: L[1/4] requires norm sizes N^{1/d²} instead of N^{1/d}. Verified by optimization: plugging N^{1/d²} into the NFS balance shifts cubic→quartic. Nested number field norms don't achieve this (composition = full norm, an invariant). This is the precise mathematical barrier *within the current framework*.

**Important**: No one has proved L[1/3] is a hard lower bound. The barrier is empirical and heuristic — every known approach hits it, but that doesn't mean every *possible* approach does. The historical progression L[1] → L[1/2] → L[1/3] took decades at each step, and each breakthrough came from an unexpected direction (polynomial structure for QS, number field norms for NFS). The next breakthrough will likely come from a similarly unexpected connection.

## Explored directions (no improvement found yet)

**Smoothness-based approaches explored:**
- **Schnorr lattice factoring**: Lattice dimension grows with smoothness bound. Ducas (CWI): 0 relations in 1000 trials.
- **Lattice crypto techniques (BDD/dual/hybrid/primal) for factoring**: Kannan embedding already in NFS. No useful BDD formulation (Coppersmith needs partial info). Dual attack degenerates to Gaussian elimination (no noise unlike LWE). Hybrid BKZ+MITM adds cost at GNFS optimal saddle point. Primal BKZ: LLL already near-optimal in 2D; high-dim has no useful relation mapping. Fundamental: p is not a short vector in any known poly-dim lattice, and factoring needs π(B) relations not one vector.
- **LLL polynomial selection**: Degree-2 already optimal (shortest vector). Degree-3+ needs NFS.
- **Multi-image NFS**: Improves L[1/3] constants only, not the exponent.
- **Degree-4+ number field norms**: LLL polynomial selection gives norms below N^{2/3} (e.g., d=4 LLL: α=0.45, d=6 LLL: α=0.31). But the N^{2/3} figure is specific to degree-3 base-m; it's not a fundamental barrier. The L[1/3, c] form is preserved because the NFS parameter optimization converges to it regardless of d or polynomial selection. LLL improves c by ~5-10% (constant factor). Cyclotomic/thin fields (||f||=O(1)) give much smaller norms but require N to have special algebraic form (SNFS regime).
- **Multi-multiplier sieve (MMS)**: ~60% more relations (constant factor) but doesn't shrink polynomial values.
- **Multiplicative Lattice Descent (MLD)**: LP expansion from B to B^α multiplies matrix dimension by ≥ R for R× more relations. L-exponent cannot improve — this is why LP variations improve constants but not the exponent.
- **Best-of-K polynomial selection**: Smoothness improvement scales as K^{u/ln(B)} ≈ K^{0.5}, always sublinear in K. Never beats sequential.
- **Smoothness correlations across polynomials**: No significant correlation found — sieve conditions at different x or different multipliers are independent.
- **CRT-structured candidates**: Higher raw smoothness rate is illusory — CRT forces large x, cofactor ≈ 2√N regardless.
- **Algebraic curve parametric families for smoothness**: Tested Pell CF convergents, generalized Pell, Cornacchia near-representations, Pellian cross-terms, and Lehman multi-k. Cornacchia shows 50-500x raw smoothness rate vs QS baseline BUT rank of relation matrix is only 20-190 out of 1229 needed — relations cluster on same few primes. Pell values are individually small (|x²-Ny²|<2√N) but sparse (only O(log N) convergents). Lehman x²-kN is equivalent to MPQS. No family achieves both high smoothness rate AND full-rank relation matrix. Value size remains the dominant factor via Dickman ρ.
- **Batch GCD collision rate**: Drops from 11.3% (40d) to 0.93% (60d). Constant-factor improvement only.
- **Multiplicative lattice relations (MLR)**: LLL finds small P(e) mod N, but CRT entanglement means small mod N ≠ small mod p. Confirms Schnorr is theoretically limited.
- **Correlated norms across number fields**: Measured ~15% positive correlation at 40 digits, but entirely SIZE BIAS (both norms are functions of the same (a,b) coordinates), NOT algebraic correlation. Effect weakens at larger N. Already fully exploited by lattice sieving. Cannot be amplified beyond current NFS techniques.
- **Additive structure of power residues via CRT**: Full enumeration with tiny primes shows power residue sets decompose via CRT but reconstructing the decomposition from the combined set requires evaluating Jacobi sums mod N, which is as hard as factoring.
- **Eigenvalue structure of matrices over Z/NZ**: Eigenvalues DO leak factoring information (they decompose via CRT). But computing eigenvalues over Z/NZ requires polynomial root-finding mod N, which reduces to known methods (Berlekamp → p-1 condition). No new attack beyond p-1/ECM.
- **QS cofactor product smoothness (random self-reducibility)**: For standard 1-LP QS, cofactor products cannot provide additional relations — cofactors are prime by necessity, and their products are not smooth.

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
- **Quaternion norms for smoothness**: Quadratic per variable. Quaternion norm (a²+b²+c²+d²) smoothness only 6% above random baseline for matched magnitude.
- **Non-commutative ring norms (M_2(Z), quaternions, Z[S_3])**: Matrix determinant has ~25% higher smooth fraction than random integers of same size, but this is SIZE BIAS (determinant distribution skewed toward small values), not algebraic advantage. Z[S_3] standard rep norm ~37% higher (same bias). Key theoretical result: any multiplicative norm N:R→Z projects factoring in R to factoring in Z. Each matrix factorization M=M₁M₂ projects to det(M)=det(M₁)·det(M₂), the SAME integer factorization. The fiber structure doesn't reveal new integer factorizations. Non-commutativity provides more PATHS to the same factorizations but doesn't create new ones.
- **Brandt matrices mod N**: For primes m < N, the quaternion norm form a²+b²+Nc²+Nd² forces c=d=0, making representation counts N-independent. The Brandt matrix is identical for all N — no factoring information leaks through spectral structure.
- **Character sum / autocorrelation**: Detecting period p requires Ω(√N) samples. Jacobi symbol autocorrelation tested extensively: semiprime vs prime statistics indistinguishable for M << p. The only "loophole" (small divisors of p-1) is exactly Pollard's p-1 method.
- **Jacobi symbol deconvolution / ICA**: (a/N) = (a/p)(a/q) is a binary ICA problem. Dead end: binary ICA on ±1 signals is information-theoretically impossible (mutual information I(s;x) = 0 when s,t independent Rademacher). Legendre symbol period = factor size (~√N), requiring Ω(√N) evaluations. Pólya-Vinogradov bounds limit correlation detection to O(p) samples. DFT requires resolving frequency 1/p. All paths require Ω(√N) work.
- **Polynomial splitting mod N**: Probability O(1/√N) per trial for monic polynomials.
- **Berlekamp over Z/NZ**: Degenerates to p-1 condition for balanced semiprimes — requires ord(r) | (q-1) in F_p, probability ≈ log(p)/p ≈ 0.
- **Approximate DL relations**: P(gcd reveals factor) ≈ 2/p. Exponential.
- **Bilinear smoothness decomposition**: Can't split norms without number fields.
- **Galois action on NFS relations**: Norm is Galois-invariant. Same relation, not new ones.
- **Sum-product phenomena over Z/NZ**: Zero divisors are too sparse (only ~2√N out of N). Detection requires k ~ √N samples — exponential in factor bit-length. E*(A)/E+(A) ratio indistinguishable from prime for N > 10^6. Extended to algebraically structured orbits (multiplicative, polynomial x²+c, affine mx+b): all orbit types show zero separation between semiprimes and primes. CRT isomorphism preserves both addition and multiplication, so sum-product ratios are ring-isomorphism invariants — fundamentally cannot distinguish Z/NZ from Z/pZ × Z/qZ.
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

Sieving: O(1) amortized per candidate but requires sequential memory access. Bernstein's batch GCD: O(log²B) per candidate but works on ARBITRARY candidate sets. For large B, batch GCD is cheaper and enables non-sequential candidate generation. However, testing showed sequential sieving wins because it naturally produces the smallest candidates (near √N). Batch GCD is only advantageous for external candidate sources, not single-N factoring.

## Research survey

- **No classical sub-L[1/3] algorithm exists** for general integer factoring (confirmed by 2020-2025 literature survey).
- **Henry Cohn (MIT)**: Argues no known barrier prevents progress below L[1/3], noting historical 1→1/2→1/3 progression. No concrete algorithm proposed.
- **Regev (2023)**: Quantum O~(n^{3/2}) gates (from Shor's O(n²)). Follow-ups by Ragavan-Vaikuntanathan, Pilatte, Ekerå-Gärtner. Hhan (EUROCRYPT 2025) proved matching lower bound in generic ring model. No classical dequantization — factoring's HSP structure resists classical simulation, the gap O(√N) vs O(polylog N) is robust.
- **Tower NFS (Barbulescu-Kim 2016)**: Sub-L[1/3] for DLP in GF(p^n) only. No factoring analog.
- **Harvey & Hittmeir (2020-2022)**: Best deterministic factoring N^{1/5+o(1)}. Exponential.
- **Schnorr lattice (2021 claim)**: Refuted experimentally by Ducas (CWI), 0/1000 relations. Never published at peer-reviewed venue.
- **GNFS L[1/3, 1.923] unchanged since 1990s**: RSA-250 (Feb 2020) remains factoring record. All progress is in practical constants (better sieving, batch smoothness).
- **NFS L[1/3] is HEURISTIC, not proven**: Best rigorous bound is L[1/2] (Lenstra-Pomerance 1992). Lee-Venkatesan (2017) proved L[1/3] only for randomized variant finding square congruences. Key unproven assumptions: smoothness heuristic ("norms behave like random integers"), independence of smoothness events, monogenic field, square root step.

## Open directions

~330 approaches have been investigated. Five structural barriers have been identified (see library/insights.txt). A breakthrough would need to violate at least one. These directions remain open:

**Violating the GF(2) bottleneck (barrier #3):**
- Every known method extracts 1 bit per relation via exponent-parity. GF(3) / ternary relations now thoroughly analyzed: (a) standard QS produces QUADRATIC relations (a² ≡ Q(x) mod N), so GF(3) kernel gives ∏Q_i = b³, yielding a² ≡ b³ — a MIXED degree relation, not c³ ≡ d³, so the cubic factorization (c-d)(c²+cd+d²) doesn't apply; (b) a dedicated cube sieve needs N^{2/3}-size values, L-exponent c ≈ 1.155 vs QS c ≈ 1.0; (c) GF(3) rank is typically HIGHER than GF(2) rank, so more relations needed; (d) information entropy per entry ≠ useful information (kernel dimension depends on rank, not per-entry entropy); (e) Z/6Z (mixed) is strictly harder (intersection of kernels). **Closed: GF(k) for k > 2 strictly worse due to fundamental mismatch between quadratic sieve structure and higher-order algebra.**
- The "super-smoothness" question (open_problems.txt, Problem 3): can a single smooth relation be made to yield >1 independent constraint? This remains open but now more precisely constrained: the algebraic object must be richer than "exponent vector mod k" for ANY k, since changing the target field doesn't help. Would need a fundamentally different algebraic coincidence, not just a different modulus for exponents.

**Violating Z-rigidity (barrier #5):**
- Z has no non-trivial ring endomorphism. But what about near-endomorphisms, or endomorphisms of structures that are not rings but faithfully encode factoring? The Deuring correspondence was explored but only through three specific channels — endomorphism ring computation is circular, BUT the correspondence also connects to modular polynomials, class polynomials, and j-invariant arithmetic, which were not deeply explored computationally.
- Quaternion algebras over Q ramified at p have computable structure (Brandt matrices were tested and found N-independent for small primes). Pizer's Ramanujan graphs now tested via supersingular j-invariant arithmetic over Z/NZ: 6 approaches (polynomial GCD with x^N-x, discriminant GCD, direct modular polynomial evaluation, iterated Frobenius Y^{N^k}, cross-GCD of Φ_l, and Hecke trace walk). Approaches 1-5 all fail with O(1/√N) success probability — the CRT opacity barrier blocks detection of isogeny graph structure differences. Approach 6 (Hecke trace walk) factors successfully but IS Pollard rho with specific polynomial T(x) = x² - 1488x + 162000 — generic, not isogeny-specific. **Key insight: isogeny graph differences ARE real mod p vs mod q, but computationally inaccessible through polynomial operations. No "smoothness amplification" mechanism exists for graph walks (unlike ECM which has group order smoothness).** Closed.
- **Hilbert class polynomials**: Now tested for h(D) = 1,2,3,4,5. Result: 0% success across all 125+ (D, N) pairs. Even when Kronecker symbols differ ((D/p) ≠ (D/q)), gcd(H_D(x), x^N - x) mod N = gcd(H_D(x), x^{N} - x) gives degree 0. The reason: x^N mod H_D(x) mod p computes x^{q mod (p-1)} (via Fermat), NOT the Frobenius x^p. So the test checks whether roots of H_D satisfy x^{q-1} ≡ 1 mod p, which requires specific divisibility conditions on p±1 — reducing to Pollard p+1 method. **Closed: reduces to p+1 method variant.**

**Violating CRT opacity (barrier #4):**
- CRT composition is opaque to all tested approaches. But the tensor decomposition investigation found the rank-2 structure IS detectable — it's just destroyed by mod-N wrapping for samples > √N. What if there's a way to work with structured samples that avoid wrapping? This would require a non-random sampling strategy correlated with the CRT decomposition — which sounds circular, but iterative/adaptive approaches might gradually accumulate partial CRT information.
- The sum-product investigation found E*(A)/E+(A) is indistinguishable for N > 10^6. This now extends to algebraically structured sets (multiplicative orbits, polynomial orbits, affine orbits) — all indistinguishable. CRT ring isomorphism preserves sum-product ratios.

**Approximate / relaxed approaches (tested, no improvement):**
- **LLL on cofactor log-lattice**: LLL minimizes sum(e_i * log(r_i)) in the reals, but smoothness is a discrete/arithmetic property. Products of distinct large primes never cancel to smooth values — LLL found 0 smooth combinations across all test cases. The continuous relaxation doesn't capture integer factorization structure. Known LP matching (exact collision) remains optimal.
- **Cofactor correlation in QS sieve**: Tested whether nearby sieve positions (x, x+k) have correlated cofactors. Result: cofactors are completely independent random integers. GCD rate = 0 (stochastic flukes only). Mod-l residue bias fully explained by "coprime to factor base" selection effect. Log-size correlation r < 0.02 at all offsets. Nearby products NOT smoother than random (ratio 0.97). Cofactor uniqueness 99.99%+. Confirms the standard independence assumption used in large-prime analysis.
- Algebraic curve families (Pell, Cornacchia, Lehman): high raw smoothness rates possible but relation rank is low (clustered on few primes). No family achieves both abundance and full-rank relations.

**Cross-disciplinary connections (mostly explored, see below):**
- Additive combinatorics: Freiman/Plünnecke-Ruzsa/BSG applied to QS smooth numbers. Result: smooth x-positions have Freiman dimension ~7 (the sieve lattice structure), dc/pred ≈ 0.37 vs 0.55 random. BUT this IS the sieve — QS already exploits it. No semiprime-vs-prime distinction. Q(x) values themselves have zero additive structure (dc/pred = 1.0000). Sumset smoothness enrichment is 1.2-1.4x (trivial multiplicative bias). PR/BSG bounds vacuous (K ≈ 93 too large). **Closed: additive combinatorics adds nothing beyond QS sieve.**
- Algebraic topology: nerve complex of factor base coverage IS non-trivial (β_1 ~ 10-30, β_2 ~ 100-300). BUT Betti numbers do not distinguish semiprimes from primes — they depend on smooth-value density and factor base size, which scale identically. CRT decomposition visible in sub-complex structure but equivalent to standard sieve data. **Closed: topology adds no information beyond the relation matrix.**
- Proof complexity: factoring certificates are O(n) bits (Pratt), but proving NON-factorability (of primes) requires AKS-class work. The asymmetry now analyzed: MR witnesses carry O(log log N) bits about the 2-adic structure of p-1, q-1. The extended MR squaring chain IS Pollard's p-1 restricted to the 2-part — factors N in poly time when v_2(p-1) ≠ v_2(q-1) (~2/3 probability). But general factoring requires probing ALL prime components of the group order, not just 2. Multi-witness combinations reduce to finding nontrivial square roots of 1. **The gap persists because compositeness detects a thin slice (CRT mismatch at ONE point) while factoring needs full resolution of the multiplicative group structure.**

## Compact catalog (~330 approaches)


See above sections for details. Summary of ~310 investigated approaches across categories:

**Smoothness-based (all L[1/2] or L[1/3], no improvement):** Schnorr lattice, LLL polynomial selection, multi-image NFS, 3-field NFS, MMS, MLD, Best-of-K, CRT candidates, batch GCD, MLR, K-large-prime, CF-base descent, multi-d CF, sparse polynomial NFS, group algebra Z[C_k], lattice tower/BKZ, lattice on O_K, batch GCD + non-sequential, higher-order residues GF(k), smooth bit-patterns, algebraic curve families (Pell/Cornacchia/Lehman).

**Group-order-based (all L[1/2] in p or worse):** ECM, p-1/p+1, CM-ECM, iterated Frobenius, class group 2-Sylow, division polynomials, Schoof-like, algebraic tori T_n, genus-2 HECM, T_6/XTR compression, Weil/Tate pairings, dynamical systems (Chebyshev/Lattès), curve point-counting.

**Algebraic/structural (all blocked):** Spectral/Cayley, Deuring/isogeny, Hecke operators, Brandt matrices, character sums, polynomial splitting, Berlekamp, DL relations, bilinear decomposition, Galois action, sum-product (random + structured orbits), tensor decomposition, p-adic lifting, Groebner, batched p-1, Newton/Hensel, isogeny walks, non-abelian HSP, subfield lifting, virtual Frobenius, higher residue symbols, exterior algebra/SNF, counting arguments, modular forms/theta, carry propagation, algebraic varieties, RS codes over Z/NZ, synthetic Frobenius, tropical geometry, higher reciprocity laws, Cayley graph word metric, Kloosterman/exponential sums, Newton identities/power sums (s_k=p^k+q^k needs s_1=p+q = circular), zeta function of Z/NZ (evaluation at any s encodes factorization), nerve complex of factor base (non-trivial topology but no factoring info), additive combinatorics/Freiman/BSG (rediscovers sieve structure), Hilbert class polynomials (reduces to p+1 method), supersingular j-invariant arithmetic (O(1/√N) barrier, Hecke walk = Pollard rho), non-commutative ring norms / M_2(Z) / Z[S_3] (multiplicative norm functor projects back to Z), GF(3)/ternary relations (degree mismatch, worse L-exponent), Jacobi symbol deconvolution/ICA (binary ICA impossible, Ω(√N) barrier), MR witness fingerprinting (= p-1 restricted to 2-part), degree-4+ norm optimization (constant factor, preserves L[1/3]), smooth value structure mining (gaps/APs/co-occurrence all match random, sieve already optimal).

**Geometric/cohomological (all reduce to classical):** Arakelov geometry (product formula = same norms as NFS, arithmetic RR gives same counts), condensed mathematics (Z/NZ is discrete, all frameworks trivialize), prismatic cohomology (requires choosing p first = circular), derived algebraic geometry, Shimura varieties, motivic integration, weight 3/2 modular forms/Shimura correspondence (coefficients encode factoring info but extraction has interpretation/computation/decomposition barriers).

**Analytic/approximation (all vacuous or blocked):** abc conjecture (N=pq squarefree makes all bounds trivial), Baker theory (|log(p/q)| bound gives 0 search-space bits), stochastic resonance (smoothness has no threshold dynamics; noise hurts via Jensen), simulated annealing (cofactor landscape has no spatial correlation), Gauss sum spectra (|G(χ)|² reveals factors but requires O(N) character enumeration), wavelet/multiresolution analysis (detects structure but signal weak, equivalent to trial division), Dedekind zeta of Q(√N) (requires O(√N) terms).

**Information-theoretic/complexity (barriers confirmed):** 1-bit/relation bottleneck, small-subgroup DL, nearby factorizations, oracle classification, interactive proofs, ML prediction, Coppersmith + partial info, precomputation S·T bound, analytic function impossibility, classical QFT analog, period-finding dequantization, non-standard computation, hybrid quantum, Z-descent barrier, L[1/3] robustness, Kolmogorov complexity (K(p|N)=O(1) but K^poly(p|N) large iff factoring hard — restates problem), Pell/class groups (genus theory gives only O(1) bits, CF period indistinguishable SP vs prime), communication complexity (CC bounds O(log log N) only, Coppersmith n/4 certificate tight for lattice methods), entropy/MI analysis (each trial yields ~10^{-3} bits, MI/cost optimum matches NFS parameter selection), SDP/LP relaxations (weak due to global carry structure — Lasserre level O(1) insufficient), descriptive complexity (factoring definable in ESO∩USO but FO+LFP definition = P algorithm = open).

**Dynamical/ergodic (all equivalent to known methods):** Ergodic theory of multiplication maps (mixing/orbit properties indistinguishable SP vs prime for generic observables, symbolic dynamics also indistinguishable), proof/circuit complexity (factoring likely not in AC0/TC0 but natural proofs barrier blocks proving it, OBDD exponential for mult bits, circuit inversion hits treewidth barrier).

**Parameterized/structural (no improvement):** FPT with k=log(p) via ECM at exp(√(k log k)), Coppersmith FPT with gap parameter, NFS NOT FPT in factor size. Matroid theory of relations (standard GF(2) linear matroid, no useful second matroid). Additive combinatorics / sumsets (QR doubling constant encodes CRT but information flows factorization→structure only). AG codes (need base field, Z/NZ not a field). Waring / r_4(N) (computing divisor sums = factoring). Compressed sensing (2-sparse divisor indicator, O(log N) measurements suffice info-theoretically but computable measurements circular). Fixed-point / PPAD (factoring likely in PPAD/PPP, finding fixed points still hard). Topos theory (functoriality barrier, effective topos respects Turing complexity). FI-modules / representation stability (wrong algebraic side). Derandomization (NFS randomness not the barrier, smoothness heuristic is). Model theory / ultraproducts (zero divisors not FO-definable without naming). Diophantine geometry / rational points (genus 0, Faltings inapplicable, Chabauty needs genus>rank). Tensor networks (Shor has O(n) entanglement, MPS/PEPS insufficient). FHE analogy (CRT is natural encryption, bootstrapping=factoring). Spin glass (glassy landscape, BP diverges, SA fails above ~15 bits). Incidence geometry / ST over Z/NZ (no statistical difference SP vs prime). PIT (universal vs existential statement mismatch). Bio-inspired (flat fitness landscape, no gradient). NCG/Connes (spectral action = character theory, reformulation only). HoTT (all homotopy trivial, complexity-oblivious). Dequantization (exponential gap robust, Tang-style inapplicable). Reverse math (FTA in RCA_0, no large cardinals needed). Sheaf theory on divisibility poset (descriptive not prescriptive). Game theory (query complexity n/2 optimal). Analytic interpolation / Carlson (power sums need s_1=p+q). Abstract interpretation (interval/octagon/polyhedra all give known bounds). ML on binary representation (no exploitable statistical regularity, carry chain destroys locality).
