# Expert Knowledge

## The smoothness bottleneck

All known sub-exponential factoring algorithms rely on finding **smooth numbers**. The L-exponent is determined by the size of values being tested for smoothness:
- QS: tests values of size ~√N → L[1/2]
- NFS: uses algebraic number fields to reduce values to ~N^(2/3) → L[1/3]
- Hypothetical L[1/4]: would need values of size ~N^(1/4) or equivalent

Any approach that still tests numbers of size √N for smoothness will remain L[1/2]. To do better, you must either generate smaller candidates, or avoid smoothness entirely.

## Why degree-2 is special (Thue's theorem)

The CF of √N produces convergents with **bounded residues**: |p_k² - N·q_k²| < 2√N regardless of k. This infinite family of small residues is what makes QS possible.

For degree ≥ 3, Thue's theorem proves |x^d - N·y^d| < C has **finitely many** solutions. No infinite family of small higher-degree residues exists over Z. The ONLY way to use higher-degree polynomials is via **number fields** (NFS), where the norm provides a different notion of size.

## Why higher-dimensional algebras don't help

Quadratic extensions, quaternion algebras, and higher — all have norms that are **degree 2 in each variable**. More variables don't reduce norm size. NFS achieves L[1/3] not by having a better norm, but by splitting across two different norm maps and requiring simultaneous smoothness. Multi-image NFS only improves L[1/3] constants.

## Why L[1/3] is the fundamental classical barrier

**The L-exponent formula**: α = 1/(k+1) where k = number of independent reduction stages.
- QS: k=1 (smoothness only) → α = 1/2
- NFS: k=2 (polynomial degree + smoothness) → α = 1/3
- FFS: k=3+ (+ tower descent via Frobenius) → α → 0

**Why NFS is stuck at k=2**: A 3rd stage requires an iterable compression map that reduces element size while preserving smoothness structure. Integer size is **archimedean** — cannot be iteratively reduced by algebraic maps. Polynomial degree is **non-archimedean** — Frobenius + quotient reduction decreases it freely. This asymmetry is the deep reason FFS reaches quasi-polynomial while NFS cannot break L[1/3].

**Constants**: GNFS c = (64/9)^{1/3} ≈ 1.923 (tight). MNFS floor: c = (32/9)^{1/3} ≈ 1.526 (tight). All improvements affect only sub-leading terms.

Every known approach falls into three categories, each with a barrier:
1. **Smoothness-based** (QS, NFS): L-exponent locked by polynomial degree and Thue's theorem.
2. **Group-order-based** (ECM, p-1/p+1, rho): For balanced semiprimes, matches QS at L[1/2].
3. **Algebraic structure** (Hecke, quaternions, isogenies, lattices): Circular, O(N^{1/4}), or reduces to smoothness.

The function field sieve breaks L[1/3] for DLP via two properties absent over Z: (a) Frobenius endomorphism for degree reduction, (b) efficient polynomial factoring (Berlekamp).

## Information-theoretic structure

Each smooth relation provides exactly **1 useful bit** (GF(2) parity vector) despite ~7.7 bits of "surprise." The 1-bit-per-relation bottleneck is fundamental to exponent-parity linear algebra. GF(2) is optimal — higher-order residues (GF(3), GF(5)) shrink the factor base, making smooth values exponentially rarer.

## What would be needed to beat L[1/3]

- A new algebraic structure over Z where "norms" are smaller than N^{2/3} per component
- OR a non-smoothness-based sub-exponential approach (none known)
- OR a framework beyond GF(2) exponent-parity that extracts >1 bit per relation
- OR a "descent" mechanism for Z analogous to Frobenius. Z has no non-trivial ring endomorphism (initial object in Ring).

**Concrete target**: L[1/4] requires norm sizes N^{1/d²} instead of N^{1/d}. Nested number field norms don't achieve this (composition = full norm, an invariant).

## Research survey

- **No classical sub-L[1/3] algorithm exists** for general integer factoring (2020-2025 literature survey).
- **GNFS L[1/3, 1.923] unchanged since 1990s**. RSA-250 (Feb 2020) remains factoring record.
- **NFS L[1/3] is heuristic, not proven**: Best rigorous bound is L[1/2] (Lenstra-Pomerance 1992).
- **Henry Cohn (MIT)**: No known barrier prevents progress below L[1/3]. No concrete algorithm proposed.
- **Regev (2023)**: Quantum O~(n^{3/2}) gates. Hhan (EUROCRYPT 2025) proved matching lower bound. No classical dequantization — the gap is genuinely exponential.
- **Tower NFS (Barbulescu-Kim 2016)**: Sub-L[1/3] for DLP in GF(p^n) only. No factoring analog.
- **Harvey & Hittmeir (2020-2022)**: Best deterministic factoring N^{1/5+o(1)}. Exponential.
- **FFS quasi-polynomial DLP** (Barbulescu-Gaudry-Joux-Thomé 2013): Exploits Frobenius + efficient polynomial factoring. No analog over Z.

## Open directions

After ~290 investigations, the remaining open space is extremely narrow:

- **Non-smoothness-based sub-exponential relations**: A fundamentally different relation type extracting >1 bit per relation could change the game. No candidate is known.
- **Approximate/relaxed algebraic structures**: Most approaches seek exact algebraic properties. Relaxed versions with controlled error might open new avenues. No concrete proposal exists.
- **Unpredictable mathematical breakthroughs**: New connections to algebraic geometry (motives, étale cohomology), additive combinatorics, or the Langlands program.

## Explored directions (~290 approaches, all concluded)

**Smoothness-based (all L[1/2] or L[1/3]):** Schnorr lattice, LLL polynomial selection, multi-image NFS, 3-field NFS, MMS, MLD, Best-of-K, CRT candidates, batch GCD, MLR, CF-base descent, lattice tower/BKZ, lattice on O_K, higher-order residues GF(k), smooth bit-patterns, correlated norms (size bias only).

**Group-order-based (all L[1/2] in p):** ECM, p-1/p+1, CM-ECM, iterated Frobenius, class group 2-Sylow, division polynomials, Schoof-like, algebraic tori T_n, genus-2 HECM, T_6/XTR, Weil/Tate pairings, dynamical systems.

**Algebraic/structural (all blocked):** Spectral/Cayley, Deuring/isogeny (3 channels), Hecke operators (circular), Brandt matrices (N-independent), character sums (Ω(√N) samples), polynomial splitting (O(1/√N)), Berlekamp over Z/NZ (degenerates to p-1), sum-product (√N samples), tensor decomposition (1/√N SNR), p-adic/Newton (chaotic), Groebner (2/√N rate), isogeny walks (reduces to Berlekamp), non-abelian HSP (order 4, trivial), subfield lifting (CRT is the secret), virtual Frobenius (O(√N)), L-functions/class numbers (equivalent to factoring), modular forms/theta (O(N) to compute), RS codes (random baseline), approximate DL (exponential).

**Geometric/cohomological:** Arakelov (same norms as NFS), condensed math (trivializes), prismatic cohomology (circular), derived AG, Shimura varieties, motivic integration, weight 3/2 forms (extraction barriers).

**Analytic/approximation:** abc conjecture (trivial for N=pq), Baker theory (0 bits), stochastic resonance (no threshold), simulated annealing (no spatial correlation), Gauss sums (O(N) enumeration), wavelets (equivalent to trial division).

**Information-theoretic/complexity:** 1-bit/relation bottleneck, precomputation S·T bound, period-finding dequantization (exponential gap), Kolmogorov complexity (restates problem), communication complexity (O(log log N) bounds only), SDP/LP relaxations (carry structure defeats).

**Parameterized/structural:** FPT via ECM, matroid theory (standard GF(2)), additive combinatorics (one-way information), Waring/r₄(N) (=factoring), compressed sensing (circular measurements), PPAD, topos theory, FI-modules, model theory, Diophantine geometry, tensor networks, FHE analogy, spin glass, NCG/Connes, HoTT (trivial), reverse math (RCA₀), sheaf theory, game theory, abstract interpretation, ML on bits (no pattern).
