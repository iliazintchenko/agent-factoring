# Expert Knowledge

## The smoothness bottleneck

All known sub-exponential factoring algorithms rely on finding **smooth numbers**. The L-exponent is determined by the size of values being tested for smoothness:
- QS: tests values of size ~√N → L[1/2]
- NFS: uses algebraic number fields to reduce values to ~N^(2/3) → L[1/3]
- Hypothetical L[1/4]: would need values of size ~N^(1/4) or equivalent

Any approach that still tests numbers of size √N for smoothness will remain L[1/2]. To do better, you must either generate smaller candidates, or avoid smoothness entirely.

The quasi-polynomial DLP breakthrough (Barbulescu-Gaudry-Joux-Thomé 2013) shows L[1/3] barriers CAN be broken for related problems in function fields. The technique exploits systematic factorization of polynomials over finite fields — no known analog over Z. The "translation problem" from function fields to number fields is a major open question.

## Schnorr lattice factoring (reviewed, does not scale)

Schnorr (2021) proposed reducing factoring to SVP in a lattice built from smooth number factorizations. Claims polynomial time, but the lattice dimension grows with the smoothness bound (~50K dimensions for 90-digit numbers). LLL/BKZ is infeasible at that scale. Confirmed by Ducas (CWI) — 0 relations in 1000 trials with Schnorr's claimed parameters.

**Lesson**: Any lattice approach that depends on a large smooth factor base inherits the same sub-exponential dimension growth. Haven't tested yet — maybe worth exploring to understand the failure mode empirically.

## Smooth Subsum Search — Hittmeir 2023 (implemented, BETTER SCALING than MPQS beyond 50 digits)

Uses CRT to construct values guaranteed divisible by several factor base primes, then tests the smaller cofactor for smoothness. Replaces sieving with structured candidate generation. L[1/2] complexity class.

**Benchmark results** (worst of 5 semiprimes, single core, seed=42):

| Digits | MPQS | SSS | Ratio |
|--------|------|-----|-------|
| 30 | 0.055s | 0.06s | 1.1x |
| 40 | 0.288s | 0.78s | 2.7x |
| 50 | 1.79s | 17.4s | 9.7x |
| 55 | 25.3s | 55.0s | 2.2x |
| 60 | 72.9s | 147.1s | 2.0x |

**Key finding**: SSS has dramatically better scaling from 50 to 60 digits:
- MPQS 50→60: 40.7x increase
- SSS 50→60: 8.4x increase

L[1/2] fit extrapolation predicts SSS overtakes MPQS around **75-80 digits**. This is because SSS generates candidates with guaranteed divisibility by ~6 factor base primes, making the remaining cofactor smaller. At larger sizes, this structural advantage compensates for the per-candidate overhead.

**Optimization**: Replaced hash map collision counting with sorted array — 33% speedup. Further optimization possible: batch product/remainder tree for smoothness testing (currently using direct trial division).

## Spectral methods on Cayley graphs (dead end classically)

The Cayley graph of (Z/NZ)* with small prime generators encodes the group structure. Eigenvalues would reveal the factorization. **Problem**: the group has φ(N) elements, so computing the spectrum is exponential classically. This is exactly what Shor's QFT solves. No known way to extract useful spectral information without quantum superposition.

## Research survey

- **No classical sub-L[1/3] algorithm exists** for general integer factoring.
- **Regev (2023)**: Improved quantum factoring — O~(n^{3/2}) gates vs Shor's O~(n^2). Fundamentally quantum, no classical dequantization known.
- **Stange (2022)**: Multiplicative relation framework — reduces factoring to index calculus in (Z/NZ)*. Clean formulation but achieves same L[1/2] or L[1/3] with NFS techniques. Not an improvement.
- **Umans & Wang (2025)**: Deterministic N^{1/6+o(1)} factoring conditional on a number-theoretic conjecture. Still exponential.
- **Harvey-Hittmeir**: Best deterministic bound N^{1/5+o(1)}. Still exponential.

## Failed approaches

- **Dixon's random squares with batch smoothness**: Random x² mod N values are ~N in size, far too large for smoothness. Cannot work without polynomial structure to reduce value size.
- **CFRAC**: L[1/2] but single expansion limits relation generation rate. Our CFRAC implementation works but is ~10-50x slower than SIQS due to sequential nature. Not competitive with multi-polynomial QS above 30 digits.
- **ECM for balanced semiprimes**: Complexity depends on smallest factor. For balanced semiprimes, factors are ~N^(1/2), so ECM is L[1/2] in N and not competitive above ~55 digits. Our scaling data confirms this: ECM at 55 digits takes 43s.
- **GNFS for small numbers (< 60 digits)**: NFS algebraic norms are dominated by polynomial coefficients (~N^{1/d}). For degree-3 on a 30-digit number, norms are ~70 bits vs factor base of 400, making algebraic smoothness probability essentially zero. NFS only helps for numbers > ~70-80 digits where the algebraic norm advantage outweighs the overhead.

## Open directions

These are starting points, not an exhaustive list.

- **Algebraic group structure**: Z_N* ≅ Z_{p-1} × Z_{q-1} but we can't see this decomposition. Can random walks, character sums, or higher-dimensional algebraic groups reveal it without smooth numbers?
- **Third image / multi-image NFS**: NFS uses two images (rational and algebraic). Could a third image reduce norms further, yielding L[1/4]? Literature (Barbulescu et al.) shows MNFS with k fields achieves L[1/3, c_k] with improving constants but cannot break the L[1/3] barrier for general integers. The reason: simultaneous smoothness across multiple fields constrains the sieve region proportionally.
- **Smooth polynomial families**: Do certain polynomial constructions yield systematically smoother values than random?
- **Batch smoothness (Bernstein)**: Product/remainder trees for smoothness detection. Better I/O complexity than sieving for large bounds. Could improve constants even if exponent stays the same.
- **Function field analogies**: Can techniques from the quasi-polynomial DLP breakthrough be adapted to number fields?
