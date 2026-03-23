# Expert Knowledge

## The smoothness bottleneck

All known sub-exponential factoring algorithms rely on finding **smooth numbers**. The L-exponent is determined by the size of values being tested for smoothness:
- QS: tests values of size ~√N → L[1/2]
- NFS: uses algebraic number fields to reduce values to ~N^(2/3) → L[1/3]
- Hypothetical L[1/4]: would need values of size ~N^(1/4) or equivalent

Any approach that still tests numbers of size √N for smoothness will remain L[1/2]. To do better, you must either generate smaller candidates, or avoid smoothness entirely.

The quasi-polynomial DLP breakthrough (Barbulescu-Gaudry-Joux-Thomé 2013) shows L[1/3] barriers CAN be broken for related problems in function fields. The technique exploits systematic factorization of polynomials over finite fields — no known analog over Z. The "translation problem" from function fields to number fields is a major open question.

## Tested approaches

### Smooth Subsum Search — Hittmeir 2023 (better scaling than MPQS beyond 50 digits)

Uses CRT to construct values guaranteed divisible by several factor base primes, then tests the smaller cofactor. Still L[1/2] but with better scaling constants.

| Digits | MPQS | SSS | Ratio |
|--------|------|-----|-------|
| 30 | 0.055s | 0.06s | 1.1x |
| 40 | 0.288s | 0.78s | 2.7x |
| 50 | 1.79s | 17.4s | 9.7x |
| 55 | 25.3s | 55.0s | 2.2x |
| 60 | 72.9s | 147.1s | 2.0x |

SSS scales ~8.4x per +10 digits (50→60) vs MPQS's ~40.7x. Extrapolation suggests SSS overtakes MPQS around 75-80 digits. Still needs batch product/remainder tree optimization from the paper — current implementation uses direct trial division for cofactor smoothness.

## Dead ends

- **Schnorr lattice factoring (2021)**: Lattice dimension grows with smoothness bound (~50K for 90d). Ducas (CWI) confirmed: 0 relations in 1000 trials. Any lattice approach depending on a large factor base inherits the same dimension blowup.
- **Spectral methods on Cayley graphs**: Cayley graph of (Z/NZ)* encodes group structure, but computing the spectrum is exponential classically — this is exactly what Shor's QFT solves.
- **Multi-image NFS**: Adding more number fields (Barbulescu et al.) improves L[1/3] constants but cannot break the L[1/3] barrier. Simultaneous smoothness across fields constrains the sieve region proportionally.
- **Dixon's random squares**: Values are ~N in size, far too large for smoothness without polynomial structure.
- **CFRAC**: Single expansion limits relation rate. Not competitive above 30 digits.
- **ECM for balanced semiprimes**: L[1/2] in N for balanced factors. Not competitive above ~55 digits.

## Research survey

- **No classical sub-L[1/3] algorithm exists** for general integer factoring.
- **Regev (2023)**: Quantum factoring in O~(n^{3/2}) gates. Fundamentally quantum, no classical dequantization known.
- **Stange (2022)**: Index calculus in (Z/NZ)*. Clean but same L[1/2] or L[1/3]. Not an improvement.
- **Umans & Wang (2025)**: Conditional deterministic N^{1/6+o(1)}. Still exponential.
- **Harvey-Hittmeir**: Deterministic N^{1/5+o(1)}. Still exponential.

## Open directions

These are starting points, not an exhaustive list.

- **Algebraic group structure**: Z_N* ≅ Z_{p-1} × Z_{q-1} but we can't see this decomposition. Can random walks, character sums, or higher-dimensional algebraic groups reveal it without smooth numbers?
- **Smooth polynomial families**: Do certain polynomial constructions yield systematically smoother values than random?
- **Batch smoothness (Bernstein)**: Product/remainder trees for smoothness detection. Better I/O complexity than sieving for large bounds.
- **Function field analogies**: Can techniques from the quasi-polynomial DLP breakthrough be adapted to number fields?
