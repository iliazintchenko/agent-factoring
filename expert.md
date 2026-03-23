# Expert Knowledge

## Why factoring is hard: the smoothness bottleneck

All known sub-exponential factoring algorithms (QS, NFS, CFRAC) rely on finding **smooth numbers** — numbers that factor completely over a small prime base. The fundamental tension:

1. **Finding smooth numbers**: The probability that a random number near N^(1/2) is B-smooth is roughly ρ(u) where u = log(N^(1/2))/log(B) and ρ is the Dickman function. Larger B makes smoothness more likely but...
2. **Linear algebra over GF(2)**: ...a larger factor base means a larger matrix to reduce.

The optimal B balances these: QS gets L[1/2, 1+o(1)], NFS gets L[1/3, (64/9)^(1/3)] by using algebraic number fields to generate numbers that are "smaller" relative to their smoothness bound.

**Key insight**: The smoothness bottleneck comes from the _size_ of the numbers being tested. If we could generate candidate relations where the numbers to factor are much smaller than sqrt(N), we could get a better complexity. NFS achieves this by working in number fields; is there a way to do even better?

**Possible escape routes:**
1. Find smooth numbers without sieving (batch smoothness, algebraic construction)
2. Use a third "image" to split values further (could give L[1/4]?)
3. Avoid smooth numbers entirely (period-finding, lattice methods, spectral methods)
4. Exploit special structure of balanced semiprimes specifically

## Schnorr lattice factoring (reviewed, does not scale)

Schnorr (2021) proposed reducing factoring to finding short vectors in a lattice constructed from the prime factorizations of smooth numbers near sqrt(N). The idea: if you can find a short enough vector, it encodes a multiplicative relation that yields a factor.

The paper claims polynomial time, but analysis of the approach shows the lattice dimension grows roughly as the number of primes in the factor base, which for 90-digit numbers means ~50K dimensions. LLL/BKZ reduction in 50K dimensions is completely infeasible.

**Conclusion**: The lattice dimension is tied to the smoothness bound, which grows sub-exponentially with N. So this doesn't escape the smoothness bottleneck — it just reformulates it as a lattice problem that's equally hard. Any lattice-based approach that depends on a large smooth factor base will have the same problem. Tested by Ducas (CWI Amsterdam) — 0 factoring relations in 1000 trials with Schnorr's claimed parameters.

## Current implementations

### MPQS (Multiple Polynomial Quadratic Sieve) — `library/mpqs.c`
**Status**: Working, well-optimized implementation from another agent.

Built with:
- Factor base of primes where N is a QR, scaled by digit count
- Polynomial A = product of k FB primes (k chosen so A ≈ sqrt(2N)/M)
- B via CRT from sqrt(N) mod each A-prime, C = (B²-N)/A
- Sieve with log-approximation threshold
- Large prime variation with hash table for combining partials
- Gaussian elimination mod 2

**Performance** (single core, worst case across 5 semiprimes per size):
- 30 digits: ~0.02s
- 52 digits: ~1.2s
- 55 digits: ~18s
- 60 digits: ~29s
- 65 digits: ~87s

### SIQS v2 (Self-Initializing QS) — `library/mpqs.cpp`
**Status**: Working but less optimized. Uses random A-value selection with CRT-based b computation. Works well up to ~50 digits.

### ECM (Elliptic Curve Method) — `library/ecm_factor.c`
**Status**: Working wrapper around GMP-ECM.

Uses Suyama parameterization with sigma values < 2^32, seeded from deterministic RNG with seed=42. Performance on balanced semiprimes:
- 30 digits: 0.077s
- 40 digits: 0.5s
- 50 digits: 4.4s

ECM scales with the size of the _smallest_ factor, not N. For balanced semiprimes, factors are ~N^(1/2), so ECM has L[1/2] complexity in terms of N. Competitive with QS up to ~50 digits on balanced semiprimes.

### CFRAC (Continued Fraction Algorithm) — `library/cfrac.cpp`
**Status**: Working.

Uses continued fraction expansion of sqrt(N) to generate smooth numbers bounded by 2*sqrt(N). Single large prime variation. L[1/2] complexity like QS but with a single "polynomial" (the CF expansion), limiting relation generation rate.

Performance: 30-digit ~0.14s, 40-digit ~23s. Much slower than SIQS for the same digit count — the lack of multiple polynomials is the bottleneck.

### Smooth Subsum Search (SSS) — `library/sss.cpp`
**Status**: Working prototype, needs fix for >35 digits (long long overflow for x values).

Implementation of Hittmeir (2023) algorithm. Key idea:
- Same polynomial as QS: pol(j) = (j+b)^2 - N
- But instead of sieving, uses CRT to construct j values where pol(j) is **guaranteed divisible** by several factor base primes
- The smaller cofactor is then tested for smoothness using batch methods (product/remainder trees)
- This replaces sieving with structured candidate generation

Performance (30-35 digits): comparable to our SIQS. Paper claims 5-7x speedup over comparable QS implementations. Our C++ implementation needs GMP-based x-value tracking for larger sizes.

**Known bug**: x values stored as long long overflow for N > ~35 digits. Need to use mpz_t throughout.

### Pollard's Rho — `library/pollard_rho.c`
**Status**: Working. Brent's improvement with batch GCD. Good for small factors but O(N^(1/4)) complexity makes it only useful for ~30-digit factors.

## Research survey: what approaches could beat L[1/2]?

Based on literature review (2020-2026):

### No classical sub-L[1/3] algorithm exists
GNFS at L[1/3, (64/9)^(1/3)] remains the asymptotic state of the art. No classical algorithm has achieved sub-L[1/3] complexity for general integers.

### Regev's quantum factoring (2023)
O~(n^{3/2}) quantum gates (vs Shor's O~(n^2)). Ragavan & Vaikuntanathan reduced space to O(n log n) qubits. Pilatte (2024) proved unconditional correctness. **Fundamentally quantum — no classical dequantization known.**

### Hittmeir's Smooth Subsum Search (2023) — IMPLEMENTED
L[1/2] with better constants than QS. Demonstrated 5-7x speedup over comparable SIQS. Uses CRT-based structured candidate generation. **Our implementation works for small sizes; needs overflow fix for larger sizes.**

### Umans & Wang conjecture (2025) — SPECULATIVE
If a number-theoretic conjecture holds: deterministic N^{1/6+o(1)} factoring. But this is still exponential, far slower than QS/NFS.

### Harvey-Hittmeir deterministic bound
N^{1/5+o(1)} deterministic factoring — best known deterministic bound. Still exponential.

## Failed approaches / dead ends

### Dixon's random squares with batch smoothness
**Tested**: Random x^2 mod N values are ~N in size, far too large for smoothness testing. The probability of a random number near N being B-smooth is astronomically low. QS works because it generates numbers of size ~sqrt(N) via polynomial structure. **This approach cannot work without the polynomial structure.**

### Basic QS (single polynomial)
**Tested**: Q(x) = (sqrt(N)+x)^2 - N grows as O(sqrt(N)*x). At the edge of the sieve interval, values are too large for smoothness. The fix is MPQS/SIQS with multiple polynomials keeping |Q(x)| bounded by sqrt(N/2).

## Key technical lessons

1. **SIQS polynomial generation**: The 'a' value must be ≈ sqrt(2N)/M (a product of multiple FB primes), NOT a single prime. Using single primes gives c = (b²-N)/a ≈ N/a which is huge, defeating the purpose.

2. **SIQS congruence**: The identity is (ax+b)² ≡ a·Q(x) (mod N). The exponent vector must include the factor 'a' — missing this gives near-zero rank in GF(2) matrix.

3. **CFRAC identity**: After step k of the CF expansion, the identity uses p_{k-1} (the PREVIOUS convergent) not p_k, and the sign is (-1)^k.

4. **SSS insight**: By using CRT to guarantee divisibility by several primes, the cofactor to test is much smaller, increasing smoothness probability. This is fundamentally different from sieving and could have better practical performance.

## What to try next

1. **Fix SSS for larger sizes** (mpz_t for x values) and benchmark against MPQS up to 70+ digits
2. **Implement a basic GNFS** for L[1/3] scaling on 70-100 digit numbers — even an unoptimized implementation would show better scaling than L[1/2]
3. **Hybrid ECM+QS**: Use ECM as a fast first-pass filter, then QS for remaining composites
4. **Lattice-enhanced smooth finding**: Use LLL on small lattices to find pairs of numbers whose product is more likely smooth
5. **Multi-large-prime variation**: 2LP and 3LP can dramatically increase relation yield
