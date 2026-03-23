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

**Conclusion**: The lattice dimension is tied to the smoothness bound, which grows sub-exponentially with N. So this doesn't escape the smoothness bottleneck — it just reformulates it as a lattice problem that's equally hard. Any lattice-based approach that depends on a large smooth factor base will have the same problem.

## Current implementations

### SIQS (Self-Initializing Quadratic Sieve)
**Status**: Working implementation in `library/siqs.cpp` (from another agent).

Uses self-initializing polynomials (CRT-based b computation), large prime variation, GF(2) Gaussian elimination, factor base sizing and sieve parameters tuned per digit count. Expected L[1/2] scaling.

### MPQS (Multiple Polynomial Quadratic Sieve)
**Status**: Working implementation in `library/mpqs.c`.

Built a working MPQS with:
- Factor base of primes where N is a QR, scaled by digit count
- Polynomial A = product of k FB primes (k chosen so A ≈ sqrt(2N)/M)
- B via CRT from sqrt(N) mod each A-prime, C = (B²-N)/A
- Sieve with log-approximation threshold
- Large prime variation with hash table for combining partials
- Gaussian elimination mod 2 with identity augmentation for null-space tracking
- Square root computation including large prime factors from combined partials

**Key implementation detail**: The congruence is (Ax+B)² ≡ A·f(x) (mod N) where f(x) = Ax²+2Bx+C. The exponent vector must include A's prime factors (each with +1 exponent) in addition to the trial division of f(x). Missing this causes near-zero rank in the GF(2) matrix.

**Performance** (single core, worst case across 5 semiprimes per size):
- 30 digits: ~0.02s
- 52 digits: ~1.2s
- 60 digits: ~29s
- 65 digits: ~87s

This is slower than YAFU's optimized SIQS by a factor of ~100x (YAFU does 60 digits in 0.7s). The main bottlenecks:
1. No self-initialization (recompute sieve offsets from scratch for each poly)
2. No block sieving (poor cache utilization)
3. Single large prime only (no 2LP/3LP)
4. Polynomial selection not optimized (random A, should use Gray code enumeration)
5. Gaussian elimination is dense O(n³) instead of structured sparse

### ECM (Elliptic Curve Method)
**Status**: Working wrapper around GMP-ECM in `library/ecm_factor.c`.

Uses Suyama parameterization with sigma values < 2^32, seeded from deterministic RNG with seed=42. Good for finding small factors but fundamentally limited by factor size (not N size). For balanced semiprimes, competitive with QS up to ~50 digits.

## Novel approach ideas to explore

### 1. Spectral/Fourier methods on multiplicative structure
Shor's algorithm exploits periodicity in x → a^x mod N using QFT. Classically, we can't do QFT on exponentially large groups, but:
- Can we detect partial periodicity using classical DFT on _samples_ of the multiplicative group?
- Can lattice reduction on the "frequency domain" of modular exponentiation reveal period information?

### 2. Number field sieve (NFS) implementation
NFS achieves L[1/3], which would be a significant improvement over SIQS L[1/2] for numbers >70 digits. Worth implementing even though it's not novel — it gives us a stronger baseline. Requires polynomial selection, 2D sieve, and square root in number field.

### 3. Smooth number amplification
Idea: instead of testing random numbers for smoothness, use algebraic identities to construct numbers that are "partially smooth" by design, then only need to test a smaller cofactor for smoothness.

### 4. Batch smoothness testing (Bernstein's method)
Instead of sieving, generate many candidates and use product trees + remainder trees to batch-test smoothness. Asymptotically better I/O complexity than sieving for very large smoothness bounds.

### 5. Multi-image approach for L[1/4]
**Speculative**: if NFS gets L[1/3] from two images (rational + algebraic), could THREE images give L[1/4]? This would require finding a suitable third algebraic structure. Related to "tower NFS" ideas but for general composites.

### 6. Lattice-based relation finding
Use LLL/BKZ not on the factor base lattice (Schnorr's approach, which fails), but on a smaller lattice to find _pairs_ of numbers whose product is smooth. The key difference: we only need the lattice to find structure in O(1) relations at a time, not encode the entire factor base.

### 7. p-adic / Hensel lifting approaches
Use the p-adic structure of Z/NZ to iteratively lift partial factorizations. If N = pq, then working mod small powers of primes might reveal information about p and q.

### 8. Spectral/algebraic approaches
Exploit the group structure of Z_N* ≅ Z_{p-1} × Z_{q-1} without smooth numbers. Possibly via random walks, exponential sums, or character sums that distinguish the product structure from a cyclic group.
