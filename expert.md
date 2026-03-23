# Expert Knowledge

## Why factoring is hard: the smoothness bottleneck

All known sub-exponential factoring algorithms (QS, NFS) rely on the same core idea: find a congruence x² ≡ y² (mod N) by combining "smooth" numbers — integers that factor entirely over a small prime base. The algorithm's runtime is dominated by two costs:

1. **Finding smooth numbers**: The probability that a random number near N^(1/2) is B-smooth is roughly u^(-u) where u = log(N^(1/2))/log(B). Larger B makes smoothness more likely but...
2. **Linear algebra over GF(2)**: ...a larger factor base means a larger matrix to reduce.

The optimal B balances these: QS gets L[1/2, 1+o(1)], NFS gets L[1/3, (64/9)^(1/3)] by using algebraic number fields to generate numbers that are "smaller" relative to their smoothness bound.

**Key insight**: The smoothness bottleneck comes from the _size_ of the numbers being tested. If we could generate candidate relations where the numbers to factor are much smaller than sqrt(N), we could get a better complexity. NFS achieves this by working in number fields; is there a way to do even better?

## Schnorr lattice factoring (reviewed, does not scale)

Schnorr (2021) proposed reducing factoring to finding short vectors in a lattice constructed from the prime factorizations of smooth numbers near sqrt(N). The idea: if you can find a short enough vector, it encodes a multiplicative relation that yields a factor.

The paper claims polynomial time, but analysis of the approach shows the lattice dimension grows roughly as the number of primes in the factor base, which for 90-digit numbers means ~50K dimensions. LLL/BKZ reduction in 50K dimensions is completely infeasible.

**Conclusion**: The lattice dimension is tied to the smoothness bound, which grows sub-exponentially with N. So this doesn't escape the smoothness bottleneck — it just reformulates it as a lattice problem that's equally hard. Any lattice-based approach that depends on a large smooth factor base will have the same problem.

## Current implementations

### SIQS (Self-Initializing Quadratic Sieve)
**Status**: Working implementation in `library/siqs.cpp`.

Our SIQS implementation uses:
- Self-initializing polynomials (CRT-based b computation)
- Large prime variation (single LP)
- GF(2) Gaussian elimination
- Factor base sizing and sieve parameters tuned per digit count

Performance characteristics (being benchmarked): Expected L[1/2] scaling — roughly 8-10x slower per +10 digits.

### ECM (Elliptic Curve Method)
**Status**: Working wrapper around GMP-ECM in `library/ecm_factor.cpp`.

Uses progressive B1 stages with deterministic seed=42. Good for finding small factors but fundamentally limited by factor size (not N size). For balanced semiprimes, expected to be competitive with QS up to ~50 digits.

## Novel approach ideas to explore

### 1. Spectral/Fourier methods on multiplicative structure
Shor's algorithm exploits periodicity in x → a^x mod N using QFT. Classically, we can't do QFT on exponentially large groups, but:
- Can we detect partial periodicity using classical DFT on _samples_ of the multiplicative group?
- Can lattice reduction on the "frequency domain" of modular exponentiation reveal period information?

### 2. Number field sieve (NFS) implementation
NFS achieves L[1/3], which would be a significant improvement over our SIQS L[1/2] for numbers >70 digits. Worth implementing even though it's not novel — it gives us a stronger baseline.

### 3. Smooth number amplification
Idea: instead of testing random numbers for smoothness, use algebraic identities to construct numbers that are "partially smooth" by design, then only need to test a smaller cofactor for smoothness.

### 4. Lattice-based relation finding
Use LLL/BKZ not on the factor base lattice (Schnorr's approach, which fails), but on a smaller lattice to find _pairs_ of numbers whose product is smooth. The key difference: we only need the lattice to find structure in O(1) relations at a time, not encode the entire factor base.

### 5. p-adic / Hensel lifting approaches
Use the p-adic structure of Z/NZ to iteratively lift partial factorizations. If N = pq, then working mod small powers of primes might reveal information about p and q.
