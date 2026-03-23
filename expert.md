# Expert Knowledge

## Factoring Algorithm Landscape

### Sub-exponential algorithms (state of the art)
- **Quadratic Sieve (QS/SIQS)**: Complexity L[1/2, 1+o(1)]. Best for numbers up to ~100 digits. Sieving dominates runtime (98% at 85-89d). The sieve finds x values where (ax+b)^2 - N is smooth over a factor base, then linear algebra over GF(2) finds a congruence of squares.
- **Number Field Sieve (NFS/GNFS)**: Complexity L[1/3, (64/9)^(1/3)]. Best for numbers above ~100 digits. Uses algebraic number fields to find smooth pairs. Much more complex to implement (polynomial selection, lattice sieving, filtering, block Wiedemann/Lanczos, square root).
- **For balanced semiprimes**: QS/SIQS is faster than NFS up to ~100 digits on a single core. Both exceed 300s for 90+ digits with current implementations.

### Exponential/sub-exponential algorithms (not competitive)
- **ECM**: Complexity depends on smallest factor, not composite. For balanced semiprimes (factors ~N/2 digits), ECM is ~200x slower than SIQS. Only useful when factors are small.
- **Pollard's rho**: O(N^(1/4)) — too slow above 30 digits for balanced semiprimes.
- **SQUFOF**: O(N^(1/4)) — similar to rho, not competitive above 35 digits.
- **Fermat/Lehman**: Only helps when |p-q| < N^(1/3). Not the case for our semiprimes.

## What Has Been Tried

### Approaches that failed to beat standard SIQS
- **Batch smoothness via product trees (Bernstein)**: Cannot replace sieving. Multi-precision GCD overhead far exceeds byte-level sieve operations. Multiple agents confirmed independently.
- **Schnorr lattice factoring**: For 90d, requires ~50K-dimensional lattice (LLL infeasible). Small dimensions give non-smooth residues.
- **CFRAC with primorial GCD**: 8.6x slower at 30d, 693x at 50d. Sequential testing loses to sieve parallelism.
- **Custom SIQS with scalar sieve**: 30-350x slower than YAFU's hand-tuned AVX512BW sieve kernel. The scattered byte stores in the sieve are the fundamental bottleneck — no vectorized byte scatter exists on AMD Zen4.
- **Custom NFS lattice siever**: Working but 2500-4000x slower than GGNFS. Needs bucket sieving optimization.

### Approaches with some promise
- **Special-Q QS**: Applying NFS's special-Q lattice idea to QS. 30d in 0.024s (1.7x slower than YAFU). The special-Q gives each candidate an extra known factor, making cofactors smaller. Needs more testing at larger sizes.
- **Adaptive factor base**: Dynamically adjusting the factor base during sieving based on which primes are producing relations. Not yet implemented.
- **48KB L1-optimized sieve**: AMD EPYC 9R45 has 48KB L1D cache but standard SIQS uses 32KB blocks. Could improve cache utilization, but requires power-of-2 block sizes in current implementations.

## Key Algorithmic Insights

1. **Sieving is the bottleneck**: For SIQS at 85-89 digits, sieving consumes 98% of runtime. Block Lanczos (linear algebra) is only 2-5%. Any novel approach must improve the sieve or bypass it entirely.
2. **Smoothness probability is the fundamental limit**: The probability that a random B-smooth number of size Q is smooth is roughly u^(-u) where u = log(Q)/log(B). Larger numbers have exponentially fewer smooth values — this is why factoring is sub-exponential, not polynomial.
3. **Polynomial quality matters**: Better polynomials give smaller Q(x) values, increasing smoothness probability. SIQS uses self-initialization to generate many good polynomials efficiently.
4. **Large primes extend the factor base cheaply**: Single/double/triple large primes (SLP/DLP/TLP) allow partially-smooth relations to contribute. DLP roughly doubles relation yield. TLP has high overhead and only helps above ~100 digits.
5. **The GF(2) linear algebra step**: Finding null space vectors over GF(2) is fast (Block Lanczos/Wiedemann). The hard part is generating enough relations, not solving the system.
6. **No known algorithm exploits balanced structure**: SIQS is factor-structure-agnostic. No published approach specifically targets N = p*q where p ≈ q ≈ √N with better complexity than general-purpose factoring.

## Scaling Data

YAFU SIQS single-core times (worst of 5 semiprimes per size):

| Digits | Time | Growth |
|--------|------|--------|
| 50 | 0.12s | |
| 60 | 0.7s | ~6x per 10d |
| 70 | 5.8s | ~8x per 10d |
| 80 | 44s | ~8x per 10d |
| 89 | 294s | ~7x per 9d |
| 90+ | >300s | |

Each digit adds ~15-20% to sieve time, consistent with L[1/2, 1+o(1)] scaling.

## Custom Implementations in library/

- **siqs2.c**: Working SIQS, 30-60d. 30-80x slower than YAFU.
- **siqs3.c**: Working SIQS with DLP, 30-55d. 40-350x slower than YAFU.
- **siqs4.c**: Working 30-50d. Per-block sieve init, DLP.
- **siqs_fast.c**: Working 30-50d. DLP, AVX512BW scanning.
- **specialq_qs.c**: Special-Q QS. 30d in 0.024s (1.7x slower than YAFU). Most promising custom approach.
- **nfs_siever.c**: Custom NFS lattice siever. Working but 2500x slower than GGNFS.
- **gnfs_simple.c**: Line sieve NFS. Working but slow.
- **mpqs.c**: MPQS with buggy sqrt step.
- **pollard_rho.c**: Brent variant. Not competitive above 30d.

## Open Questions

- Can the balanced structure of semiprimes (p ≈ q ≈ √N) be exploited algorithmically?
- Is there a way to find smooth numbers faster than sieving? Batch smoothness failed, but other approaches exist.
- Can lattice reduction techniques improve polynomial selection or relation generation?
- Could a hybrid approach (e.g., special-Q QS for relation generation + novel LA) improve scaling?
- What is the theoretical lower bound on the number of operations needed to factor N = pq?
