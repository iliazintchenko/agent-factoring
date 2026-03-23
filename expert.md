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
- **Special-Q QS**: Tested more thoroughly — fails to collect enough relations for 40d+ because trial division can't replace sieving. The sparse candidate positions (only M/q per special prime) don't compensate for the lack of sieve amortization. Not viable without a full sieve within each special-Q class.
- **DLP-SIQS (Double Large Prime)**: Implemented with Pollard rho cofactor splitting and union-find DLP graph. 5-100x slower than YAFU (30d: 0.07s, 50d: 12.7s, 55d: 64.5s). DLP matching produces few combined relations — the bottleneck is the scalar sieve, not the LP strategy.

### Approaches with some promise
- **SPQS (Smooth Polynomial QS / Multi-Polynomial Batch Sieve)**: Novel approach that sieves MULTIPLE SIQS polynomials simultaneously over the same sieve block. For each 'a' value, generates BATCH_POLYS (=4) b-values and processes them in a single pass over the factor base. This amortizes the outer loop over FB primes and produces ~4x more smooth candidates per sieve iteration. Results: 1.4x slower than YAFU at 30d, 8x at 50d, 44x at 60d. **Best custom implementation.** See scaling data below.
- **CADO-NFS**: Successfully built from source. NFS implementation with L[1/3] scaling. Factored 60d in ~28s wallclock (139s CPU multi-threaded). Running benchmarks for 70-90d. The key question is whether NFS shows better scaling than QS for 80-100d numbers despite higher overhead.

## Key Algorithmic Insights

1. **Sieving is the bottleneck**: For SIQS at 85-89 digits, sieving consumes 98% of runtime. Block Lanczos (linear algebra) is only 2-5%. Any novel approach must improve the sieve or bypass it entirely.
2. **Smoothness probability is the fundamental limit**: The probability that a random B-smooth number of size Q is smooth is roughly u^(-u) where u = log(Q)/log(B). Larger numbers have exponentially fewer smooth values — this is why factoring is sub-exponential, not polynomial.
3. **Polynomial quality matters**: Better polynomials give smaller Q(x) values, increasing smoothness probability. SIQS uses self-initialization to generate many good polynomials efficiently.
4. **Large primes extend the factor base cheaply**: Single/double/triple large primes (SLP/DLP/TLP) allow partially-smooth relations to contribute. DLP roughly doubles relation yield. TLP has high overhead and only helps above ~100 digits.
5. **The GF(2) linear algebra step**: Finding null space vectors over GF(2) is fast (Block Lanczos/Wiedemann). The hard part is generating enough relations, not solving the system.
6. **No known algorithm exploits balanced structure**: SIQS is factor-structure-agnostic. No published approach specifically targets N = p*q where p ≈ q ≈ √N with better complexity than general-purpose factoring.
7. **Multi-polynomial batching helps**: SPQS shows that processing multiple polynomials per sieve block amortizes the FB prime iteration overhead. At 4 polynomials per batch, smooth candidate yield increases ~3-4x per sieve pass with only ~15% memory overhead. The improvement narrows as FB size grows (larger sizes have more primes, and the inner loop dominates over the outer loop).
8. **DLP graph matching is hard**: Even with aggressive LP bounds (300x FB max), the DLP graph has too few edges for cycle formation at moderate sizes. The birthday paradox works against us — need ~sqrt(LP space) DLP relations for matches, but LP space grows faster than relation yield.
9. **Trial division after sieve is the main custom implementation bottleneck**: YAFU uses sieve-informed trial division (only tests primes whose sieve roots match x), plus bucket sieving for large primes. Custom implementations doing naive trial division lose 10-50x here.

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

SPQS Batch Sieve (worst of 5 per size):

| Digits | Time | Growth | vs YAFU |
|--------|------|--------|---------|
| 30 | 0.02s | | 1.4x |
| 35 | 0.042s | ~2x/5d | |
| 40 | 0.148s | 3.5x/5d | 8.7x |
| 45 | 0.48s | 3.2x/5d | |
| 50 | 0.998s | 2.1x/5d | 8.3x |
| 55 | 5.18s | 5.2x/5d | |
| 60 | 31.1s | 6.0x/5d | 44x |
| 65 | running | | |

DLP-SIQS (worst of 5 per size):

| Digits | Time | vs YAFU |
|--------|------|---------|
| 30 | 0.07s | 5x |
| 35 | 0.31s | |
| 40 | 1.27s | 75x |
| 45 | 6.65s | |
| 50 | 12.7s | 106x |
| 55 | 64.5s | |

Each digit adds ~15-20% to sieve time, consistent with L[1/2, 1+o(1)] scaling.

## Custom Implementations in library/

### Working implementations (best to worst)
- **spqs.c**: **Best custom.** Multi-polynomial batch sieve SIQS. 1.4-44x slower than YAFU. 30-65d+. `gcc -O3 -march=native -o spqs library/spqs.c -lgmp -lm`
- **dlp_siqs.c**: SIQS with DLP (Pollard rho splitting), union-find graph. 5-100x slower than YAFU. 30-55d. `gcc -O3 -march=native -mavx512bw -o dlp_siqs library/dlp_siqs.c -lgmp -lm`
- **siqs2.c**: Working SIQS, SLP, Gray code. 30-80x slower than YAFU. `gcc -O2 -march=native -o siqs2 library/siqs2.c -lgmp -lm`
- **siqs3.c**: SIQS with DLP, inline Block Lanczos. 40-350x slower. `gcc -O3 -march=native -mavx512bw -o siqs3 library/siqs3.c -lgmp -lm`
- **lattice_siqs.c**: Clean SIQS for scaling measurement. `gcc -O3 -march=native -o lattice_siqs library/lattice_siqs.c -lgmp -lm`

### Experimental / Novel
- **specialq_qs.c**: Special-Q QS. NEGATIVE RESULT: can't collect enough relations without sieving.
- **nfs_siever.c**: Custom NFS lattice siever. Working but 2500x slower than GGNFS.
- **gnfs_simple.c**: Line sieve NFS. Working but slow.
- **batch_smooth.c**, **batch_qs.c**, **batch_siqs.c**: Batch smoothness approaches. NEGATIVE RESULT.
- **lattice_factor.c**, **lattice_factor_v2.c**, **lattice_factor_batch.c**: Lattice-based approaches. NEGATIVE RESULT.

### Other
- **mpqs.c**: MPQS with buggy sqrt step.
- **pollard_rho.c**: Brent variant. Not competitive above 30d.
- **factor_oracle.c**: Multi-strategy oracle. Useful for small factors.
- **special_factor.c**: Pollard p-1, Williams p+1, ECM.

## Open Questions

- **Can NFS beat QS at 85-100 digits on single core?** CADO-NFS is now built. Need to benchmark with single-thread constraint. NFS overhead (poly selection, filtering, LA, square root) may exceed QS sieve savings at these sizes.
- **Can the batch polynomial approach be pushed further?** SPQS uses 4 polynomials per batch. More polynomials = more memory but also more candidates. The optimal batch size depends on L1 cache pressure vs relation yield.
- **Can SIMD improve the sieve inner loop?** The key bottleneck is scattered byte stores to random positions. AVX512 vscatterdq exists but is slow on Zen4. Alternative: use vectorized sieve scanning (already done in some implementations) but the sieve update loop itself remains scalar.
- **Can we bypass sieving entirely?** All QS variants sieve. The theoretical alternative is batch smoothness testing (product trees), but it's proven slower in practice. What about number-theoretic approaches (class groups, ideal theory)?
- **Can TLP (triple large primes) actually help at 80-90 digits?** The literature says overhead dominates below 100d. But with aggressive ECM cofactorization and efficient hypergraph cycle finding, the crossover point might be lower.
