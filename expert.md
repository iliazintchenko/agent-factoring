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
8. **DLP graph matching is hard**: Even with aggressive LP bounds (300x FB max), the DLP graph has too few edges for cycle formation at moderate sizes. **Update**: DLP with SQUFOF cofactorization tested through 65d — zero composite DLP cofactors found. After trial division by FB, residues > LP_bound are almost always prime at 30-65d. DLP only becomes useful at 80+ digits.
9. **Trial division after sieve is the main custom implementation bottleneck**: YAFU uses sieve-informed trial division (only tests primes whose sieve roots match x), plus bucket sieving for large primes. Custom implementations doing naive trial division lose 10-50x here.
10. **Bucket sieve is critical for large factor bases**: Without bucket sieve, sieve hits for primes > BLOCK_SIZE cause random memory access (cache misses). Bucket sieve pre-sorts hits by block, converting to sequential access. Key implementation detail: buckets MUST create new slices when they fill up (>75% of BUCKET_ALLOC). Without this, entries are silently dropped, causing 2.3x fewer relations per polynomial.
11. **Parameter tuning: sieve interval size is critical**: YAFU uses tiny sieve intervals (1-8 blocks per side) because its AVX512 sieve is fast enough that polynomial overhead is small. For scalar implementations, larger intervals (10-50 blocks) are better because they amortize the per-polynomial overhead (mpz operations for b,c computation, ainv computation, etc.).
12. **DLP with LP columns vs matching**: SLP (single large prime) is best handled by hash-matching (combine two relations with same LP). DLP (double large prime) can use LP columns in the GF(2) matrix instead of graph cycle detection. The LP-column approach is simpler but requires more relations. For SLP-only, matching is always better (no extra matrix columns needed).
13. **Adaptive sieve threshold improves SLP matching**: Lowering the threshold by ~6 bits for larger sizes dramatically increases SLP partials, improving matching rate by ~33% at 55d. The tradeoff is more TD on false positives, but SLP gains compensate. Best scaling improvement found so far.
14. **Negative results**: MCFRAC (multi-CF, 70-300x slower, sequential), PairQS (paired smooth, LESS smooth than individual), batch=8/16 (no improvement, inner loop dominates), bucket sieve without batch poly (slower than SPQS).

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

SIQS-Opt Bucket Sieve (worst of 5 per size):

| Digits | Time | Growth | vs YAFU |
|--------|------|--------|---------|
| 30 | 0.028s | | 2.0x |
| 35 | 0.050s | 1.8x/5d | |
| 40 | 0.179s | 3.6x/5d | 10.5x |
| 45 | 0.811s | 4.5x/5d | ~16x |
| 50 | 2.227s | 2.7x/5d | 18.6x |
| 55 | 4.844s | 2.2x/5d | ~12x |
| 60 | 19.203s | 4.0x/5d | 27x |
| 65 | 106.5s | 5.5x/5d | ~43x |
| 70 | 265.2s | 2.5x/5d | 46x |
| 75 | >295s | | |

SPQS-DLP Adaptive Threshold (worst of 5 per size):

| Digits | Time | Growth | vs YAFU | vs SPQS |
|--------|------|--------|---------|---------|
| 30 | 0.024s | | 1.7x | 1.2x |
| 40 | 0.163s | ~6.8x/10d | 9.6x | 1.1x |
| 50 | 1.06s | ~6.5x/10d | 8.8x | 1.1x |
| 55 | 3.49s | 3.3x/5d | - | **0.67x** |
| 60 | 29.6s | 8.5x/5d | 42x | **0.95x** |
| 61 | 40.9s | 1.4x/1d | - | - |

Key finding: adaptive threshold gives best improvement at 55d (33% faster than SPQS).
At 60d, still competitive. DLP itself contributes zero relations below 65d.

HyperSIQS with DLP+TLP (worst of 5 per size):

| Digits | Time | Growth | vs YAFU |
|--------|------|--------|---------|
| 30 | 0.071s | | 5x |
| 35 | 0.066s | ~1x/5d | |
| 40 | 0.166s | 2.5x/5d | 10x |
| 45 | 0.772s | 4.6x/5d | |
| 50 | 1.578s | 2.0x/5d | 13x |
| 55 | 5.688s | 3.6x/5d | |
| 60 | 20.72s | 3.6x/5d | 30x |
| 65 | 70.23s | 3.4x/5d | |

Each digit adds ~15-20% to sieve time, consistent with L[1/2, 1+o(1)] scaling.

## Custom Implementations in library/

### Working implementations (best to worst)
- **spqs_dlp.c**: SPQS + DLP (SQUFOF) + adaptive threshold. **Best at 55d** (3.5s vs 5.2s SPQS). DLP doesn't contribute below 65d but adaptive threshold improves SLP matching. `gcc -O3 -march=native -o spqs_dlp library/spqs_dlp.c -lgmp -lm`
- **spqs2.c / spqs.c**: Best custom at 30-65d. SPQS2 adds bucket sieve to batch approach. 1.4-44x slower than YAFU. `gcc -O3 -march=native -o spqs library/spqs.c -lgmp -lm`
- **fast_siqs.c**: Bucket sieve SIQS with __int128 trial division, sieve-informed TD, Gray code self-init, SLP. 15-25x slower than YAFU. 30-65d. `gcc -O3 -march=native -o fast_siqs library/fast_siqs.c -lgmp -lm`
- **siqs_opt.c**: Bucket sieve SIQS, Gray code poly, SLP matching, tracked offsets and 64-bit scanning. 2-46x slower. 30-70d. `gcc -O3 -march=native -o siqs_opt library/siqs_opt.c -lgmp -lm`
- **dlp_opt.c**: DLP SIQS with LP columns in GF(2) matrix. Hybrid SLP matching + DLP. Testing at 70d+.
- **siqs_bucket.c**: Similar optimized SIQS with bucket sieve. Comparable to fast_siqs. `gcc -O3 -march=native -o siqs_bucket library/siqs_bucket.c -lgmp -lm`
- **sqqs.c**: Special-Q QS (novel, uses special primes to reduce Q(x) size). NEGATIVE RESULT: overhead exceeds benefit.
- **dlp_siqs.c**: SIQS with DLP (Pollard rho splitting). 5-100x slower. 30-55d.
- **siqs2.c**: Working SIQS, SLP, Gray code. 30-80x slower.
- **siqs3.c**: SIQS with DLP, inline Block Lanczos. 40-350x slower.

### Experimental / Novel
- **hyper_siqs.c**: SIQS with Triple Large Primes (TLP), bucket sieve, Contini Gray code with incremental root updates, 64-bit Pollard rho for DLP/TLP cofactor splitting. Novel TLP variation accepts cofactors with up to 3 large primes. Results: 5x slower than YAFU at 30d, 30x at 60d. Best at 65d among custom implementations (70.2s vs spqs2's 94.1s). TLP not yet contributing useful relations at these sizes—needs larger LP bounds.
- **ecm_siqs.c**: SIQS with ECM cofactorization. DLP matching was wrong; LP-column approach in dlp_opt.c is correct.
- **latsieve_qs.c**: SIQS with ECM cofactor splitting. NEEDS more work.
- **nfs_factor.c**: NFS skeleton with base-m poly selection. Algebraic sqrt NOT implemented.
- **specialq_qs.c**: Special-Q QS. NEGATIVE RESULT: can't collect relations without sieving.
- **nfs_siever.c**: Custom NFS lattice siever. Working but 2500x slower than GGNFS.
- **batch_smooth.c**, **batch_qs.c**, **batch_siqs.c**: Batch smoothness approaches. NEGATIVE RESULT.
- **lattice_factor.c**, **lattice_factor_v2.c**, **lattice_factor_batch.c**: Lattice-based approaches. NEGATIVE RESULT.

### Other
- **pollard_rho.c**: Brent variant. Not competitive above 30d.
- **factor_oracle.c**: Multi-strategy oracle. Useful for small factors.
- **special_factor.c**: Pollard p-1, Williams p+1, ECM.

## Key Bottleneck Analysis

The gap between custom SIQS (15-25x slower) and YAFU comes from:
1. **Scalar sieve inner loop** (~60% of gap): YAFU uses AVX512BW vectorized sieve; custom uses scalar byte stores. No portable C workaround exists.
2. **Trial division efficiency** (~20%): YAFU uses sieve-root-informed TD with multiplication-by-inverse; custom uses mpz_divisible_ui_p fallback for some primes.
3. **Polynomial overhead** (~10%): YAFU's self-init is highly optimized with incremental root updates; custom recomputes roots from scratch.
4. **Parameter tuning** (~10%): YAFU has decades of tuning; custom parameters are close but not optimal.

The 15-25x gap is NOT closeable with pure C without SIMD. The sieve inner loop (scattered byte stores at stride p) is inherently memory-bound and benefits enormously from wide SIMD gather/scatter.

## Open Questions

- **Can NFS beat QS at 85-100 digits on single core?** NFS has L[1/3] vs QS's L[1/2] scaling. But NFS implementation complexity is enormous (polynomial selection, lattice sieve, filtering, algebraic square root). The crossover depends on constant factors.
- **Can ECM cofactorization improve SIQS?** ECM can split cofactors into DLP relations, but overhead (~10μs per ECM call) is high vs sieve amortization. Needs DLP graph (union-find) and proper exponent tracking. Potentially useful above 70d where sieve cost dominates.
- **Can a simplified NFS be implemented?** The main barrier is the algebraic square root (Couveignes' algorithm). Line sieving + GF(2) LA + base-m polynomial selection are all implementable. CRT-based algebraic sqrt requires: root finding mod large primes (Cantor-Zassenhaus), square root mod each prime (Tonelli-Shanks), polynomial interpolation, CRT accumulation.
- **Can we exploit balanced semiprime structure?** No known algorithm specifically targets N = p*q with p ≈ q. The Fermat/Lehman approaches only help when |p-q| is small relative to N^(1/3).
- **Can TLP (triple large primes) help at 80-90 digits?** Literature says overhead dominates below 100d. With ECM cofactorization and hypergraph cycle finding, crossover might be lower.
