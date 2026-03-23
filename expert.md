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
- **SIQS-Bucket with Gray Code + DLP→SLP Pipeline** (siqs_bucket.c): Key optimizations:
  1. YAFU-calibrated parameters (reduced sieve interval from 14x too large)
  2. Gray code self-initialization (O(FB) additions per poly switch instead of O(FB) mod_inverse calls)
  3. DLP→SLP pipeline matching: DLP cofactors (p1,p2) checked against SLP hash; if either matches, creates new SLP partial with remaining LP
  4. Bucket sieve for large FB primes
  Results: 50d=0.86s (7x YAFU), 60d=8.9s (13x YAFU), 65d=70.7s (best custom by 25%). At 65d, DLP pipeline contributes 24% of all relations.
- **SPQS2 (SPQS + Bucket Sieve + AVX512 LA)**: Extends SPQS with bucket sieve for large primes (p > 32768), 4x small prime loop unrolling, native 64-bit trial division fast path, and AVX512 GF(2) linear algebra. **Best custom implementation.** 1.4x faster than SPQS at 60d (19s vs 31s). Extends to 70d (139s) and 75d (~250s). See scaling data below.
- **SPQS (Smooth Polynomial QS / Multi-Polynomial Batch Sieve)**: Sieves 4 SIQS polynomials simultaneously per sieve block. Amortizes FB prime iteration overhead. Superseded by SPQS2.
- **CADO-NFS**: Successfully built from source. NFS implementation with L[1/3] scaling. Factored 60d in ~28s wallclock (139s CPU multi-threaded). Running benchmarks for 70-90d. The key question is whether NFS shows better scaling than QS for 80-100d numbers despite higher overhead.

## Key Algorithmic Insights

1. **Sieving is the bottleneck**: For SIQS at 85-89 digits, sieving consumes 98% of runtime. Block Lanczos (linear algebra) is only 2-5%. Any novel approach must improve the sieve or bypass it entirely.
2. **Smoothness probability is the fundamental limit**: The probability that a random B-smooth number of size Q is smooth is roughly u^(-u) where u = log(Q)/log(B). Larger numbers have exponentially fewer smooth values — this is why factoring is sub-exponential, not polynomial.
3. **Polynomial quality matters**: Better polynomials give smaller Q(x) values, increasing smoothness probability. SIQS uses self-initialization to generate many good polynomials efficiently.
4. **Large primes extend the factor base cheaply**: Single/double/triple large primes (SLP/DLP/TLP) allow partially-smooth relations to contribute. DLP roughly doubles relation yield. TLP has high overhead and only helps above ~100 digits.
5. **The GF(2) linear algebra step**: Finding null space vectors over GF(2) is fast (Block Lanczos/Wiedemann). The hard part is generating enough relations, not solving the system.
6. **No known algorithm exploits balanced structure**: SIQS is factor-structure-agnostic. No published approach specifically targets N = p*q where p ≈ q ≈ √N with better complexity than general-purpose factoring.
7. **Multi-polynomial batching helps**: SPQS shows that processing multiple polynomials per sieve block amortizes the FB prime iteration overhead. At 4 polynomials per batch, smooth candidate yield increases ~3-4x per sieve pass with only ~15% memory overhead. The improvement narrows as FB size grows (larger sizes have more primes, and the inner loop dominates over the outer loop).
8. **DLP→SLP pipeline is the key DLP strategy**: Direct DLP-DLP matching (graph cycles) is useless below 80d (birthday paradox). Instead, DLP→SLP pipeline: for each DLP cofactor (p1,p2), check if p1 or p2 is in the SLP hash. If matched, combine DLP+SLP → new SLP partial with remaining LP. At 65d, this contributes 24% of all relations. Previous implementations used DLP graph matching which produced ~0 combined relations.
9. **Trial division after sieve is the main custom implementation bottleneck**: YAFU uses sieve-informed trial division (only tests primes whose sieve roots match x), plus bucket sieving for large primes. Custom implementations doing naive trial division lose 10-50x here.
10. **Bucket sieve is critical for large factor bases**: Without bucket sieve, sieve hits for primes > BLOCK_SIZE cause random memory access (cache misses). Bucket sieve pre-sorts hits by block, converting to sequential access. Key implementation detail: buckets MUST create new slices when they fill up (>75% of BUCKET_ALLOC). Without this, entries are silently dropped, causing 2.3x fewer relations per polynomial.
11. **Parameter tuning: sieve interval + Gray code**: YAFU uses tiny sieve intervals (1-3 blocks per side) because Gray code self-initialization makes polynomial switching O(FB) simple additions. With Gray code, even scalar implementations should use small M (2-4 blocks). Reducing M from 28 blocks to 2 blocks AND adding Gray code cut 60d time from 21s to 7.8s. However, without Gray code, larger intervals (10-50 blocks) are better to amortize mpz-based polynomial generation.
12. **DLP→SLP positive feedback loop at 65-72d**: The DLP→SLP pipeline creates a feedback loop: DLP splits generate new SLP partials, which increase the SLP hash, which enables more DLP matches. At 65d: 24% DLP contribution. At 70d: 15% steady-state. At 73d (borderline): 40% of rels from DLP pipeline. The 65→70d growth is only 1.5x/5d (vs 7.4x without DLP), demonstrating genuine scaling improvement. Key threshold: need >5K SLP partials for the feedback loop to activate.
13. **Adaptive sieve threshold improves SLP matching**: Lowering the threshold by ~6 bits for larger sizes dramatically increases SLP partials, improving matching rate by ~33% at 55d. The tradeoff is more TD on false positives, but SLP gains compensate. Best scaling improvement found so far.
14. **Negative results**: MCFRAC (multi-CF, 70-300x slower, sequential), PairQS (paired smooth, LESS smooth than individual), batch=8/16 (no improvement, inner loop dominates), bucket sieve without batch poly (slower than SPQS).
15. **Structured GE enables larger FB**: Singleton removal (iteratively remove weight-1 columns and their rows) reduces matrix 10-20% before dense Gaussian elimination. For 12K x 12K matrix, reduces to ~10K x 10K, cutting GE time by ~30%. Combined with the 48KB L1-cache-optimized sieve blocks, this is the turbo_siqs approach that first reached 72d.
16. **Sieve width has diminishing returns**: Experiments with nblocks=32, 56, 80 for 73d show the relation rate PER SECOND is roughly constant (28-32 rels/s) regardless of sieve width. This is because wider sieve means larger Q(x) at edges, exactly offsetting the extra positions. The optimal nblocks is where polynomial switching overhead is minimal.
17. **Structured GE with doubleton merging enables 73-74d**: Doubleton merging (XOR rows sharing a weight-2 column, remove both column and one row) reduces 20K matrices to ~15K. Combined with singleton removal, a 4-5 pass structured GE takes a 20K x 20K matrix to 14K x 14K. Dense GE on the reduced matrix takes ~25-30s. This enabled FB=18-20K which gives 3x more smooth values per polynomial, pushing the reliable factoring limit from 72d to 74d.
18. **74d is the reliable limit; 75d needs ~20% speedup**: At 74d, worst case is 290s (sieve 260s + LA 30s). At 75d, sieve needs ~340s, exceeding the 295s timeout. Batch polynomial sieving (~15% speedup) combined with a working Block Lanczos (~saves 25s LA) would bring 75d into range (~289+1 = 290s).

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

**Turbo SIQS** (worst of 5 per size) — **current best custom implementation**:

| Digits | Time | Growth | vs YAFU |
|--------|------|--------|---------|
| 30 | 0.04s | | 2.9x |
| 40 | 0.14s | | 8.2x |
| 50 | 0.97s | ~7x/10d | 8.1x |
| 55 | 3.79s | 3.9x/5d | |
| 60 | 26.1s | 6.9x/5d | 37x |
| 65 | 59.7s | 2.3x/5d | |
| 67 | 89.0s | | |
| 68 | 78.8s | | |
| 70 | 97.6s | 1.6x/5d | 16.8x |
| 71 | 172.6s | | |
| 72 | 247.8s | ~2.5x/2d | |
| 73 | 294.2s | | |
| 74 | 289.9s | | |

Key features: bucket sieving, Gray code self-init, SLP matching (LP_mult=150-200), structured Gaussian elimination with doubleton merging (SGE reduces 20K matrix to ~14K). 48KB L1-cache-optimized sieve blocks. FB=18-20K for 73-74d (enabled by SGE). First custom SIQS to reliably factor 74-digit balanced semiprimes within 300s single-core.

HyperSIQS2 with TLP + optimizations (worst of 5 per size):

| Digits | Time | Growth | vs YAFU |
|--------|------|--------|---------|
| 70 | 69.1s | | **12x** |
| 72 | 144.9s | ~2.1x/2d | |
| 73 | 139.9s | ~1x/1d | |
| 74 | 201s | 1.4x/1d | |
| 75 | 286s | 1.4x/1d | |

Key innovation: TLP (triple large primes) + incremental root offsets (no per-block mod div) + AVX512 GF(2) LA + interleaved two-root sieve. **First implementation to reliably factor all 5 semiprimes at 75 digits.** TLP gives fundamentally better relation yield per polynomial than SLP/DLP alone.

SIQS-Bucket with Gray Code + DLP→SLP Pipeline (worst of 5 per size):

| Digits | Time | Growth | vs YAFU |
|--------|------|--------|---------|
| 30 | 0.033s | | 2.4x |
| 40 | 0.174s | 5.3x/10d | 10.2x |
| 50 | 0.849s | 4.9x/10d | 7.1x |
| 55 | 4.97s | 5.9x/5d | |
| 60 | 9.86s | 2.0x/5d | 14.1x |
| 65 | 72.5s | 7.4x/5d | |
| 70 | 88.7s | 1.3x/5d | **15.3x** |
| 72 | 188.8s | ~2.1x/2d | |
| 73 | 295s (4/5) | ~1.6x/1d | |

**Key innovations**: (1) DLP→SLP pipeline: positive feedback loop at 65-73d. (2) TLP with primality pre-check. (3) 48KB L1 sieve blocks. (4) Gray code self-init. (5) Incremental sieve offsets. **Best custom at 70d (88.7s) and 72d (188.8s).**

SPQS2 Bucket Sieve (worst of 5 per size):

| Digits | Time | Growth | vs YAFU |
|--------|------|--------|---------|
| 30 | 0.022s | | 1.6x |
| 35 | 0.041s | ~2x/5d | |
| 40 | 0.139s | 3.4x/5d | 8.2x |
| 45 | 0.446s | 3.2x/5d | |
| 50 | 0.956s | 2.1x/5d | 8.0x |
| 55 | 4.65s | 4.9x/5d | |
| 60 | 19s | 4.1x/5d | 27x |
| 65 | 90s | 4.7x/5d | |
| 70 | 139s | 1.5x/5d | 24x |
| 75 | ~250s | ~1.8x/5d | |

Turbo SIQS2 (AVX512 LA + DLP→SLP pipeline, worst of 5 per size):

| Digits | Time | Growth | vs YAFU |
|--------|------|--------|---------|
| 30 | 0.058s | | 4.1x |
| 40 | 0.144s | ~2.5x/10d | 8.5x |
| 50 | 0.986s | ~6.8x/10d | 8.2x |
| 55 | 3.80s | 3.9x/5d | |
| 60 | 26.3s | 6.9x/5d | 38x |
| 65 | 60.0s | 2.3x/5d | |
| 70 | 97.7s | 1.6x/5d | **17x** |

Note: turbo_siqs2 at 70d (97.7s) has the best scaling 60→70d (3.7x/10d vs YAFU's 8.3x/10d).
The 48KB L1-cache sieve blocks + structured GE + DLP graph give better scaling at large sizes.

SPQS Batch Sieve (worst of 5 per size):

| Digits | Time | Growth | vs YAFU |
|--------|------|--------|---------|
| 30 | 0.02s | | 1.4x |
| 35 | 0.042s | ~2x/5d | |
| 40 | 0.148s | 3.5x/5d | 8.7x |
| 50 | 0.998s | | 8.3x |
| 55 | 5.18s | 5.2x/5d | |
| 60 | 31.1s | 6.0x/5d | 44x |

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

HyperSIQS with DLP+TLP + Structured GE (worst of 5 per size):

| Digits | Time | Growth | vs YAFU |
|--------|------|--------|---------|
| 30 | 0.074s | | 5x |
| 50 | 1.57s | | 13x |
| 60 | 18.6s | ~12x/10d | 27x |
| 65 | 64.3s | 3.5x/5d | |
| 70 | **74.9s** | **1.2x/5d** | **13x** |
| 72 | 164s | | |
| 73 | 168s | | |
| 75 | 262s | 3.5x/5d | |

**Best custom at 70d** (beats turbo_siqs 97.6s). Key insight: smaller FB (10-15K) with much larger LP multiplier (200-500x) and structured GE for fast LA. The growth rate drops dramatically at 65-70d due to aggressive SLP matching with larger LP bounds.

Each digit adds ~15-20% to sieve time, consistent with L[1/2, 1+o(1)] scaling.

## Custom Implementations in library/

### Working implementations (best to worst)
- **turbo_siqs.c**: **Best custom. Fastest at 70d+ (97.6s). First to reliably factor 72d (247.8s).** Bucket sieve, Gray code, SLP (LP_mult=200), structured GE (singleton removal + dense GE). 48KB L1-optimized blocks. 3-17x slower than YAFU. 30-72d. `gcc -O3 -march=native -o turbo_siqs library/turbo_siqs.c -lgmp -lm`
- **spqs2.c**: SPQS with bucket sieve + AVX512 LA. 1.6-24x slower than YAFU. 30-75d. `gcc -O3 -march=native -mavx512f -o spqs2 library/spqs2.c -lgmp -lm`
- **siqs_bucket.c**: Gray code + DLP→SLP pipeline + bucket sieve. 50d=0.86s (7x YAFU), 60d=8.9s (13x YAFU). `gcc -O3 -march=native -o siqs_bucket library/siqs_bucket.c -lgmp -lm`
- **hyper_siqs.c**: **Best at 60d** (single number). TLP SIQS with bucket sieve, Gray code, Pollard rho cofactor splitting. 9s on one 60d number. But at 70d (218s), LA dominates due to huge 22000x22000 matrix. `gcc -O3 -march=native -o hyper_siqs library/hyper_siqs.c -lgmp -lm`
- **hybrid_siqs.c**: SPQS multi-poly batch (4 polys) combined with bucket sieve for large primes. Slower at small sizes (overhead) but competitive at 65d (86.2s). `gcc -O3 -march=native -o hybrid_siqs library/hybrid_siqs.c -lgmp -lm`
- **spqs_dlp.c**: SPQS + DLP (SQUFOF) + adaptive threshold. **Best at 55d** (3.5s). `gcc -O3 -march=native -o spqs_dlp library/spqs_dlp.c -lgmp -lm`
- **fast_siqs.c**: Bucket sieve SIQS with __int128 TD, sieve-informed TD, Gray code, SLP. 18s at 60d. `gcc -O3 -march=native -o fast_siqs library/fast_siqs.c -lgmp -lm`
- **siqs_native.c**: Batch polynomials (4) + Gray code + SLP + native 128-bit TD. 19s at 60d. `gcc -O3 -march=native -o siqs_native library/siqs_native.c -lgmp -lm`
- **siqs_opt.c**: Bucket sieve SIQS, tracked offsets, 64-bit scanning. 2-46x slower. 30-70d. `gcc -O3 -march=native -o siqs_opt library/siqs_opt.c -lgmp -lm`
- **dlp_opt.c**: DLP SIQS with LP columns in GF(2) matrix. Testing at 70d+.
- **spqs.c**: Multi-polynomial batch sieve SIQS. 1.4-44x slower than YAFU. 30-65d. `gcc -O3 -march=native -o spqs library/spqs.c -lgmp -lm`
- **dlp_siqs.c**: SIQS with DLP (Pollard rho splitting), union-find graph. 5-100x slower than YAFU. 30-55d. `gcc -O3 -march=native -mavx512bw -o dlp_siqs library/dlp_siqs.c -lgmp -lm`
- **siqs2.c**: Working SIQS, SLP, Gray code. 30-80x slower than YAFU. `gcc -O2 -march=native -o siqs2 library/siqs2.c -lgmp -lm`
- **siqs3.c**: SIQS with DLP, inline Block Lanczos. 40-350x slower. `gcc -O3 -march=native -mavx512bw -o siqs3 library/siqs3.c -lgmp -lm`
- **lattice_siqs.c**: Clean SIQS for scaling measurement. `gcc -O3 -march=native -o lattice_siqs library/lattice_siqs.c -lgmp -lm`

### Experimental / Novel
- **hyper_siqs.c**: SIQS with Triple Large Primes (TLP), bucket sieve, Contini Gray code with incremental root updates, 64-bit Pollard rho for DLP/TLP cofactor splitting. Novel TLP variation accepts cofactors with up to 3 large primes. Results: 5x slower than YAFU at 30d, 30x at 60d. Best at 65d among custom implementations (70.2s vs spqs2's 94.1s). TLP not yet contributing useful relations at these sizes—needs larger LP bounds.
- **ecm_siqs.c**: SIQS with ECM cofactorization. DLP matching was wrong; LP-column approach in dlp_opt.c is correct.
- **latsieve_qs.c**: SIQS with ECM cofactor splitting. NEEDS more work.
- **nfs_factor.c**: NFS skeleton with base-m poly selection. Algebraic sqrt NOT implemented.
- **specialq_qs.c**: Special-Q QS. NEGATIVE RESULT: can't collect relations without sieving.
- **nfs_siever.c**: Custom NFS lattice siever. Working but 2500x slower than GGNFS.
- **batch_smooth.c**, **batch_qs.c**, **batch_siqs.c**: Batch smoothness approaches. NEGATIVE RESULT.
- **lattice_factor.c**, **lattice_factor_v2.c**, **lattice_factor_batch.c**: Lattice-based approaches. NEGATIVE RESULT.

### NFS (Novel)
- **nfs_complete.c**: Complete GNFS with degree-3 polynomial, line sieving, GF(2) LA with 100 QC columns, AND two algebraic sqrt approaches (Hensel lift + Couveignes CRT). Sieve works (8000 rels in 73s for 30d). LA works (64 deps). Both sqrt methods verified: T^2=S confirmed at full precision. **Remaining issue**: all dependencies give trivial gcd(X±Y,N) = N or 1. Both Hensel and CRT approaches produce the same sign. Root cause likely in QC column computation or the LA matrix setup, not the sqrt itself. Reference implementation: stubbscroll/nfs on GitHub. `gcc -O3 -march=native -o nfs_complete library/nfs_complete.c -lgmp -lm`
- **siqs_adaptive.c**: SIQS with Gray code self-init, YAFU-calibrated parameters with linear interpolation, bucket sieving. Reaches 65d (271s). Not competitive with best custom SIQS. `gcc -O3 -march=native -o siqs_adaptive library/siqs_adaptive.c -lgmp -lm`

### Other
- **pollard_rho.c**: Brent variant. Not competitive above 30d.
- **factor_oracle.c**: Multi-strategy oracle. Useful for small factors.
- **special_factor.c**: Pollard p-1, Williams p+1, ECM.

## Key Bottleneck Analysis

### Profiling Results (siqs_native, 60d)
- **Sieve inner loop: 73%** of total sieve time. The scattered byte stores `sieve[j] += logp` for each FB prime dominate.
- **Candidate scan + trial division: 25%**. Scanning sieve for candidates above threshold, then GMP trial division.
- **Polynomial setup: 2%**. Gray code solution updates are fast.

### Gap vs YAFU (15-30x slower)
1. **Scalar sieve inner loop** (~60% of gap): YAFU uses AVX512BW vectorized sieve; custom uses scalar byte stores. No portable C workaround exists.
2. **Trial division efficiency** (~20%): YAFU uses sieve-root-informed TD with multiplication-by-inverse; custom uses mpz_divisible_ui_p fallback for primes where Q(x) exceeds 128 bits. Native 128-bit TD helps for smaller numbers.
3. **Polynomial overhead** (~10%): YAFU's self-init is highly optimized with incremental root updates; custom recomputes roots from scratch. Gray code helps but isn't as fast as YAFU's approach.
4. **Parameter tuning** (~10%): YAFU has decades of tuning; custom parameters are close but not optimal.

### TLP (Triple Large Primes) Trade-off
- At 60d: TLP gives 2x speedup (9s vs 19s) by dramatically increasing relation yield per sieve pass
- At 70d: TLP HURTS (218s vs 165s) because the large FB (22000 primes) creates a 22000x22000 GF(2) matrix, and Gaussian elimination takes 164s (75% of total time)
- **Key insight**: TLP requires Block Lanczos or Block Wiedemann to handle the larger matrices efficiently. Without it, moderate FB + SLP/DLP is faster for 70d+

### SPV (Small Prime Variation) - NEGATIVE RESULT
- Skip sieving with primes < 256 and adjust threshold by expected contribution
- Saves ~40% of sieve work but introduces many false positives
- With native 128-bit TD, false positives are cheap to reject, but the imprecise threshold adjustment causes more misses too
- Net effect: approximately zero improvement or slight regression
- Would need a more sophisticated SPV implementation (e.g., sieve initialization with precomputed small prime contributions)

The 15-25x gap is NOT closeable with pure C without SIMD. The sieve inner loop (scattered byte stores at stride p) is inherently memory-bound and benefits enormously from wide SIMD gather/scatter.

## Open Questions

- **Can NFS beat QS at 85-100 digits on single core?** NFS has L[1/3] vs QS's L[1/2] scaling, but L-notation analysis shows NFS constants are WORSE than QS until ~100-130d. For 80d: QS L-complexity ≈ exp(31), NFS ≈ exp(33). NFS only wins above 100d in theory, and a simple NFS implementation would need 120-150d to beat optimized QS. Not viable for our 30-100d range. The existing nfs_factor.c skeleton has a broken polynomial selection (c_0 coefficient is ~N instead of ~m) and no algebraic square root.
- **Can ECM cofactorization improve SIQS?** ECM can split cofactors into DLP relations, but overhead (~10μs per ECM call) is high vs sieve amortization. Needs DLP graph (union-find) and proper exponent tracking. Potentially useful above 70d where sieve cost dominates.
- **NFS algebraic sqrt progress**: nfs_complete.c has working Hensel-lift AND CRT algebraic sqrt:
  1. **Hensel approach**: Compute S(x) at full precision, take sqrt via Newton iteration (T and (2T)^{-1}). Verified T^2 = S at full precision. But always gives X = ±Y (trivial gcd).
  2. **CRT approach (Couveignes)**: Use large primes where f is irreducible. Multiply S by f'(α)^2 for sign consistency. Evaluate T(m) mod p at each CRT prime. CRT reconstruction gives Y=T(m) mod N since Y < N < Πp. Sign fixed via norm comparison (from stubbscroll/nfs reference implementation).
  3. **Key insight**: CRT of T(m) evaluations DOES work: Y = T(m) mod N is small (< N), so CRT with Πp > N recovers it. My earlier analysis that "CRT won't work" was wrong.
  4. **Remaining issue**: ALL dependencies give trivial X = ±Y even with combined deps of size 20-60, 100 QC columns, and f'(α)^2 trick. The sign correlation is systematic. Root cause likely: incorrect QC computation or missing algebraic conditions in the LA matrix.
  5. **References**: stubbscroll/nfs on GitHub has working Couveignes implementation for degree-3 NFS. Thomé's paper (LORIA) covers all NFS sqrt algorithms.
- **Can we exploit balanced semiprime structure?** No known algorithm specifically targets N = p*q with p ≈ q. The Fermat/Lehman approaches only help when |p-q| is small relative to N^(1/3).
- **Can TLP (triple large primes) help at 80-90 digits?** YES! TLP is now the dominant approach at 70-74d. hyper_siqs with TLP provides 56-61% of relations via SLP combination, dramatically reducing sieve time. Results: 70d=75s (vs 90s without TLP), 73d=170s (all 5), 74d=232s (4/5). TLP + Block Lanczos pushed the frontier from 72d → 74d. The literature claim that TLP doesn't help below 100d is WRONG for our implementations - the DLP→SLP pipeline makes TLP effective at 70d+.

## Agent-4 Specific Findings

### Bucket Sieve Implementation
- Bucket sieve for large primes (p >= BLOCK_SIZE) is CRITICAL at 60d+ where >60% of FB primes exceed 32768
- **KEY BUG**: Bucket overflow handling must create new slices when buckets fill (>75% of BUCKET_ALLOC). Without this, sieve hits are silently dropped, causing 2.3x fewer relations per polynomial
- The YAFU-style bucket format (32-bit entries: 16-bit FB index + 16-bit sieve offset) is efficient and SIMD-friendly

### DLP (Double Large Prime) at 65-70d
- DLP with LP-column approach (adding columns to GF(2) matrix for each unique LP) is NOT viable at 65-70d
- Root cause: LP space (~10^9) is vastly larger than DLP relation count (~500). Birthday paradox requires ~37K DLP relations for 50% collision but we only get ~50 per 1000 polys
- Extending SLP matching to DLP LP range also hurts (lower match rate with larger LP space)
- DLP→SLP pipeline (matching DLP cofactors against SLP hash) is the only effective DLP strategy, contributing 24% of relations at 65d

### Multi-Factor-Base Sieve (Novel, NEGATIVE)
- Concept: sieve with small FB (cheap), trial-divide candidates with extended FB (free for confirmed candidates)
- Result: 2.2x SLOWER than standard SIQS. The sieve must iterate over all FB primes to be effective
- The reduced sieve sensitivity from fewer sieve primes loses more smooth candidates than extended TD recovers
- This is fundamentally because the sieve is an ADDITIVE accumulator - every omitted prime reduces the score, and lowering the threshold only partially compensates

### Block Lanczos vs Gaussian Elimination
- At 70d with FB=10000: Gaussian elimination takes ~74s on 10150x10001 matrix (88% of total time!)
- Block Lanczos on same matrix: 2-4s (20x improvement)
- Key: QS matrices are SPARSE (~20-30 non-zeros per row). BL exploits this; Gauss doesn't.
- BL implemented in lanczos.h: proper Montgomery-style block iteration with 64-wide block vectors
- At 75d: BL eliminates LA as bottleneck entirely. The sieve is the ONLY bottleneck.
- turbo_siqs already uses BL + "Structured GE" for optimal LA performance

### 75d Feasibility Analysis (Updated with TLP)
- Without TLP: 75d needs ~390s for sieve, infeasible
- With TLP (hyper_siqs): 75d[1] still times out. TLP provides ~56% combined rels but the total sieve+combination time still exceeds 295s at 75d
- The TLP approach pushed the frontier from 72d → 74d (4/5), a 2-digit improvement
- 75d may be achievable with further TLP optimization (e.g., better cycle finding, more aggressive LP bounds)
- YAFU can do 75d in ~16-20s, gap is ~15x with TLP (improved from ~20x without)

### Key Constant Factor Analysis
- Custom SIQS implementations are 13-46x slower than YAFU across 30-70d
- ~60% of the gap is from the scalar sieve inner loop (vs YAFU's AVX512BW)
- ~20% from trial division efficiency
- ~10% from polynomial overhead
- ~10% from parameter tuning
- No algorithmic improvement can overcome the AVX512 advantage; need SIMD for competitive performance
