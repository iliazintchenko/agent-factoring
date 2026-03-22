# Agent Expert Knowledge

## Key Finding: YAFU SIQS with AVX512BW is Fastest for Single-Core Balanced Semiprimes

YAFU's SIQS (Self-Initializing Quadratic Sieve), rebuilt with AVX512BW sieve kernels, is the fastest single-core approach for balanced semiprimes from 30-88 digits. For 89+ digits, single-core SIQS exceeds 300s.

### Why SIQS beats ECM for balanced semiprimes
- ECM's complexity depends on the **smallest factor size**, not the composite
- For balanced semiprimes, factors are ~N/2 digits — ECM treats these as hard
- SIQS complexity depends on the **composite size** — sub-exponential in N
- The crossover where SIQS beats ECM for balanced semiprimes is below 30 digits

### Single-Core Performance (YAFU SIQS AVX512BW rebuild, -threads 1 -seed 42)

| Digits | Worst-case time | Notes |
|--------|----------------|-------|
| 30-44  | 0.014-0.018s   | **smallmpqs** command (3-7x faster than siqs) |
| 45-50  | 0.08-0.12s     | siqs (smallmpqs segfaults for 45+d) |
| 51-57  | 0.14-0.38s     | |
| 58-62  | 0.58-1.2s      | |
| 63-68  | 1.7-3.7s       | |
| 69-70  | 5.8-6.1s       | |
| 71-74  | 8.0-13.3s      | |
| 75-78  | 18.4-29.2s     | |
| 79-80  | 40.1-47.1s     | |
| 81-84  | 54.8-95.7s     | With -siqsNB 16 |
| 85-86  | 138.8-151.6s   | With -siqsNB 16 |
| 87     | 193.7s         | With -siqsNB 16 |
| 88     | 228.5s         | With -siqsNB 16 |
| 89     | 311s (worst)   | 4/5 under 285s, 89d[4] needs 303s sieving + 8s BL |
| 90+    | >300s          | Not achievable single-core |

### YAFU Build Configuration (CRITICAL)
```bash
make -f Makefile.gcc clean && make -f Makefile.gcc yafu NO_ZLIB=1 ECM=1 USE_AVX2=1 SKYLAKEX=1 VBITS=256 -j48
```
- `SKYLAKEX=1`: enables `USE_AVX512F`, `USE_AVX512BW`, `-march=skylake-avx512`
- **USE_AVX512BW is critical**: Enables hand-written AVX512BW sieve and resieve kernels (`med_sieveblock_32k_avx512bw`, `resieve_medprimes_32k_avx512bw`). These give **14% improvement** over AVX2-only build across all sizes and **3-7x on small sizes** (reduced startup overhead).
- `VBITS=256`: 256-bit Block Lanczos vectors. VBITS=512 is NOT supported (build fails).
- `NO_ZLIB=1`: Required if zlib not installed; avoids link errors.
- Binary: `/tmp/agent-factoring-2/yafu/yafu` (SKYLAKEX rebuild)
- Agent-1 binary also available: `/tmp/agent-factoring-1/yafu/yafu`
- **Previous build bug**: SKYLAKEX was set but `USE_AVX512BW` guards in `med_sieve_32k_avx2.c` and `tdiv_resieve_32k_avx2.c` weren't being triggered, leaving AVX512BW sieve functions undefined. Fixed by ensuring proper `#ifdef USE_AVX512BW` compilation.

### Parameter Tuning Results (89d as test case)
All tested on 89d[3] (hardest number, 375s with old binary, 293s with AVX512BW):

| Parameter | Values Tested | Result |
|-----------|--------------|--------|
| siqsB (FB size) | 40000, 50000, 55000, 60000, 90000 vs default 69888 | Default is optimal. Smaller FB = slower sieve, larger FB = slower LA. |
| siqsNB (sieve blocks) | 1-64 tested | **NB=10 for 73-79d, NB=12 for 80-81d, NB=14 for 82-84d, NB=18 for 85-89d**. Default optimal for <73d. |
| siqsB (FB override) | 50K-120K | **B=70K for 85d, B=80K for 86/88d, B=90K for 87d, B=100K for 89d**. Only helps 85+d. |
| siqsM (LP multiplier) | 50 vs default | No improvement (324s vs 318s) |
| forceDLP | Yes vs default | No effect (DLP already used by default at this size) |
| forceTLP | Yes | Works (see TLP section below), but 44% slower due to batch overhead |
| siqsLPB (LP bound, bits) | 30 vs default ~28 | Slightly slower (380s vs 375s) |
| siqsMFBD (DLP cofactor bound) | 1.7, 1.8, 1.85, 1.9(default), 1.95, 2.0 | Default 1.9 is optimal. MFBD=1.8 is 1% better on individual numbers but 2% worse on worst-case across all sizes. |
| siqsTF (closnuf override) | 89, 90, 93, 95 on 88d | All within noise of default (159-166s). No improvement. |
| NB (wider range) | 16, 20, 22, 24 on 89d[0] | All ~215s, no difference from NB=18. |
| NB on 71-76d | 8-12 tested | Default is already optimal for 71-76d. NB tuning only helps 73+d. |

### Why 89d+ Exceeds 300s Single-Core
Timing breakdown for 89d (system GMP build, -siqsNB 18):
- **With NB=18, B=100K**: 89d[0-3] pass (204-274s). 89d[4] = ~295s (borderline).
- **Sieving dominates**: 222s of 228s total for 88d (98% sieve, 2% Block Lanczos)
- 89d[4] needs ~100K relations at ~4800 rels/sec, requiring ~290+s
- **Branch misprediction**: perf shows 5.96% branch miss rate, ~20% cycle waste (inherent to sieve algorithm)
- **IPC**: 2.67 instructions/cycle (good for AVX512 workload)
- All parameters exhaustively tested:
  - siqsB (40K-120K), siqsNB (8-64), siqsM (50-200), forceDLP, forceTLP, forceQLP
  - siqsLPB, siqsMFBD (1.7-2.0), siqsTF (89-95), combined NB+B, NB+M
  - NB=16-32 with B=80K-100K on 89d — all ~215s on 89d[0]
- GNFS via YAFU: -xover 85 works for 90d but needs ~365s total (47s polyselect + 318s sieve + LA)
- CADO-NFS single-core: server/client overhead prevents efficient single-core operation
- PGO/LTO/O3/march=znver4: no improvement (hand-written AVX512 intrinsics)
- Each additional digit adds ~35-50% more sieving time

### TLP (Triple Large Primes) Analysis
YAFU supports TLP (use_dlp=2) with `-forceTLP`. Key findings:
- TLP **does work** on CPU single-threaded (batch GCD + micro-ECM cofactoring)
- The batch processing code in tdiv.c handles single-threaded mode when `sconf->rb[0].num_relations > target_relations`
- For 85d: TLP = 134.6s vs DLP = 93.3s (**44% slower**)
- Only 1.5% of cycles use TLP relations (689 of 46690 partial cycles)
- TLP overhead (batch init ~4s, batch GCD processing, 3LP filtering/graph building) outweighs the benefit of extra relations
- TLP is designed for 100+ digit numbers where sieve dominance is extreme
- The `-siqsBT` parameter controls batch size (default 1M, but works with smaller values for testing)
- Code location: batch processing in `tdiv.c:574-434`, TLP filtering in `SIQS.c:325` (was disabled with `if(0)` but the tdiv.c code path works)

### GNFS for 90d Analysis
YAFU GNFS with `-xover 85` on 90d:
- Polynomial selection: ~47s (degree 4 poly)
- Sieve: needs 1,460,000 relations at ~4500 rels/sec = ~325s
- Total: ~365-380s — **65-80s over budget**
- GGNFS lattice sievers are functional but not fast enough single-threaded
- CADO-NFS single-thread: server/client architecture adds overhead; stuck in polling

### Alternatives Explored
| Approach | Result |
|----------|--------|
| GMP-ECM | Much slower for balanced semiprimes |
| msieve SIQS | 2.3x slower than YAFU single-threaded (80d: 142s vs 62s) |
| CADO-NFS NFS | Uses multiprocessing (1561 CPU-sec for 89d via multiple worker processes). Violates single-core rule. |
| CADO-NFS SIQS | Not a standalone SIQS; the `sieve/siqs` binary is actually lattice sieve |
| Custom SIQS | library/siqs.c exists but broken |
| VBITS=512 | Not supported by msieve's lanczos.h (limited to 64/128/256) |
| C-Quadratic-Sieve (Michel Leonard) | 10x slower than YAFU on 60d, fails on 70d+. Not competitive. |
| TLP (forceTLP) | 44% slower than DLP on 85d. Overhead outweighs benefit at <100d. |
| GNFS (YAFU -xover 85) | ~365s for 90d. Too slow for 300s budget. |
| Balanced semiprime exploitation | No known algorithm exploits balance. SIQS is factor-structure-agnostic. Fermat/Lehman only help when |p-q| < N^(1/3). |

### Sieve Architecture (for future optimization attempts)
- **Hot function**: `med_sieveblock_32k_avx512bw()` — 32-way SIMD, 64 scattered byte subtractions per iteration
- **Sieve array**: 32KB (fits L1D cache). Bottleneck is read-modify-write throughput, not cache misses.
- **Root update**: AVX512BW masked add with boundary check via `_mm512_cmplt_epu16_mask`
- **Polynomial generation**: Gray code self-initialization, O(1) amortized per polynomial
- **Parameter auto-tuning**: Table in `siqs_aux.c:325-376`, interpolated by bit count
- **AMD EPYC 9R45**: L1D=48KB (sieve fits), no efficient byte scatter/gather, manual unroll is optimal

### Resume Feature
YAFU can save/resume via siqs.dat. Key findings:
- Use `-siqsT <seconds>` for clean shutdown (avoids corrupted save files)
- **DO NOT** use `timeout` command — SIGTERM corrupts siqs.dat
- Resume adds ~25-30s overhead — slower than continuous run for borderline numbers
- Only useful if factoring must span multiple invocations

## Tools
- `yafu/yafu`: YAFU binary (AVX512BW+VBITS=256). Use `siqs(N)` with `-threads 1 -seed 42`.
- `library/bench_single.py`: Benchmark script (run 5 semiprimes per size in parallel)
- `library/update_best.py`: Update best-algos.json from benchmark results
- `library/yafu_resume.sh`: YAFU wrapper with save/resume support
