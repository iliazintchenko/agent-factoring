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
| 30-50  | 0.04-0.12s     | Startup-dominated |
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
| 89     | 307s (worst)   | 4/5 under 300s, 89d[4] needs 301s sieving |
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
| siqsNB (sieve blocks) | 1, 8 (default), 12, 16, 24 | **NB=16 is best for 80+d (2-10% improvement)**. NB=1 much slower. NB=16 hurts <80d. |
| siqsM (LP multiplier) | 50 vs default | No improvement (324s vs 318s) |
| forceDLP | Yes vs default | No effect (DLP already used by default at this size) |
| forceTLP | Yes | Crashes (unsupported) |
| siqsLPB (LP bound, bits) | 30 vs default ~28 | Slightly slower (380s vs 375s) |

### Why 89d+ Exceeds 300s Single-Core
Timing breakdown for 89d (system GMP build, -siqsNB 16):
- **With NB=16**: 89d[0-3] = 207-282s (all pass). 89d[4] = >300s (timeout).
- **Sieving dominates**: 222s of 228s total for 88d (98% sieve, 2% Block Lanczos)
- 89d[4] needs ~70K relations at ~3200 sieve ops/sec, requiring ~300+s
- **Branch misprediction**: perf shows 5.96% branch miss rate, ~20% cycle waste (inherent to sieve algorithm)
- **IPC**: 2.67 instructions/cycle (good for AVX512 workload)
- All parameters exhaustively tested: siqsB, siqsNB (8-64), siqsM, siqsMFBD, siqsSSidx, forceDLP, forceTLP, siqsLPB, -noopt, -siqsNobat, -inmem — none help
- GNFS via YAFU refuses (gnfs_xover starts at 91d)
- CADO-NFS single-core: ~720 CPU-sec (much slower)
- PGO/LTO/O3/march=znver4: no improvement (hand-written AVX512 intrinsics)
- Each additional digit adds ~35-50% more sieving time

### Alternatives Explored
| Approach | Result |
|----------|--------|
| GMP-ECM | Much slower for balanced semiprimes |
| msieve SIQS | 2.3x slower than YAFU single-threaded (80d: 142s vs 62s) |
| CADO-NFS NFS | Uses multiprocessing (1561 CPU-sec for 89d via multiple worker processes). Violates single-core rule. |
| CADO-NFS SIQS | Not a standalone SIQS; the `sieve/siqs` binary is actually lattice sieve |
| Custom SIQS | library/siqs.c exists but broken |
| VBITS=512 | Not supported by msieve's lanczos.h (limited to 64/128/256) |

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
