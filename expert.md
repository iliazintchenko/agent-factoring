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
| 81-84  | 56.3-98.9s     | |
| 85-86  | 142.2-155.4s   | |
| 87     | 199.7s         | |
| 88     | 252.6s         | Close to 280s limit |
| 89     | 511.7s (worst) | 4/5 numbers under 280s, worst=512s. IMPOSSIBLE. |
| 90+    | >600s          | Not achievable single-core |

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
| siqsNB (sieve blocks) | 1, 8 (default) | NB=1 much slower (555s). Default NB=8 is optimal. |
| siqsM (LP multiplier) | 50 vs default | No improvement (324s vs 318s) |
| forceDLP | Yes vs default | No effect (DLP already used by default at this size) |
| forceTLP | Yes | Crashes (unsupported) |
| siqsLPB (LP bound, bits) | 30 vs default ~28 | Slightly slower (380s vs 375s) |

### Why 89d+ Exceeds 280s Single-Core
Timing breakdown for 89d[4] (hardest number, AVX512BW rebuild):
- **Sieving**: ~340s (70K relations needed, ~3200 sieve ops/sec, ~190 useful rels/sec)
- **Block Lanczos**: ~170s
- **Total**: ~510s — neither phase alone fits in 280s
- The other 4 numbers of 89d complete in 222-372s. Only the worst-case fails.
- All parameter tuning exhaustively tested: siqsB (40K-100K), siqsM (100-200), siqsNB (12-20), forceDLP, forceTLP, siqsLPB — ALL time out on this number.
- GNFS via YAFU refuses (gnfs_xover starts at 91d), GNFS via cado-nfs also times out.
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
