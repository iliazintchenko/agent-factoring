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
| 88     | 240s           | Close to 300s limit (system GMP build) |
| 89     | 213-300+s      | 4/5 pass, 1 timeout. Hardest number exceeds 300s. |
| 90+    | >300s          | Not achievable single-core |

### YAFU Build Configuration (CRITICAL)
```bash
# IMPORTANT: Use system GMP instead of bundled GMP for ~5-13% speedup
# Modify Makefile.gcc: INC += -I/usr/include, LIBS += -L/usr/lib64, -L/usr/local/lib
make -f Makefile.gcc clean && make -f Makefile.gcc yafu ECM=1 USE_AVX2=1 SKYLAKEX=1 VBITS=256 -j48
```
- `SKYLAKEX=1`: enables `USE_AVX512F`, `USE_AVX512BW`, `-march=skylake-avx512`
- `VBITS=256`: 256-bit Block Lanczos vectors (2x faster LA than VBITS=64)
- **System GMP 6.2.1** vs bundled GMP 6.2.0: **5-13% faster** (88d: 240s vs 252s)
- PGO and LTO: **no measurable improvement** (sieve uses hand-written AVX512 intrinsics)
- `-O3` vs `-O2`: **no improvement** for same reason

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

### Why 89d+ Exceeds 300s Single-Core
Timing breakdown for 89d (AVX512BW rebuild):
- **Sieving dominates**: 287-312s for sieving, only 5s for Block Lanczos
- 89d needs ~70K relations at ~3100-3200 sieve ops/sec
- The hardest 89d numbers require 312s of sieving alone — no parameter tuning helps
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
