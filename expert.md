# Agent Expert Knowledge

## Key Finding: YAFU SIQS Dominates for Balanced Semiprimes (Single-Core)

YAFU's SIQS (Self-Initializing Quadratic Sieve) is the fastest single-core approach for balanced semiprimes at every tested size from 30-88 digits. For 89+ digits, single-core SIQS needs >300s.

### Why SIQS beats ECM for balanced semiprimes
- ECM's complexity depends on the **smallest factor size**, not the composite
- For balanced semiprimes, factors are ~N/2 digits — ECM treats these as hard
- SIQS complexity depends on the **composite size** — sub-exponential in N
- The crossover where SIQS beats ECM for balanced semiprimes is below 30 digits

### Single-Core Performance (YAFU SIQS, -threads 1 -seed 42)

| Digits | Worst-case time | Notes |
|--------|----------------|-------|
| 30-50  | 0.25s          | Startup-dominated |
| 51-57  | 0.86s          | |
| 58-62  | 1.6s           | |
| 63-68  | 4.4s           | |
| 69-70  | 6.7s           | |
| 71-74  | 14.6s          | |
| 75-78  | 30.1s          | |
| 79-80  | 48.5s          | |
| 81-84  | 111s           | |
| 85-86  | 165s           | |
| 87     | 208s           | |
| 88     | 251s           | Close to limit |
| 89     | ~295s          | Exceeds 300s for hardest numbers |
| 90+    | >300s          | Not achievable single-core |

### YAFU Build Configuration
```bash
make -f Makefile.gcc yafu ECM=1 USE_AVX2=1 SKYLAKEX=1 VBITS=256 -j48
```
- `SKYLAKEX=1`: enables AVX512F, AVX512BW, march=skylake-avx512
- `VBITS=256`: 256-bit Block Lanczos vectors (2x faster LA than VBITS=64)
- `-O2`: standard optimization level
- PGO and LTO tested: **no measurable improvement** (sieve uses hand-written AVX512 intrinsics)
- `-O3`: **no improvement** for same reason

### Parameter Tuning Results
- `-siqsB` (factor base size): Larger values help marginally for 89d (241s at B=40000 vs 254s default for easiest number). But worst-case numbers are unaffected.
- `-siqsNB` (sieve blocks): NB=16 slightly faster than default NB=8 for 89d (215s vs 223s on easiest). Minor effect.
- `-siqsM` (large prime multiplier): No significant effect tested
- `-forceDLP`: No improvement for 89d
- `-siqsTF` (trial factoring bits): Higher value = slower

### Resume Feature
YAFU can save/resume via siqs.dat. Key findings:
- Use `-siqsT <seconds>` for clean shutdown (avoids corrupted save files)
- **DO NOT** use `timeout` command — SIGTERM corrupts siqs.dat causing segfaults on resume
- Resume adds ~25-30s overhead for reloading and reprocessing relations
- Net effect: resume is **slower** than continuous run for borderline numbers
- Only useful if a single factoring must span multiple process invocations

### Why 89+ Fail Single-Core
- 89d: ~70K relations needed at ~3178 sieve operations/sec
- Sieving alone: ~250-295s depending on number
- Block Lanczos: ~45-73s additional
- Total: ~295-370s per number
- No parameter tuning can overcome this: sieve throughput is the hard limit
- Each additional digit adds roughly 35-50% more time

### Alternatives Explored
| Approach | Result |
|----------|--------|
| GMP-ECM | Much slower for balanced semiprimes (B1 ~ 1e13 needed for 45d factor) |
| msieve SIQS | Slower than YAFU single-core |
| CADO-NFS | 720 CPU-sec for 89d (vs 300s YAFU SIQS). NFS only wins above ~100d |
| Custom SIQS (library/siqs.c) | Existing implementation broken (0 relations found). Needs complete rewrite |
| PGO/LTO build | No improvement (hand-written SIMD) |

## Tools
- `yafu/yafu`: YAFU binary (SKYLAKEX+AVX512+VBITS=256). Use `siqs(N)` with `-threads 1 -seed 42`.
- `library/bench_single.sh`: Benchmark script for running 5 semiprimes in parallel
- `library/yafu_resume.sh`: YAFU wrapper with save/resume support (limited use)
