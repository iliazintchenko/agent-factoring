# Agent Expert Knowledge

## Key Finding: YAFU SIQS with AVX512BW is Fastest for Single-Core Balanced Semiprimes

YAFU's SIQS (Self-Initializing Quadratic Sieve), rebuilt with AVX512BW sieve kernels, is the fastest single-core approach for balanced semiprimes from 30-89 digits. For 90+ digits, neither SIQS nor NFS can complete within 300s single-core.

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
| 89     | 294.4s         | With -siqsNB 18 -siqsB 100000. All 5 pass. |
| 90+    | >300s          | **Definitively infeasible** single-core (see NFS analysis) |

### YAFU Build Configuration (CRITICAL)
```bash
# Step 1: Fix monty.h - add 'static' to ALL __inline functions
# In include/monty.h, change ALL bare '__inline' to 'static __inline':
# submod, addmod, mulredc, sqrredc, mulredc63, sqrredc63
sed -i 's/^__inline uint64_t/static __inline uint64_t/' include/monty.h

# Step 2: Build with AVX512BW
make -f Makefile.gcc clean && make -f Makefile.gcc yafu NO_ZLIB=1 ECM=1 USE_AVX2=1 SKYLAKEX=1 VBITS=256 -j48
```
- `SKYLAKEX=1`: enables `USE_AVX512F`, `USE_AVX512BW`, `-march=skylake-avx512`
- **USE_AVX512BW is critical**: Enables hand-written AVX512BW sieve and resieve kernels (`med_sieveblock_32k_avx512bw`, `resieve_medprimes_32k_avx512bw`). These give **14% improvement** over AVX2-only build across all sizes and **3-7x on small sizes** (reduced startup overhead).
- **monty.h static inline fix**: Add `static` to `__inline` in monty.h (submod/addmod/mulredc etc.). Required for PGO builds. Performance impact: **marginal (1-2%) at best** — initial measurement of 16% was a load-noise artifact. Confirmed by concurrent A/B test: fixed and unfixed binaries produce identical times (44.4s vs 44.7s on 80d[0], 42.0s vs 42.1s on 80d[2]).
- `VBITS=256`: 256-bit Block Lanczos vectors. VBITS=512 is NOT supported (build fails).
- `NO_ZLIB=1`: Required if zlib not installed; avoids link errors.
- Binary: `/tmp/agent-factoring-10/yafu_build/yafu` (monty.h fix + AVX512BW, best known build)
- Agent-4 binary (no monty.h fix): `/tmp/agent-factoring-4/yafu/yafu`
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
- GNFS via YAFU: works with `-xover 85` + GGNFS sievers. 90d: 404s total (55s poly + 300s sieve + 50s filter/LA).
- CADO-NFS single-core (`-t 1 --slaves 1`): 90d: ~570s est. las siever 1500 rels/sec (3x slower than GGNFS)
- msieve NFS: 90d: poly select 167s alone, sieve >240s additional. Total >400s.
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

### GNFS for 89-90d Analysis
YAFU GNFS with `-xover 85` and default GGNFS sievers:
- **90d[0] completed in 280s wallclock (259s user)** on low-load machine (load ~8.75)
  - Poly select: ~47s (degree 4, score 4.512e-08)
  - Sieve: ~220s (1.46M rels needed, GGNFS I=12, lpb=25/25, ~23K rels per 5s batch)
  - Filter+LA+sqrt: ~13s
- On loaded machine (load ~10+), same config times out at >295s due to cache contention
- **AVX512 GGNFS sievers are SLOWER** than default sievers (timed out vs 280s)
- `-psearch min` reduces poly select from 47s to 15s but worse poly increases sieve time
- Custom params (lpb=24, I=11 from CADO-NFS) produce 3x fewer rels/q = much slower sieving
- Sieve rate highly sensitive to L3 cache contention from parallel processes
- CADO-NFS las siever crashes with C++ exception (AVX512BW not in build flags)
- msieve standalone NFS: poly select alone takes 148-167s for 89-90d, not viable
- **Key insight**: GNFS is feasible for 89-90d only on idle/low-load machine (user time ~250-260s)

### YAFU SIQS on 90d — Closest Attempt
NB=20 B=120K on 90d (all 5 semiprimes):
- 90d[0]: timeout (>295s)
- 90d[1]: **248.6s** ✓
- 90d[2]: timeout (>295s, sieve done but BL didn't finish)
- 90d[3]: **282.3s** ✓
- 90d[4]: **281.2s** ✓
- 3/5 pass, 2/5 exceed 300s. 90d is definitively not achievable single-core.

### Alternatives Explored
| Approach | Result |
|----------|--------|
| GMP-ECM | Much slower for balanced semiprimes |
| msieve SIQS | 2.3x slower than YAFU single-threaded (80d: 142s vs 62s) |
| msieve NFS (-n flag) | Poly selection alone takes 167s for 90d, then sieving times out. Total >400s. |
| YAFU GNFS (GGNFS sievers) | 90d: 404s total. Poly select ~55s, sieve ~300s (1.46M rels at 4800/sec), filter+LA ~50s. |
| CADO-NFS single-thread | 90d: 245s CPU sieve for 43% of relations, ETA ~570s total. las sieve rate 1497 rels/sec (3x slower than GGNFS) but needs fewer rels (852K vs 1.46M). |
| CADO-NFS NFS | Uses multiprocessing. Violates single-core rule. |
| AVX512 gather-scatter sieve | **20% slower** than scalar on AMD Zen4. AMD's scatter is ~10 cycles/element vs 1 cycle/element for scalar stores. |
| Custom SIQS | library/siqs.c exists but broken |
| VBITS=512 | Not supported by msieve's lanczos.h (limited to 64/128/256) |
| C-Quadratic-Sieve (Michel Leonard) | 10x slower than YAFU on 60d, fails on 70d+. Not competitive. |
| yamaquasi (Rust SIQS) | 2.3x slower than YAFU (70d: 13.6s vs 5.8s, 85d: 264s vs 136s). |
| TLP (forceTLP) | 44% slower than DLP on 85d. Overhead outweighs benefit at <100d. |
| GNFS (YAFU -xover 85) | ~365s for 90d. Too slow for 300s budget. |
| USE_BATCHPOLY build | 11% slower (80d: 52.9s vs 47.0s). Batch polynomial root updates hurt performance. |
| inmem=100 (all in-memory) | No improvement on 80d or 88d vs default inmem cutoff of 70d. |
| CPU pinning (taskset) | No improvement. Single-threaded YAFU doesn't benefit from core affinity. |
| PGO (profile-guided optimization) | **WORKS after fix**: Add `static` to `__inline` in `monty.h` for mulredc/sqrredc/mulredc63. PGO binary at `/tmp/agent-factoring-3/yafu/yafu.pgo`. Solo 89d[4]=294.1s (from 297s with O3, 300s with O2). ~1-2% improvement. |
| **monty.h static inline fix (ALL functions)** | Marginal improvement (~1-2%). Add `static` to ALL `__inline` functions in monty.h. Initial "16% improvement" was **load-noise artifact** — A/B test with both binaries running concurrently shows identical times (44.4s vs 44.7s on 80d[0]). The fix is still useful for enabling PGO builds. |
| Balanced semiprime exploitation | No known algorithm exploits balance. SIQS is factor-structure-agnostic. Fermat/Lehman only help when |p-q| < N^(1/3). |
| Custom SIQS (library/cqs/siqs_gmp.c) | Working SIQS in C+GMP. 30d: 0.04s, 40d: 2.9s, 50d: 4.1s. Factors 30-50d but ~50x slower than YAFU for 50d. Missing: Gray code self-init (generates new 'a' per poly instead of 2^(s-1) B values), AVX512 sieve, optimized trial division. The gap is fundamental engineering, not algorithmic. |
| YAFU 48KB sieve blocks | BLOCKSIZE is hardcoded to 32768 with values scattered throughout AVX512BW assembly. Changing to 48KB (AMD L1D size) would require extensive surgery. Not viable. |

### Sieve Architecture (for future optimization attempts)
- **Hot function**: `med_sieveblock_32k_avx512bw()` — 32-way SIMD, 64 scattered byte subtractions per iteration
- **Sieve array**: 32KB (fits L1D cache). Bottleneck is read-modify-write throughput, not cache misses.
- **Root update**: AVX512BW masked add with boundary check via `_mm512_cmplt_epu16_mask`
- **Polynomial generation**: Gray code self-initialization, O(1) amortized per polynomial
- **Parameter auto-tuning**: Table in `siqs_aux.c:325-376`, interpolated by bit count
- **AMD EPYC 9R45**: L1D=48KB (sieve fits), no efficient byte scatter/gather, manual unroll is optimal

### NFS vs SIQS Crossover Analysis (90d single-core)
- **SIQS (YAFU)**: 90d easiest = 245s, hardest = >300s (timeout). Needs ~72K relations at ~200 useful/sec.
- **GNFS (YAFU+GGNFS)**: 90d = 404s. Sieving dominates (~300s for 1.46M relations). Poly select ~55s.
- **NFS (msieve built-in)**: 90d = >400s. Poly select ~167s + sieve >240s.
- **NFS (CADO-NFS las)**: 90d estimated ~570s. Better filtering (852K vs 1.46M rels) but slower siever.
- **Conclusion**: SIQS is faster than NFS for 90d. Crossover is around 100-110d for single-core. 90-100d is a dead zone where both algorithms exceed 300s.

### 90d Parameter Exhaustion
Tested on 90d[0] (hardest number): ALL combos fail under 295s:
- NB: 10, 12, 14, 16, 18, 20, 22, 24 — all timeout
- B: 80K, 90K, 100K, 110K, 120K, 130K, 150K — all timeout
- MFBD: 2.0, 2.2 — no help
- LPB: 30, 32 — no help
- M: 150, 200 — no help
90d confirmed infeasible with any SIQS parameter combination.

### YAFU Source-Level Analysis (compile-time options)
- **USE_BATCHPOLY**: Commented out in `qs_impl.h`. Would batch bucket root updates every 4 polys. But sets `FORCE_GENERIC=1`, disabling AVX512BW sieve. Author's comment: "unable to get this to run faster."
- **Parameter table** (`siqs_aux.c:327`): AVX512F table has NB=8 for 281-298 bits (85-90d). Our experiments show NB=14-18 + B=70-100K are optimal. But modifying the table hurts because the interpolation code averages NB from bounding rows instead of interpolating. Always use `-siqsNB` and `-siqsB` command-line overrides.
- **Sieve kernel** (`med_sieve_32k_avx2.c:689`): Store-to-load forwarding stall in inner loop (512-bit store then 16-bit loads) is unavoidable on AMD Zen4. No vpscatterb instruction exists.
- **Block Lanczos**: Only 2-8% of total time at 85-89d. VBITS=512 would halve BL time but isn't supported.
- **DLP cofactoring**: Batch GCD path disabled (`&& 0`). Online microECM per candidate is faster for DLP.

### Resume Feature
YAFU can save/resume via siqs.dat. Key findings:
- Use `-siqsT <seconds>` for clean shutdown (avoids corrupted save files)
- **DO NOT** use `timeout` command — SIGTERM corrupts siqs.dat
- Resume adds ~25-30s overhead — slower than continuous run for borderline numbers
- Only useful if factoring must span multiple invocations

## Custom SIQS Implementations

### siqs2.c — Working Custom SIQS (agent-6 rewrite)
- **Status**: Working, 30-50d factoring confirmed, correct results
- **Performance**: 30d: 0.044s, 40d: 0.8s, 50d: 10s
- **Comparison to YAFU**: 3-80x slower (YAFU: 30d=0.014s, 40d=0.017s, 50d=0.12s)
- **Features**: Knuth-Schroeppel multiplier, Gray code self-init, 32KB block sieve, SLP, a*Q(x) tracking
- **Key bugs fixed**: sieve offset, GF(2) matrix a-primes, loop termination

### siqs_fast.c — Working Custom SIQS (agent-1, with DLP)
- **Status**: Working, 30-50d confirmed, DLP support
- **Performance**: 30d: 0.6s, 40d: 0.95s, 50d: 21.7s
- **Features**: AVX512BW sieve scanning, Gray-code poly switching, SLP+DLP, GF(2) GE
- **Known issue**: Multiplier (k>1) causes trivial congruences in sqrt step — disabled for now
- **Key bug**: Mirror positions Q(x)=Q(-x-2b/a) must be deduplicated (skip Y<0)

### Key Insights for Custom SIQS
1. **Multiplier handling**: With kN (k>1), sqrt step systematically produces X ≡ ±Y (mod N). May relate to LP products interacting with multiplier.
2. **Mirror positions**: Q(x) = Q(-x-2b/a), so both sides of sieve give identical Q. Must skip negative Y to avoid trivial SLP combinations.
3. **Scalar sieve bottleneck**: Without AVX512 scattered byte ops, custom SIQS is 40-180x slower than YAFU.
4. **SQUFOF/Rho don't work**: O(N^{1/4}) is dramatically slower than SIQS's sub-exponential for 30+d balanced semiprimes.
5. **The a*g(x) factor**: For SIQS, exponent matrix must track a*g(x), not g(x).
6. **Sieve threshold**: log2(M * sqrt(N)) minus small-prime correction. Too high = few candidates, too low = false positives.

## PGO Build
- PGO for YAFU: ~1-2% improvement at best
- Fix needed: `__inline` → `static __inline` in monty.h
- Not impactful enough to change best-algos.json entries

## Tools
- `yafu/yafu`: YAFU binary (AVX512BW+VBITS=256). Use `siqs(N)` with `-threads 1 -seed 42`.
- `library/siqs_fast.c`: Custom SIQS (DLP, working 30-50d). Compile: `gcc -O3 -march=native -mavx512bw -o siqs_fast library/siqs_fast.c -lgmp -lm`
- `library/siqs2.c`: Custom SIQS (SLP only, working 30-60d)
- `library/mpqs.c`: Custom MPQS (sieve works, sqrt step buggy)
- `library/pollard_rho.c`: Pollard rho (too slow for balanced semiprimes)
- `library/bench_single.py`: Benchmark script (run 5 semiprimes per size in parallel)
- `library/update_best.py`: Update best-algos.json from benchmark results
- `library/yafu_resume.sh`: YAFU wrapper with save/resume support
