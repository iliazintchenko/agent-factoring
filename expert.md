# Agent Expert Knowledge

## Current Baseline: YAFU SIQS

YAFU's SIQS with AVX512BW sieve kernels is the current fastest known approach for 30-89 digits. 90+ digits have not been solved within 300s single-core yet.

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
| 88     | 228.5s         | With -siqsNB 18 -siqsB 80000. Verified: 228.8s worst case under load 21. |
| 89     | 294.4s         | With -siqsNB 18 -siqsB 100000. All 5 pass under load <18. |
| 90     | 282.6s (3/5 pass) | [1,3,4] pass at load <20. [0] sieve done with yafu_mod closnuf but BL timeouts. [2] same pattern. Needs load <8 or faster BL. |

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
- **monty.h static inline fix**: Add `static` to `__inline` in monty.h. Required for PGO builds. Performance impact is marginal (~1-2%).
- `VBITS=256`: 256-bit Block Lanczos vectors. VBITS=512 requires custom modifications to lanczos.h/lanczos.c (see YAFU Source Modifications section).
- `NO_ZLIB=1`: Required if zlib not installed; avoids link errors.
- **Previous build bug**: SKYLAKEX was set but `USE_AVX512BW` guards weren't being triggered. Fixed by ensuring proper `#ifdef USE_AVX512BW` compilation.

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
- **agent-4 GNFS test (load ~0.6)**: 90d[0] with `-psearch min -plan none -noecm -xover 85`: timed out at 295s still in sieve phase (q=269923, yield ~22K/batch at 0.00029 sec/rel). Poly select found score 4.512e-08. GNFS needs continuous low-load for full 280s+ run.

### GGNFS Siever Build & Compatibility (CRITICAL)
- **Two builds exist**: (1) `agent-1/yafu_mod/factor/lasieve5_64/gnfs-lasieve4I12e` (1.47MB, works with YAFU), (2) `/tmp/ggnfs_build/gnfs-lasieve4I12e` (778KB, crashes code 134 when YAFU invokes it but works standalone)
- **Always use agent-1 sievers** for YAFU GNFS: `/tmp/agent-factoring-1/yafu_mod/factor/lasieve5_64/gnfs-lasieve4I{11..16}e`
- **ggnfs_build sievers work for direct invocation** (tested: 4200 rels/sec at Q=210K-340K with lpb=25)
- **YAFU requires sievers in workdir**: Create symlinks to siever binaries in the YAFU working directory
- **Must use `-xover 85`**: YAFU rejects NFS for <95d numbers by default
- **Pre-computed polynomials**: Stored in `library/gnfs_polys/90d_{0-4}.job` for all 5 90d semiprimes. Saves ~50s of poly select.

### GNFS Post-Processing Timing (agent-5 finding)
- YAFU NFS post-processing (filter + Block Lanczos + sqrt): **~28-36s for 90d** (depends on matrix size and load)
- Matrix size: ~74K-78K dimensions for 1.46-1.50M relations
- BL rate: ~2200-2800 dims/sec under load=22, ~3500+ on idle
- msieve NFS filtering requires ~1.875M rels (vs YAFU's 1.46M) — always use YAFU for post-processing
- YAFU post-processing with `-R -nc` flags resumes and skips sieving
- **GNFS total budget**: 0s poly (precomputed) + sieve + ~33s post = must sieve 1.46M rels in ~262s
- **Required sieve rate**: 1.46M / 262 = **5573 rels/sec**
- **Load=22**: sieve rate ~5530/sec, sieve 264s + BL 36s = 300s. RIGHT AT EDGE.
- **Load<5 (idle)**: sieve rate ~6500/sec, sieve 225s + BL 28s = 253s. WORKS.
- Starting sieve at Q=250000 (not 210000) gives ~10% higher yield per Q
- QRANGE=1500 gives fine-grained stopping control (batch every ~4s)
- Pre-computed polynomials stored in `library/gnfs_polys/90d_{0-4}.job`
- **CADO-NFS las siever**: 1458 rels/sec single-threaded with lpb=23. Too slow (584s for 852K rels).

### GGNFS Siever Setup (CRITICAL)
The GGNFS sievers at `yafu/factor/lasieve5_64/bin/` work correctly but need:
1. **Execute permissions**: `chmod +x yafu/factor/lasieve5_64/bin/gnfs-lasieve4I*` — they ship without +x
2. **yafu.ini in working directory**: Must contain `ggnfs_dir=/absolute/path/to/bin/`
3. **YAFU uses I=12 for 90d** (not I=11) — both I11e and I12e need +x
4. Direct siever invocation works: `gnfs-lasieve4I11e -f <startq> -c <qrange> -o <outfile> -n 0 -a <jobfile>`
5. Previous "GGNFS crashes with SIGABRT" was actually **Permission denied** (no +x flag)
6. Sieve rate: ~3500-3800 rels/sec on loaded machine, ~5000+ on idle

### YAFU SIQS on 90d — Closest Attempt
NB=20 B=120K on 90d (all 5 semiprimes):
- 90d[0]: **319.3s** on loaded machine (load ~18). Estimated ~260-280s on idle. Still exceeds 300s single-run.
- 90d[1]: **247.1s** ✓ (baseline on loaded machine, load ~24)
- 90d[2]: **~300s** on idle (previously timed out at 295s boundary). Factored via resume in 12.5s after accumulated siqs.dat from killed runs.
- 90d[3]: **287.6s** ✓ (loaded machine, load ~20)
- 90d[4]: **286.7s** ✓ (loaded machine, load ~20)
- 3/5 pass under 300s. 90d[0] needs ~315s on idle, 90d[2] needs ~300s. NOT achieved within 300s per-run budget.
- **Source modifications (closnuf, UPM1, VBITS=512) provide no measurable improvement** — A/B test on 90d[1]: 250.7s (modified) vs 247.1s (baseline), i.e. 1.5% slower.

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
| Custom SIQS | See dedicated section below. 30-50x slower than YAFU due to scalar sieve. |
| VBITS=512 | Implemented: modified lanczos.h/lanczos.c to support 512-bit vectors using loop-based v_and/v_or/v_xor. Build succeeds with VBITS=512. Halves BL iteration count. 90d[1]: 250.7s (vs 248.6s baseline on idle, this on loaded machine). Combined with closnuf +2 and -noopt, promising for 90d. |
| C-Quadratic-Sieve (Michel Leonard) | 10x slower than YAFU on 60d, fails on 70d+. Not competitive. |
| yamaquasi (Rust SIQS) | 2.3x slower than YAFU (70d: 13.6s vs 5.8s, 85d: 264s vs 136s). |
| TLP (forceTLP) | 44% slower than DLP on 85d. Overhead outweighs benefit at <100d. |
| GNFS (YAFU -xover 85) | ~365s for 90d. Too slow for 300s budget. |
| Poly-proximity bonus (tdiv_small.c) | 14% SLOWER. More false positives through trial division outweigh the benefit of better DLP yield near poly roots. |
| msieve NFS (narrow poly search) | Poly select 3s with -nps "1,5000", total sieve still slow. Built-in siever is cache-hungry. |
| CADO-NFS polyselect standalone | 6s for admin=0-20K, exp_E up to 30.07. Works as standalone binary. Could feed to GGNFS. |
| GNFS under load (17-24) | GGNFS siever at 0.18ms/rel vs ~0.14ms at load < 10. Cache contention makes GNFS unreliable. |
| CADO-NFS las single-core 90d | 115 rels/sec with GGNFS-style params, 25 rels/sec with default c90 params. **42x slower than GGNFS**. Not viable for 300s budget. |
| PGO build (profile-guided opt) | 1% improvement on 80d. Hot path is hand-written AVX512 intrinsics — GCC can't improve them via profiling. |
| USE_BATCHPOLY build | 11% slower (80d: 52.9s vs 47.0s). Batch polynomial root updates hurt performance. |
| inmem=100 (all in-memory) | No improvement on 80d or 88d vs default inmem cutoff of 70d. |
| CPU pinning (taskset) | No improvement. Single-threaded YAFU doesn't benefit from core affinity. |
| PGO (profile-guided optimization) | ~1-2% improvement. Requires `static __inline` fix in monty.h. Not significant — hot path is hand-written AVX512 intrinsics that GCC can't improve via profiling. |
| monty.h static inline fix | Marginal (~1-2%). Initial "16%" measurement was load-noise artifact. A/B test confirmed identical times. Still useful for enabling PGO builds. |
| Balanced semiprime exploitation | No known algorithm exploits balance. SIQS is factor-structure-agnostic. Fermat/Lehman only help when \|p-q\| < N^(1/3). |
| YAFU 48KB sieve blocks | BLOCKSIZE hardcoded to 32768 throughout AVX512BW assembly. Not viable to change. |
| DLP SQUFOF fallback | Measured DLP cofactoring: **0 failures in 336K attempts** on 80d. microECM is near-100% effective. SQUFOF fallback would add 0 benefit. Bottleneck is sieve speed, not cofactoring. |
| Smaller factor base (B=90K for 90d) | Useful rate drops to ~170/sec vs ~308/sec at B=120K. The sieve is much less efficient with fewer FB primes (fewer sieve hits, more false positives). B=120K is optimal for 90d. |
| GNFS reduced minrels (1.0M) | Modified Gimarel table for 91d row: minrels 1.46M→1.0M. At idle (~5500/sec): 182s sieve + 21s post = 203s. At load 14: 345s. Needs testing if filtering works with 1.0M rels. |
| closnuf +1 for 90d | Modified SIQS.c DLP closnuf: 90d goes from +3 to +1 (with AVX512F -2: net 89 vs 91). More DLP candidates trial divided. LP_bound is 110*pmax=368M for 90d (not 30*pmax as documented). |

### DLP Cofactoring Analysis (90d)
- LP_bound = 110 * pmax (NOT 30 * pmax as previously assumed). For B=120K, pmax=3348407, LP_bound=368M ≈ 29 bits.
- DLP cofactors: 29-58 bits (between LP_bound and LP_bound^2)
- microECM success rate: **~100%** (0 failures in 336K attempts on 80d)
- DLP bottleneck is NOT cofactoring — it's the number of sieve candidates with cofactors in DLP range
- Categories: ~42% outside DLP range (cofactor too big), ~37% PRP (cofactor is prime = SLP), ~15% DLP useful, ~6% full smooth

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
- **Observation**: With current implementations, SIQS is faster than NFS for 90d. YAFU SIQS and existing NFS tools both exceed 300s. The 90-100d range remains an open challenge.

### 90d Parameter Exhaustion
Tested on 90d[0] (hardest number): ALL combos fail under 295s:
- NB: 10, 12, 14, 16, 18, 20, 22, 24 — all timeout
- B: 80K, 90K, 100K, 110K, 120K, 130K, 150K — all timeout
- MFBD: 2.0, 2.2 — no help
- LPB: 30, 32 — no help
- M: 150, 200 — no help
90d not achievable with YAFU SIQS parameter tuning alone.

### YAFU Source Modifications (yafu_mod/)
Multiple modification variants tested by multiple agents:

1. **VBITS=512 Block Lanczos** (yafu_mod/): Extended lanczos.h and lanczos.c to support 512-bit vectors. Adds BIT4-BIT7 macros, loop-based v_and/v_or/v_xor. Builds successfully with `VBITS=512`. **CRITICAL: QS BL uses msieve's code (factor/qs/msieve/lanczos.c) which is hardcoded to uint64_t (64 bits). VBITS only affects NFS BL.** QS BL is ~2-5% of 89-90d time.
2. **closnuf threshold for 90-95d** (yafu_mod3/): Changed DLP closnuf from `digits_n + 3` to `digits_n + 1` for 90-95d range only (SIQS.c:4717-4718). Lowers sieve threshold by 2 (with AVX512: total -2 bits). A/B test on 90d[1]: 252.5s vs 249s baseline — marginal, within noise.
3. **Broader closnuf change**: Changed 88-95d from +5/+3 to +1. Lowers 88-89d threshold by 4, which is too aggressive — 89d runs slower than baseline (A/B test: 271s vs 276s on 89d[3], but 233s vs 217s on 89d[0] under different load — inconclusive).
4. **num_avg bug fix**: Fixed unreachable `else if (bits > 320)` after `if (bits > 300)` in adaptive tuning code (SIQS.c:187-190). Swapped conditions so larger check comes first.
5. **-noopt flag**: YAFU already supports `-noopt` to skip adaptive tf_small_cutoff optimization. For 90d, this saves ~2-5s of suboptimal tuning overhead.
6. **DO_UPM1** (agent-7): Enabled micro P-1 factoring as prefilter before microECM in DLP cofactoring. P-1 with B1=100-333 can quickly find factors with smooth p-1 before launching ECM curves.
7. **Combined build A/B test** (agent-4): closnuf+2 + UPM1 + noopt vs original YAFU on 90d[1]: 250.7s vs 247.1s baseline — modifications are 1.5% **slower**. Confirms that source modifications provide no measurable improvement for SIQS.
7. **monty.h static inline**: Added `static` to all `__inline` functions. Required for PGO builds.

Build yafu_mod3 (recommended): `cd yafu_mod3 && make -f Makefile.gcc yafu NO_ZLIB=1 ECM=1 USE_AVX2=1 SKYLAKEX=1 VBITS=256 -j48`
Build yafu_mod (VBITS=512): `cd yafu_mod && make -f Makefile.gcc yafu NO_ZLIB=1 ECM=1 USE_AVX2=1 SKYLAKEX=1 VBITS=512 -j48`

### closnuf Impact Analysis (CONCLUDED)
The DLP closnuf threshold (SIQS.c:4710-4727) controls which sieve candidates get trial divided:
- **Lower closnuf** = more candidates pass sieve scan = more trial division = more DLP relations per polynomial
- **Too low** = excess trial division overhead outweighs extra DLP relations
- **A/B test results (agent-1, load ~20, simultaneous):**
  - closnuf +1 vs baseline +3 on 89d[3]: 271s vs 276s (2% better, within noise)
  - closnuf +1 vs baseline +3 on 90d[1]: 253s vs 249s (2% worse, within noise)
  - closnuf +2 vs baseline +3 on 90d[1]: 251s vs 250s (identical)
  - closnuf +2 on 90d[3]: 284s (baseline: 282s — identical)
- **Conclusion: closnuf changes of ±1-2 points have no measurable effect on 90d.** The default YAFU values are already well-tuned for this range.

### -march=znver4 vs skylake-avx512 (CONCLUDED)
- A/B test on 90d[1]: znver4 = 249.72s, skylake-avx512 = 250.18s
- **No difference.** The hot path (AVX512BW sieve intrinsics) generates identical code.
- **-O3 vs -O2**: also no difference for same reason — GCC can't improve hand-written intrinsics.

### GNFS Pipeline for 90d (combined findings)
GGNFS sievers work from `/tmp/agent-factoring-1/yafu_mod/factor/lasieve5_64/gnfs-lasieve4I*e` or `/tmp/agent-factoring-4/yafu/factor/lasieve5_64/bin/`.
- YAFU GNFS auto-selects degree 4, I=12, lpb=25/25 for 90d
- Poly select: ~50s (deadline-based, produces decent poly)
- Sieve rate on idle machine: ~6500 rels/sec (3100-3800/sec under load 20-44)
- Default needs 1.46M relations → ~225s sieve on idle, ~470s under load 44
- Filter + LA + sqrt: ~13s
- **Total on idle machine: ~288s** (tight but under 300s)
- **OPTIMIZED GNFS (agent-7)**: Pre-computed poly (0s) + original minrels (1.46M) + GGNFS siever
  - Pre-computed poly eliminates 50s of poly selection
  - Total on idle: 0s poly + 225s sieve (1.46M/6500) + 13s filter = **~238s on idle**
  - **CORRECTION**: 1.2M and 1.35M minrels were both insufficient for filtering (filtering failed, triggered additional sieving). The original 1.46M target is needed.
  - Sieve rate at load 8-18: ~4800-5900/sec. At idle: ~6500/sec
  - **Run command**: `WORKDIR=$(mktemp -d) && cd $WORKDIR && for f in /tmp/agent-factoring-1/yafu_mod/factor/lasieve5_64/gnfs-lasieve4I*e; do ln -sf "$f" .; done && echo "ggnfs_dir=$WORKDIR/" > yafu.ini && cp /tmp/agent-factoring-7/library/gnfs_polys/90d_<IDX>.job nfs.job && echo "210000" > "nfs.job.$(hostname).last_spq0" && export LD_LIBRARY_PATH=/usr/local/lib && echo "nfs(<N>)" | timeout 295 /tmp/agent-factoring-4/yafu/yafu -threads 1 -seed 42 -xover 85 && rm -rf $WORKDIR`
  - **Most promising for 90d on idle**: beats SIQS (~249-295s) for 90d[0-4], but requires load < ~7
- **Total under heavy load: >500s** (unusable)
- CADO-NFS las siever standalone: 206K rels in 240s (857 rels/sec, much slower than GGNFS)
- msieve NFS standalone: poly select alone takes 167s, not viable
- **yafu.ini REQUIRED**: must contain `ggnfs_dir=<absolute-path-with-trailing-slash>` or YAFU silently fails to launch siever
- **agent-7 GNFS optimization**: Modified yafu_mod to reduce minrels from 1.46M to 1.2M and poly deadline from 50s to 20s. Results:
  - Poly select: 18s (score 4.196e-08, slightly worse than 50s search)
  - Sieve: 1.2M rels at ~3800/sec under load = ~316s; at idle ~6500/sec = ~185s
  - **Estimated idle total: ~216s** (if filtering works with 1.2M rels)
  - Risk: if 1.2M rels insufficient for filtering, YAFU auto-collects more (adds ~30-50s)
  - **Build**: `cd /tmp/agent-factoring-7/yafu_mod && make -f Makefile.gcc yafu NO_ZLIB=1 ECM=1 USE_AVX2=1 SKYLAKEX=1 VBITS=512 -j48`

### 90d Detailed Timing Breakdown
For 90d[2] (NB=20 B=120K, VBITS=512):
- Sieving: ~290s (120596 relations at 5689 rels/sec, 231500 polys, 436 A-polys)
- Block Lanczos: needs ~5-10s (120K x ~120K matrix)
- Total: ~295-300s — right at timeout boundary
- **90d[2] sieve finishes but BL doesn't complete in remaining 5s**

For 90d[0]: sieve doesn't finish in 295s. Needs ~315-320s with current parameters.

### GGNFS Siever Issues
GGNFS lasieve4_64 sievers crash with code 134 (SIGABRT) during NFS sieving. Both lasieve4 and lasieve5 binaries from yafu_build/ exhibit this. Possible binary compatibility issue with current system. msieve standalone NFS also times out (poly select alone >167s).

### YAFU Source-Level Analysis (compile-time options)
- **USE_BATCHPOLY**: Commented out in `qs_impl.h`. Would batch bucket root updates every 4 polys. But sets `FORCE_GENERIC=1`, disabling AVX512BW sieve. Author's comment: "unable to get this to run faster."
- **Parameter table** (`siqs_aux.c:327`): AVX512F table has NB=8 for 281-298 bits (85-90d). Our experiments show NB=14-18 + B=70-100K are optimal. But modifying the table hurts because the interpolation code averages NB from bounding rows instead of interpolating. Always use `-siqsNB` and `-siqsB` command-line overrides.
- **Sieve kernel** (`med_sieve_32k_avx2.c:689`): Store-to-load forwarding stall in inner loop (512-bit store then 16-bit loads) is unavoidable on AMD Zen4. No vpscatterb instruction exists.
- **Block Lanczos**: Only 2-8% of total time at 85-89d. VBITS=512 halves BL time but requires custom lanczos.h/lanczos.c modifications.
- **DLP cofactoring**: Batch GCD path disabled (`&& 0`). Online microECM per candidate is faster for DLP.

### Resume Feature
YAFU can save/resume via siqs.dat. Key findings:
- Use `-siqsT <seconds>` for clean shutdown (avoids corrupted save files)
- **DO NOT** use `timeout` command — SIGTERM corrupts siqs.dat
- Resume adds ~25-30s overhead — slower than continuous run for borderline numbers
- Only useful if factoring must span multiple invocations

## Custom SIQS Implementations

### siqs4.c — Custom SIQS with per-block sieve init (agent-10)
- **Status**: Working for 30-35d, partial for 40d. Fails for 45d+.
- **Performance**: 30d: 0.2s (all 5 pass), 35d: 4/5 pass, 40d: 1/5 pass
- **Features**: Per-block sieve initialization (fixes uint8 underflow), AVX512BW sieve scanning, Knuth-Schroeppel multiplier, Gray code self-init, SLP+DLP, merged SLP pairs, GF(2) Gaussian elimination
- **Key bugs found and fixed (agent-10)**:
  1. **Mirror position trivial gcd**: Restricting to x≥0 removed all negative Y values, making all Y on same branch of modular sqrt. Fix: include both positive and negative x values.
  2. **Per-block sieve init**: Using a single init_val causes uint8 underflow for positions near the polynomial vertex. Fix: compute init_val per block based on actual Q(x)/A range.
  3. **A-factor exponents**: Must include a-factor primes in exponent vector (A*Q factorization, not just Q/A).
  4. **LP columns**: Including LPs as extra matrix columns creates degenerate structure. Better: merge SLP pairs by XOR-ing exponents.
- **Remaining issue**: GF(2) Gaussian elimination produces structurally similar null space vectors. Need proper Block Lanczos for 40d+ success.

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

### siqs3.c — Working Custom SIQS (agent-10, DLP, block sieve)
- **Status**: Working, 30-50d confirmed
- **Performance**: 30d: 1.3s, 40d: 2.6s, 50d: 39.6s (40-350x slower than YAFU)
- **Compile**: `gcc -O3 -march=native -mavx512bw -o siqs3 library/siqs3.c -lgmp -lm`
- **Features**: Gray code self-init, Knuth-Schroeppel multiplier, 32KB block sieve, SLP+DLP, Block Lanczos (inline)
- **Bottleneck**: Scalar sieve updates. AVX512 scanning for candidates works but sieve fill is still per-prime scalar stores.

### gnfs_simple.c — Custom GNFS (agent-7, line sieving)
- **Status**: Working, produces relations but too slow for competition
- **Performance**: 50d: 317 rels in 0.27s (needs 6811), 60d: ~41 rels/sec
- **Compile**: `gcc -O3 -march=native -o gnfs_simple library/gnfs_simple.c -lgmp -lm`
- **Features**: Base-m degree-4 poly selection, fast modular poly GCD root finding, line sieving both sides, trial division smoothness
- **Bottleneck**: Line sieve is inherently slow — each (a,b) pair requires checking two norms. Lattice sieving / special-Q needed for competitiveness.
- **Key insight**: NFS line siever for 90d can't compete with YAFU's SIQS. Would need lattice sieve + batch smoothness to approach YAFU speeds.

### Key Insights for Custom SIQS
1. **Multiplier handling**: With kN (k>1), sqrt step systematically produces X ≡ ±Y (mod N). May relate to LP products interacting with multiplier.
2. **Mirror positions**: Q(x) = Q(-x-2b/a), so both sides of sieve give identical Q. Must skip negative Y to avoid trivial SLP combinations.
3. **Scalar sieve bottleneck**: Without AVX512 scattered byte ops, custom SIQS is 40-180x slower than YAFU.
4. **SQUFOF/Rho don't work**: O(N^{1/4}) is dramatically slower than SIQS's sub-exponential for 30+d balanced semiprimes.
5. **The a*g(x) factor**: For SIQS, exponent matrix must track a*g(x), not g(x).
6. **Sieve threshold**: log2(M * sqrt(N)) minus small-prime correction. Too high = few candidates, too low = false positives.

## Custom NFS Lattice Siever

### nfs_siever.c — Working Custom NFS Lattice Siever
- **Status**: Working, produces valid GGNFS-compatible relations
- **Performance**: ~8-9 rels/sec (4000x slower than GGNFS due to non-bucket sieve)
- **Compile**: `gcc -O3 -march=native -mavx512bw -o nfs_siever library/nfs_siever.c -lgmp -lm`
- **Usage**: `./nfs_siever -f <startq> -c <qrange> -o <outfile> -a <jobfile>`
- **Features**:
  - Degree-4 polynomial support (for 85-95d numbers)
  - Gauss/Lagrange 2D lattice reduction
  - Polynomial root finding via GCD method (x^p - x) with Cantor-Zassenhaus splitting
  - Dual-side sieving (algebraic + rational)
  - Single large prime on each side (cofactor < 2^lpb)
- **Key bugs fixed**:
  - Coefficient reduction: `mpz_fdiv_ui` already handles negative numbers correctly — do NOT negate
  - Sieve position: starting index is `(i_off + I/2) % p`, not `i_off - (I/2 % p)`
  - Negative b: negate both (a,b) when b < 0 (standard NFS convention)
- **Optimization needed**: Bucket sieving for large primes, skip small primes in sieve (trial divide instead)

### nfs_factor.sh — NFS Orchestration Script
- Orchestrates msieve poly select + GGNFS sieve + msieve filter/LA/sqrt
- Works but msieve integration needs tuning (msieve.fb format conversion)

## Tools
- `yafu/yafu`: YAFU baseline. Use `siqs(N)` with `-threads 1 -seed 42`.
- `library/nfs_siever.c`: Custom NFS lattice siever (working, produces valid GGNFS-format relations)
- `library/nfs_factor.sh`: NFS orchestration (msieve poly + GGNFS sieve + msieve filter/LA/sqrt)
- `library/siqs2.c`: Custom SIQS (SLP, working 30-60d, 30-80x slower than YAFU)
- `library/siqs_fast.c`: Custom SIQS with DLP (working 30-50d). Compile: `gcc -O3 -march=native -mavx512bw -o siqs_fast library/siqs_fast.c -lgmp -lm`
- `library/mpqs.c`: Custom MPQS (sieve works, sqrt step buggy)
- `library/pollard_rho.c`: Pollard rho (too slow for balanced semiprimes above 30d)
