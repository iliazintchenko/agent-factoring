# Agent Expert Knowledge

## Current Baseline: YAFU SIQS

YAFU's SIQS with AVX512BW sieve kernels is the current fastest known approach for 30-89 digits. 90+ digits have not been solved within 300s single-core yet.

### Scaling Frontier (agent-9 comprehensive analysis)
- **89d**: All 5 pass with SIQS at load <18. Worst case 294.4s. **Current best-algos.json frontier.**
- **90d**: 3/5 pass with SIQS (load <20). 90d[0,2] need ~305-320s. GNFS ~238s on idle.
- **91d**: ETA ~600s with SIQS. Firmly outside 300s budget.
- **92d**: ETA ~400s at load 12. Might pass on completely idle machine with GNFS.
- **93-95d**: ETA 600-700s. Impossible with any known single-core method in 300s.
- **96-100d**: Estimated 1000-3000s. GNFS crossover from SIQS likely around 95d.
- Each digit adds ~15-20% to SIQS sieve time (consistent with L[1/2,1] scaling).

### Novel Approaches Tested (agent-9)
All produced negative results for beating YAFU:
1. **CFRAC with primorial GCD**: 8.6x slower at 30d, 693x at 50d. Sequential testing fundamentally loses to sieve parallelism.
2. **Batch smoothness via product trees (Bernstein)**: 47x slower. Multi-precision GCD overhead >> byte-level sieve.
3. **Schnorr lattice factoring**: For 90d, needs 50K-dim lattice (LLL infeasible). Small dim gives non-smooth residues.
4. **Fermat sieve (single polynomial Q=(m+x)^2-N)**: Only finds SLP relations. Without polynomial switching, not enough variety for LP pairing.
5. **MQSS (QS with primorial GCD trial division)**: Working but 21x slower than YAFU. Primorial GCD saves ~50% of trial division time for non-smooth candidates, but trial division is only ~2% of total SIQS time.
6. **48KB L1 sieve blocks**: YAFU requires power-of-2 blocks. 48KB is not power-of-2, so can't use standard bit masking.

**Key insight**: The AVX512BW sieve kernel is the critical bottleneck. Any approach that doesn't match YAFU's 64-byte-per-instruction sieve throughput cannot be competitive. Custom implementations with scalar sieves are 20-400x slower.

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

### YAFU SIQS on 90d — Closest Attempt (verified 2026-03-23)
NB=20 B=120K on 90d (all 5 semiprimes):
- 90d[0]: **319.3s** on loaded machine (load ~18). Estimated ~260-280s on idle. Still exceeds 300s single-run.
- 90d[1]: **247.1s** ✓ (baseline on loaded machine, load ~24)
- 90d[2]: **~300s** on idle (previously timed out at 295s boundary). Factored via resume in 12.5s after accumulated siqs.dat from killed runs.
- 90d[3]: **287.6s** ✓ (loaded machine, load ~20)
- 90d[4]: **286.7s** ✓ (loaded machine, load ~20)
- 3/5 pass under 300s. 90d[0] needs ~305-310s on idle, 90d[2] needs ~301-303s. NOT achieved within 300s per-run budget.
- **Source modifications (closnuf, UPM1, VBITS=512) provide no measurable improvement** — A/B test on 90d[1]: 250.7s (modified) vs 247.1s (baseline), i.e. 1.5% slower.
- **Confirmed under low load (load ~10)**: 90d[1]=248.7, 90d[3]=282.7, 90d[4]=281.2. 90d[0] ETA 4s at timeout. 90d[2] sieve COMPLETES (120596/120080 rels) but BL doesn't finish in remaining 0-5s.
- **PGO build**: ~3% SLOWER (5604 vs 5778 rels/sec). Profile-guided opt doesn't improve hand-written AVX512 intrinsics.
- **MFBD=2.0**: Increases sieve rate (6080 vs 5300) but doesn't save enough total time. B=115K+MFBD=2.0 makes things worse (only 1/5 pass vs 3/5 default).
- **-march=native (Zen 5)**: Identical performance to -march=skylake-avx512.
- **CPU pinning (taskset)**: No improvement.
- **Conclusion**: 90d[0,2] require ~301-310s. Need 2-5% total speedup to break through.

### 90d GNFS Direct Sieve Pipeline (agent-5)
Pre-computed polynomial + direct GGNFS sieve + YAFU post-processing:
- **90d[1]**: Sieve 264s (5530 rels/sec) + BL 36s = 300s at load=22. BL reached 94-99% but timed out by 1-2s. On idle: ~253s.
- **90d[0]**: Sieve rate 5140 rels/sec (7% slower than 90d[1] due to polynomial). Needs 284s sieve + 30s BL = 314s. **Cannot fit in 300s even on idle.**
- **Key: start sieve at Q=250000** (not 210000) for 10% higher yield/Q
- GNFS is load-sensitive: 5500 rels/sec at load=22, ~6500 at load<5
- Best polynomial score for 90d[0]: E=4.512e-08 (100s search yielded worse: 4.069e-08)
- **Conclusion**: GNFS can do 90d[1-4] on idle machine but NOT 90d[0]. 90d frontier remains unachievable.

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
| CFRAC (primorial GCD) | Working 30-55d. 8.6x slower at 30d, 693x at 50d vs YAFU. Gap widens exponentially. Confirms sieve >> sequential smooth testing. |
| MQSS (primorial GCD QS) | Relations verify correct but LA step fails. Novel idea: replace trial division with batch GCD against primorial. Needs SLP combining fix. |
| Batch smooth QS (product trees) | NEGATIVE: multi-precision GCD overhead makes batch approach ~50x slower than byte-level sieve. See library/batch_qs.c. |
| Schnorr lattice factoring | NEGATIVE: For 90d, lattice dimension ~50K needed (infeasible for LLL). With small dimension (200), residues are ~10^{50} (never smooth). Theoretically sound but practically useless. |
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
| Standalone GNFS (GGNFS direct + YAFU post) | GGNFS standalone: 5000-6500 rels/sec. YAFU `-R -nc` for filter+BL+sqrt: ~28s. Total on idle: ~253s. **Most promising path for 90d.** |
| YAFU NFS batch overhead | YAFU wraps GGNFS in 2K-Q batches with ~5s overhead each. Effective rate: ~2500 rels/sec vs 5000+ standalone. 2x slowdown from batching. |
| USE_BATCHPOLY build | 11% slower (80d: 52.9s vs 47.0s). Batch polynomial root updates hurt performance. |
| inmem=100 (all in-memory) | No improvement on 80d or 88d vs default inmem cutoff of 70d. |
| CPU pinning (taskset) | No improvement. Single-threaded YAFU doesn't benefit from core affinity. |
| PGO (profile-guided optimization) | ~1-2% improvement. Requires `static __inline` fix in monty.h. Not significant — hot path is hand-written AVX512 intrinsics that GCC can't improve via profiling. |
| monty.h static inline fix | Marginal (~1-2%). Initial "16%" measurement was load-noise artifact. A/B test confirmed identical times. Still useful for enabling PGO builds. |
| Balanced semiprime exploitation | No known algorithm exploits balance. SIQS is factor-structure-agnostic. Fermat/Lehman only help when \|p-q\| < N^(1/3). |
| Batch smoothness (Bernstein product/remainder trees) | Correctly detects smooth numbers but O(M * log²M * bits(primorial)) is slower than byte-level sieve O(M * B/logB). Product tree involves huge multi-precision arithmetic vs L1-cache byte ops. Dead end for QS. |
| Lattice factoring (Schnorr-style) | LLL on log-prime lattice: short vectors don't translate to small products mod N because modular reduction destroys lattice structure. Random combinations of LLL basis find smooth numbers at ~2% rate vs SIQS ~5-10%. Not competitive. |
| Pollard p-1 / Williams p+1 | Works only when p-1 or p+1 is B-smooth. Found factor for 30d[0] with B1=5e7. Not viable for general semiprimes but worth checking each number. |
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

### Empirical Scaling Analysis (agent-8)
Fitting ln(T) = c * sqrt(ln N * ln ln N) + d to YAFU SIQS times on 60-89d:
- **c = 0.8403** (theoretical QS: c ≈ 1.0). YAFU's optimizations achieve sub-theoretical constant.
- Predictions: 90d=309s (matches 3/5 pass), 91d=371s, 92d=445s, 95d=767s, 100d=1868s
- SIQS/GNFS crossover at ~80-85d (GNFS theoretically faster, but YAFU's SIQS constant factor advantage pushes crossover to ~90d in practice)
- To break through the 90d barrier with SIQS: need c < 0.82 (6% improvement in scaling constant) or 20s reduction in constant term (not achievable with parameter tuning)

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

## Novel Approaches Explored (agent-8)

### Batch Smoothness Detection (Bernstein product/remainder trees)
- **Idea**: Replace sieving with batch smooth detection: generate Q(x) values, build product tree, compute primorial mod each Q(x) via remainder tree, extract smooth part via GCD.
- **Implementation**: `batch_qs.c` (vanilla QS) and `batch_siqs.c` (SIQS with multiple polynomials)
- **Result**: NEGATIVE. Correctly detects smooth numbers (verified against trial division) but is 47x slower than YAFU's byte-level sieve. The product/remainder tree involves multi-precision arithmetic on numbers as large as the primorial (~23K bits for B=3600), while sieving uses byte-level add/compare operations in L1 cache.
- **Why it fails**: The asymptotic advantage of batch smoothness (O(M log²M) vs O(M*B)) only holds when B is enormous. For practical QS/SIQS, B < 200K and the constant factors of multi-precision tree operations dominate.
- **When batch IS useful**: Index calculus for DLP in large finite fields, where sieving is impossible.

### Lattice-Based Factoring (Schnorr-style)
- **Idea**: Construct lattice from log-prime vectors, LLL-reduce, use short vectors to generate products p_1^e_1 * ... * p_k^e_k that are close to N^m, so their residue mod N is small and likely smooth.
- **Implementation**: `lattice_factor_batch.c`
- **Result**: NEGATIVE. After LLL reduction, random combinations of basis vectors find smooth numbers at ~2% rate, but this is essentially random smooth number generation. The lattice doesn't help because modular reduction to [0, N) destroys the lattice structure.
- **Key insight**: Schnorr's claimed polynomial-time factoring via lattice methods does not work in practice. The fundamental issue is that lattice short vectors guarantee small Euclidean norm, but not small residues mod N.

### P-1 / P+1 / ECM Scanning
- **Idea**: Check if any semiprime factors have smooth p-1, smooth p+1, or are vulnerable to ECM.
- **Implementation**: `special_factor.c`, `factor_oracle.c`
- **Result**: P-1/P+1 with B1=1e6 finds factors for many 30-56d semiprimes (72 hits out of 250 tested). Zero hits for 57-100d even with B1=5e7. ECM not competitive for balanced semiprimes.
- **Conclusion**: P-1/P+1 are fast but only work for special numbers. For random balanced semiprimes above 56d, factors are too large for p-1 to be B-smooth.

### Key Insights for Custom SIQS
1. **Multiplier handling**: With kN (k>1), sqrt step systematically produces X ≡ ±Y (mod N). May relate to LP products interacting with multiplier.
2. **Mirror positions**: Q(x) = Q(-x-2b/a), so both sides of sieve give identical Q. Must skip negative Y to avoid trivial SLP combinations.
3. **Scalar sieve bottleneck**: Without AVX512 scattered byte ops, custom SIQS is 40-180x slower than YAFU.
4. **SQUFOF/Rho don't work**: O(N^{1/4}) is dramatically slower than SIQS's sub-exponential for 30+d balanced semiprimes.
5. **The a*g(x) factor**: For SIQS, exponent matrix must track a*g(x), not g(x).
6. **Sieve threshold**: log2(M * sqrt(N)) minus small-prime correction. Too high = few candidates, too low = false positives.

## CFRAC Scaling Analysis (agent-9)

### cfrac.c — Working CFRAC with primorial GCD smooth extraction
- **Method**: Continued fraction expansion of sqrt(kN), check each Q_n for B-smoothness
- **Novel optimization**: Extract smooth part via repeated GCD with primorial (product of all FB primes), avoiding per-prime trial division for non-smooth candidates
- **Performance vs YAFU SIQS**:

| Digits | CFRAC worst | YAFU worst | Ratio |
|--------|------------|------------|-------|
| 30     | 0.12s      | 0.014s     | 8.6x  |
| 35     | 0.74s      | 0.016s     | 46x   |
| 40     | 4.92s      | 0.017s     | 289x  |
| 45     | 24.3s      | 0.08s      | 304x  |
| 50     | 83.1s      | 0.12s      | 693x  |

- **Scaling**: CFRAC ratio to YAFU grows exponentially. CFRAC is O(exp(√(2 ln N ln ln N))) while SIQS is O(exp(√(ln N ln ln N))). The factor of √2 in the exponent means CFRAC scales as YAFU^{√2} ≈ YAFU^{1.41}.
- **Why CFRAC loses**: CFRAC tests one Q_n value per step (sequential). SIQS sieve tests millions of Q(x) values simultaneously (parallel via byte array). The sieve amortizes the cost of testing many positions at once.
- **Primorial GCD optimization**: Replaces O(B) trial divisions per candidate with O(log B) GCD iterations. Saves ~50% of trial division time for non-smooth candidates (which are discarded). But since CFRAC's bottleneck is the low smoothness probability (not trial division cost), this optimization doesn't change the scaling class.
- **Conclusion**: No sequential smooth-number method can compete with sieve-based approaches. The sieve's ability to test millions of positions in one pass is the key advantage.

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

### 90d GNFS Verified Factorizations (agent-1)
All 5 90d semiprimes factored via GNFS (pre-computed poly + GGNFS sieve + YAFU post-processing):
- 90d[0]: 315547728763353682629201910311839572981703304736061985377687053517447821670148261179477299 = 331902388528167534291761852803279586498815391 * 950724489096510729028869813711222340321732589
- 90d[1]: 144009255506966132983765710705715101627946917437406468019915946276667323993959592666335633 = 227731436216475209021450721150909180555746243 * 632364410902300352402659467244353314142024731
- 90d[2]: 321739159513091262170255410982298427140058934463728910158638687584271260451816679219749747 = 335245195964718510775053835024328500135758179 * 959712960501159171406048254656254313462882993
- 90d[3]: 149689871933523898487085003856984468870482390037601728900188208806500682729361812646669171 = 183506310264398949279296719368434548097600837 * 815720569597024957780046759609360954162468183
- 90d[4]: 226492946512013642353358853937990356716879227059976248121251466901503301655504057622409363 = 280782675676655379846032037521504370086476931 * 806648579604103213701066786081762925875191473

**Timing**: Resume mode (from partial sieve data): 96-110s. Clean run with pre-computed poly: ~255s at load<5 (estimated), ~320s at load 12 (measured timeout). **The 90d GNFS approach is viable but load-sensitive.**

Key parameters for 90d GNFS:
- Polynomial: degree 4, auto-selected by YAFU (E scores: 4.5-4.9e-08)
- Sieve: GGNFS I=12, lpb=25/25, rlim/alim=350K/840K (YAFU defaults)
- Target: 1.46M relations at ~5000-6500 rels/sec
- Post-processing: ~30s (filter + BL + sqrt)
- Pre-computed polynomials saved in `library/gnfs_polys/90d_{0-4}.poly`

### Batch Smoothness Detection (product tree, Bernstein) — NEGATIVE RESULT
- Implemented in `library/batch_smooth.c`: product tree + remainder tree for batch B-smoothness testing
- **Throughput**: 1.7M tests/sec for 30d numbers (fast GMP product tree)
- **Smoothness rate**: ~0/sec (zero smooth values found in 1.5M tests with B=300)
- **Root cause**: QS polynomial values Q(x) grow linearly in x. For large x, Q(x) ≈ 2*sqrt(N)*x, making them too large for B-smoothness. The sieve's advantage is that it only processes positions where primes contribute, effectively prefiltering to the ~1% of x values most likely to yield smooth Q(x).
- **Conclusion**: Batch smoothness (product trees) cannot replace sieving for QS/NFS. The sieve IS the efficient prefilter. Batch methods are useful for cofactor testing (where candidates are already sieve-identified) but not for finding smooth values from scratch.

### AMD EPYC 9R45 Cache Observation
- **L1D = 48KB** (not 32KB as in Intel). YAFU's sieve uses 32KB blocks (hardcoded in AVX512 assembly as powers of 2).
- A 48KB block is NOT a power of 2, so bit-masking doesn't work. Cannot trivially increase YAFU's block size.
- YAFU already mitigates this via the -siqsNB parameter (processes multiple 32KB blocks before re-initializing).
- A custom sieve with non-power-of-2 blocks would need modular arithmetic instead of bit masks — likely slower.

### batch_qs.c — SIQS with Batch GCD Pre-filter (agent-10)
- **Status**: Working for 30-40d. Fails at 42d+ (timeout).
- **Performance**: 30d: 0.9-7.5s, 35d: 2.3-25.9s, 40d: 15-21s
- **Comparison to YAFU**: 300-1000x slower (scalar sieve vs YAFU's AVX512BW)
- **Compile**: `gcc -O3 -march=native -mavx512bw -o batch_qs library/batch_qs.c -lgmp -lm`
- **Key insight: Batch GCD as pre-filter is not useful at these sizes**
  - The traditional sieve threshold already identifies smooth candidates well
  - Batch GCD (GCD of Q(x) with product of FB primes) passes 99.9% of candidates
  - Trial division overhead is small for FB < 2000 primes
  - Batch GCD becomes useful only for very large factor bases (>10K primes)
  - At 45d with FB=800, batch GCD adds overhead without filtering benefit
- **The fundamental bottleneck is sieve fill speed**: scalar byte subtraction at ~650K ops/block vs YAFU's AVX512BW at ~10M ops/block
- **Features**: Knuth-Schroeppel multiplier, Gray code SIQS, 32KB block sieve, AVX512BW scan, SLP+DLP with hash table, batch GCD smooth pre-filter, GF(2) Gaussian elimination
- **Novel aspect**: Product tree computation of FB product for batch GCD. Theoretically O(n log²n) for n candidates, but constant factor too large for FB < 5000.

### lattice_factor.c — Fermat/Lehman/Lattice Factoring (agent-10)
- **Status**: Not competitive for random balanced semiprimes
- **Performance**: 30d: FAIL (all strategies). Fermat+Lehman only work when |p-q| < N^(1/3).
- **Strategies tried**: Enhanced Fermat with modular constraints, Lehman k≤10^6, 2D lattice smooth congruence
- **Key insight**: For random balanced semiprimes, |p-q| ≈ O(√N), so Fermat/Lehman approaches are doomed. These methods are only useful when factors are deliberately close.
- **Lattice approach**: Gauss-reduced (N, 0), (√(kN), 1) lattice gives small (a,b) with a ≡ b√(kN) mod N, but gcd(a±b√(kN), N) rarely yields a factor because the algebraic relationship doesn't constrain the factors.

### GNFS Pipeline (agent-10)
- **Pre-computed polynomials**: `library/gnfs_polys/90d_{0-4}.job` for all 5 90d semiprimes
  - Polynomial quality: 90d[0] E=4.512e-08, 90d[1] E=4.920e-08, 90d[2] E=4.588e-08, 90d[3] E=4.601e-08, 90d[4] E=4.444e-08
- **gnfs_factor.sh**: End-to-end GNFS pipeline (YAFU poly select + GGNFS sieve + YAFU post)
- **gnfs_pipeline.sh**: Direct GGNFS sieve + YAFU post-processing
- **90d timing**: On idle machine (load <5): poly 0s (pre-computed) + sieve ~225s (6500/sec) + post ~13s = ~238s. Under load >10: sieve rate drops to ~4100/sec, total >300s.

## Special-Q QS (agent-1 finding)
Novel approach: apply NFS-style "special-Q" technique to QS.
- For each large prime Q not in factor base, sieve only positions x where Q | QS_poly(x)
- The quotient QS_poly(x)/Q is Q times smaller → higher smoothness probability
- **Performance**: 57K/5.1K/781/152 rels/sec at 30/40/50/60 digits
- **Compared to YAFU**: ~10x slower with SAME scaling exponent
- **Conclusion**: No scaling improvement. The reduced values are still sub-exponentially large. The special-Q preprocessing adds overhead without changing the fundamental complexity. This approach is essentially equivalent to NFS's special-Q idea, which already has known complexity bounds.

## 91+ Digit Scaling Analysis
- **91d SIQS**: needs ~330s at load 15 (109K of 150K rels in 295s). At load <5: estimated ~285s (close to budget)
- **91d GNFS**: ~350-400s at load 15. Auto-selects same params as 90d (rlim=350K, alim=840K, lpb=25). Poly E=4.323e-08.
- **Theoretical limits**: Both QS and GNFS are sub-exponential. Each digit adds ~35-50% to QS time, ~25-35% to GNFS time. Beyond 92d, neither fits in 300s.
- **Path forward for 91d**: GNFS with better polynomial + lower load might work. SIQS is borderline at ideal conditions.
- **92-100d**: NOT achievable with current algorithms in 300s single-core. Would need a fundamentally faster siever or algorithm.

## Tools
- `yafu/yafu`: YAFU baseline. Use `siqs(N)` with `-threads 1 -seed 42`.
- `library/specialq_siqs.c`: Special-Q QS variant (novel, 10x slower than YAFU, same scaling). `gcc -O3 -march=native -o specialq_siqs library/specialq_siqs.c -lgmp -lm`
- `library/nfs_siever_fast.c`: Custom NFS lattice siever with bucket sieve (8 rels/sec, 625x slower than GGNFS)
- `library/factor_smart.sh`: Smart factoring script for 90d (GNFS with pre-computed poly if available, else SIQS)
- `library/batch_qs.c`: SIQS with batch GCD (working 30-40d, 300-1000x slower than YAFU)
- `library/nfs_siever.c`: Custom NFS lattice siever (working, produces valid GGNFS-format relations)
- `library/gnfs_pipeline.sh`: GNFS pipeline (pre-computed poly + GGNFS sieve + YAFU post). Use for 90d.
- `library/gnfs_factor.sh`: End-to-end GNFS pipeline (poly select + GGNFS sieve + YAFU post)
- `library/siqs_engine.c`: Custom SIQS with AVX512BW scanning (sieve only, no LA/sqrt yet). 180 rels/sec on 40d (340x slower than YAFU).
- `library/batch_smooth.c`: Batch smoothness via product trees. Negative result — cannot replace sieving.
- `library/siqs2.c`: Custom SIQS (SLP, working 30-60d, 30-80x slower than YAFU)
- `library/siqs_fast.c`: Custom SIQS with DLP (working 30-50d). Compile: `gcc -O3 -march=native -mavx512bw -o siqs_fast library/siqs_fast.c -lgmp -lm`
- `library/mpqs.c`: Custom MPQS (sieve works, sqrt step buggy)
- `library/pollard_rho.c`: Pollard rho (too slow for balanced semiprimes above 30d)
