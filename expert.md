# Agent Expert Knowledge

## Algorithm Performance for Balanced Semiprimes

### YAFU multi-threaded SIQS - FASTEST for 60-85 digits
- Compiled as full yafu binary in yafu/yafu
- Usage: `echo "siqs(<N>)" | LD_LIBRARY_PATH=/usr/local/lib yafu/yafu -threads 48`
- Uses multi-threaded sieving (all 48 cores), dramatically faster than msieve
- Performance (worst case across 5 balanced semiprimes):
  - 60d: 1.8s | 65d: 2.6s | 70d: 2.6s | 75d: 4.1s | 77d: 4.3s
- **BUG WARNING**: YAFU hangs on some inputs during linear algebra phase
  - Affected: ~20% of numbers at 67+d (67,68,69,74,78,80-86+ seen)
  - Symptom: "rels found" message repeats infinitely
  - Must use a timeout and fall back to msieve

### msieve SIQS (compiled from yafu/) - Reliable fallback
- Usage: `./msieve -q -s /tmp/session.dat <N>`
- Single-threaded sieving, much slower than YAFU for 70+d
- Performance: 30d: 0.05s | 50d: 0.23s | 60d: 1.87s | 70d: 18.7s | 75d: ~100s
- Reliable: no hang bugs seen

### Best Strategy: YAFU first, msieve fallback
- library/fast_factor.sh: tries YAFU with short timeout, falls back to msieve
- YAFU timeout = digits/3 seconds (scales with problem size)
- For 30-60d: msieve alone is fine (<2s)
- For 65-85d: YAFU when it works (2-12s), msieve fallback (~100s) when it hangs
- For 86-100d: Need investigation - YAFU mostly hangs, msieve very slow

### Parallel ECM (factor.c) - GMP-ECM
- 48-core parallel ECM, good for unbalanced semiprimes
- Fails for balanced semiprimes above ~60 digits
- Not competitive vs SIQS for this benchmark

## Key Performance Insights
1. Multi-threaded SIQS (YAFU) is 10-50x faster than single-threaded (msieve) at 70+d
2. YAFU's block Lanczos has a bug causing hangs on ~20% of inputs
3. msieve is reliable but wastes 47/48 cores during sieving
4. For 80-100d balanced semiprimes, need either fixed YAFU or a custom multi-threaded SIQS

## SKYLAKEX YAFU Build (from /tmp/agent-factoring-1/yafu/yafu)
- Built with: SKYLAKEX=1 USE_AVX2=1 ECM=1 VBITS=256
- Uses AVX-512 SIMD for sieving: 16-lane vectorization, 2x faster than SSE2
- Performance with all 48 threads on single number:
  - 80d: 15s | 85d: ~35s | 90d: ~90s | 91d: 59s | 93d: 182s | 95d: 137s | 100d: 302s
- No hang bugs observed with `siqs()` command
- CRITICAL: must use `siqs()` not `factor()` (latter may hang)
- Must run in own temp directory (working dir conflicts)

## 91-100 Digit Challenge
- 5 numbers must be factored in parallel within 280s wallclock
- With 48 cores split 5 ways (~10 threads each), 91d numbers timeout
- Single number with 48 threads: 91d=59s, 93d=182s, 95d=137s, 100d=302s
- Key bottleneck: thread scaling is sublinear (7x speedup from 1 to 48 threads)
- For 95-100d: even single number barely fits in 300s
- Possible improvements: better polynomial selection, NFS for 95+

## SIQS Architecture (from YAFU source study)
- Sieve: 32KB blocks (L1 cache fit), byte array with log-prime subtraction
- SIMD: AVX2 processes 8 primes/cycle, AVX512 does 16 primes/cycle
- Bucket sieving for large primes with AVX512 gather/scatter
- Linear algebra: Block Lanczos (not Gaussian) - multi-threaded, ~64 dependencies
- Large prime: single, double, and triple LP variations for better relation yield
- Multiplier selection: optimizes kN mod 8 for best sieve efficiency
- Factor base: 175 primes at 30d, scales to 115,500 at 100d

## Reference Implementations
- /tmp/agent-factoring-1/yafu/yafu: SKYLAKEX YAFU - fastest, use this
- yafu/yafu: Our YAFU build - crashes on 80+d (missing SKYLAKEX flags)
- yafu/msieve: msieve SIQS, reliable single-threaded fallback
- cado-nfs/: GNFS implementation for very large numbers (100+d)
