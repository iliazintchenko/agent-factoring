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

## Reference Implementations
- yafu/yafu: Full YAFU with multi-threaded SIQS (best when it doesn't hang)
- yafu/msieve: msieve SIQS, reliable single-threaded fallback
- cado-nfs/: GNFS implementation for very large numbers (100+d)
