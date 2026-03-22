# Agent Expert Knowledge

## Algorithm Performance for Balanced Semiprimes

### msieve SIQS (compiled from yafu/) - FASTEST for all sizes 30-78
- Compiled via: `make -f Makefile.gcc msieve NO_ZLIB=1` in yafu/
- Usage: `./msieve -q -s /tmp/session.dat <N>` (unique session file for parallel runs!)
- Performance (worst case across 5 balanced semiprimes per size):
  - 30d: 0.05s | 40d: 0.11s | 50d: 0.23s | 60d: 1.87s
  - 65d: 10.8s | 70d: 18.7s | 75d: 55.7s | 76d: 102s
  - 77d: 155s | 78d: 154s
  - 79-80d: EXCEEDS 280s deadline (single 80d number: 289s)
- SIQS sieving is SINGLE-THREADED. -t flag only helps linear algebra.
- Threading actually SLOWS DOWN SIQS (cache contention) for <85 digit numbers
- msieve has AVX2 sieve variants in yafu but build issues prevented use

### Parallel ECM (factor.c) - 48-core parallel GMP-ECM
- Fork-based: spawns 48 workers with escalating B1 levels
- Performance: 30-38d <0.15s, 39-48d <4.2s, 49-50d ~20s
- Fails for balanced semiprimes above ~60 digits (factors >30 digits)
- Useful as concurrent racer alongside msieve for 60-80d numbers

### Performance Tiers
1. **30-65 digits**: msieve SIQS, well under 11s worst case
2. **66-78 digits**: msieve SIQS, under 155s worst case, within deadline
3. **79-85 digits**: PROBLEM ZONE - msieve SIQS exceeds 280s deadline
   - Need: hybrid approach (msieve + parallel ECM racing)
   - Or: optimized SIQS with AVX2 sieve
4. **85-100 digits**: HARD - SIQS gets exponentially slower
   - msieve GNFS (-n flag) is slower than SIQS for 80d, likely better at 95+
   - 90-digit test still running after 15+ minutes
   - May need cado-nfs for these sizes

### Critical Insights
- msieve SIQS is single-threaded for sieving (the bottleneck phase)
- With 48 cores available, 43 cores are wasted when running msieve SIQS
- The hybrid approach (race msieve + parallel ECM) uses all cores
- GNFS crossover with SIQS is typically 95-110 digits
- For 80-90 digit balanced semiprimes, SIQS is still the best algorithm

## SIQS Algorithm Notes (from studying yafu/msieve code)
- Polynomial selection: a = product of s factor base primes, b via CRT
- Sieve uses 32K or 64K blocks with SIMD optimization (SSE2 default, AVX2 available)
- Factor base size scales with input: ~150 at 30d, ~5000 at 80d, ~15000 at 95d
- Large prime variation: single and double large primes supported
- Square root extraction: works mod kN, then removes multiplier k afterward
- Linear algebra: Block Lanczos (multi-threaded)

## Reference Implementations
- yafu/: msieve SIQS + GNFS source, compiled as ./msieve
- cado-nfs/: Best-in-class GNFS implementation for very large numbers
