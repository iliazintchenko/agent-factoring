# Agent Expert Knowledge

## Algorithm Performance for Balanced Semiprimes

### msieve SIQS (compiled from yafu/) - FASTEST for all sizes tested
- Compiled via: `make -f Makefile.gcc msieve NO_ZLIB=1` in yafu/
- Usage: `./msieve -q -s /tmp/session.dat <N>`
- Performance (worst case across 5 balanced semiprimes per size):
  - 30 digits: 0.05s | 40 digits: 0.11s | 50 digits: 0.23s
  - 60 digits: 1.87s | 65 digits: 10.8s | 70 digits: 18.7s
  - 75-100 digits: benchmarking in progress
- Dramatically faster than ECM for balanced semiprimes at ALL sizes

### Parallel ECM (factor.c) - 48-core parallel GMP-ECM
- Fork-based: spawns 48 workers with escalating B1 levels
- B1 table: 15d→2K, 20d→11K, 25d→50K, 30d→250K, 35d→1M, 40d→3M, 45d→11M, 50d→43M
- Performance: 30-38d <0.15s, 39-48d <4.2s, 49-50d ~20s
- Fails for balanced semiprimes above ~60 digits (factors >30 digits)
- Useful as fallback or for unbalanced semiprimes

### Key Insights
1. msieve SIQS beats parallel ECM by 10-100x for balanced semiprimes
2. ECM is only competitive when factors are small (<25 digits)
3. For this benchmark (balanced semiprimes), msieve is the clear winner at all sizes
4. Parallel ECM useful as fallback when msieve is unavailable

## SIQS Algorithm Notes (from studying yafu/msieve code)
- Polynomial selection: a = product of s factor base primes, b via CRT
- Sieve uses 32K or 64K blocks with SIMD optimization (SSE4.1/AVX2)
- Factor base size scales with input: ~150 primes at 30d, ~5000 at 80d
- Large prime variation: single and double large primes supported
- Square root extraction: works mod kN, then removes multiplier k
- Linear algebra: Block Lanczos

## Reference Implementations
- yafu/: msieve SIQS source, compiled as ./msieve
- cado-nfs/: GNFS implementation for very large numbers (>100 digits)
