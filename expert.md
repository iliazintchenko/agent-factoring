# Agent Expert Knowledge

## Algorithm Performance for Balanced Semiprimes

### YAFU multi-threaded SIQS - Primary tool for all sizes
- Binary: `yafu/yafu` (built with `make -f Makefile.gcc yafu NO_ZLIB=1 ECM=1 USE_AVX2=1`)
- Usage: `echo "siqs(<N>)" | LD_LIBRARY_PATH=/usr/local/lib yafu/yafu -threads 48 -seed 42 [-siqsB <fb_size>]`

**CRITICAL: YAFU has a parameter selection bug for 78-90 digit numbers.**
- Default chooses factor base ~45000 for 80 digits (should be ~15000)
- This causes "matrix is corrupt" error and wastes time re-sieving
- Fix: use `-siqsB <size>` flag with these values:
  - 70d: `-siqsB 8000`
  - 80d: `-siqsB 15000`
  - 85d: `-siqsB 25000`
  - 90d: `-siqsB 40000`
  - 95d: `-siqsB 60000`
  - 100d: `-siqsB 90000`

**YAFU also has a block Lanczos hang bug.**
- On some matrices, the Lanczos iteration diverges (infinite loop)
- Fix applied in `yafu/factor/qs/msieve/lanczos.c`:
  - Added 10-second wall-clock timeout on Lanczos iteration
  - On timeout, bail out and retry with different random seed
  - Seed is rotated using LCG on each retry attempt

### Performance (when running alone on 48 cores)
- 30-50d: 0.1-1.0s (msieve faster for <50d at ~0.05-0.2s)
- 55-65d: 1-3s
- 70-77d: 2-7s
- 78-85d: 4-18s (requires -siqsB tuning)
- 90d: ~23s (requires -siqsB 40000)
- 95d: ~100-190s (from agent-1 benchmarks with 24 threads)

### msieve SIQS - Reliable single-threaded fallback
- Usage: `./msieve -q -s /tmp/msieve_$$.dat <N>`
- Single-threaded sieving, reliable but slow at 78+d
- Performance: 30-60d <2s, 70d ~19s, 75d ~100s, 80d+ exceeds deadline

### Key Insights
1. YAFU multi-threaded SIQS is 10-50x faster than msieve at 70+d
2. Must use `-siqsB` to override YAFU's buggy FB size selection at 78+d
3. Block Lanczos can hang - need timeout+retry with different seed
4. Resource contention from multiple agents severely impacts performance
5. Use workdir approach: `WORKDIR=$(mktemp -d) && cd $WORKDIR && ... && rm -rf $WORKDIR`

## Reference Implementations
- yafu/yafu: Multi-threaded SIQS with AVX2 sieve (primary tool)
- yafu/msieve: Single-threaded SIQS fallback (reliable)
- cado-nfs/: GNFS for very large numbers (100+d)
