# Agent Expert Knowledge

## Algorithm Selection by Size
- **30-45 digits**: Parallel ECM (B1=500-11000) handles these in <4s with 48 cores
- **45-60 digits**: ECM with B1=50000-250000, worst case ~20s at 50 digits
- **60-70 digits**: ECM with B1=250000-3000000, gets slow for balanced semiprimes
- **70-100 digits**: Need SIQS - ECM alone is too slow for balanced semiprimes

## ECM Strategy (current implementation)
- Use fork() to spawn NCORES (48) parallel ECM workers
- Each worker runs GMP-ECM with different random sigma
- Escalating B1 levels: 500 → 2000 → 11000 → 50000 → 250000 → 1000000 → 3000000 → ...
- First worker to find a factor kills all others via pipe+SIGKILL
- GMP-ECM B1 table for factor digit counts: 15d:2e3, 20d:11e3, 25d:5e4, 30d:25e4, 35d:1e6, 40d:3e6, 45d:11e6, 50d:43e6

## SIQS (in development)
- Working implementation but needs optimization
- Key insight from YAFU: use SIQS with self-initializing polynomials
- Factor base sizes from YAFU: 50d:~50, 70d:~12000, 100d:~115000 (in bits not digits)
- Need: proper large prime variation, better parameter tuning, possibly multi-threaded sieving

## Performance Observations
- Trial division to 1M: negligible cost
- Pollard's rho (hand-coded): ~1.4s for 15-digit factors - MUCH slower than ECM
- GMP-ECM: highly optimized, dominates rho for all sizes when parallelized
- Parallel speedup nearly linear with 48 cores for ECM

## Reference Code
- yafu: automated factoring tool that picks the best algorithm (trial division → ECM → SIQS → NFS) based on input size. Current SOTA SIQS implementation. Best single tool for numbers up to ~160 digits.
- cado-nfs: NFS-only implementation with the best polynomial selection, sieving, and linear algebra. Used for all recent factoring records.
