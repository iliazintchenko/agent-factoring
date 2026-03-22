# Agent Expert Knowledge

## Key Finding: YAFU SIQS Dominates for Balanced Semiprimes

YAFU's SIQS (Self-Initializing Quadratic Sieve) is the fastest approach for balanced semiprimes at every tested size from 30-100 digits. This holds against parallel ECM, msieve SIQS (single-threaded), and Pollard's rho.

### Why SIQS beats ECM for balanced semiprimes
- ECM's complexity depends on the **smallest factor size**, not the composite
- For balanced semiprimes, factors are ~N/2 digits — ECM treats these as hard
- SIQS complexity depends on the **composite size** — sub-exponential in N
- The crossover where SIQS beats ECM for balanced semiprimes is below 30 digits

### YAFU Build Configuration (CRITICAL)
The Makefile.gcc ignores config.mk — flags MUST be passed on the command line:
```bash
make -f Makefile.gcc yafu ECM=1 USE_AVX2=1 SKYLAKEX=1 VBITS=256 -j48
```
- `SKYLAKEX=1`: enables AVX512F, AVX512BW, march=skylake-avx512
- `VBITS=256`: 256-bit Block Lanczos vectors (2x faster LA than VBITS=64)
- IFMA/ICELAKE: SLOWER on this CPU — do NOT use
- VBITS=128: also slower than VBITS=256

### Performance Data (YAFU SIQS, 5 parallel, 48-core AMD EPYC 9R45)
| Digits | Worst-case time | Threads | Notes |
|--------|----------------|---------|-------|
| 30-57  | ~1.0-1.7s      | 9 each  | Startup-dominated |
| 58-65  | 1.7-3.3s       | 9 each  | |
| 66-75  | 3.4-20.2s      | 9 each  | |
| 76-85  | 13-75.5s       | 9 each  | |
| 86-92  | 80-226s        | 9 each  | |
| 93-96  | 140-236s       | 48 each | 5*48=240 threads oversubscribed |
| 97-100 | FAILS          | any     | Cannot complete 5 parallel in 300s |

### Thread Scaling
- 9 threads: 75% CPU efficiency (sweet spot for 5 parallel instances)
- 48 threads: only 17% efficiency (terrible, sieving doesn't scale)
- For 93-96d: oversubscription (5*48 threads) works due to shorter total CPU time
- For 97+: total CPU time exceeds what 48 cores can deliver in 300s

### Why 97-100 Fail
- 97d: ~1700 CPU-sec/number. 5 * 1700 = 8500 needed vs 48*300*0.6 = 8640 budget (marginal)
- 100d: 373s even for a single instance — impossible within 300s
- Root cause: L3 cache thrashing (192MB / 5 instances) + poor SIQS thread scaling

## Tools
- `yafu/yafu`: YAFU binary (SKYLAKEX+AVX512+VBITS=256). Use `siqs(N)` command.
- `library/factor.c`, `pfactor.c`: C-based trial div + rho + ECM
- `library/run_yafu.sh`: YAFU wrapper (temp dir, factor extraction)
- `yafu/msieve`: msieve binary (single-threaded SIQS)
