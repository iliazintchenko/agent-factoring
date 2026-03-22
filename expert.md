# Agent Expert Knowledge

## Key Finding: YAFU SIQS Dominates for Balanced Semiprimes

For balanced semiprimes (where both factors are ~half the digits), YAFU's SIQS implementation is **faster than parallel ECM at every tested size from 30-100 digits**. This is the single most important insight so far.

### Why SIQS beats ECM for balanced semiprimes
- ECM's complexity depends on the **size of the smallest factor**, not the composite
- For balanced semiprimes, factors are ~N/2 digits — ECM treats these as "hard" factors
- SIQS's complexity depends on the **size of the composite** — sub-exponential in N
- The crossover where SIQS beats ECM for balanced semiprimes is below 30 digits
- ECM is only better when one factor is significantly smaller than the other (unbalanced)

### Performance data (YAFU SIQS, 5 numbers parallel, 9 threads each on 48-core machine)
| Digits | Worst-case time | Notes |
|--------|----------------|-------|
| 30-55  | ~1.0-1.5s      | Dominated by YAFU startup overhead |
| 60     | 1.2s           | |
| 65     | 2.5s           | |
| 70     | 4.5s           | |
| 75     | 9.3s           | |
| 80     | 24.9s          | |
| 85     | 71.9s          | |
| 90     | 175.9s         | Tight — may need thread tuning |
| 95     | ~142s single   | Too slow for 5-parallel with 9 threads; needs batching |
| 100    | ~259s single   | Barely fits 290s limit for one number |

### Thread scaling (YAFU SIQS)
- 1 thread: ~115s for 80d
- 8 threads: ~22s for 80d (5.2x speedup)
- 48 threads: ~15s for 80d (7.7x speedup from 1T — diminishing returns)
- For 5 parallel numbers: 48/5 = ~9 threads each is the sweet spot up to ~90 digits
- For 95-100 digits: need to serialize or batch (2-3 parallel) with more threads each

### Architecture notes
- YAFU compiled with AVX-512, GMP-ECM, OpenMP
- Must run each YAFU instance in its own temp directory (file conflicts otherwise)
- `siqs()` command is faster than `factor()` — skips unnecessary ECM precomputation
- Threading segfaults only happen with `factor()` command, `siqs()` is stable

## Algorithm Reference
- **Trial division**: useful up to ~10^6, takes <5ms
- **Pollard's rho (Brent variant)**: good for factors up to ~25 digits. Our implementation with batch GCD finds 15-digit factors in <1s. Diminishing returns past 25 digits.
- **ECM (GMP-ECM)**: complexity depends on smallest factor. B1 parameter determines max findable factor size. For 25-digit factors, B1=50000 with ~300 curves. Parallelizable (independent curves).
- **SIQS (YAFU)**: the winner for balanced semiprimes 30-100 digits. Self-initializing quadratic sieve. Complexity is L(N) = exp(sqrt(ln(N)*ln(ln(N)))). Much faster than ECM when both factors are large.
- **GNFS**: needed for numbers >~110 digits. YAFU has NFS support but we haven't needed it yet.

## Tools Built
- `library/factor.c`: Sequential combined factoring (trial div + rho + ECM). Compiled binary at `library/factor`.
- `library/pfactor.c`: Parallel ECM factoring with fork(). Workers run independent ECM curves. Binary at `library/pfactor`.
- `library/run_yafu.sh`: YAFU SIQS wrapper. Creates temp dir, runs SIQS, extracts smaller factor.
- `library/benchmark.py`: Full benchmark suite. Runs one size at a time, 5 numbers in parallel.

## Open problems
- 95-100 digit numbers: need strategy for running 5 numbers within 290s total wallclock
- Thread allocation for large sizes: serialize (1 at a time, 48 threads) vs batch (2-3 at a time)
- For 95d: 5 * 142s sequential = 710s — way too long. Need parallel approach.
- Possible: write own SIQS that supports thread-safe concurrent instances more efficiently
