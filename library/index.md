# Library Index

## Factoring Tools (Best to Worst)

### race_factor.sh - BEST for 65-100 digits
- **Usage**: `bash library/race_factor.sh <N> [timeout]`
- **Method**: Races 6 YAFU SIQS instances with different random seeds (8 threads each), falls back to msieve
- **Why**: YAFU's block Lanczos has a hang bug on ~20% of inputs. Different seeds avoid the hang.
- **Performance**: 65d: 2-5s, 75d: 3-8s, 85d: 10-30s, needs benchmarking for 90-100d

### fast_factor.sh - Simple YAFU+msieve fallback
- **Usage**: `bash library/fast_factor.sh <N> [timeout]`
- **Method**: Single YAFU attempt with timeout, fallback to msieve
- **Best for**: 60-77 digits where YAFU usually doesn't hang

### yafu/yafu - Multi-threaded SIQS
- **Usage**: `echo "siqs(<N>)" | LD_LIBRARY_PATH=/usr/local/lib yafu/yafu -threads 48`
- **Performance**: 10-50x faster than msieve for 65-85d when it doesn't hang
- **Bug**: Hangs in block Lanczos on ~20-30% of inputs at 80+d

### msieve (compiled from yafu/) - Reliable single-threaded SIQS
- **Binary**: `./msieve -q -s <session_file> <N>`
- **Performance**: 30-60d: <2s, 70d: 19s, 75d: ~100s, 80d: ~250s, 85d+: >300s
- **Limitation**: Single-threaded sieving, too slow for 80+ digits

### factor.c / pfactor.c - Parallel ECM
- **Usage**: `./factor <N>` or `./pfactor <N>`
- **Best for**: Unbalanced semiprimes (one factor much smaller)
- **Not competitive** for balanced semiprimes above 60 digits

## Benchmark Scripts
- bench_race.py - Benchmark race_factor.sh
- bench_fast.py - Benchmark fast_factor.sh
- bench_yafu.py - Benchmark YAFU SIQS directly
- bench_msieve.py - Benchmark msieve
