# Library Index

## Factoring Tools

### msieve (compiled from yafu/)
- **Binary**: `./msieve -q -s <session_file> <N>`
- **Method**: Auto-selects trial division → ECM → SIQS → GNFS based on input size
- **Status**: Working. Primary factoring tool for all sizes.
- **Performance**: See best-algos.json for per-size benchmarks

### factor.c - ECM-based factoring
- **Compile**: `gcc -O2 -o factor library/factor.c -lgmp -lecm -lm`
- **Usage**: `./factor <N> [deadline_seconds]`
- **Method**: Trial division + Pollard rho + GMP-ECM
- **Best for**: 30-60 digit balanced semiprimes (factors up to ~30 digits)
- **Status**: Working

### factorize.sh - Best-method wrapper
- **Usage**: `./library/factorize.sh <N> [deadline]`
- **Method**: Tries msieve first, falls back to ECM
- **Status**: Working

## Benchmark Tools

### run_all.sh - Full benchmark suite
- **Usage**: `bash library/run_all.sh [tool] [start_digits] [end_digits] [deadline]`
- **Status**: Working

## In Development

### siqs.c - Custom SIQS implementation
- **Status**: Buggy (square root extraction fails). Using msieve instead.
- **Issue**: Exponent vectors for combined partial relations produce trivial null vectors.

### qs.c - Basic QS (single polynomial)
- **Status**: Too slow for practical use. Single polynomial yields insufficient smooth values.
