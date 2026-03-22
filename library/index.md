# Library Index

## Best Approach: YAFU SIQS with AVX512BW (Single-Core)

### yafu/yafu - Single-threaded SIQS (AVX512BW)
- **Usage**: `echo "siqs(<N>)" | LD_LIBRARY_PATH=/usr/local/lib /tmp/agent-factoring-3/yafu/yafu -threads 1 -seed 42`
- **Build**: `make -f Makefile.gcc clean && make -f Makefile.gcc yafu NO_ZLIB=1 ECM=1 USE_AVX2=1 SKYLAKEX=1 VBITS=256 -j48`
- **Performance**: 30d: 0.04s, 60d: 0.7s, 70d: 5.8s, 80d: 47s, 85d: 142s, 88d: 252s
- **Limit**: 88 digits (89d needs >300s for hardest numbers)
- **Critical**: Binary MUST be rebuilt with SKYLAKEX=1 to enable AVX512BW sieve kernels. Default build misses them.

## Other Tools

### msieve (compiled from yafu/) - Standalone single-threaded SIQS
- **Binary**: `yafu/msieve -q -s <session_file> <N>`
- **Performance**: 2.3x slower than YAFU. 80d: ~142s, 89d: ~648s
- **Use**: Fallback if YAFU hangs on a specific number

### factor.c / pfactor.c - ECM-based factoring
- **Usage**: `gcc -O3 library/factor.c -lgmp -lecm -o factor && ./factor <N>`
- **Not competitive** for balanced semiprimes above 40 digits

### siqs.c - Custom SIQS (BROKEN)
- Does not produce relations correctly. Needs complete rewrite.

## Benchmark Scripts
- **bench_single.py**: Main benchmark (runs 5 semiprimes per size in parallel, single-threaded each)
- **update_best.py**: Updates best-algos.json with benchmark results
- **run_all_sizes.py**: Run all sizes sequentially (from agent-1)
- **find_optimal_siqsB.py**: Parameter search for siqsB (from agent-1)
