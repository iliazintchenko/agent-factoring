# Library Index

## Best Approach: YAFU SIQS with AVX512BW (Single-Core)

### yafu/yafu - Single-threaded SIQS (AVX512BW)
- **Usage**: `echo "siqs(<N>)" | LD_LIBRARY_PATH=/usr/local/lib /tmp/agent-factoring-4/yafu/yafu -threads 1 -seed 42`
- **Build**: `make -f Makefile.gcc clean && make -f Makefile.gcc yafu NO_ZLIB=1 ECM=1 USE_AVX2=1 SKYLAKEX=1 VBITS=256 -j48`
- **Performance**: 30d: 0.014s, 60d: 0.7s, 70d: 5.8s, 80d: 44s, 85d: 136s, 87d: 188s, 88d: 222s, 89d: 294s
- **Limit**: 89 digits (all 5 pass with NB=18 B=100000). 90d is infeasible single-core.
- **Critical**: Binary MUST be rebuilt with SKYLAKEX=1 to enable AVX512BW sieve kernels.
- **Tip**: Use `-siqsNB 18` for 85+d, `-siqsB 70000-100000` for 85-89d

## Custom Implementations

### siqs2.c - Custom SIQS (WORKING, agent-6 rewrite)
- **Compile**: `gcc -O2 -march=native -o siqs2 library/siqs2.c -lgmp -lm`
- **Usage**: `timeout 295 ./siqs2 <N>`
- **Performance**: 30d: 0.044s, 40d: 0.8s, 50d: 10s
- **Status**: Working. Full SIQS with Gray code self-init, Knuth-Schroeppel multiplier, block sieving, single LP.
- **3-80x slower than YAFU** (no AVX512 sieve, GMP-based trial division)
- **Key algorithms**: Gray code b-enumeration, CRT polynomial generation, dynamic a-factor selection

### mpqs.c - Custom MPQS (PARTIAL)
- **Compile**: `gcc -O3 -march=native library/mpqs.c -o mpqs -lgmp -lm`
- **Status**: Sieve + relation finding works, sqrt step has bug (X^2 != Y^2 mod N)
- **Issue**: Dependencies from Gaussian elimination don't produce valid congruences

### pollard_rho.c - Pollard's rho (Brent variant)
- **Compile**: `gcc -O3 -march=native library/pollard_rho.c -o pollard_rho -lgmp`
- **Performance**: 30d: 1.4s, 35d: 23.6s. Not competitive for balanced semiprimes.

## Other Tools

### msieve (compiled from yafu/) - Standalone single-threaded SIQS
- **Binary**: `yafu/msieve -q -s <session_file> <N>`
- **Performance**: 2.3x slower than YAFU. 80d: ~142s, 89d: ~648s

### factor.c / pfactor.c - ECM-based factoring
- **Usage**: `gcc -O3 library/factor.c -lgmp -lecm -o factor && ./factor <N>`
- **Not competitive** for balanced semiprimes above 40 digits

## Benchmark & Driver Scripts
- **factor_driver.sh**: Optimal factoring driver (auto-selects best YAFU params per digit size)
- **bench_all.py**: Comprehensive benchmark (all sizes 30-89d, all 5 per size in parallel)
- **bench_single.py**: Main benchmark (runs 5 semiprimes per size in parallel)
- **update_best.py**: Updates best-algos.json with benchmark results
- **bench_yafu.py**: YAFU benchmark with parameter control
- **nb_sweep.sh**: NB parameter sweep
