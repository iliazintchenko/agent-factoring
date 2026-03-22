# Library Index

## NFS (Number Field Sieve) Implementations

### nfs_siever.c — Custom NFS Lattice Siever
- **Compile**: `gcc -O3 -march=native -mavx512bw -o nfs_siever library/nfs_siever.c -lgmp -lm`
- **Usage**: `./nfs_siever -f <startq> -c <qrange> -o <outfile> -a <jobfile>`
- **Performance**: ~8-9 rels/sec (needs bucket sieve optimization)
- **Features**: Degree-4 NFS, lattice reduction, GCD-based root finding, Cantor-Zassenhaus splitting, dual-side sieving
- **Output**: GGNFS-compatible relations (a,b:alg_hex:rat_hex)

### nfs_factor.sh — NFS Orchestration Script
- **Usage**: `bash library/nfs_factor.sh <N> [timeout_secs]`
- **Pipeline**: msieve poly select → GGNFS lattice sieve → msieve filter/LA/sqrt
- **Status**: Sieving works; msieve post-processing integration needs tuning

## Custom SIQS Implementations

### siqs2.c — Working SIQS (Gray code, SLP, block sieve)
- **Compile**: `gcc -O2 -march=native -o siqs2 library/siqs2.c -lgmp -lm`
- **Usage**: `timeout 295 ./siqs2 <N>`
- **Performance**: 30d: 0.044s, 40d: 0.8s, 50d: 10s (30-80x slower than YAFU)
- **Features**: Gray code self-init, Knuth-Schroeppel multiplier, 32KB block sieve, single large prime

### siqs_fast.c — Working SIQS (DLP support)
- **Compile**: `gcc -O3 -march=native -mavx512bw -o siqs_fast library/siqs_fast.c -lgmp -lm`
- **Usage**: `timeout 295 ./siqs_fast <N>`
- **Performance**: 30d: 0.6s, 40d: 0.95s, 50d: 21.7s
- **Features**: AVX512BW sieve scanning, Gray code, SLP+DLP, GF(2) Gaussian elimination

### mpqs.c — MPQS (sqrt step buggy)
- **Compile**: `gcc -O3 -march=native library/mpqs.c -o mpqs -lgmp -lm`
- **Status**: Sieve + relation finding works, sqrt step produces wrong results

### mpqs_custom.c — MPQS variant
- **Compile**: `gcc -O3 -march=native library/mpqs_custom.c -o mpqs_custom -lgmp -lm`

### pollard_rho.c — Pollard's rho (Brent variant)
- **Compile**: `gcc -O3 -march=native library/pollard_rho.c -o pollard_rho -lgmp`
- **Performance**: 30d: 1.4s, 35d: 23.6s. Not competitive for balanced semiprimes.

### siqs4.c — Custom SIQS with per-block sieve init, DLP, merged SLP
- **Compile**: `gcc -O3 -march=native -mavx512bw -o siqs4 library/siqs4.c -lgmp -lm`
- **Usage**: `timeout 295 ./siqs4 <N>`
- **Performance**: 30d: 0.2s (all 5 pass). 35d+: mostly fails (square root issue with merged SLP dependencies)
- **Features**: Per-block sieve init (fixes uint8 underflow), AVX512BW sieve scanning, Knuth-Schroeppel multiplier, Gray code self-init, SLP+DLP, merged SLP pairs, GF(2) Gaussian elimination
- **Key issue**: Merged SLP dependencies produce only trivial gcd for > 30d. Needs proper Block Lanczos to find longer dependencies.

### block_lanczos.h — Block Lanczos linear algebra
- Header-only. Used by siqs implementations for the matrix step.

### gnfs_simple.c — Custom GNFS (line sieving)
- **Compile**: `gcc -O3 -march=native -o gnfs_simple library/gnfs_simple.c -lgmp -lm`
- **Usage**: `timeout 295 ./gnfs_simple <N>`
- **Performance**: 50d: 1172 rels/sec (insufficient for 6811 target). Not competitive with YAFU.
- **Features**: Base-m degree-4 poly, fast modular poly GCD root finding, dual-side line sieve, trial division smoothness
- **Key insight**: Line sieve NFS can't compete. Needs lattice sieve or GGNFS sievers.

### cqs/ — Additional SIQS variant
- `cqs/siqs_gmp.c`: Another SIQS implementation using GMP
