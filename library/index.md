# Library Index

## Custom Implementations

### SIQS
- **siqs2.c** — Working. Gray code, SLP, 32KB block sieve. 30-60d. `gcc -O2 -march=native -o siqs2 library/siqs2.c -lgmp -lm`
- **siqs3.c** — Working. DLP, inline Block Lanczos. 30-55d. `gcc -O3 -march=native -mavx512bw -o siqs3 library/siqs3.c -lgmp -lm`
- **siqs4.c** — Working 30-50d, sqrt fails 35d+. Per-block sieve init, DLP. `gcc -O3 -march=native -mavx512bw -o siqs4 library/siqs4.c -lgmp -lm`
- **siqs_fast.c** — Working 30-50d. DLP, AVX512BW scanning. `gcc -O3 -march=native -mavx512bw -o siqs_fast library/siqs_fast.c -lgmp -lm`
- **siqs.c**, **siqs_avx.c**, **siqs_optimized.c**, **cqs/siqs_gmp.c** — Other variants.

### MPQS
- **mpqs.c** — Sieve works, sqrt step buggy. `gcc -O3 -march=native library/mpqs.c -o mpqs -lgmp -lm`
- **mpqs_custom.c** — Variant.

### NFS
- **nfs_siever.c** — Custom lattice siever, GGNFS-compatible output. ~8-9 rels/sec. `gcc -O3 -march=native -mavx512bw -o nfs_siever library/nfs_siever.c -lgmp -lm`
- **gnfs_simple.c** — Line sieve NFS. Working but slow. `gcc -O3 -march=native -o gnfs_simple library/gnfs_simple.c -lgmp -lm`
- **gnfs_factor.c** — NFS orchestration in C.

### Other
- **pollard_rho.c** — Brent variant. Not competitive above 30d.
- **block_lanczos.h** — Block Lanczos linear algebra (header-only).
