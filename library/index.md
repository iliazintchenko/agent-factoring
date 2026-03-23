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

### Experimental / Novel Approaches
- **siqs_engine.c** — SIQS with AVX512BW scan, DLP. Sieve-only (no LA/sqrt). 180 rels/sec on 40d (340x slower than YAFU). `gcc -O2 -march=native -mavx512bw -o siqs_engine library/siqs_engine.c -lgmp -lm`
- **batch_smooth.c** — Batch B-smoothness via product trees (Bernstein). NEGATIVE RESULT: cannot replace sieving. `gcc -O3 -march=native -o batch_smooth library/batch_smooth.c -lgmp -lm`
- **batch_qs.c** — SIQS with batch GCD pre-filter / Bernstein batch smooth. Working 30-40d. NEGATIVE RESULT: 47x slower than sieve. `gcc -O3 -march=native -mavx512bw -o batch_qs library/batch_qs.c -lgmp -lm`
- **batch_siqs.c** — SIQS with batch smooth detection (product/remainder trees). Working but slower than sieve. `gcc -O3 -march=native -o batch_siqs library/batch_siqs.c -lgmp -lm`
- **batch_debug.c** — Debug tool: compares batch smooth detection vs trial division. Confirms correctness. `gcc -O3 -march=native -o batch_debug library/batch_debug.c -lgmp -lm`
- **lattice_factor.c** — Enhanced Fermat + Lehman + lattice smooth congruence search (agent-7). `gcc -O3 -march=native -o lattice_factor library/lattice_factor.c -lgmp -lm`
- **lattice_factor_v2.c** — Lattice-based factoring combining LLL with QS. `gcc -O3 -march=native -o lattice_factor_v2 library/lattice_factor_v2.c -lgmp -lm`
- **lattice_factor_batch.c** — Schnorr-style LLL on log-prime lattice. NEGATIVE RESULT: 2% smooth rate from random combinations, not better than random search. `gcc -O3 -march=native -o lattice_factor_batch library/lattice_factor_batch.c -lgmp -lm`
- **lattice_siqs.c** — Clean SIQS for scaling measurement. `gcc -O3 -march=native -o lattice_siqs library/lattice_siqs.c -lgmp -lm`
- **special_factor.c** — Pollard p-1, Williams p+1, ECM via GMP-ECM library. Works for numbers with smooth p-1/p+1. No hits above 56d with B1=1e6. `gcc -O3 -march=native -o special_factor library/special_factor.c -lgmp -lecm -lm`
- **factor_oracle.c** — Multi-strategy oracle: trial div → Fermat → rho → p-1 → p+1 → ECM → fail. `gcc -O3 -march=native -o factor_oracle library/factor_oracle.c -lgmp -lecm -lm`

### CFRAC
- **cfrac.c** — **Working.** Continued Fraction factoring with primorial GCD smooth extraction. Factors 30-55d. Novel: replaces trial division with batch GCD against primorial. 8.6x slower than YAFU at 30d, 693x at 50d. Good for scaling analysis. `gcc -O3 -march=native -o cfrac library/cfrac.c -lgmp -lm`

### MQSS (WIP)
- **mqss.c** — Novel QS variant with primorial GCD smooth extraction. Relations verify correctly but GF(2) solver doesn't produce non-trivial factors. Needs SLP combining fix. `gcc -O2 -march=native -o mqss library/mqss.c -lgmp -lm`

### Other
- **pollard_rho.c** — Brent variant. Not competitive above 30d.
- **lattice_factor.c** — Fermat/Lehman/lattice. Not competitive for random balanced semiprimes.
- **block_lanczos.h** — Block Lanczos linear algebra (header-only).

## GNFS Pipeline
- **gnfs_factor.sh** — End-to-end GNFS pipeline (YAFU poly select + GGNFS sieve + YAFU post-processing)
- **gnfs_pipeline.sh** — Direct GGNFS sieve + YAFU post-processing with monitoring
- **gnfs_polys/** — Pre-computed GNFS polynomials for 90d semiprimes (90d_{0-4}.job)
