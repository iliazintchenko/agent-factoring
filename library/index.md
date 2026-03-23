# Library Index

## Best Implementations (ranked by 70d performance)

1. **turbo_siqs.c** — **Best at 70d** (97.8s worst-of-5). 48KB L1-optimized blocks, BL + structured GE, DLP graph. `gcc -O3 -march=native -o turbo_siqs library/turbo_siqs.c -lgmp -lm`
2. **siqs_bucket.c** — **Best at 65d** (70.7s). Gray code + DLP→SLP pipeline + bucket sieve + Block Lanczos. `gcc -O3 -march=native -o siqs_bucket library/siqs_bucket.c -lgmp -lm`
3. **spqs2.c** — SPQS + bucket sieve. 139s at 70d. `gcc -O3 -march=native -o spqs2 library/spqs2.c -lgmp -lm`
4. **spqs_dlp.c** — SPQS + DLP (SQUFOF) + adaptive threshold. **Best at 55d** (3.5s). `gcc -O3 -march=native -o spqs_dlp library/spqs_dlp.c -lgmp -lm`
5. **hyper_siqs.c** — TLP SIQS. Fast at 60d (15s) but LA bottleneck at 70d. `gcc -O3 -march=native -o hyper_siqs library/hyper_siqs.c -lgmp -lm`
6. **siqs_native.c** — Batch polys (4) + Gray code + SLP. 256s at 70d. `gcc -O3 -march=native -o siqs_native library/siqs_native.c -lgmp -lm`

## Key LA Implementations
- **lanczos.h** — Block Lanczos for sparse GF(2) matrices. 20x faster than Gaussian elimination for 10000x10000 matrices.
- **block_lanczos.h** — Simplified BL (power iteration, doesn't converge properly).

## Novel/Experimental

- **gnfs_work.c** — Complete GNFS with base-m poly selection, line sieve, dual-side smoothness, GF(2) LA, quadratic characters, Hensel lifting algebraic sqrt. Sieve+LA works; algebraic sqrt produces T_j^2=S_j but Y^2≠X^2 mod N (QC/norm tracking issue). `gcc -O3 -march=native -o gnfs_work library/gnfs_work.c -lgmp -lm`
- **nfs_hensel_sqrt.c** — Standalone Hensel lifting module for NFS algebraic sqrt. Root lifting + sqrt lifting + Lagrange interpolation.
- **latsieve_qs.c** — SIQS with ECM cofactor splitting. ECM finds DLP relations but needs proper DLP graph. `gcc -O3 -march=native -o latsieve_qs library/latsieve_qs.c -lgmp -lm -lecm`
- **dlp_opt.c** — DLP with LP columns in GF(2) matrix + Pollard rho cofactor splitting. NEGATIVE RESULT at 65-70d (LP space too large). `gcc -O3 -march=native -o dlp_opt library/dlp_opt.c -lgmp -lm`
- **mfbs_siqs.c** — Multi-Factor-Base Sieve. NEGATIVE RESULT (sieve reduction hurts more than extended TD helps). `gcc -O3 -march=native -o mfbs_siqs library/mfbs_siqs.c -lgmp -lm`
- **ecm_siqs.c** — SIQS with ECM cofactorization. DLP matching was wrong.
- **sqqs.c** — Special-Q QS. NEGATIVE RESULT: overhead exceeds benefit.
- **nfs_factor.c** — NFS skeleton (algebraic sqrt NOT implemented).
- **turbo_siqs.c** — SIQS variant.

## Older SIQS Variants
- **dlp_siqs.c**, **siqs2.c**, **siqs3.c**, **siqs4.c** — Older, slower implementations.
- **siqs_bucket.c (agent-4)** — Initial bucket sieve (precursor to siqs_opt.c).
- **siqs_fast.c**, **siqs_hybrid.c**, **siqs_engine.c** — Various optimizations.

## Other Algorithms
- **pollard_rho.c** — Brent variant. Not competitive above 30d.
- **factor_oracle.c** — Multi-strategy oracle.
- **special_factor.c** — Pollard p-1, Williams p+1, ECM.
- **batch_smooth.c**, **batch_qs.c**, **batch_siqs.c** — Batch smoothness. NEGATIVE RESULT.
- **lattice_factor*.c** — Lattice-based. NEGATIVE RESULT.
- **nfs_siever.c**, **gnfs_simple.c** — NFS implementations (slow).

## Negative Results Summary
1. **DLP LP-column approach** (dlp_opt.c): LP space too large at 65-70d for birthday collisions
2. **MFBS** (mfbs_siqs.c): Reducing sieve FB and adding extended TD hurts more than helps
3. **Batch smoothness** (batch_*.c): Product tree GCD overhead exceeds sieve
4. **Lattice factoring** (lattice_*.c): LLL infeasible for needed dimensions
5. **Special-Q QS** (sqqs.c, specialq_qs.c): Can't collect relations without full sieve
6. **MCFRAC** (mcfrac.c): Sequential CF expansion can't compete with parallel sieve
7. **PairQS** (pairqs.c): Product of two Q(x) is LESS smooth than individual
8. **Batch poly >4** (spqs batch=8,16): No improvement, inner loop dominates
