# Library Index

## Best Implementations (ranked by 70d performance)

1. **hyper_siqs2.c** — **BEST OVERALL. 70d=39.6s (6.8x YAFU), 75d=256s (all 5/5).** TLP + incremental root offsets + AVX512 structured GE + optimized sqrt (exponent tracking) + interleaved sieve. `gcc -O3 -march=native -mavx512f -mavx512bw -o hyper_siqs2 library/hyper_siqs2.c -lgmp -lm`
2. **hyper_siqs.c** — 74.9s at 70d. TLP SIQS with DLP graph, structured GE, Gray code. `gcc -O3 -march=native -o hyper_siqs library/hyper_siqs.c -lgmp -lm`
3. **turbo_pp.c** — 80s at 70d. turbo_siqs + interleaved sieve + incremental offsets. `gcc -O3 -march=native -mavx512f -mavx512bw -o turbo_pp library/turbo_pp.c -lgmp -lm`
4. **turbo_siqs.c** — 97.8s at 70d. 48KB blocks, structured GE, DLP graph. `gcc -O3 -march=native -o turbo_siqs library/turbo_siqs.c -lgmp -lm`
5. **siqs_bucket2.c** — Best at 50-60d (0.73s/9.6s). siqs_bucket + AVX512 LA. `gcc -O3 -march=native -mavx512f -o siqs_bucket2 library/siqs_bucket2.c -lgmp -lm`
6. **siqs_bucket.c** — Gray code + DLP→SLP pipeline + bucket sieve. `gcc -O3 -march=native -o siqs_bucket library/siqs_bucket.c -lgmp -lm`
7. **spqs2.c** — SPQS + bucket sieve + AVX512 LA. `gcc -O3 -march=native -mavx512f -o spqs2 library/spqs2.c -lgmp -lm`

## Key LA Implementations
- **structured_gauss.h** — Structured GE with singleton removal + AVX512 GF(2) XOR. Used by hyper_siqs2.
- **block_lanczos.h** — Block Lanczos (included but falls back to structured GE).

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
