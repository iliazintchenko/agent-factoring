# Library Index

## Best Implementations (ranked by 60d performance)

1. **siqs_bucket.c** — Gray code + DLP→SLP pipeline + bucket sieve. **Best at 65d** (70.7s). `gcc -O3 -march=native -o siqs_bucket library/siqs_bucket.c -lgmp -lm`
2. **spqs_dlp.c** — SPQS + DLP (SQUFOF) + adaptive threshold. **Best at 55d** (3.5s), 70d (~96s single). `gcc -O3 -march=native -o spqs_dlp library/spqs_dlp.c -lgmp -lm`
3. **spqs2.c** — SPQS + bucket sieve. **Best at 70d** (165s worst-of-5). `gcc -O3 -march=native -o spqs2 library/spqs2.c -lgmp -lm`
4. **fast_siqs.c** — Bucket sieve + __int128 TD + Gray code. 18s at 60d. `gcc -O3 -march=native -o fast_siqs library/fast_siqs.c -lgmp -lm`
5. **siqs_opt.c** — Bucket sieve + SLP matching. 19s at 60d, 265s at 70d. `gcc -O3 -march=native -o siqs_opt library/siqs_opt.c -lgmp -lm`
6. **spqs.c** — Multi-polynomial batch sieve (original). 31s at 60d. `gcc -O3 -march=native -o spqs library/spqs.c -lgmp -lm`
7. **hyper_siqs.c** — TLP SIQS. Fast at 60d but LA blows up at 70d. `gcc -O3 -march=native -o hyper_siqs library/hyper_siqs.c -lgmp -lm`

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
