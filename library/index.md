# Library Index

## Implementations

- **lgsh.c** — MPQS factoring with Knuth multiplier fallback. Self-initializing multi-polynomial quadratic sieve with LP matching. The main working factorer. Compile: `gcc -O3 -o lgsh lgsh.c -lgmp -lm`
- **siqs.c** — SIQS with batch GCD cofactor matching (experimental). Novel approach using Bernstein product-tree batch GCD to find shared factors among cofactors. Compile: `gcc -O3 -o siqs siqs.c -lgmp -lm`
- **hsd.c** — Hierarchical Smooth Decomposition (experimental). Allows large cofactors and uses sort-based + triple LP matching. Compile: `gcc -O3 -o hsd hsd.c -lgmp -lecm -lm -I/usr/local/include -L/usr/local/lib`
- **qs_simple.c** — Minimal single-polynomial QS for debugging/verification.

## Utilities

- **bench.py** — Benchmark a binary across semiprimes. Usage: `python3 bench.py <binary> <approach_name> [sizes...]`. Records results to `experiments.log`.

- **ccd_factor.c** — Cofactor Collision Descent: QS-style sieve + batch cofactor matching + GF(2) linear algebra. Uses aggressive large prime collection and proper dedup matching. Compile: `gcc -O3 -o ccd_factor ccd_factor.c -lgmp -lecm -lm`
- **ecm_factor.c** — Minimal ECM wrapper using GMP-ECM with escalating bounds. Utility for testing/verification. Compile: `gcc -O3 -o ecm_factor ecm_factor.c -lgmp -lecm -lm`

## Experimental (dead ends — documented in expert.md)

- **mld.c** — Multiplicative Lattice Descent: tested LP expansion with alpha=3. Result: 2.16x worse than standard QS.
- **corr_smooth.c** — Smoothness correlation measurement across polynomial pairs. Result: no significant correlation.
- **spectral_walk.c** — Multi-base p-1, multi-polynomial rho, ECM accumulation. All fail for balanced semiprimes.
- **schoof_factor2.c** — Schoof-like torsion factoring via x^N mod (x³+ax+b, N). Result: O(√N) per curve, same as trial division.
- **divpoly_factor.c** — Division polynomial evaluation for factoring. Result: ECM with single-prime stage 1, strictly worse.
- **batch_gcd.c/h** — Bernstein product-tree batch GCD infrastructure. Used by siqs.c.

## Archived/Debug

- **factor_core.h** — Common utilities header (partially used).
- **run_bench.sh** — Shell script for benchmarking across semiprime sizes.
