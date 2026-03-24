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

## Archived/Debug

- **factor_core.h** — Common utilities header (partially used).
- **verify_mpqs.c** — MPQS relation verification tool.
