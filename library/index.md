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

## Experimental/Analysis

- **cofactor_analysis.c** — Analyze cofactor size distribution after B-smooth extraction. Measures how cofactor ratio scales with N size. Key finding: cofactor ratio ~0.76-0.92 depending on B, confirming smoothness bottleneck. Compile: `gcc -O3 -o cofactor_analysis cofactor_analysis.c -lgmp -lm`
- **poly_split.c** — Polynomial splitting factoring attempt (DEAD END). Tests whether polynomial root-finding mod composite N can reveal factors. Probability of success is O(1/√N) per trial. Compile: `gcc -O3 -o poly_split poly_split.c -lgmp -lm`
- **batch_ecm.c** — Custom Montgomery curve ECM with Suyama parametrization. Useful for controlled experiments but slower than GMP-ECM (no stage 2). Compile: `gcc -O3 -o batch_ecm batch_ecm.c -lgmp -lm`
- **ecm_baseline.c** — GMP-ECM wrapper with proper Suyama parametrization for baseline timing. Factors 30-50 digit semiprimes reliably. Compile: `gcc -O3 -o ecm_baseline ecm_baseline.c -lgmp -lecm -lm -I/usr/local/include -L/usr/local/lib`
- **qs_engine.c** — Single-polynomial QS with LP pair matching and GF(2) linear algebra. Working for 30-36 digits, good for testing variations. Compile: `gcc -O3 -o qs_engine qs_engine.c -lgmp -lm`

## Archived/Debug

- **factor_core.h** — Common utilities header (partially used).
- **verify_mpqs.c** — MPQS relation verification tool.
