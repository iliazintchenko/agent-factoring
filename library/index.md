# Library Index

## Research Documents

- **factoring_theory_summary.txt** — Comprehensive theoretical summary of the factoring landscape, synthesizing ~75 investigated approaches. Covers: three classical families, why L[1/3] is the barrier (α=1/(k+1), archimedean vs non-archimedean), what doesn't work and why, what would be needed.
- **open_problems.txt** — 15 specific open mathematical problems whose resolution would advance factoring research. Covers: algebraic/structural, complexity theory, quantum-classical interface, effective algebraic NT, linear algebra.

## Implementations

- **nfs_poly.c** — LLL-improved NFS polynomial selection. Achieves 64.5% more smooth values than standard base-m method by finding polynomials with 2-7x smaller max coefficient. Uses GMP.
- **siqs.c** — Self-Initializing Quadratic Sieve. Factors 30-46 digit semiprimes (100% success). Log-based sieving, bit-packed Gaussian elimination mod 2, square root extraction. Uses GMP.

## Utilities

- **bench.py** — Benchmark a binary across semiprimes. Usage: `python3 bench.py <binary> <approach_name> [sizes...]`. Records results to `experiments.log`.
