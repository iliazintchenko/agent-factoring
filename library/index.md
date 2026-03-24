# Library Index

## Research Documents

- **top5_insights.txt** — The definitive top 5 theoretical insights from ~142 investigated approaches: (1) L[1/3] structural barrier from α=1/(k+1), (2) archimedean/non-archimedean size gap, (3) 1-bit-per-relation GF(2) bottleneck, (4) CRT opacity, (5) Z-rigidity (no endomorphism). Includes synthesis showing how all five interrelate.
- **factoring_theory_summary.txt** — Comprehensive theoretical summary synthesizing ~75 approaches. Covers: three classical families, why L[1/3] is the barrier, what doesn't work and why, what would be needed.
- **open_problems.txt** — 15 specific open mathematical problems whose resolution would advance factoring research.

## Implementations

- **nfs_poly.c** — LLL-improved NFS polynomial selection. Achieves 64.5% more smooth values than standard base-m method by finding polynomials with 2-7x smaller max coefficient. Uses GMP.
- **siqs.c** — Self-Initializing Quadratic Sieve. Factors 30-46 digit semiprimes (100% success). Log-based sieving, bit-packed Gaussian elimination mod 2, square root extraction. Uses GMP.

## Utilities

- **bench.py** — Benchmark a binary across semiprimes. Usage: `python3 bench.py <binary> <approach_name> [sizes...]`. Records results to `experiments.log`.
