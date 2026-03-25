# Library Index

## Research Documents

- **insights.txt** — The definitive top 5 theoretical insights from ~330 investigated approaches: (1) L[1/3] structural barrier from α=1/(k+1), (2) archimedean/non-archimedean size gap, (3) 1-bit-per-relation GF(2) bottleneck, (4) CRT opacity, (5) Z-rigidity (no endomorphism). Includes synthesis showing how all five interrelate.
- **open_problems.txt** — 15 specific open mathematical problems whose resolution would advance factoring research.

## Implementations

- **nfs_poly.c** — LLL-improved NFS polynomial selection. Achieves 64.5% more smooth values than standard base-m method by finding polynomials with 2-7x smaller max coefficient. Uses GMP.
- **siqs.c** — Basic Quadratic Sieve (not true SIQS). L[1/2] baseline that factors 30-46 digit semiprimes. Log-based sieving, bit-packed Gaussian elimination mod 2, square root extraction. Uses GMP.

## Utilities

- **bench.py** — Benchmark a binary across semiprimes. Usage: `python3 bench.py <binary> <approach_name> [sizes...]`.
