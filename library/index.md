# Library Index

## Research Documents

- **barrier_synthesis.txt** — Formal analysis of the L[1/3] barrier: requirements (R1-R5) a barrier-breaking structure must satisfy, 20 structural families analyzed against them, gap analysis, candidate object types, meta-theorem, and 8 prioritized research directions.
- **open_problems.txt** — 15 specific open mathematical problems whose resolution would advance factoring research.

## Implementations

- **nfs_poly.c** — LLL-improved NFS polynomial selection. Achieves 64.5% more smooth values than standard base-m method by finding polynomials with 2-7x smaller max coefficient. Uses GMP.
- **siqs.c** — Basic Quadratic Sieve (not true SIQS). L[1/2] baseline that factors 30-46 digit semiprimes. Log-based sieving, bit-packed Gaussian elimination mod 2, square root extraction. Uses GMP.

## Utilities

- **bench.py** — Benchmark a binary across semiprimes. Usage: `python3 bench.py <binary> <approach_name> [sizes...]`.
