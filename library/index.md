# Library Index

## Implementations

- **srg.c** — Smooth Residue Graph factoring. Multi-stage: trial division → Pollard's rho → ECM → sieve with multi-large-prime relations + extended GF(2) linear algebra. Compile: `gcc -O2 -fno-lto -o srg srg.c -lgmp -L/usr/local/lib -lecm -lm`

## Utilities

- **bench.py** — Benchmark a binary across semiprimes. Usage: `python3 bench.py <binary> <approach_name> [sizes...]`
- **run_scaling.py** — Full scaling test. Records worst-case times to `algo-scaling.json` and results to `experiments.log`.
