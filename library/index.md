# Library Index

## Factoring Tools

- **ecm_factor.c** — ECM factoring using GMP-ECM library. Suyama parametrization with deterministic sigma from seed=42 RNG. Tries increasing B1 bounds from 10K to 100B. Compile: `gcc -O3 -o ecm_factor ecm_factor.c -lgmp -lecm -L/usr/local/lib -I/usr/local/include -lm`

- **mpqs.cpp** — Self-Initializing Quadratic Sieve (SIQS). Multi-polynomial QS using products of fb primes for 'a' values. Includes single large prime variation and GF(2) Gaussian elimination. Seed=42. Compile: `g++ -O3 -o mpqs mpqs.cpp -lgmp -lm`

- **pollard_rho.c** — Pollard's rho with Brent's cycle detection and batch GCD. Good for small factors. Seed=42. Compile: `gcc -O3 -o pollard_rho pollard_rho.c -lgmp -lm`

## Utilities

- **bench.py** — Python benchmarking script. Runs a factoring binary against semiprimes.json, records to experiments.log, prints scaling data. Usage: `python3 bench.py <binary> <approach_name> [sizes...]`

- **run_bench.sh** — Shell benchmarking script (deprecated, use bench.py instead).
