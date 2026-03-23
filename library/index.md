# Library Index

## Factoring Tools

- **mpqs.c** — Well-optimized MPQS. The fastest general-purpose factoring tool for 30-85 digit numbers. Compile: `gcc -O3 -o mpqs mpqs.c -lgmp -lm`

- **mpqs.cpp** — Self-Initializing Quadratic Sieve (SIQS). Multi-polynomial QS using products of fb primes for 'a' values. Includes single large prime variation and GF(2) Gaussian elimination. Seed=42. Compile: `g++ -O3 -o mpqs mpqs.cpp -lgmp -lm`

- **ecm_factor.c** — ECM factoring using GMP-ECM. Suyama parametrization, seed=42. Compile: `gcc -O3 -o ecm_factor ecm_factor.c -lgmp -lecm -L/usr/local/lib -I/usr/local/include -lm`

- **ecm_factor.cpp** — ECM factoring using GMP-ECM library (C++ version). Usage: `./ecm_factor <number>`. Outputs `factor1 factor2` or returns 1 on failure. Seed=42, single core, internal 290s timeout.

- **siqs_factor.cpp** — Self-Initializing Quadratic Sieve. Usage: `./siqs_factor <number>`. Features: multi-polynomial SIQS, single large prime variation, GF(2) Gaussian elimination. Seed=42, single core.

- **cfrac.cpp** — CFRAC using CF expansion of sqrt(N). Single large prime variation. Compile: `g++ -O3 -o cfrac cfrac.cpp -lgmp -lm`

- **sss.cpp** — Smooth Subsum Search (Hittmeir 2023). CRT-based structured candidate generation. Works 30-35 digits; needs mpz_t fix for larger sizes. Compile: `g++ -O3 -o sss sss.cpp -lgmp -lm`

- **pollard_rho.c** — Pollard's rho with Brent's improvement and batch GCD. Good for small factors. Seed=42. Compile: `gcc -O3 -o pollard_rho pollard_rho.c -lgmp -lm`

- **siqs.cpp** — Earlier basic QS. Superseded by mpqs.c.

- **dixon_batch.cpp** — Dixon's method with batch smoothness. Does NOT work (dead end).

## Utilities

- **bench.py** — Benchmarking script. Usage: `python3 bench.py <binary> <approach_name> [sizes...]`

- **run_scaling.py** — Runs a factoring binary across all semiprimes of each digit size, records worst-case times to `algo-scaling.json` and individual results to `experiments.log`. Usage: `python3 run_scaling.py <binary> <approach_name> [min_digits] [max_digits]`

- **run_bench.sh** — Shell benchmarking script (deprecated, use bench.py instead).
