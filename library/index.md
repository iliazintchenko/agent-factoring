# Library Index

## Factoring Tools

- **mpqs.c** — Well-optimized MPQS. The fastest general-purpose factoring tool for 30-85 digit numbers. Compile: `gcc -O3 -o mpqs mpqs.c -lgmp -lm`

- **mpqs.cpp** — Alternative SIQS implementation. Random A-value selection, CRT b. Works up to ~50 digits. Compile: `g++ -O3 -o siqs_v2 mpqs.cpp -lgmp -lm`

- **ecm_factor.c** — ECM factoring using GMP-ECM. Suyama parametrization, seed=42. Compile: `gcc -O3 -o ecm_factor ecm_factor.c -lgmp -lecm -L/usr/local/lib -I/usr/local/include -lm`

- **cfrac.cpp** — CFRAC using CF expansion of sqrt(N). Single large prime variation. Compile: `g++ -O3 -o cfrac cfrac.cpp -lgmp -lm`

- **sss.cpp** — Smooth Subsum Search (Hittmeir 2023). CRT-based structured candidate generation. Works 30-35 digits; needs mpz_t fix for larger sizes. Compile: `g++ -O3 -o sss sss.cpp -lgmp -lm`

- **pollard_rho.c** — Pollard's rho with Brent's improvement. Compile: `gcc -O3 -o pollard_rho pollard_rho.c -lgmp -lm`

- **siqs.cpp** — Earlier basic QS. Superseded by mpqs.c.

- **dixon_batch.cpp** — Dixon's method with batch smoothness. Does NOT work (dead end).

## Utilities

- **bench.py** — Benchmarking script. Usage: `python3 bench.py <binary> <approach_name> [sizes...]`

- **run_bench.sh** — Shell benchmarking script for single semiprimes.
