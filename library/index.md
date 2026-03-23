# Library Index

## Factoring implementations
- `ecm_factor.c` — ECM factoring using GMP-ECM library. Usage: `./ecm_factor <N> [B1] [curves]`. Compiles with `-I/usr/local/include -L/usr/local/lib -lecm -lgmp -lm`.
- `mpqs.c` — Multiple Polynomial Quadratic Sieve. Usage: `./mpqs <N>`. Compiles with `-lgmp -lm`. Full MPQS with CRT-based polynomial generation and large prime variation.
- `siqs.c` — Basic QS (window-shifting). Not competitive — exists for reference only. MPQS is strictly better.

## Benchmarking
- `bench_ecm.sh` — Run ECM across all semiprime sizes, parallel execution.
- `bench_mpqs.sh` — Run MPQS across semiprime sizes 30-70, parallel execution.
- `run_bench.sh` — Generic benchmark runner for any factoring binary.
