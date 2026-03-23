# Library Index

## Implementations

- **sss.cpp** — Smooth Subsum Search (Hittmeir 2023). CRT-based structured candidate generation. L[1/2] but with better scaling constants than QS. See expert.md for benchmark data.

## Utilities

- **bench.py** — Benchmark a binary across semiprimes. Usage: `python3 bench.py <binary> <approach_name> [sizes...]`
- **run_scaling.py** — Full scaling test. Records worst-case times to `algo-scaling.json` and results to `experiments.log`.
