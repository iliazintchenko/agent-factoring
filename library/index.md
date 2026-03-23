# Library Index

## Implementations

- **srg.c** — Smooth Residue Graph factoring. Multi-stage: trial division → Pollard's rho → ECM → sieve with multi-large-prime relations + extended GF(2) linear algebra. Compile: `gcc -O2 -fno-lto -o srg srg.c -lgmp -L/usr/local/lib -lecm -lm`
- **hybrid.c** — Hybrid P-1/P+1/ECM engine using GMP-ECM library with progressive bounds. Good for 30-55+ digit semiprimes. Compile: `gcc -O3 -o hybrid hybrid.c -I/usr/local/include -L/usr/local/lib -lecm -lgmp -lm`
- **ifm.c** — Iterated Frobenius Map (novel but dead end — O(N^{1/4}) with 100x overhead vs rho). Compile: `gcc -O3 -o ifm ifm.c -lgmp -lm`
- **cad.c** — Cascaded Algebraic Descent (WIP). Sieve + aggressive cofactor splitting via rho/ECM. Compile: `gcc -O3 -o cad cad.c -I/usr/local/include -L/usr/local/lib -lecm -lgmp -lm`

## Utilities

- **bench.py** — Benchmark a binary across semiprimes. Usage: `python3 bench.py <binary> <approach_name> [sizes...]`
- **run_scaling.py** — Full scaling test. Records worst-case times to `algo-scaling.json` and results to `experiments.log`.
