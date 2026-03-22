GOAL:

The goal is twofold:

1. **Become the world's best semiprime factoring expert** — accumulate deep knowledge in `expert.md` and build powerful tools in `library/`. This knowledge compounds: what you learn on one number size should make you better at all future sizes.
2. **Find the fastest possible algorithm for each digit size (30–100)** — by any means necessary. For each size there are 5 semiprimes; your score is the **worst-case (longest) wallclock time** across those 5. That's the number to minimize. You might implement classic algorithms (trial division, Pollard's rho, ECM, quadratic sieve, GNFS), combine them, tune parameters, exploit structure, use SIMD intrinsics, whatever works. The wallclock time is what matters, not the method.

ENVIRONMENT:

- expert.md: the core knowledge base. Contains everything you know about factoring semiprimes — strategies, tricks, rules of thumb, algorithm comparisons, parameter tuning insights, failure modes, and so on. This is the most important file in the project. It should be a live reflection of your current understanding — not a polished report you write at the end. Update it after every significant experiment or discovery. Write down what you tried, what happened, and what it suggests, even if you're not sure yet. Partial insights are valuable. If you've run 3+ experiments without updating expert.md, you're falling behind. And updates are not append-only: as your understanding evolves, restructure, revise, or remove things that turned out to be wrong. The document should always reflect your current best view, not a chronological log.
- library/: your codebase. All code lives here — factoring implementations, analysis tools, utilities, everything. `library/index.md` has a full overview. Write code directly in the library as reusable, well-structured C/C++ source files. Every technique described in `expert.md` should have a corresponding implementation here. Like `expert.md`, the library is a living thing — update, improve, or remove code as your understanding evolves.
- experiments.log: append-only log of every experiment you run. Format per line:
  ```
  [YYYY-MM-DD HH:MM:SS] size: <digits> | approach: <short description> | time: <seconds or FAIL> | notes: <what was learned>
  ```
- best-algos.json: tracks the best (fastest) algorithm found for each number size. Format:
  ```json
  {
    "<digit_count>": {
      "time": <longest wallclock time in seconds among the 5 semiprimes of this size>,
      "run": "<instructions for how to compile and run, pointing to library/ code>"
    }
  }
  ```
  `time` is the worst-case (maximum) wallclock time across the 5 semiprimes of a given digit count — this is the number to minimize. `run` is a self-contained shell command that compiles and runs the relevant library code on a given number, e.g. `"gcc -O2 library/pollard_rho.c -o pollard_rho -lgmp && ./pollard_rho <N>"`. Keep it up to date whenever you find a faster approach.
- semiprimes.json: the frozen test suite. Contains balanced semiprimes from 30 to 100 digits, 5 random ones per size. Format: `{"30": ["num1", "num2", ...], "31": [...], ...}`. **Do not modify this file.**
- yafu/: source code of YAFU (Yet Another Factoring Utility). Strong SIQS implementation. Useful reference for algorithms, data structures, and parameter tuning.
- cado-nfs/: source code of CADO-NFS, the best-in-class NFS implementation. Useful reference for polynomial selection, sieving strategies, and linear algebra.

RULES:

- Max runtime for any single script/process is 300 seconds — this means the entire execution, not each step within it. Put appropriate time guards so this is never exceeded.
- If something gets terminated because of the timeout, make sure to at least have logs to learn from.
- **Single core only.** All factoring code must run on a single CPU core. No multithreading, no multiprocessing, no OpenMP, no pthreads. The benchmark measures single-core performance.
- **Seed is always 42.** Any random seed used anywhere (RNG initialization, ECM curves, etc.) must be 42. No seed hacking — you may not search over seeds to find ones that happen to work well on specific inputs.
- Never stop. Only the user can stop you. Nobody else.
- You can use the browser, read papers, or any other tool at your disposal.
- All factoring code must be in C or C++, compiled and run locally, CPU only. Use GMP (`-lgmp`) for big integer arithmetic and GMP-ECM (`-lecm`) for elliptic curve factoring. You may also use other C/C++ libraries if they help. Use yafu/ and cado-nfs/ as references to learn from.
- Make sure any code you write is fast — use `-O2` or `-O3`, consider SIMD, avoid unnecessary allocations, profile hotspots.
- There are multiple agents working on this repo simultaneously. Pull before starting work and commit+push frequently — every improvement to `best-algos.json`, every update to `expert.md`, every new or changed file in `library/`, every batch of `experiments.log` entries. Other agents depend on your commits to avoid duplicating work and to build on your findings. If you haven't committed in 10 minutes, you're falling behind.

TIPS:

- Start by understanding the benchmark landscape — look at semiprimes.json to see the range (30–100 digits) and sample sizes.
- Different approaches may be most suitable at different sizes. Find the crossover points.
- After any progress or learning, update `expert.md` and commit. Knowledge that isn't written down is knowledge lost.
- If you're not making progress on a size, move on to a different one or a different approach. Revisit later with fresh knowledge.
- When timing, always measure the worst case across all 5 semiprimes of a given size — that's the number that goes into best-algos.json.
- Always run the 5 semiprimes of a given size in parallel. 
