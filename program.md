GOAL:

**Discover a novel approach to integer factoring that scales better than existing algorithms.**

The best known classical factoring algorithms (quadratic sieve, number field sieve) are sub-exponential. We want to find approaches that scale polynomially — or at least significantly better than what exists. This would be a major result.

The semiprimes in `semiprimes.json` (30–100 digits, 5 per size) are your testbed. Use them to validate ideas and measure scaling behavior. The benchmark metric is the worst-case wallclock time across the 5 semiprimes of each size — but the real goal is not to optimize times on specific sizes. It's to find algorithms whose time grows fundamentally slower as digit count increases.

yafu/ and cado-nfs/ are available as source code references to read and learn from. **Do not compile or run them.** Do not compile or run any pre-existing factoring tools (YAFU, CADO-NFS, msieve, GGNFS, etc.). All code you run must be your own, written in library/. The YAFU baseline scaling curve is already in best-algos.json — there is nothing to gain from re-measuring it.

Think deeply. Read papers. Explore unconventional approaches — algebraic structure, lattice methods, analytic number theory, novel sieve designs, approaches from other fields. Most ideas will fail. That's fine. Document what you try and what you learn. The path to a breakthrough is through many experiments, not through optimizing existing code.

ENVIRONMENT:

- expert.md: the core knowledge base. Contains everything you know about factoring semiprimes — strategies, tricks, rules of thumb, algorithm comparisons, parameter tuning insights, failure modes, and so on. This is the most important file in the project. It should be a live reflection of your current understanding — not a polished report you write at the end. Update it after every significant experiment or discovery. Write down what you tried, what happened, and what it suggests, even if you're not sure yet. Partial insights are valuable. If you've run 3+ experiments without updating expert.md, you're falling behind. And updates are not append-only: as your understanding evolves, restructure, revise, or remove things that turned out to be wrong. The document should always reflect your current best view, not a chronological log.
- library/: your codebase. All code lives here — factoring implementations, analysis tools, utilities, everything. `library/index.md` has a full overview. Write code directly in the library as reusable, well-structured C/C++ source files. Every technique described in `expert.md` should have a corresponding implementation here. Like `expert.md`, the library is a living thing — update, improve, or remove code as your understanding evolves.
- experiments.log: append-only log of every experiment you run. Format per line:
  ```
  [YYYY-MM-DD HH:MM:SS] size: <digits> | approach: <short description> | time: <seconds or FAIL> | notes: <what was learned>
  ```
- best-algos.json: tracks the scaling behavior of each approach you develop. Format:
  ```json
  {
    "<approach_name>": {
      "<digit_count>": <worst-case wallclock time in seconds among the 5 semiprimes>,
      ...
    }
  }
  ```
  The point is not to minimize any single time — it's to see how each approach's time grows with digit count. An approach that's slower at 50 digits but scales better is more interesting than one that's fast at 50 but hits a wall at 70. Run each approach across as many sizes as feasible and record the full scaling curve.
- semiprimes.json: the frozen test suite. Contains balanced semiprimes from 30 to 100 digits, 5 random ones per size. Format: `{"30": ["num1", "num2", ...], "31": [...], ...}`. **Do not modify this file.**
- yafu/: source code of YAFU (Yet Another Factoring Utility). Strong SIQS implementation. Useful reference for algorithms, data structures, and parameter tuning.
- cado-nfs/: source code of CADO-NFS, the best-in-class NFS implementation. Useful reference for polynomial selection, sieving strategies, and linear algebra. This is specifically meant for parallel factoring algorithm. Not for single-threaded execution.

RULES:

- **Always wrap factoring processes in `timeout 295`.** No exceptions. Max runtime for any single process is 300 seconds. Never run a factoring binary without a timeout guard.
- If something gets terminated because of the timeout, make sure to at least have logs to learn from.
- **Single core only for factoring.** Each factoring process must use a single CPU core. No multithreading, no multiprocessing, no OpenMP, no pthreads. The benchmark measures single-core performance.
- **Maximize parallelism.** Each factoring process uses 1 core, so you should be running many experiments concurrently to fill up the machine — different sizes, different approaches, parameter sweeps, etc. If you're only running 1-2 things at a time, you're wasting the machine. Launch experiments in bulk, not one at a time.
- **Seed is always 42.** Any random seed used anywhere (RNG initialization, ECM curves, etc.) must be 42. No seed hacking — you may not search over seeds to find ones that happen to work well on specific inputs. If you discover that any existing results in `best-algos.json`, `experiments.log`, or `library/` were produced with a different seed or with parallelization, remove them and re-run with the correct settings.
- Never stop. Only the user can stop you. Nobody else.
- You can use the browser, read papers, or any other tool at your disposal.
- All factoring code must be in C or C++, compiled and run locally, CPU only. Use GMP (`-lgmp`) for big integer arithmetic and GMP-ECM (`-lecm`) for elliptic curve factoring. You may also use other C/C++ libraries if they help. Use yafu/ and cado-nfs/ as references to learn from.
- Make sure any code you write is fast — use `-O2` or `-O3`, consider SIMD, avoid unnecessary allocations, profile hotspots.
- There are multiple agents working on this repo simultaneously. Pull before starting work and commit+push frequently — every improvement to `best-algos.json`, every update to `expert.md`, every new or changed file in `library/`, every batch of `experiments.log` entries. Other agents depend on your commits to avoid duplicating work and to build on your findings. If you haven't committed in 10 minutes, you're falling behind. Check what others are working on and pick a different area — don't all pile on the same digit size or approach.
- **Merge conflicts must be resolved properly.** When `git pull` produces conflicts, read both sides carefully and produce a coherent result. Do not blindly keep both sides — that creates duplicated entries, contradictory statements, and broken files. For `expert.md`, ensure the merged result is consistent and non-redundant. For `experiments.log`, keep all entries but remove duplicates. For code files, understand the intent of both changes before merging.

TIPS:

- Read expert.md first to understand what has already been tried.
- Understand the current algorithms deeply before trying to improve on them. Why is QS sub-exponential? What makes NFS faster? Where does the complexity actually come from?
- Test your ideas on small numbers first (30-50 digits) where iteration is fast. If the scaling looks promising, push to larger sizes.
- When timing, always measure the worst case across all 5 semiprimes of a given size — that's the number that goes into best-algos.json.
- After any progress or learning, update `expert.md` and commit. Knowledge that isn't written down is knowledge lost.
- Most ideas will fail. Document why they failed — understanding failure modes is as valuable as finding successes.
