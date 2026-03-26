An autonomous AI agent that teaches itself to become the world's top expert on semiprime factoring. Given a list of balanced semiprimes from 30 to 100 digits, 5 random ones of each size, it discovers novel strategies and iteratively refines its toolbox to reduce runtime on a single CPU core. The runtime for a given size is taken as the longest wallclock time across the set of 5 semiprimes of this size.

## How it works

1. An AI agent (e.g. Claude Code) reads `program.md` for instructions
2. It reads `expert.md` for accumulated knowledge from prior runs
3. It reads `code/` and `docs/` for tools and research documents
4. It explores novel factoring approaches, discovers what works, updates everything
5. It pushes its findings to this repo so other agents can build on its findings

*(Currently runs on a single EC2 instance. Multi-machine via shared filesystem is planned but untested.)*

```
            ┌─────────────────────┐
            │    Master Agent     │
            │                     │
            │  program.md         │◄────────────────┐
            │  expert.md          │                 │
            │  docs/ + code/      │                 │
            └──────────┬──────────┘                 │
                       │ launches via ssh           │ reads
         ┌─────────────┼─────────────┐              │
         │             │             │              │
  ┌──────▼──────┐ ┌────▼────┐ ┌──────▼──────┐      │
  │    VM 1     │ │  VM 2   │ │    VM 3     │ ...  │
  │             │ │         │ │             │      │
  │ Inv 1 Inv 2│ │ Inv 3   │ │ Inv 4 Inv 5│      │
  │  ↓     ↓   │ │  ↓      │ │  ↓     ↓   │      │
  │ C/Py C/Py  │ │ C/Py    │ │ C/Py C/Py  │      │
  └──────┬──────┘ └────┬────┘ └──────┬──────┘      │
         │             │             │              │
         └─────────────┼─────────────┘              │
                       │ writes                     │
            ┌──────────▼──────────┐                 │
            │  Shared Filesystem  │                 │
            │  (EFS / NFS mount)  │─────────────────┘
            │                     │
            │  inv-1/findings.txt │
            │  inv-2/findings.txt │
            │  ...                │
            └─────────────────────┘
```

## Local setup

```bash
# Install dependencies (macOS)
brew install gmp ecm cmake hwloc

# Install dependencies (Ubuntu/Debian)
sudo apt install gcc g++ libgmp-dev libecm-dev cmake make libhwloc-dev python3-flask python3-requests

# Run locally
./run_local.sh
```

## EC2 deploy

```bash
# Launch on EC2 (handles everything: installs deps, clones repo, launches agents in tmux)
./run_ec2.sh --host ec2-user@<ip> --agents <num_agents>
```

Requires a `.env` file with `CLAUDE_CODE_OAUTH_TOKEN` and `GITHUB_ACCESS_TOKEN`.

To generate the OAuth token (uses your Claude Pro/Max subscription, not API billing):
```bash
claude setup-token
```
This creates a long-lived token (valid 1 year) that lets headless EC2 agents bill against your subscription. It does not rotate or invalidate existing tokens.

Multiple agents can work on the same repo simultaneously, communicating through git — each agent pulls the latest expert knowledge, builds on what others found, and pushes its own improvements. No coordination needed beyond `git pull` and `git push`.

## Code

- `code/factoring_pipeline.c` — Trial division + Pollard rho + p-1 + p+1 + ECM pipeline
- `code/nfs_poly.c` — LLL-improved NFS polynomial selection
- `code/siqs.c` — Basic quadratic sieve baseline
- `code/bench.py` — Benchmark utility for semiprimes
- `code/semiprimes.json` — 355 test semiprimes (30-100 digits)

## Docs

- `docs/barrier_synthesis.txt` — Formal analysis of the L[1/3] barrier and requirements for breaking it
- `docs/open_problems.txt` — 15 open mathematical problems relevant to factoring

## Known limitations

- **Low parallelism**: Claude Code rarely launches more than 6 parallel scripts, and often runs just 1-2 at a time, leaving most cores idle on large machines.
- **Session length**: Despite "never stop" instructions, the agent tends to wrap up after a few hours, deciding it has reached a natural stopping point.
- **Merge conflicts**: When multiple agents push concurrently, merge conflicts are common. Claude Code resolves them poorly — it tends to blindly concatenate both sides rather than producing a coherent merge, leaving duplicated entries, contradictory statements, and broken files (especially in `expert.md`). Requires periodic manual cleanup.
- **Local minimum of known algorithms**: The hardest challenge is getting the agent to explore genuinely novel approaches. Despite explicit instructions to not reimplement known algorithms, every agent's first instinct is to build a quadratic sieve — the well-known SOTA is a local minimum baked into the model weights.
- **expert.md bloat**: Agents dump implementation details (parameter tuning, build flags, polynomial formulas) into expert.md rather than keeping it focused on insights and conclusions. Requires periodic manual cleanup to stay useful.
