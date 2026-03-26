An autonomous AI agent that teaches itself to become the world's top expert on semiprime factoring. Given a list of balanced semiprimes from 30 to 100 digits, 5 random ones of each size, it discovers novel strategies and iteratively refines its toolbox to reduce runtime on a single CPU core. The runtime for a given size is taken as the longest wallclock time across the set of 5 semiprimes of this size.

## How it works

1. An AI agent (e.g. Claude Code) reads `program.md` for instructions
2. It reads `expert.md` for accumulated knowledge from prior runs
3. It reads `code/` and `docs/` for tools and research documents
4. It explores novel factoring approaches, discovers what works, updates everything
5. It pushes its findings to this repo so other agents can build on its findings

```
              ┌─────────────────────┐
              │    Master Agent     │
              │                     │
              │  program.md         │◄──────────────────────┐
              │  expert.md          │                       │
              │  docs/ + code/      │                       │
              └──────────┬──────────┘                       │
                         │ launches via ssh                 │ 
         ┌───────────────┼───────────────┬────────────┐     │
         │               │               │            │     │
  ┌──────▼──────┐ ┌──────▼──────┐ ┌──────▼──────┐     ▼     │
  │    VM 1     │ │    VM 2     │ │    VM 3     │    ...    │  
  │             │ │             │ │             │           │  reads
  │ Inv 1 Inv 2 │ │ Inv 3 Inv 4 │ │ Inv 5 Inv 6 │           │
  │  ↓     ↓    │ │  ↓     ↓    │ │  ↓     ↓    │           │
  │ C/Py  C/Py  │ │ C/Py  C/Py  │ │ C/Py  C/Py  │           │
  └──────┬──────┘ └──────┬──────┘ └──────┬──────┘           │
         │               │               │                  │
         └───────────────┼───────────────┘                  │
                         │ writes                           │
              ┌──────────▼──────────┐                       │
              │  Shared Filesystem  │                       │
              │  (EFS / NFS mount)  │───────────────────────┘
              │                     │
              │  inv-1/findings.txt │
              │  inv-2/findings.txt │
              │  ...                │
              └─────────────────────┘
```

*(So far been tested with only a single EC2 instance. Extension to multi-machine via shared filesystem for higher-throughput runs is tritival.)*

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

- **Session length**: Despite "never stop" instructions, the master agent tends to wrap up after a few hours, deciding it has reached a natural stopping point.
- **Orphan processes**: When the master is killed or an investigator times out, child processes (Python scripts, compiled programs) can survive and consume CPU indefinitely. Requires manual cleanup (`ps aux` filtered by CPU).
- **Diminishing returns**: The master enters a grinding mode after ~10 batches, mechanically launching variations on exhausted themes rather than thinking strategically. The later investigations produce well-written dead-end proofs but no new leads.
- **Local minimum of known algorithms**: Investigators' first instinct is to build a quadratic sieve. The master kills them when caught, but it wastes investigator slots.
- **expert.md bloat**: The master appends findings without checking for duplicates or cleaning up. Entries from different runs end up in wrong sections. Requires periodic manual cleanup.
- **Merge conflicts**: If multiple masters push to the same repo concurrently, merge conflicts on expert.md are common and resolved poorly (blind concatenation). Single-master setup avoids this.
