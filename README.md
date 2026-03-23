An autonomous AI agent that teaches itself to become the world's top expert on semiprime factoring. Given a list of balanced semiprimes from 30 to 100 digits, 5 random ones of each size, it discovers novel strategies and iteratively refines its toolbox to reduce runtime on a single CPU core. The runtime for a given size is taken as the longest wallclock time across the set of 5 semiprimes of this size. 

## How it works

1. An AI agent (e.g. Claude Code) reads `program.md` for instructions
2. It reads `expert.md` for accumulated knowledge from prior runs
3. It reads `library/` and `algo-scaling.json` for the tools and best algos it has so far
4. It explores novel factoring approaches, discovers what works, updates everything
5. It pushes its findings to this repo so other agents can build on its findings

```
                              ┌─────────────────┐
                              │   Agent Brain   │
                              │                 │
                              │ expert.md       │
                              │ library/        │
                              │ algo-scaling.json │
                              │ experiments.log │
                              └───────┬─────────┘
                        git pull/push │
                 ┌────────────┬───────┴─────────┬─────────────┐
                 │            │                 │             │
          ┌──────▼──────┐ ┌───▼────────┐ ┌──────▼──────┐     ...
          │   VM  1     │ │   VM  2    │ │   VM  3     │
          │             │ │            │ │             │
          │ ┌─────────┐ │ │ ┌────────┐ │ │ ┌─────────┐ │
          │ │ Agent 1 │ │ │ │Agent 3 │ │ │ │ Agent 5 │ │
          │ │ ┌─┬─┬─┐ │ │ │ │┌─┬─┬─┐ │ │ │ │ ┌─┬─┬─┐ │ │
          │ │ │F│F│F│ │ │ │ ││F│F│F│ │ │ │ │ │F│F│F│ │ │
          │ │ └─┴─┴─┘ │ │ │ │└─┴─┴─┘ │ │ │ │ └─┴─┴─┘ │ │
          │ ├─────────┤ │ │ ├────────┤ │ │ ├─────────┤ │
          │ │ Agent 2 │ │ │ │Agent 4 │ │ │ │ Agent 6 │ │
          │ │ ┌─┬─┬─┐ │ │ │ │┌─┬─┬─┐ │ │ │ │ ┌─┬─┬─┐ │ │
          │ │ │F│F│F│ │ │ │ ││F│F│F│ │ │ │ │ │F│F│F│ │ │
          │ │ └─┴─┴─┘ │ │ │ │└─┴─┴─┘ │ │ │ │ └─┴─┴─┘ │ │
          │ └─────────┘ │ │ └────────┘ │ │ └─────────┘ │
          └─────────────┘ └────────────┘ └─────────────┘

          F = factoring process
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

## Library

All code the agent writes lives in `library/`:

## Known limitations

- **Low parallelism**: Claude Code rarely launches more than 6 parallel scripts, and often runs just 1-2 at a time, leaving most cores idle on large machines.
- **Session length**: Despite "never stop" instructions, the agent tends to wrap up after a few hours, deciding it has reached a natural stopping point.
- **Merge conflicts**: When multiple agents push concurrently, merge conflicts are common. Claude Code resolves them poorly — it tends to blindly concatenate both sides rather than producing a coherent merge, leaving duplicated entries, contradictory statements, and broken files (especially in `expert.md`). Requires periodic manual cleanup.
