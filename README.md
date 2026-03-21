An autonomous AI agent that teaches itself to become the world's top expert on semiprime factoring. Given a list of balanced semiprimes from 30 to 100 digits, 5 random ones of each size, it discovers novel strategies and iteratively refines its toolbox to reduce runtime on a single CPU core. The runtime for a given size is taken as the longest wallclock time across the set of 5 semiprimes of this size. 

## How it works

1. An AI agent (e.g. Claude Code) reads `program.md` for instructions
2. It reads `expert.md` for accumulated knowledge from prior runs
3. It reads `library/` and `best-algos.json` for the tools and best algos it has so far
4. It explores novel factoring approaches, discovers what works, updates everything
5. It pushes its findings to this repo so other agents can build on its findings

```
                              ┌─────────────────┐
                              │   Agent Brain   │
                              │                 │
                              │ expert.md       │
                              │ library/        │
                              │ best-algos.json │
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
brew install gmp cmake

# Install dependencies (Ubuntu/Debian)
sudo apt install gcc g++ libgmp-dev cmake make

# Clone reference factoring tools (vendored without git history)
git clone --depth 1 https://github.com/bbuhrow/yafu.git && rm -rf yafu/.git
git clone --depth 1 https://gitlab.inria.fr/cado-nfs/cado-nfs.git && rm -rf cado-nfs/.git

# Run locally
./run_local.sh
```

## EC2 deploy

```bash
# Launch on EC2 (handles everything: installs deps, clones repo, launches agents in tmux)
./run_ec2.sh --host ec2-user@<ip> --agents 3
```

Requires a `.env` file with `CLAUDE_CODE_API_KEY` and `GITHUB_ACCESS_TOKEN`. The API key is auto-refreshed from your local Claude Code login on each deploy.

Multiple agents can work on the same repo simultaneously, communicating through git — each agent pulls the latest expert knowledge, builds on what others found, and pushes its own improvements. No coordination needed beyond `git pull` and `git push`.

## Library

All code the agent writes lives in `library/`:

## Known limitations

- **Low parallelism**: Claude Code rarely launches more than 6 parallel scripts, and often runs just 1-2 at a time, leaving most cores idle on large machines.
- **Session length**: Despite "never stop" instructions, the agent tends to wrap up after a few hours, deciding it has reached a natural stopping point.

The agent maintains `expert.md` as a living knowledge base and improves the library as it learns.
