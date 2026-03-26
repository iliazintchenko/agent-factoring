# Can an AI break the L[1/3] barrier in integer factoring?

AI agent attempts to solve one of the deepest open problems in computational mathematics: **can integer factoring be done faster than L[1/3]?** The Number Field Sieve has been the state of the art since 1993 — over 30 years with no improvement to the exponent. We let an AI systematically explore whether a better approach exists.

**Live probability estimate**: P(classical algorithm ever beats L[1/3]) ~ **2-3%**.

## Discovered so far

**530 approaches investigated over ~27 hours of autonomous research across smoothness-based methods, group-order algorithms, algebraic geometry, lattice methods, cohomological invariants, dynamical systems, information theory, and 20+ other mathematical areas.** No scaling breakthrough yet, but has produced arguably the most comprehensive survey of factoring approaches ever compiled. The structural insights about *why* L[1/3] is so hard are genuinely new.

**Identified five structural barriers** explain why every approach fails:
1. **L[1/3] wall** — any sieve with k parameters achieves L[1/(k+1)]; NFS has k=2
2. **Archimedean/non-archimedean gap** — integer size can't be iteratively reduced like polynomial degree
3. **GF(2) bottleneck** — each smooth relation yields exactly 1 bit; proved tight via universality of order-2 elements
4. **CRT opacity** — every computable invariant of Z/NZ either decomposes via CRT or costs as much as factoring
5. **Z-rigidity** — Z has no non-trivial ring endomorphism, which is why Frobenius descent (the key to function field breakthroughs) can't work over the integers

These form a **causal chain**: Z is rigid (#5) -> CRT is opaque (#4) -> search is binarized (#3) -> no metric shortcut (#2) -> L[1/3] wall (#1).

**The key meta-theorem** (informal): any algorithm using (i) algebraic structures over Z/NZ, (ii) smooth number detection, and (iii) GF(2) linear algebra achieves at best L[1/3]. Breaking through requires violating at least one of these three conditions.

## How it works

A master agent (Claude Code) reads `program.md` for its research directive, `expert.md` for accumulated knowledge, and `docs/` for formal barrier analysis. It spawns concurrent investigator subprocesses, each exploring a specific mathematical question. When an investigator finishes, the master harvests its findings into `expert.md` and launches a replacement. This runs continuously.

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

*(So far tested on a single EC2 instance. Extension to multi-machine via shared filesystem is trivial.)*

## What's in the repo

**Knowledge base:**
- `expert.md` — 530 investigated approaches with conclusions, organized by category. The core reference.
- `docs/barrier_synthesis.txt` — Formal analysis of the L[1/3] barrier: five requirements a barrier-breaking structure must satisfy, 20 structural families analyzed, a meta-theorem, and 8 prioritized research directions.
- `docs/open_problems.txt` — 15 formally stated open mathematical problems whose resolution would advance factoring.

**Agent instructions:**
- `program.md` — The master agent's directive: how to launch investigators, what to look for, when to kill, how to harvest.

**Code:**
- `code/factoring_pipeline.c` — Trial division + Pollard rho + p-1 + p+1 + ECM pipeline (GMP)
- `code/nfs_poly.c` — LLL-improved NFS polynomial selection (GMP)
- `code/siqs.c` — Basic quadratic sieve baseline (GMP)
- `code/bench.py` — Benchmark utility across 355 semiprimes (30-100 digits)

## Running it yourself

```bash
# macOS
brew install gmp ecm cmake hwloc

# Ubuntu/Debian
sudo apt install gcc g++ libgmp-dev libecm-dev cmake make libhwloc-dev

# Run locally (launches in tmux)
./run_local.sh

# Or deploy to EC2
./run_ec2.sh --host ec2-user@<ip>
```

Requires a `.env` with `CLAUDE_CODE_OAUTH_TOKEN` (from `claude setup-token`) and `GITHUB_ACCESS_TOKEN`.

## Selected highlights from 530 investigations

- **Shor = the perfectoid tilt**: Shor's quantum algorithm is precisely "efficient computation of the perfectoid tilt" — quantum parallelism evaluates O(N) Frobenius iterates in O(log N) time. This is the cleanest known explanation of the quantum advantage for factoring.

- **Third channel systematically closed**: All 5 candidates for a third NFS smoothness source (second number field, p-adic, function field, adelic, modular forms) are either circular or reduce to MNFS constant improvements. Z provably cannot support a third independent smoothness channel.

- **The cohomological trichotomy**: Every cohomological invariant of Z/NZ falls into exactly one of three categories: (1) computable but CRT-decomposed, (2) encodes factoring but costs >= factoring, or (3) circular. No sweet spot exists.

- **1-bit-per-relation is a theorem**: Proved tight within congruence-of-squares via universality of order-2 elements. GF(k) for any k > 2 is strictly worse.

- **Tower descent axiomatics**: Precisely characterized the 4 requirements for FFS-style descent over Z — all impossible. The deep reason: in function fields, "degree" and "complexity" are independent parameters; in number fields they are the same.

## Known limitations

- **Session length**: The master wraps up after a few hours despite "never stop" instructions.
- **Orphan processes**: Investigator child processes can survive timeouts. Requires manual cleanup.
- **Diminishing returns**: After ~10 batches, the master enters grinding mode — mechanically launching variations rather than thinking strategically.
- **Local minimum**: Investigators' first instinct is to build a quadratic sieve. The master kills them when caught, but it wastes slots.
- **expert.md bloat**: The master appends without deduplicating. Requires periodic manual cleanup.
