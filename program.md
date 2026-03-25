GOAL:

**Discover a novel approach to integer factoring that scales better than existing algorithms.**

The best known classical factoring algorithms are sub-exponential: quadratic sieve achieves L[1/2], number field sieve achieves L[1/3]. NFS has been the state of the art since 1993 — over 30 years with no improvement in the scaling exponent. We want to do better.

The ultimate goal is polynomial scaling. But any improvement in the L exponent — even L[1/4] — would be a significant result that nobody has achieved classically.

A key observation: Shor's quantum factoring algorithm achieves polynomial time by reducing factoring to period-finding in the multiplicative group mod N. The quantum Fourier transform finds this period efficiently — but *why*? What is it about the algebraic structure of Z_N* that makes factoring reducible to period-finding? Is there a way to exploit that same structure classically? Nobody has proven there isn't.

**Any approach that is L[1/2] or L[1/3] is not progress.** If you can show that your approach has the same asymptotic scaling as QS or NFS — even if it has a different name or different constants — stop working on it and move on to a different idea. Optimizing constants within L[1/2] is not the goal. The goal is to improve the exponent.

HOW YOU WORK:

You manage this research project by launching investigator subprocesses — separate Claude Code instances that each explore a specific question. You do not write factoring code yourself. Instead, you:

1. **Think** about what directions to explore, informed by expert.md (what's been tried) and the current theoretical landscape
2. **Launch investigators** to test specific hypotheses or explore specific questions
3. **Monitor** their progress by reading their logs
4. **Harvest** useful findings when they finish, updating expert.md and library/
5. **Repeat** — always keeping 5 investigators running

LAUNCHING INVESTIGATORS:

For each investigator, create a working directory, copy in whatever files they need, and launch:

```bash
mkdir -p /tmp/inv-N
cp semiprimes.json /tmp/inv-N/   # plus any other files they need
TIMEOUT_CMD=$(command -v timeout || command -v gtimeout)
cd /tmp/inv-N && nohup $TIMEOUT_CMD --kill-after=30 900 \
  claude -p '<task description>' \
  --dangerously-skip-permissions --verbose --output-format text \
  > log.txt 2>&1 &
echo $! > pid.txt
```

Note: the 15-minute timeout (900s) prevents investigators from running forever. `--kill-after=30` sends SIGKILL if SIGTERM doesn't work. The timeout command kills the entire process group, including any child processes (compiled programs, Python scripts). Use `timeout` on Linux, `gtimeout` on macOS (from `brew install coreutils`).

The task description should be specific enough that the investigator knows exactly what to do, but open enough that they can discover unexpected things. Include:
- What question to investigate
- What tools are available (GMP, GMP-ECM, etc.)
- What concrete first step to take
- What to write to `findings.txt` when done
- A time expectation ("this should take about 10-15 minutes")

Example good tasks:
- "Investigate whether the Deuring correspondence between supersingular elliptic curves and maximal orders in quaternion algebras ramified at p can leak information about p when working over Z/NZ. Start by implementing point-counting on random curves over Z/NZ for 20-digit semiprimes. Use GMP. Write findings to findings.txt."
- "Test whether the distribution of element orders in Z_N* for N=pq differs detectably from the distribution for N=p'q' where p'q'≈pq. Sample 10000 random elements, compute order mod small primes (2,3,5,7,11), compare distributions for 5 different 30-digit semiprimes. Use GMP. Write findings to findings.txt."
- "Read the paper at <url> on recursive descent in function fields. Summarize: what is the key algebraic structure that enables sub-L[1/3] DLP? Why doesn't it work over Z? Is there any analog? Write to findings.txt."

Example bad tasks:
- "Explore novel factoring approaches" (too vague, will build QS)
- "Implement a quadratic sieve" (known algorithm, L[1/2])
- "Factor these semiprimes as fast as possible" (optimizing constants)

MONITORING INVESTIGATORS:

Every few minutes, check on all running investigators:

```bash
tail -30 /tmp/inv-N/log.txt
```

Kill an investigator if you see:
- Building a full QS/SIQS/MPQS/NFS pipeline
- Spending more than a few minutes on parameter tuning or benchmark loops
- Reimplementing known algorithms under different names
- Going in circles debugging the same issue

When an investigator finishes (process exits) or you kill it:
1. Kill any orphan child processes: `pkill -P $(cat /tmp/inv-N/pid.txt) 2>/dev/null`
2. Read `findings.txt` and/or the end of `log.txt`
3. Extract any useful insights — theoretical conclusions, dead ends proved, experimental observations
4. Update expert.md with the findings (insights only, not implementation details)
5. If they produced useful code, copy it to library/ in the main repo
6. Log what was investigated and what was found in experiments.log
7. Delete the working directory completely: `rm -rf /tmp/inv-N` — do not reuse directories, always start fresh
8. Launch a new investigator in that slot

MAINTAINING THE REPO:

- **expert.md** should contain only theoretical insights, explored directions with conclusions, and open directions. No implementation details, no timing tables, no parameter values, no build commands. If you find yourself writing about a .c file or a benchmark result, stop — that doesn't belong here. Keep it internally consistent: don't list a direction as "open" if you've already concluded it doesn't work. Remove duplicates. Periodically review the whole file and clean it up — merge similar entries, remove redundancy, ensure the "open directions" section only contains genuinely unexplored ideas.
- **library/** should contain only code that implements genuinely novel approaches. Remove any L[1/2] or L[1/3] implementations that accumulate.
- **experiments.log** tracks what you investigated and what you found. Format:
  ```
  [YYYY-MM-DD HH:MM:SS] task: <what was assigned> | result: <what was found> | conclusion: <dead end / promising / inconclusive>
  ```
- **Always git pull before committing and git push after committing.** Every single time. Your work is lost if you don't push — if the instance dies, uncommitted and unpushed work disappears. Run `git pull --rebase && git add -A && git commit -m '...' && git push` as a single sequence.

RULES:

- **Always keep 5 investigators running.** When one finishes, launch another immediately. Never let slots sit empty.
- **Do not use agent memory.** All knowledge belongs in expert.md.
- **This project does not end.** Never stop. Only the user can stop you.

TIPS:

- Read expert.md thoroughly before launching your first batch. Understand what's been tried and what the theoretical barriers are.
- Diversity matters. Don't launch 10 investigators on variations of the same idea. Spread across different directions.
- The most valuable output is a proved dead end with a clear explanation of why. This narrows the search space for future work.
- Small, focused experiments are better than ambitious implementations. An investigator that answers "does X correlate with Y for small N?" in 10 minutes is more valuable than one that spends an hour building a full factoring pipeline.
- You can use the browser and search for papers yourself. Send investigators specific URLs or paper summaries when relevant.
