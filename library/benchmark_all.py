#!/usr/bin/env python3
"""Run comprehensive benchmark for all sizes, one size at a time."""
import json, subprocess, sys, time, os
from concurrent.futures import ProcessPoolExecutor, as_completed

YAFU = "/tmp/agent-factoring-1/yafu/yafu"
SEMIPRIMES = json.load(open("semiprimes.json"))
TIMEOUT = 280

def run_yafu_siqs(args):
    size, idx, n, threads, timeout = args
    workdir = f"/tmp/yafu_{os.getpid()}_{idx}"
    os.makedirs(workdir, exist_ok=True)
    t0 = time.monotonic()
    try:
        r = subprocess.run(
            f"cd {workdir} && echo 'siqs({n})' | timeout {timeout} {YAFU} -threads {threads}",
            shell=True, capture_output=True, text=True, timeout=timeout + 10
        )
        elapsed = time.monotonic() - t0
        # Extract factors
        factors = []
        for line in r.stdout.split('\n'):
            line = line.strip()
            if line.startswith('P') and ' = ' in line:
                factors.append(line.split(' = ')[1])
        subprocess.run(f"rm -rf {workdir}", shell=True)
        if factors:
            return (size, idx, elapsed, min(factors, key=lambda x: int(x)), "")
        return (size, idx, elapsed, "FAIL", r.stderr[-200:] if r.stderr else "no factors found")
    except subprocess.TimeoutExpired:
        subprocess.run(f"rm -rf {workdir}", shell=True)
        return (size, idx, time.monotonic() - t0, "TIMEOUT", "")

def main():
    if len(sys.argv) > 1:
        sizes = [int(x) for x in sys.argv[1:]]
    else:
        sizes = list(range(30, 101))

    results_all = {}
    for s in sizes:
        nums = SEMIPRIMES[str(s)]
        # Thread allocation: 5 parallel * T threads = 48 cores
        if s <= 90:
            threads = 9
            max_parallel = 5
        elif s <= 95:
            threads = 24
            max_parallel = 2
        else:
            threads = 48
            max_parallel = 1

        tasks = [(s, idx, n, threads, TIMEOUT) for idx, n in enumerate(nums)]
        results = []

        with ProcessPoolExecutor(max_workers=max_parallel) as ex:
            futures = {ex.submit(run_yafu_siqs, t): t for t in tasks}
            for f in as_completed(futures):
                r = f.result()
                results.append(r)
                status = f"OK {r[2]:.3f}s" if r[3] not in ("FAIL", "TIMEOUT") else r[3]
                print(f"  {s}d #{r[1]}: {status} (T={threads})", flush=True)

        times = [r[2] for r in results if r[3] not in ("FAIL", "TIMEOUT")]
        fails = sum(1 for r in results if r[3] in ("FAIL", "TIMEOUT"))
        if fails > 0:
            print(f"  => {s}d: FAIL ({fails}/5)", flush=True)
        else:
            max_t = max(times)
            print(f"  => {s}d: worst={max_t:.3f}s", flush=True)
            results_all[str(s)] = round(max_t, 3)

    print("\n=== SUMMARY ===")
    for k in sorted(results_all.keys(), key=int):
        print(f"  {k:>3}d: {results_all[k]:>10.3f}s")

    # Update best-algos.json
    try:
        existing = json.load(open("best-algos.json"))
    except:
        existing = {}

    for k, v in results_all.items():
        s = int(k)
        if s <= 90:
            threads = 9
        elif s <= 95:
            threads = 24
        else:
            threads = 48
        run_cmd = f"WORKDIR=$(mktemp -d /tmp/yafu_XXXXXX) && cd $WORKDIR && echo 'siqs(<N>)' | timeout 280 {YAFU} -threads {threads} && rm -rf $WORKDIR"
        if k not in existing or v < existing[k]["time"]:
            existing[k] = {"time": v, "run": run_cmd}

    with open("best-algos.json", "w") as f:
        json.dump(existing, f, indent=2, sort_keys=lambda x: int(x) if x.isdigit() else x)

    print("\nbest-algos.json updated")

if __name__ == "__main__":
    main()
