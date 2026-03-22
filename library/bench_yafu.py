#!/usr/bin/env python3
"""Benchmark YAFU SIQS on semiprimes, running all 5 per size in parallel."""
import json, subprocess, sys, time, os

YAFU = sys.argv[1] if len(sys.argv) > 1 else "/tmp/agent-factoring-4/yafu/yafu"
SIZES = sys.argv[2:] if len(sys.argv) > 2 else ["85", "86", "87", "88", "89"]

with open("/tmp/agent-factoring-9/semiprimes.json") as f:
    semiprimes = json.load(f)

# NB/B params per size (from expert.md / best-algos.json)
PARAMS = {
    "73": "-siqsNB 10", "74": "", "75": "", "76": "",
    "77": "-siqsNB 11", "78": "-siqsNB 11", "79": "-siqsNB 11",
    "80": "-siqsNB 12", "81": "-siqsNB 12",
    "82": "-siqsNB 14", "83": "-siqsNB 14", "84": "-siqsNB 14",
    "85": "-siqsNB 18 -siqsB 70000", "86": "-siqsNB 18 -siqsB 80000",
    "87": "-siqsNB 18 -siqsB 90000", "88": "-siqsNB 18 -siqsB 80000",
    "89": "-siqsNB 18 -siqsB 100000",
}

for size in SIZES:
    if size not in semiprimes:
        print(f"Size {size} not in semiprimes.json")
        continue

    nums = semiprimes[size]
    params = PARAMS.get(size, "")
    procs = []

    for i, N in enumerate(nums):
        workdir = f"/tmp/yafu_bench_{size}_{i}"
        os.makedirs(workdir, exist_ok=True)
        cmd = f'cd {workdir} && echo "siqs({N})" | timeout 295 {YAFU} -threads 1 -seed 42 {params} 2>&1'
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        procs.append((i, N, p, time.time()))

    results = []
    for i, N, p, t0 in procs:
        out = p.communicate()[0].decode()
        elapsed = time.time() - t0
        # Find factor in output
        factor = ""
        for line in out.split('\n'):
            if line.startswith('P') and '=' in line:
                factor = line.split('=')[1].strip()
                break
        if p.returncode == 124:
            results.append((i, "TIMEOUT", elapsed))
        elif factor:
            results.append((i, factor, elapsed))
        else:
            results.append((i, "FAIL", elapsed))
        # Cleanup
        workdir = f"/tmp/yafu_bench_{size}_{i}"
        subprocess.run(f"rm -rf {workdir}", shell=True)

    worst = max(r[2] for r in results)
    print(f"{size}d: worst={worst:.1f}s | " + " | ".join(f"[{r[0]}]={r[2]:.1f}s{'✓' if r[1] not in ('TIMEOUT','FAIL') else '✗'}" for r in results))
