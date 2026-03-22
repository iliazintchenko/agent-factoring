#!/usr/bin/env python3
"""Benchmark YAFU SIQS single-threaded across all sizes."""
import json, subprocess, time, sys, os, tempfile

YAFU = "/tmp/agent-factoring-3/yafu/yafu"

with open("semiprimes.json") as f:
    semiprimes = json.load(f)

sizes = sorted(semiprimes.keys(), key=int)
if len(sys.argv) > 1:
    sizes = sys.argv[1:]

# Extra args for YAFU
extra_args = "-siqsNB 16"  # NB=16 gives 2-10% improvement over default NB=8
if "--siqsB" in sys.argv:
    idx = sys.argv.index("--siqsB")
    extra_args += f" -siqsB {sys.argv[idx+1]}"
    sizes = [s for s in sizes if s not in ["--siqsB", sys.argv[idx+1]]]

for size in sizes:
    nums = semiprimes.get(size, [])
    if not nums:
        continue

    worst = 0
    times = []
    all_ok = True

    # Run all 5 in parallel as separate processes
    procs = []
    workdirs = []
    starts = []
    for n in nums:
        wd = tempfile.mkdtemp(prefix="yafu_")
        workdirs.append(wd)
        cmd = f"echo 'siqs({n})' | {YAFU} -threads 1 -seed 42 {extra_args}"
        starts.append(time.time())
        p = subprocess.Popen(cmd, shell=True, cwd=wd,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        procs.append(p)

    for i, (p, wd) in enumerate(zip(procs, workdirs)):
        try:
            out, _ = p.communicate(timeout=280)
            elapsed = time.time() - starts[i]
            out_str = out.decode(errors='replace')
            if 'P' in out_str and 'ans' in out_str:
                times.append(elapsed)
            else:
                times.append(float('inf'))
                all_ok = False
        except subprocess.TimeoutExpired:
            p.kill()
            p.communicate()
            times.append(float('inf'))
            all_ok = False
        finally:
            subprocess.run(f"rm -rf {wd}", shell=True)

    worst = max(times)
    status = "OK" if all_ok else "FAIL"
    times_str = ", ".join(f"{t:.1f}" if t < 1000 else "FAIL" for t in times)
    print(f"{size}d: worst={worst:.1f}s [{times_str}] {status} {extra_args}")
