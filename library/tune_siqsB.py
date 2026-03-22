#!/usr/bin/env python3
"""Systematically find optimal siqsB for a digit size by testing ALL 5 semiprimes."""
import json, subprocess, time, sys, os, tempfile

YAFU = "/tmp/agent-factoring-2/yafu/yafu"

with open("semiprimes.json") as f:
    semiprimes = json.load(f)

size = sys.argv[1]
nums = semiprimes[size]

# siqsB values to test
if len(sys.argv) > 2:
    b_values = [int(x) for x in sys.argv[2:]]
else:
    d = int(size)
    if d < 75:
        b_values = [0] + list(range(3000, 15000, 2000))
    elif d < 82:
        b_values = [0] + list(range(5000, 25000, 3000))
    elif d < 88:
        b_values = [0] + list(range(10000, 50000, 5000))
    elif d < 93:
        b_values = [0] + list(range(20000, 70000, 5000))
    elif d < 97:
        b_values = [0] + list(range(30000, 100000, 10000))
    else:
        b_values = [0] + list(range(50000, 150000, 10000))

best_worst = float('inf')
best_b = None

for b in b_values:
    siqsb_arg = f"-siqsB {b}" if b > 0 else ""
    label = f"siqsB={b}" if b > 0 else "default"

    procs = []
    workdirs = []
    starts = []
    for n in nums:
        wd = tempfile.mkdtemp(prefix="yafu_")
        workdirs.append(wd)
        cmd = f"echo 'siqs({n})' | {YAFU} -threads 1 -seed 42 {siqsb_arg}"
        starts.append(time.time())
        p = subprocess.Popen(cmd, shell=True, cwd=wd,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        procs.append(p)

    times = []
    all_ok = True
    timeout_val = min(280, best_worst * 1.5 + 5) if best_worst < float('inf') else 280
    for i, (p, wd) in enumerate(zip(procs, workdirs)):
        try:
            out, _ = p.communicate(timeout=timeout_val)
            elapsed = time.time() - starts[i]
            if b'ans' in out:
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
            subprocess.run(f"rm -rf {wd}", shell=True, capture_output=True)

    worst = max(times)
    status = "OK" if all_ok else "FAIL"
    times_str = ", ".join(f"{t:.1f}" if t < 1000 else "FAIL" for t in times)
    marker = " <-- BEST" if all_ok and worst < best_worst else ""
    print(f"  {label}: worst={worst:.1f}s [{times_str}] {status}{marker}")

    if all_ok and worst < best_worst:
        best_worst = worst
        best_b = b

b_str = f"siqsB={best_b}" if best_b and best_b > 0 else "default"
print(f"\nBest for {size}d: {b_str} at {best_worst:.1f}s")
