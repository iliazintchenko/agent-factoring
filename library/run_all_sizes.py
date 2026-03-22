#!/usr/bin/env python3
"""Run all 5 semiprimes for each size, single-threaded, and update best-algos.json."""
import json, subprocess, time, sys, os, tempfile, datetime

YAFU = "/tmp/agent-factoring-2/yafu/yafu"  # SKYLAKEX+AVX512+VBITS=256 build

with open("semiprimes.json") as f:
    semiprimes = json.load(f)

# Optimal siqsB values per size (empty = use default)
SIQSB = {}  # Will be populated as we discover optimal values

# Parse args
sizes = sorted(semiprimes.keys(), key=int)
if len(sys.argv) > 1:
    sizes = sys.argv[1:]

update_json = "--update" in sys.argv
sizes = [s for s in sizes if s != "--update"]

try:
    with open("best-algos.json") as f:
        best_algos = json.load(f)
except:
    best_algos = {}

for size in sizes:
    nums = semiprimes.get(size, [])
    if not nums:
        continue

    siqsb_arg = f"-siqsB {SIQSB[int(size)]}" if int(size) in SIQSB else ""

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
    for i, (p, wd) in enumerate(zip(procs, workdirs)):
        try:
            out, _ = p.communicate(timeout=280)
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
    print(f"{size}d: worst={worst:.1f}s [{times_str}] {status}")

    if all_ok and update_json:
        run_cmd = f"WORKDIR=$(mktemp -d /tmp/yafu_XXXXXX) && cd $WORKDIR && echo 'siqs(<N>)' | timeout 280 {YAFU} -threads 1 -seed 42 {siqsb_arg} && rm -rf $WORKDIR"
        if size not in best_algos or worst < best_algos[size].get("time", float('inf')):
            best_algos[size] = {"time": round(worst, 1), "run": run_cmd}
            print(f"  -> Updated best-algos.json for {size}d")

    # Log to experiments.log
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open("experiments.log", "a") as f:
        siqsb_note = f", siqsB={SIQSB[int(size)]}" if int(size) in SIQSB else ""
        f.write(f"[{ts}] size: {size} | approach: YAFU SIQS 1-thread{siqsb_note} | time: {worst:.1f} | notes: worst of 5, all times: [{times_str}]\n")

if update_json:
    with open("best-algos.json", "w") as f:
        json.dump(best_algos, f, indent=2, sort_keys=lambda x: int(x) if x.isdigit() else x)
        f.write("\n")
