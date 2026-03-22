#!/usr/bin/env python3
"""
Benchmark YAFU SIQS on all semiprimes 30-89 digits.
Runs all 5 numbers per size in parallel (each single-threaded).
Reports worst-case time per size.

Usage: python3 library/bench_all.py [start_size] [end_size] [yafu_binary]
"""
import json, subprocess, sys, time, os, tempfile

# Defaults
START = int(sys.argv[1]) if len(sys.argv) > 1 else 30
END = int(sys.argv[2]) if len(sys.argv) > 2 else 89
YAFU = sys.argv[3] if len(sys.argv) > 3 else "/tmp/agent-factoring-4/yafu/yafu"

# Load semiprimes
with open(os.path.join(os.path.dirname(__file__), "..", "semiprimes.json")) as f:
    semiprimes = json.load(f)

# Best known params per size (from best-algos.json / expert.md)
PARAMS = {}
for s in range(30, 45):
    PARAMS[str(s)] = {"cmd": "smallmpqs", "extra": ""}
for s in range(45, 73):
    PARAMS[str(s)] = {"cmd": "siqs", "extra": ""}
for s in range(73, 77):
    PARAMS[str(s)] = {"cmd": "siqs", "extra": "-siqsNB 10" if s == 73 else ""}
PARAMS["77"] = {"cmd": "siqs", "extra": "-siqsNB 11"}
PARAMS["78"] = {"cmd": "siqs", "extra": "-siqsNB 11"}
PARAMS["79"] = {"cmd": "siqs", "extra": "-siqsNB 11"}
PARAMS["80"] = {"cmd": "siqs", "extra": "-siqsNB 12"}
PARAMS["81"] = {"cmd": "siqs", "extra": "-siqsNB 12"}
for s in range(82, 85):
    PARAMS[str(s)] = {"cmd": "siqs", "extra": "-siqsNB 14"}
PARAMS["85"] = {"cmd": "siqs", "extra": "-siqsNB 18 -siqsB 70000"}
PARAMS["86"] = {"cmd": "siqs", "extra": "-siqsNB 18 -siqsB 80000"}
PARAMS["87"] = {"cmd": "siqs", "extra": "-siqsNB 18 -siqsB 90000"}
PARAMS["88"] = {"cmd": "siqs", "extra": "-siqsNB 18 -siqsB 80000"}
PARAMS["89"] = {"cmd": "siqs", "extra": "-siqsNB 18 -siqsB 100000"}

def run_one(yafu, N, cmd, extra_params, timeout_sec=295):
    """Run YAFU on one number. Returns (wallclock_time, factor_found)."""
    workdir = tempfile.mkdtemp(prefix="yafu_")
    expr = f"{cmd}({N})"
    full_cmd = f'cd {workdir} && echo "{expr}" | timeout {timeout_sec} {yafu} -threads 1 -seed 42 {extra_params}'

    t0 = time.monotonic()
    try:
        result = subprocess.run(full_cmd, shell=True, capture_output=True, text=True, timeout=timeout_sec + 5)
        elapsed = time.monotonic() - t0
        # Check if factor was found
        output = result.stdout + result.stderr
        found = "P" in output and "ans" in output
        return (elapsed, found)
    except subprocess.TimeoutExpired:
        return (timeout_sec, False)
    finally:
        subprocess.run(f"rm -rf {workdir}", shell=True, capture_output=True)

results = {}
for size in range(START, END + 1):
    s = str(size)
    if s not in semiprimes:
        continue

    params = PARAMS.get(s, {"cmd": "siqs", "extra": ""})
    nums = semiprimes[s]

    # Run all 5 in parallel
    import concurrent.futures
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        futures = []
        for i, N in enumerate(nums):
            futures.append(executor.submit(run_one, YAFU, N, params["cmd"], params["extra"]))

        times = []
        all_found = True
        for i, f in enumerate(futures):
            elapsed, found = f.result()
            times.append(elapsed)
            if not found:
                all_found = False
                print(f"  {size}d[{i}]: TIMEOUT or FAIL at {elapsed:.1f}s", file=sys.stderr)

    worst = max(times)
    results[s] = {
        "worst": worst,
        "all_found": all_found,
        "times": [round(t, 1) for t in times],
        "params": params
    }

    status = "OK" if all_found else "FAIL"
    print(f"{size}d: worst={worst:.1f}s [{', '.join(f'{t:.1f}' for t in times)}] {status} ({params['cmd']} {params['extra']})")

# Summary
print(f"\n--- Summary (load avg: {os.getloadavg()[0]:.1f}) ---")
for s in sorted(results.keys(), key=int):
    r = results[s]
    status = "OK" if r["all_found"] else "FAIL"
    print(f"  {s}d: {r['worst']:.1f}s {status}")
