#!/usr/bin/env python3
"""Benchmark YAFU SIQS single-threaded on all semiprimes."""
import json, subprocess, sys, time, os, tempfile

YAFU = "/tmp/agent-factoring-3/yafu/yafu"
SEMIPRIMES = "/tmp/agent-factoring-3/semiprimes.json"

def factor_one(n, siqsB=None, timeout=295):
    """Factor a single number with YAFU single-threaded. Returns (time, success, output)."""
    workdir = tempfile.mkdtemp(prefix="yafu_bench_")
    cmd = f'echo "siqs({n})" | LD_LIBRARY_PATH=/usr/local/lib {YAFU} -threads 1 -seed 42'
    if siqsB:
        cmd += f' -siqsB {siqsB}'
    try:
        start = time.time()
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True,
                              timeout=timeout, cwd=workdir)
        elapsed = time.time() - start
        output = result.stdout + result.stderr
        success = "***factors found***" in output and ("P" in output or "p" in output)
        return elapsed, success, output
    except subprocess.TimeoutExpired:
        return timeout, False, "TIMEOUT"
    finally:
        subprocess.run(f"rm -rf {workdir}", shell=True, capture_output=True)

def get_siqsB(digits):
    """Get optimal siqsB for digit size."""
    if digits < 70:
        return None  # default is fine
    params = {70: 8000, 75: 12000, 78: 10000, 80: 15000, 82: 15000,
              85: 25000, 87: 30000, 90: 40000, 93: 50000, 95: 60000,
              98: 75000, 100: 90000}
    # Find closest key <= digits
    best = None
    for k in sorted(params.keys()):
        if k <= digits:
            best = params[k]
    return best

def bench_size(digits, nums):
    """Benchmark all 5 numbers of a size in parallel using subprocess."""
    siqsB = get_siqsB(digits)
    print(f"\n=== {digits}d (siqsB={siqsB}) ===")

    # Launch all 5 in parallel
    procs = []
    workdirs = []
    start = time.time()
    for i, n in enumerate(nums):
        workdir = tempfile.mkdtemp(prefix=f"yafu_b{digits}_{i}_")
        workdirs.append(workdir)
        cmd = f'echo "siqs({n})" | LD_LIBRARY_PATH=/usr/local/lib {YAFU} -threads 1 -seed 42'
        if siqsB:
            cmd += f' -siqsB {siqsB}'
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           cwd=workdir)
        procs.append((p, time.time(), n))

    # Collect results
    results = []
    for p, pstart, n in procs:
        try:
            stdout, stderr = p.communicate(timeout=295)
            elapsed = time.time() - pstart
            output = stdout.decode() + stderr.decode()
            success = "***factors found***" in output
            results.append((elapsed, success, n))
            status = "OK" if success else "FAIL"
            print(f"  {n[:20]}... {elapsed:.2f}s {status}")
        except subprocess.TimeoutExpired:
            p.kill()
            p.communicate()
            results.append((295, False, n))
            print(f"  {n[:20]}... TIMEOUT")

    # Cleanup
    for wd in workdirs:
        subprocess.run(f"rm -rf {wd}", shell=True, capture_output=True)

    wall = time.time() - start
    times = [r[0] for r in results]
    successes = [r[1] for r in results]
    worst = max(times)
    all_ok = all(successes)

    print(f"  worst={worst:.2f}s wall={wall:.2f}s all_ok={all_ok}")
    return worst, all_ok, results

def main():
    with open(SEMIPRIMES) as f:
        sp = json.load(f)

    start_size = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    end_size = int(sys.argv[2]) if len(sys.argv) > 2 else 100

    results = {}
    for digits in range(start_size, end_size + 1):
        key = str(digits)
        if key not in sp:
            continue
        worst, all_ok, details = bench_size(digits, sp[key])
        results[digits] = {"worst": worst, "all_ok": all_ok}

        if not all_ok:
            print(f"  WARNING: {digits}d had failures!")
        if worst > 290:
            print(f"  STOPPING: {digits}d too slow ({worst:.1f}s)")
            break

    # Summary
    print("\n=== SUMMARY ===")
    for d, r in sorted(results.items()):
        status = "OK" if r["all_ok"] else "FAIL"
        print(f"  {d}d: {r['worst']:.2f}s {status}")

if __name__ == "__main__":
    main()
