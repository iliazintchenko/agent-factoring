#!/usr/bin/env python3
"""
Benchmark a factoring binary against semiprimes.
Usage: python3 bench.py <binary> <approach_name> [sizes...]
If no sizes specified, runs all sizes from 30-100.
Records to experiments.log and prints algo-scaling data.
"""
import json, subprocess, sys, os, time
from datetime import datetime
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SEMIPRIMES = REPO / "semiprimes.json"
LOG = REPO / "experiments.log"

def run_one(binary, n, timeout=295):
    """Run binary on a single number. Returns (factor, elapsed) or (None, elapsed)."""
    start = time.monotonic()
    try:
        result = subprocess.run(
            ["timeout", str(timeout), binary, n],
            capture_output=True, text=True, timeout=timeout + 5
        )
        elapsed = time.monotonic() - start
        stdout = result.stdout.strip()
        stderr = result.stderr.strip()

        for line in stdout.split('\n'):
            if line.startswith("FACTOR:"):
                factor = line.split(":")[1].strip()
                return factor, elapsed, stderr
        return None, elapsed, stderr
    except subprocess.TimeoutExpired:
        elapsed = time.monotonic() - start
        return None, elapsed, "TIMEOUT"

def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <binary> <approach_name> [sizes...]")
        sys.exit(1)

    binary = os.path.abspath(sys.argv[1])
    approach = sys.argv[2]

    with open(SEMIPRIMES) as f:
        data = json.load(f)

    if len(sys.argv) > 3:
        sizes = sys.argv[3:]
    else:
        sizes = sorted(data.keys(), key=int)

    results = {}
    for size in sizes:
        nums = data.get(str(size), [])
        if not nums:
            continue
        worst = 0
        all_ok = True
        for n in nums:
            factor, elapsed, stderr = run_one(binary, n)
            ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            if factor:
                with open(LOG, "a") as f:
                    f.write(f"[{ts}] size: {size} | approach: {approach} | time: {elapsed:.3f}s | notes: found factor {factor}\n")
                worst = max(worst, elapsed)
            else:
                with open(LOG, "a") as f:
                    f.write(f"[{ts}] size: {size} | approach: {approach} | time: FAIL | notes: no factor for {n[:20]}...\n")
                all_ok = False
                break  # no point continuing this size

        if all_ok:
            results[size] = round(worst, 3)
            print(f"  {size}-digit: {worst:.3f}s (worst of 5)")
        else:
            results[size] = "FAIL"
            print(f"  {size}-digit: FAIL")

    # Print JSON fragment for algo-scaling.json
    print(f"\nResults for '{approach}':")
    print(json.dumps({approach: results}, indent=2))

if __name__ == "__main__":
    main()
