#!/usr/bin/env python3
"""
Benchmark a factoring binary against semiprimes.
Usage: python3 bench.py <binary> <approach_name> [sizes...]
If no sizes specified, runs all sizes from 30-100.
Prints results to stdout.
"""
import json, subprocess, sys, os, time
from datetime import datetime
from pathlib import Path

SEMIPRIMES = Path(__file__).resolve().parent / "semiprimes.json"

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
        # Also try: last line has two space-separated numbers (factor cofactor)
        last_line = stdout.strip().split('\n')[-1] if stdout.strip() else ""
        parts = last_line.strip().split()
        if len(parts) == 2:
            try:
                f1, f2 = int(parts[0]), int(parts[1])
                n_val = int(n)
                if f1 > 1 and f2 > 1 and f1 * f2 == n_val:
                    return str(f1), elapsed, stderr
            except ValueError:
                pass
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
                worst = max(worst, elapsed)
            else:
                all_ok = False
                break  # no point continuing this size

        if all_ok:
            results[size] = round(worst, 3)
            print(f"  {size}-digit: {worst:.3f}s (worst of 5)")
        else:
            results[size] = "FAIL"
            print(f"  {size}-digit: FAIL")

    # Print summary
    print(f"\nResults for '{approach}':")
    for k, v in sorted(results.items(), key=lambda x: int(x[0])):
        print(f"  {k}-digit: {v}")

if __name__ == "__main__":
    main()
