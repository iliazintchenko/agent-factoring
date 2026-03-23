#!/usr/bin/env python3
"""Run scaling tests for a factoring approach across semiprime sizes."""

import json
import subprocess
import sys
import time
import os

def run_test(binary, number, timeout=295):
    """Run factoring binary on a number, return (time_seconds, factor1, factor2) or (None, None, None) on failure."""
    start = time.monotonic()
    try:
        result = subprocess.run(
            ["timeout", str(timeout), binary, number],
            capture_output=True, text=True, timeout=timeout + 10
        )
        elapsed = time.monotonic() - start
        if result.returncode == 0 and result.stdout.strip():
            parts = result.stdout.strip().split()
            if len(parts) == 2:
                f1, f2 = parts
                # Verify
                if int(f1) * int(f2) == int(number):
                    return elapsed, f1, f2
                else:
                    print(f"  WRONG: {f1} * {f2} != {number}", file=sys.stderr)
        if result.stderr:
            for line in result.stderr.strip().split('\n')[-3:]:
                print(f"  stderr: {line}", file=sys.stderr)
        return None, None, None
    except (subprocess.TimeoutExpired, Exception) as e:
        return None, None, None

def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <binary> <approach_name> [min_digits] [max_digits]")
        sys.exit(1)

    binary = sys.argv[1]
    approach = sys.argv[2]
    min_digits = int(sys.argv[3]) if len(sys.argv) > 3 else 30
    max_digits = int(sys.argv[4]) if len(sys.argv) > 4 else 100

    with open("semiprimes.json") as f:
        semiprimes = json.load(f)

    # Load existing scaling data
    scaling_file = "algo-scaling.json"
    if os.path.exists(scaling_file):
        with open(scaling_file) as f:
            scaling = json.load(f)
    else:
        scaling = {}

    if approach not in scaling:
        scaling[approach] = {}

    log_entries = []

    for size in range(min_digits, max_digits + 1):
        key = str(size)
        if key not in semiprimes:
            continue

        nums = semiprimes[key]
        times = []
        all_ok = True

        for i, num in enumerate(nums):
            print(f"[{approach}] {size} digits, #{i+1}/5: {num[:20]}...", file=sys.stderr)
            t, f1, f2 = run_test(binary, num)
            if t is not None:
                times.append(t)
                print(f"  -> {t:.3f}s", file=sys.stderr)
                log_entries.append(
                    f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] size: {size} | approach: {approach} | "
                    f"time: {t:.3f} | notes: factored {num[:15]}... = {f1[:15]}... * {f2[:15]}..."
                )
            else:
                all_ok = False
                print(f"  -> FAIL", file=sys.stderr)
                log_entries.append(
                    f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] size: {size} | approach: {approach} | "
                    f"time: FAIL | notes: failed to factor {num[:15]}..."
                )

        if all_ok and times:
            worst = max(times)
            scaling[approach][key] = round(worst, 3)
            print(f"[{approach}] {size} digits: worst={worst:.3f}s", file=sys.stderr)
        elif times:
            # Partial success - record but note it
            print(f"[{approach}] {size} digits: {len(times)}/5 succeeded", file=sys.stderr)
        else:
            print(f"[{approach}] {size} digits: all failed, stopping", file=sys.stderr)
            break

    # Save scaling data
    with open(scaling_file, 'w') as f:
        json.dump(scaling, f, indent=2)
        f.write('\n')

    # Append to experiments log
    with open("experiments.log", 'a') as f:
        for entry in log_entries:
            f.write(entry + '\n')

    print(f"\nScaling results for {approach}:", file=sys.stderr)
    for k in sorted(scaling[approach].keys(), key=int):
        print(f"  {k} digits: {scaling[approach][k]}s", file=sys.stderr)

if __name__ == "__main__":
    main()
