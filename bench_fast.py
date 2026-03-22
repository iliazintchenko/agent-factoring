#!/usr/bin/env python3
"""Benchmark fast_factor.sh (YAFU+msieve fallback) across digit sizes."""
import json, subprocess, sys, time, os, re

with open('semiprimes.json') as f:
    data = json.load(f)

sizes = sorted(data.keys(), key=int)
if len(sys.argv) > 1:
    lo = int(sys.argv[1])
    hi = int(sys.argv[2]) if len(sys.argv) > 2 else lo
    sizes = [s for s in sizes if lo <= int(s) <= hi]

results = {}
for size in sizes:
    nums = data[size]
    times = []
    all_ok = True
    for i, n in enumerate(nums):
        t0 = time.monotonic()
        try:
            r = subprocess.run(
                ['bash', 'library/fast_factor.sh', n, '280'],
                capture_output=True, timeout=290, cwd='/tmp/agent-factoring-3')
            elapsed = time.monotonic() - t0
            factor = r.stdout.decode().strip()
            if factor and r.returncode == 0:
                nn, ff = int(n), int(factor)
                if nn % ff == 0 and ff > 1 and ff < nn:
                    times.append(elapsed)
                else:
                    print(f"  WRONG: {size}d #{i}: {factor}")
                    all_ok = False; times.append(300)
            else:
                print(f"  FAIL: {size}d #{i}")
                all_ok = False; times.append(300)
        except subprocess.TimeoutExpired:
            print(f"  TIMEOUT: {size}d #{i}")
            all_ok = False; times.append(300)

    worst = max(times) if times else 300
    status = "OK" if all_ok else "PARTIAL"
    print(f"{size}d: worst={worst:.3f}s  [{', '.join(f'{t:.3f}' for t in times)}]  {status}")
    sys.stdout.flush()
    results[size] = {"worst": worst, "all_ok": all_ok}

print("\n--- Summary ---")
for size in sorted(results.keys(), key=int):
    r = results[size]
    print(f"  {size}d: {r['worst']:.3f}s {'OK' if r['all_ok'] else 'FAIL'}")
