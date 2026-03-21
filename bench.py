#!/usr/bin/env python3
"""Benchmark factoring across digit sizes. Runs one number at a time."""
import json, subprocess, sys, time, os

os.environ['LD_LIBRARY_PATH'] = '/usr/local/lib'

with open('semiprimes.json') as f:
    data = json.load(f)

binary = sys.argv[1] if len(sys.argv) > 1 else './factor'
sizes = sorted(data.keys(), key=int)
if len(sys.argv) > 2:
    lo = int(sys.argv[2])
    hi = int(sys.argv[3]) if len(sys.argv) > 3 else lo
    sizes = [s for s in sizes if lo <= int(s) <= hi]

results = {}
for size in sizes:
    nums = data[size]
    times = []
    all_ok = True
    for i, n in enumerate(nums):
        t0 = time.monotonic()
        try:
            r = subprocess.run([binary, n], capture_output=True, timeout=300,
                             env={**os.environ, 'LD_LIBRARY_PATH': '/usr/local/lib'})
            elapsed = time.monotonic() - t0
            if r.returncode == 0:
                factor = r.stdout.decode().strip()
                nn, ff = int(n), int(factor)
                if nn % ff == 0 and ff > 1 and ff < nn:
                    times.append(elapsed)
                else:
                    print(f"  WRONG: {size}d #{i}: {factor}")
                    all_ok = False; times.append(300)
            else:
                print(f"  FAIL: {size}d #{i}: {r.stderr.decode().strip()[:100]}")
                all_ok = False; times.append(300)
        except subprocess.TimeoutExpired:
            print(f"  TIMEOUT: {size}d #{i}")
            all_ok = False; times.append(300)

    worst = max(times) if times else 300
    status = "OK" if all_ok else "PARTIAL"
    print(f"{size}d: worst={worst:.3f}s  [{', '.join(f'{t:.3f}' for t in times)}]  {status}")
    sys.stdout.flush()
    results[size] = worst

print("\n--- Summary ---")
for size in sorted(results.keys(), key=int):
    t = results[size]
    print(f"  {size}d: {t:.3f}s {'OK' if t < 300 else 'FAIL'}")
