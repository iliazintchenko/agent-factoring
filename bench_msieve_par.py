#!/usr/bin/env python3
"""Benchmark msieve: run 5 semiprimes per size in parallel with -t 9 each (45 cores total)."""
import json, subprocess, sys, time, os, re

with open('semiprimes.json') as f:
    data = json.load(f)

sizes = sorted(data.keys(), key=int)
if len(sys.argv) > 1:
    lo = int(sys.argv[1])
    hi = int(sys.argv[2]) if len(sys.argv) > 2 else lo
    sizes = [s for s in sizes if lo <= int(s) <= hi]

THREADS_PER = 9  # 5 * 9 = 45 threads, fits in 48 cores

results = {}
for size in sizes:
    nums = data[size]
    t0 = time.monotonic()
    procs = []
    sessions = []
    for i, n in enumerate(nums):
        s = f"/tmp/msieve_{size}_{i}_{os.getpid()}.dat"
        sessions.append(s)
        p = subprocess.Popen(
            ['./msieve', '-q', '-t', str(THREADS_PER), '-s', s, n],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            cwd='/tmp/agent-factoring-3')
        procs.append(p)

    times = []
    all_ok = True
    for i, (p, n) in enumerate(zip(procs, nums)):
        try:
            out, err = p.communicate(timeout=300)
            elapsed = time.monotonic() - t0
            factors = re.findall(r'p\d+: (\d+)', out.decode())
            if factors:
                f = int(factors[0])
                nn = int(n)
                if nn % f == 0 and f > 1 and f < nn:
                    times.append(elapsed)
                else:
                    print(f"  WRONG: {size}d #{i}")
                    all_ok = False; times.append(300)
            else:
                print(f"  FAIL: {size}d #{i}: {out.decode()[:60]}")
                all_ok = False; times.append(300)
        except subprocess.TimeoutExpired:
            p.kill()
            print(f"  TIMEOUT: {size}d #{i}")
            all_ok = False; times.append(300)
        finally:
            try: os.unlink(sessions[i])
            except: pass

    worst = max(times) if times else 300
    status = "OK" if all_ok else "PARTIAL"
    print(f"{size}d: worst={worst:.3f}s  [{', '.join(f'{t:.3f}' for t in times)}]  {status}")
    sys.stdout.flush()
    results[size] = {"worst": worst, "all_ok": all_ok}

print("\n--- Summary ---")
for size in sorted(results.keys(), key=int):
    r = results[size]
    print(f"  {size}d: {r['worst']:.3f}s {'OK' if r['all_ok'] else 'FAIL'}")
