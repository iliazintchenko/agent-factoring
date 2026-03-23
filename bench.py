#!/usr/bin/env python3
"""Benchmark factoring implementations on semiprimes.json"""
import json, subprocess, time, sys

binary = sys.argv[1] if len(sys.argv) > 1 else './siqs_native'
sizes = sys.argv[2:] if len(sys.argv) > 2 else [str(s) for s in range(30, 76)]

with open('semiprimes.json') as f:
    sp = json.load(f)

results = {}
for size in sizes:
    if size not in sp:
        continue
    times = []
    all_ok = True
    for i, n in enumerate(sp[size]):
        t0 = time.monotonic()
        r = subprocess.run(['timeout', '295', binary, n], capture_output=True, text=True)
        t1 = time.monotonic()
        factor = r.stdout.strip()
        ok = factor != '' and 'FAIL' not in factor
        if ok:
            f_val = int(factor)
            N = int(n)
            ok = N % f_val == 0 and f_val > 1 and f_val < N
        elapsed_t = t1 - t0 if ok else -1
        times.append(elapsed_t)
        if not ok:
            all_ok = False
            print(f'  {size}d[{i}]: FAIL ({t1-t0:.1f}s)', file=sys.stderr)
        else:
            print(f'  {size}d[{i}]: {elapsed_t:.3f}s OK', file=sys.stderr)
    worst = max(times)
    if worst > 295:
        worst = -1
    results[size] = worst
    status = 'OK' if all_ok and worst > 0 else 'PARTIAL' if any(t > 0 for t in times) else 'FAIL'
    print(f'{size}d worst={worst:.3f}s [{status}]  all={[round(t,3) for t in times]}')
    if worst > 280 or worst < 0:
        print(f'Stopping at {size}d (timeout)')
        break

print('\n--- JSON ---')
print(json.dumps({k: round(v, 3) for k, v in results.items() if v > 0}, indent=2))
