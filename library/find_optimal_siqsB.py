#!/usr/bin/env python3
"""Find optimal siqsB for a given digit size by testing one semiprime."""
import json, subprocess, time, sys, os, tempfile

YAFU = "/tmp/agent-factoring-3/yafu/yafu"
if not os.path.exists(YAFU):
    YAFU = "/tmp/agent-factoring-2/yafu/yafu"

with open("semiprimes.json") as f:
    semiprimes = json.load(f)

size = sys.argv[1]
nums = semiprimes[size]
# Test on first number
n = nums[0]

# Test different siqsB values
if len(sys.argv) > 2:
    b_values = [int(x) for x in sys.argv[2:]]
else:
    # Auto range based on digit size
    d = int(size)
    if d < 70:
        b_values = list(range(2000, 12000, 2000))
    elif d < 80:
        b_values = list(range(5000, 25000, 3000))
    elif d < 85:
        b_values = list(range(10000, 40000, 5000))
    elif d < 90:
        b_values = list(range(15000, 60000, 5000))
    elif d < 95:
        b_values = list(range(25000, 80000, 5000))
    else:
        b_values = list(range(40000, 120000, 10000))

# Also test default (no siqsB)
test_configs = [("default", "")] + [(f"siqsB={b}", f"-siqsB {b}") for b in b_values]

best_time = float('inf')
best_config = None

for name, args in test_configs:
    wd = tempfile.mkdtemp(prefix="yafu_")
    cmd = f"echo 'siqs({n})' | {YAFU} -threads 1 -seed 42 {args}"
    start = time.time()
    try:
        result = subprocess.run(cmd, shell=True, cwd=wd, capture_output=True,
                              text=True, timeout=min(280, best_time * 2 + 10))
        elapsed = time.time() - start
        if 'ans' in result.stdout:
            print(f"  {name}: {elapsed:.1f}s")
            if elapsed < best_time:
                best_time = elapsed
                best_config = name
        else:
            print(f"  {name}: FAIL (no answer)")
    except subprocess.TimeoutExpired:
        elapsed = time.time() - start
        print(f"  {name}: TIMEOUT ({elapsed:.0f}s)")
    finally:
        subprocess.run(f"rm -rf {wd}", shell=True, capture_output=True)

print(f"\nBest for {size}d: {best_config} at {best_time:.1f}s")
