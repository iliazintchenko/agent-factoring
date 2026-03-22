#!/usr/bin/env python3
"""Benchmark factoring across semiprime sizes, one size at a time."""
import json, subprocess, sys, time, os
from concurrent.futures import ProcessPoolExecutor, as_completed

PFACTOR_BIN = "./library/pfactor"
YAFU_SCRIPT = "./library/run_yafu.sh"
SEMIPRIMES = json.load(open("semiprimes.json"))
TIMEOUT = 280

def factor_one(args):
    size, idx, n, method, timeout, extra = args
    t0 = time.monotonic()
    try:
        if method == "pfactor":
            workers = extra.get("workers", 8)
            cmd = [PFACTOR_BIN, n, str(timeout), str(workers)]
        elif method == "yafu":
            threads = extra.get("threads", 8)
            cmd = [YAFU_SCRIPT, n, str(timeout), str(threads)]
        else:
            return (size, idx, 0, "FAIL", f"unknown method {method}")

        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout + 10)
        elapsed = time.monotonic() - t0
        if r.returncode == 0:
            factor = r.stdout.strip()
            if factor and factor != "FAIL":
                return (size, idx, elapsed, factor, "")
        return (size, idx, elapsed, "FAIL", r.stderr.strip()[:200])
    except subprocess.TimeoutExpired:
        elapsed = time.monotonic() - t0
        return (size, idx, elapsed, "TIMEOUT", "")

def pick_method(digits):
    """Choose method and params based on digit count."""
    if digits <= 57:
        # Parallel ECM - 5 numbers * 8 workers = 40 processes (fits in 48 cores)
        return "pfactor", {"workers": 8}
    elif digits <= 70:
        # YAFU SIQS - 5 numbers * 9 threads = 45 (fits in 48 cores)
        return "yafu", {"threads": 9}
    elif digits <= 85:
        return "yafu", {"threads": 9}
    else:
        # For 90+ digits, fewer parallel runs but more threads each
        return "yafu", {"threads": 9}

def main():
    if len(sys.argv) > 1:
        sizes = [int(x) for x in sys.argv[1:]]
    else:
        sizes = sorted(int(k) for k in SEMIPRIMES.keys())

    summary = {}
    for s in sizes:
        method, extra = pick_method(s)
        tasks = []
        for idx, n in enumerate(SEMIPRIMES[str(s)]):
            tasks.append((s, idx, n, method, TIMEOUT, extra))

        results = []
        with ProcessPoolExecutor(max_workers=5) as ex:
            futures = {ex.submit(factor_one, t): t for t in tasks}
            for f in as_completed(futures):
                size, idx, elapsed, factor, stderr = f.result()
                results.append((idx, elapsed, factor, stderr))
                status = f"OK {elapsed:.3f}s" if factor not in ("FAIL", "TIMEOUT") else f"{factor}"
                print(f"  {s}d #{idx}: {status} [{method}]", flush=True)

        times = [e[1] for e in results if e[2] not in ("FAIL", "TIMEOUT")]
        fails = sum(1 for e in results if e[2] in ("FAIL", "TIMEOUT"))
        if fails > 0:
            print(f"  => {s}d: FAIL ({fails}/5 failed)")
        else:
            max_t = max(times)
            print(f"  => {s}d: worst={max_t:.3f}s")
            summary[str(s)] = {"time": round(max_t, 3), "method": method}

    print("\n=== SUMMARY ===")
    for k in sorted(summary.keys(), key=int):
        v = summary[k]
        print(f"  {k:>3}d: {v['time']:>10.3f}s  ({v['method']})")

if __name__ == "__main__":
    main()
