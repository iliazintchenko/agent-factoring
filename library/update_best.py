#!/usr/bin/env python3
"""Update best-algos.json with single-threaded YAFU results."""
import json, sys

def get_siqsB(digits):
    if digits < 70:
        return None
    params = {70: 8000, 75: 12000, 78: 10000, 80: 15000, 82: 15000,
              85: 25000, 87: 30000, 90: 40000, 93: 50000, 95: 60000,
              98: 75000, 100: 90000}
    best = None
    for k in sorted(params.keys()):
        if k <= digits:
            best = params[k]
    return best

def make_run_cmd(digits, siqsB=None):
    cmd = 'WORKDIR=$(mktemp -d /tmp/yafu_XXXXXX) && cd $WORKDIR && echo "siqs(<N>)" | LD_LIBRARY_PATH=/usr/local/lib /tmp/agent-factoring-3/yafu/yafu -threads 1 -seed 42'
    if siqsB:
        cmd += f' -siqsB {siqsB}'
    cmd += ' && rm -rf $WORKDIR'
    return cmd

def update(best_algos_path, results):
    """results is a dict of {digit_str: worst_time}"""
    with open(best_algos_path) as f:
        best = json.load(f)

    for digits_str, worst_time in results.items():
        digits = int(digits_str)
        siqsB = get_siqsB(digits)
        run_cmd = make_run_cmd(digits, siqsB)

        if digits_str not in best or worst_time < best[digits_str]["time"]:
            best[digits_str] = {"time": round(worst_time, 3), "run": run_cmd}
            print(f"Updated {digits_str}d: {worst_time:.3f}s")

    # Sort by digit count
    sorted_best = dict(sorted(best.items(), key=lambda x: int(x[0])))
    with open(best_algos_path, 'w') as f:
        json.dump(sorted_best, f, indent=2)
    print(f"Wrote {len(sorted_best)} entries to {best_algos_path}")

if __name__ == "__main__":
    # Usage: python3 update_best.py <results_json_file>
    # Or just import and call update()
    pass
