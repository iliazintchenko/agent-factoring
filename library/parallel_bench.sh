#!/bin/bash
# Parallel benchmarking across all semiprimes
# Usage: parallel_bench.sh <binary> <approach> [min_digits] [max_digits] [max_parallel]

BINARY="$1"
APPROACH="$2"
MIN_D="${3:-30}"
MAX_D="${4:-100}"
MAX_PAR="${5:-20}"
SCRIPT_DIR="$(dirname "$0")"
SEMIPRIMES="$SCRIPT_DIR/../semiprimes.json"
LOG="$SCRIPT_DIR/../experiments.log"
SCALING="$SCRIPT_DIR/../algo-scaling.json"
RESULTS_DIR=$(mktemp -d)

echo "Benchmarking $APPROACH from $MIN_D to $MAX_D digits (max $MAX_PAR parallel)"

# Launch all jobs for a given digit size, wait, collect results
for digits in $(seq $MIN_D $MAX_D); do
    nums=$(python3 -c "
import json
with open('$SEMIPRIMES') as f:
    data = json.load(f)
key = str($digits)
if key in data:
    for n in data[key]:
        print(n)
" 2>/dev/null)
    [ -z "$nums" ] && continue

    # Launch 5 jobs in parallel
    idx=0
    pids=()
    while IFS= read -r num; do
        bash "$SCRIPT_DIR/run_bench.sh" "$BINARY" "$APPROACH" "$digits" "$idx" "$num" \
            > "$RESULTS_DIR/${digits}_${idx}.txt" &
        pids+=($!)
        idx=$((idx + 1))
    done <<< "$nums"

    # Wait for all 5
    for pid in "${pids[@]}"; do
        wait "$pid"
    done

    # Collect results
    worst_time=0
    all_passed=true
    for i in $(seq 0 $((idx-1))); do
        result=$(cat "$RESULTS_DIR/${digits}_${i}.txt" 2>/dev/null)
        d=$(echo "$result" | awk '{print $1}')
        ix=$(echo "$result" | awk '{print $2}')
        t=$(echo "$result" | awk '{print $3}')
        f=$(echo "$result" | awk '{print $4}')

        timestamp=$(date '+%Y-%m-%d %H:%M:%S')
        if [ "$t" != "FAIL" ] && [ "$t" != "TIMEOUT" ]; then
            echo "[$timestamp] size: $digits | approach: $APPROACH | time: ${t}s | notes: factor=$f" >> "$LOG"
            is_worse=$(echo "$t > $worst_time" | bc 2>/dev/null)
            [ "$is_worse" -eq 1 ] 2>/dev/null && worst_time=$t
        else
            echo "[$timestamp] size: $digits | approach: $APPROACH | time: FAIL | notes: $t $f" >> "$LOG"
            all_passed=false
        fi
    done

    if [ "$all_passed" = true ]; then
        worst_rounded=$(printf "%.3f" "$worst_time")
        echo "  $digits digits: worst = ${worst_rounded}s"

        python3 -c "
import json, fcntl
path = '$SCALING'
try:
    with open(path) as f:
        data = json.load(f)
except:
    data = {}
if '$APPROACH' not in data:
    data['$APPROACH'] = {}
data['$APPROACH'][str($digits)] = float($worst_rounded)
with open(path, 'w') as f:
    json.dump(data, f, indent=2)
    f.write('\n')
"
    else
        echo "  $digits digits: FAILED"
        break
    fi
done

rm -rf "$RESULTS_DIR"
echo "Done."
