#!/bin/bash
# Run benchmarks on semiprimes for a given binary and approach name
# Usage: ./run_bench.sh <binary> <approach> [start_size] [end_size]

BINARY=$1
APPROACH=$2
START=${3:-30}
END=${4:-50}
SEMIPRIMES="../semiprimes.json"
LOG="../experiments.log"

if [ -z "$BINARY" ] || [ -z "$APPROACH" ]; then
    echo "Usage: $0 <binary> <approach> [start_size] [end_size]"
    exit 1
fi

for size in $(seq $START $END); do
    # Extract semiprimes for this size using python
    nums=$(python3 -c "
import json
with open('$SEMIPRIMES') as f:
    d = json.load(f)
key = str($size)
if key in d:
    for n in d[key]:
        print(n)
")
    if [ -z "$nums" ]; then continue; fi

    worst_time=0
    all_ok=1
    times=""

    for n in $nums; do
        t0=$(date +%s%N)
        result=$(timeout 295 $BINARY "$n" 2>/dev/null)
        status=$?
        t1=$(date +%s%N)
        elapsed=$(echo "scale=3; ($t1 - $t0) / 1000000000" | bc)

        if [ $status -ne 0 ]; then
            elapsed="FAIL"
            all_ok=0
        fi

        times="$times $elapsed"

        # Track worst time
        if [ "$elapsed" != "FAIL" ]; then
            is_worse=$(echo "$elapsed > $worst_time" | bc -l 2>/dev/null)
            if [ "$is_worse" = "1" ]; then
                worst_time=$elapsed
            fi
        fi
    done

    ts=$(date "+%Y-%m-%d %H:%M:%S")
    if [ $all_ok -eq 1 ]; then
        echo "[$ts] size: $size | approach: $APPROACH | time: ${worst_time}s | notes: worst of 5, times:$times"
        echo "[$ts] size: $size | approach: $APPROACH | time: ${worst_time}s | notes: worst of 5, times:$times" >> "$LOG"
    else
        echo "[$ts] size: $size | approach: $APPROACH | time: FAIL | notes: some failed, times:$times"
        echo "[$ts] size: $size | approach: $APPROACH | time: FAIL | notes: some failed, times:$times" >> "$LOG"
    fi
done
