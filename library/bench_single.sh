#!/bin/bash
# Benchmark YAFU SIQS single-threaded for a given digit size
# Usage: bash bench_single.sh <digit_size> [start_index] [count]
# Runs 5 semiprimes in parallel, each on 1 thread, reports worst-case time

DIGIT_SIZE=$1
START=${2:-0}
COUNT=${3:-5}
YAFU="/tmp/agent-factoring-1/yafu/yafu"
SEMIPRIMES="/tmp/agent-factoring-1/semiprimes.json"

# Extract the semiprimes for this size
NUMBERS=$(python3 -c "
import json
with open('$SEMIPRIMES') as f:
    data = json.load(f)
nums = data.get('$DIGIT_SIZE', [])
for i in range($START, min($START+$COUNT, len(nums))):
    print(nums[i])
")

if [ -z "$NUMBERS" ]; then
    echo "ERROR: No semiprimes found for size $DIGIT_SIZE"
    exit 1
fi

MAX_TIME=0
ALL_OK=true
TIMES=""
PIDS=()
TMPFILES=()

# Launch all in parallel
for N in $NUMBERS; do
    TMPFILE=$(mktemp /tmp/yafu_result_XXXXXX)
    TMPFILES+=("$TMPFILE")
    WORKDIR=$(mktemp -d /tmp/yafu_bench_XXXXXX)
    (
        cd "$WORKDIR"
        START_T=$(date +%s%N)
        echo "siqs($N)" | timeout 290 "$YAFU" -threads 1 -seed 42 > "$TMPFILE.out" 2>&1
        RC=$?
        END_T=$(date +%s%N)
        ELAPSED=$(echo "scale=3; ($END_T - $START_T) / 1000000000" | bc)
        if [ $RC -ne 0 ] && ! grep -q "factors found" "$TMPFILE.out"; then
            echo "FAIL $ELAPSED $N" > "$TMPFILE"
        else
            echo "OK $ELAPSED $N" > "$TMPFILE"
        fi
        rm -rf "$WORKDIR" "$TMPFILE.out"
    ) &
    PIDS+=($!)
done

# Wait for all
for PID in "${PIDS[@]}"; do
    wait $PID
done

# Collect results
for TMPFILE in "${TMPFILES[@]}"; do
    if [ -f "$TMPFILE" ]; then
        read STATUS ELAPSED NUM < "$TMPFILE"
        if [ "$STATUS" = "FAIL" ]; then
            echo "  FAIL: $NUM (${ELAPSED}s)"
            ALL_OK=false
        else
            echo "  OK: $NUM (${ELAPSED}s)"
            if (( $(echo "$ELAPSED > $MAX_TIME" | bc -l) )); then
                MAX_TIME=$ELAPSED
            fi
        fi
        TIMES="$TIMES $ELAPSED"
        rm -f "$TMPFILE"
    fi
done

if [ "$ALL_OK" = true ]; then
    echo "SIZE=$DIGIT_SIZE WORST=$MAX_TIME ALL_TIMES=$TIMES"
else
    echo "SIZE=$DIGIT_SIZE FAIL TIMES=$TIMES"
fi
