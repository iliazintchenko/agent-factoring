#!/bin/bash
# Benchmark YAFU SIQS across all sizes, single-threaded, seed=42
# Outputs worst-case time per size

YAFU="./yafu/yafu"
SEMIPRIMES="semiprimes.json"
LOGFILE="yafu_bench.log"

> "$LOGFILE"

for SIZE in $@; do
    # Extract semiprimes for this size using python
    NUMS=$(python3 -c "
import json
with open('$SEMIPRIMES') as f:
    d = json.load(f)
nums = d.get('$SIZE', [])
for n in nums:
    print(n)
")

    WORST=0
    IDX=0
    for N in $NUMS; do
        START=$(date +%s.%N)
        RESULT=$(echo "factor($N)" | timeout 295 $YAFU -threads 1 2>&1)
        END=$(date +%s.%N)
        ELAPSED=$(python3 -c "print(f'{$END - $START:.3f}')")

        # Check if it succeeded
        if echo "$RESULT" | grep -q "P[0-9]"; then
            STATUS="OK"
        else
            STATUS="FAIL"
            ELAPSED="FAIL"
        fi

        echo "[$SIZE] idx=$IDX time=${ELAPSED}s $STATUS" | tee -a "$LOGFILE"

        if [ "$STATUS" = "OK" ]; then
            WORSE=$(python3 -c "print(max($WORST, $ELAPSED))")
            WORST=$WORSE
        fi

        IDX=$((IDX + 1))
    done

    echo "[$SIZE] WORST: ${WORST}s" | tee -a "$LOGFILE"
done
