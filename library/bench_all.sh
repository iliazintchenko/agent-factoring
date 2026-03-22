#!/bin/bash
# Benchmark all sizes using YAFU with best known parameters
# Runs all 5 semiprimes per size, reports worst-case time

YAFU="/tmp/agent-factoring-4/yafu/yafu"
SEMIPRIMES="/tmp/agent-factoring-8/semiprimes.json"

# Parse digit size from argument
SIZE=$1
if [ -z "$SIZE" ]; then
    echo "Usage: $0 <digit_size> [command] [extra_args]"
    exit 1
fi

CMD="${2:-siqs}"
EXTRA="${3:-}"

# Extract semiprimes for this size using python
NUMS=$(python3 -c "
import json
with open('$SEMIPRIMES') as f:
    d = json.load(f)
nums = d.get('$SIZE', [])
for n in nums:
    print(n)
")

if [ -z "$NUMS" ]; then
    echo "No semiprimes for size $SIZE"
    exit 1
fi

MAX_TIME=0
ALL_TIMES=""
IDX=0

while IFS= read -r N; do
    WORKDIR=$(mktemp -d /tmp/yafu_XXXXXX)
    START=$(date +%s%N)
    cd $WORKDIR
    RESULT=$(echo "${CMD}(${N})" | timeout 295 $YAFU -threads 1 -seed 42 $EXTRA 2>&1)
    EXIT_CODE=$?
    END=$(date +%s%N)
    ELAPSED=$(echo "scale=3; ($END - $START) / 1000000000" | bc)
    rm -rf $WORKDIR

    if [ $EXIT_CODE -eq 124 ]; then
        echo "  [$IDX] TIMEOUT (>295s)"
        ELAPSED="TIMEOUT"
        MAX_TIME=999
    else
        # Check if factored
        FACTOR=$(echo "$RESULT" | grep "^P" | head -1)
        if [ -z "$FACTOR" ]; then
            echo "  [$IDX] FAIL: no factor found in ${ELAPSED}s"
            MAX_TIME=999
        else
            echo "  [$IDX] ${ELAPSED}s - $FACTOR"
            if (( $(echo "$ELAPSED > $MAX_TIME" | bc -l) )); then
                MAX_TIME=$ELAPSED
            fi
        fi
    fi
    ALL_TIMES="$ALL_TIMES $ELAPSED"
    IDX=$((IDX + 1))
done <<< "$NUMS"

echo "SIZE=$SIZE WORST=$MAX_TIME TIMES=$ALL_TIMES"
