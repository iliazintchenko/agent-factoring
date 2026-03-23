#!/bin/bash
# Run a factoring binary against all semiprimes of a given size.
# Usage: ./run_bench.sh <binary> <approach_name> [sizes...]
# Example: ./run_bench.sh ./ecm_factor ecm 30 40 50 60
# Records results to experiments.log and prints worst-case times for algo-scaling.json

BINARY="$1"
APPROACH="$2"
shift 2
SIZES="$@"

if [ -z "$BINARY" ] || [ -z "$APPROACH" ]; then
    echo "Usage: $0 <binary> <approach_name> [sizes...]"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
SEMIPRIMES="$REPO_DIR/semiprimes.json"
LOG="$REPO_DIR/experiments.log"

for SIZE in $SIZES; do
    # Extract semiprimes for this size using python
    NUMS=$(python3 -c "
import json
with open('$SEMIPRIMES') as f:
    data = json.load(f)
for n in data.get('$SIZE', []):
    print(n)
")

    WORST=0
    ALL_OK=true

    for NUM in $NUMS; do
        START=$(date +%s.%N)
        RESULT=$(timeout 295 "$BINARY" "$NUM" 2>/tmp/bench_stderr_$$ )
        EXIT_CODE=$?
        END=$(date +%s.%N)
        ELAPSED=$(python3 -c "print(f'{$END - $START:.3f}')")

        STDERR=$(cat /tmp/bench_stderr_$$ 2>/dev/null)
        rm -f /tmp/bench_stderr_$$

        TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')

        if echo "$RESULT" | grep -q "^FACTOR:"; then
            FACTOR=$(echo "$RESULT" | grep "^FACTOR:" | head -1 | awk '{print $2}')
            echo "[$TIMESTAMP] size: $SIZE | approach: $APPROACH | time: ${ELAPSED}s | notes: found factor $FACTOR" >> "$LOG"

            # Update worst case
            if python3 -c "exit(0 if float('$ELAPSED') > float('$WORST') else 1)"; then
                WORST=$ELAPSED
            fi
        else
            echo "[$TIMESTAMP] size: $SIZE | approach: $APPROACH | time: FAIL | notes: timeout or no factor found for $NUM" >> "$LOG"
            ALL_OK=false
            WORST="FAIL"
        fi
    done

    if [ "$ALL_OK" = true ]; then
        echo "SIZE=$SIZE WORST=$WORST"
    else
        echo "SIZE=$SIZE WORST=FAIL"
    fi
done
