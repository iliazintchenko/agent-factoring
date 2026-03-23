#!/bin/bash
# Run ECM benchmark across all semiprime sizes
# Runs up to 20 parallel jobs

cd "$(dirname "$0")"

BINARY="./ecm_factor"
SEMIPRIMES="../semiprimes.json"
LOGFILE="../experiments.log"
RESULTS_DIR="/tmp/ecm_results"
mkdir -p "$RESULTS_DIR"

# Get all digit sizes
SIZES=$(python3 -c "
import json
with open('$SEMIPRIMES') as f:
    d = json.load(f)
for k in sorted(d.keys(), key=int):
    print(k)
")

run_size() {
    local digits=$1
    local nums=$(python3 -c "
import json
with open('$SEMIPRIMES') as f:
    d = json.load(f)
for n in d.get('$digits', []):
    print(n)
")

    local worst_time=0
    local all_ok=1
    local results=""

    for n in $nums; do
        local start_time=$(date +%s.%N)
        local result=$(timeout 295 "$BINARY" "$n" 2>/tmp/ecm_stderr_${digits}_$$ )
        local exit_code=$?
        local stderr_out=$(cat /tmp/ecm_stderr_${digits}_$$ 2>/dev/null)
        local time_val=$(echo "$stderr_out" | grep -oP 'time=\K[0-9.]+' | tail -1)
        local ts=$(date '+%Y-%m-%d %H:%M:%S')

        if [ "$exit_code" -ne 0 ] || echo "$result" | grep -q "FAIL"; then
            echo "[$ts] size: $digits | approach: ecm | time: FAIL | notes: failed on $n" >> "$LOGFILE"
            all_ok=0
        else
            echo "[$ts] size: $digits | approach: ecm | time: ${time_val}s | notes: factored $n -> $result" >> "$LOGFILE"
            if [ "$(echo "$time_val > $worst_time" | bc -l 2>/dev/null)" = "1" ]; then
                worst_time=$time_val
            fi
        fi
        rm -f /tmp/ecm_stderr_${digits}_$$
    done

    if [ "$all_ok" = "1" ]; then
        echo "${digits}:${worst_time}" > "$RESULTS_DIR/${digits}.txt"
        echo "ECM $digits digits: worst=${worst_time}s"
    else
        echo "ECM $digits digits: FAILED"
    fi
}

# Run all sizes, up to 20 in parallel
active=0
for d in $SIZES; do
    run_size "$d" &
    active=$((active + 1))
    if [ $active -ge 20 ]; then
        wait -n 2>/dev/null || wait
        active=$((active - 1))
    fi
done
wait

echo "=== ECM Results ==="
for f in $(ls "$RESULTS_DIR"/*.txt 2>/dev/null | sort -t/ -k4 -n); do
    cat "$f"
done
