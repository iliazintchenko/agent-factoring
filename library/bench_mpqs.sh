#!/bin/bash
# Benchmark MPQS across semiprime sizes
# Runs 5 semiprimes per size, records worst case
cd "$(dirname "$0")"

BINARY="./mpqs"
SEMIPRIMES="../semiprimes.json"
LOGFILE="../experiments.log"
RESULTS_DIR="/tmp/mpqs_results"
mkdir -p "$RESULTS_DIR"

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

    for n in $nums; do
        local result=$(timeout 295 "$BINARY" "$n" 2>/tmp/mpqs_err_${digits}_$$ )
        local exit_code=$?
        local stderr_out=$(cat /tmp/mpqs_err_${digits}_$$ 2>/dev/null)
        local time_val=$(echo "$stderr_out" | grep -oP 'time=\K[0-9.]+' | tail -1)
        local ts=$(date '+%Y-%m-%d %H:%M:%S')

        if [ "$exit_code" -ne 0 ] || echo "$result" | grep -q "FAIL"; then
            echo "[$ts] size: $digits | approach: mpqs | time: FAIL | notes: failed on $n" >> "$LOGFILE"
            all_ok=0
        else
            echo "[$ts] size: $digits | approach: mpqs | time: ${time_val}s | notes: factored" >> "$LOGFILE"
            if [ -n "$time_val" ] && [ "$(echo "$time_val > $worst_time" | bc -l 2>/dev/null)" = "1" ]; then
                worst_time=$time_val
            fi
        fi
        rm -f /tmp/mpqs_err_${digits}_$$
    done

    if [ "$all_ok" = "1" ] && [ -n "$worst_time" ]; then
        echo "${digits}:${worst_time}" > "$RESULTS_DIR/${digits}.txt"
        echo "MPQS $digits digits: worst=${worst_time}s"
    else
        echo "MPQS $digits digits: FAILED (or incomplete)"
        echo "${digits}:FAIL" > "$RESULTS_DIR/${digits}.txt"
    fi
}

# Run sizes 30-70 in parallel (up to 10 at once)
active=0
for d in $(seq 30 70); do
    # Check if this size exists in semiprimes.json
    exists=$(python3 -c "
import json
with open('$SEMIPRIMES') as f:
    d = json.load(f)
print('yes' if str($d) in d else 'no')
")
    [ "$exists" != "yes" ] && continue

    run_size "$d" &
    active=$((active + 1))
    if [ $active -ge 10 ]; then
        wait -n 2>/dev/null || wait
        active=$((active - 1))
    fi
done
wait

echo "=== MPQS Results ==="
for f in $(ls "$RESULTS_DIR"/*.txt 2>/dev/null | sort -t/ -k4 -n); do
    cat "$f"
done
