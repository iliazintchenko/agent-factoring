#!/bin/bash
# Benchmark factor on all semiprimes, running all 5 of each size in parallel
# Usage: ./bench.sh [start_digits] [end_digits] [deadline_per_number]

START=${1:-30}
END=${2:-70}
DEADLINE=${3:-55}
RESULTS_DIR="/tmp/factor_bench_$$"
mkdir -p "$RESULTS_DIR"

echo "Benchmarking digits $START to $END, deadline=${DEADLINE}s per number"

for digits in $(seq $START $END); do
    nums=$(python3 -c "
import json
d = json.load(open('semiprimes.json'))
for n in d.get('$digits', []):
    print(n)
")
    if [ -z "$nums" ]; then continue; fi

    echo "--- $digits digits ---"
    i=0
    pids=()
    for n in $nums; do
        i=$((i+1))
        outfile="$RESULTS_DIR/${digits}_${i}.txt"
        (
            start_t=$(date +%s%N)
            result=$(timeout ${DEADLINE}s ./factor "$n" "$DEADLINE" 2>&1)
            end_t=$(date +%s%N)
            elapsed=$(echo "scale=3; ($end_t - $start_t) / 1000000000" | bc)
            if echo "$result" | grep -q "^[0-9]"; then
                factor=$(echo "$result" | head -1)
                echo "OK $elapsed $factor" > "$outfile"
            else
                echo "FAIL $elapsed" > "$outfile"
            fi
        ) &
        pids+=($!)
    done
    for pid in "${pids[@]}"; do wait $pid; done

    max_time=0
    all_ok=1
    for j in $(seq 1 5); do
        r=$(cat "$RESULTS_DIR/${digits}_${j}.txt" 2>/dev/null)
        status=$(echo "$r" | awk '{print $1}')
        t=$(echo "$r" | awk '{print $2}')
        if [ "$status" = "OK" ]; then
            printf "  #%d: %.3fs\n" $j $t
            if (( $(echo "$t > $max_time" | bc -l) )); then max_time=$t; fi
        else
            printf "  #%d: FAIL\n" $j
            all_ok=0
        fi
    done
    if [ "$all_ok" = "1" ]; then
        printf "  worst: %.3fs\n" $max_time
    else
        echo "  INCOMPLETE"
    fi
done

rm -rf "$RESULTS_DIR"
