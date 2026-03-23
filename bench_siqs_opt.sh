#!/bin/bash
# Benchmark siqs_opt across digit sizes
# Run all 5 semiprimes per size, record worst-case time
set -e

BINARY=./siqs_opt
SEMIPRIMES=semiprimes.json
OUTDIR=/tmp/siqs_opt_results
mkdir -p "$OUTDIR"

run_one() {
    local size=$1 idx=$2 N=$3
    local outfile="$OUTDIR/${size}d_${idx}.txt"
    local start_time=$(date +%s%N)
    local result
    result=$(timeout 295 $BINARY "$N" 2>"$OUTDIR/${size}d_${idx}.err")
    local exit_code=$?
    local end_time=$(date +%s%N)
    local elapsed=$(echo "scale=3; ($end_time - $start_time) / 1000000000" | bc)

    if [ $exit_code -eq 0 ] && [ -n "$result" ]; then
        echo "$elapsed OK $result" > "$outfile"
    elif [ $exit_code -eq 124 ]; then
        echo "TIMEOUT" > "$outfile"
    else
        echo "FAIL" > "$outfile"
    fi
}

# Run benchmarks for sizes 30-85 in batches
for size in $(seq 30 5 85); do
    echo "=== Benchmarking ${size}d ==="
    # Get the 5 semiprimes for this size
    NUMS=$(python3 -c "import json; d=json.load(open('$SEMIPRIMES')); [print(n) for n in d.get('$size', [])]")

    idx=0
    pids=()
    for N in $NUMS; do
        run_one $size $idx "$N" &
        pids+=($!)
        idx=$((idx + 1))
    done

    # Wait for all 5
    for pid in "${pids[@]}"; do
        wait $pid
    done

    # Report results
    worst=0
    all_ok=1
    for i in 0 1 2 3 4; do
        f="$OUTDIR/${size}d_${i}.txt"
        if [ -f "$f" ]; then
            content=$(cat "$f")
            if [[ "$content" == TIMEOUT* ]] || [[ "$content" == FAIL* ]]; then
                echo "  ${size}d[$i]: $content"
                all_ok=0
            else
                t=$(echo "$content" | awk '{print $1}')
                echo "  ${size}d[$i]: ${t}s"
                if (( $(echo "$t > $worst" | bc -l) )); then
                    worst=$t
                fi
            fi
        fi
    done

    if [ "$all_ok" = "1" ]; then
        echo "  WORST: ${worst}s"
    fi
    echo ""
done
