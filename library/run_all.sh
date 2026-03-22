#!/bin/bash
# Benchmark all semiprimes using msieve (SIQS) and ECM
# Runs all 5 semiprimes of each size in parallel
# Reports worst-case time per size

cd /tmp/agent-factoring-2

TOOL=${1:-msieve}  # msieve or factor
START=${2:-30}
END=${3:-100}
DEADLINE=${4:-280}
RESULTS="/tmp/bench_results_$$"
mkdir -p "$RESULTS"

echo "=== Benchmarking $TOOL, digits $START-$END, deadline=${DEADLINE}s ==="

for digits in $(seq $START $END); do
    nums=$(python3 -c "
import json
d = json.load(open('semiprimes.json'))
for n in d.get('$digits', []):
    print(n)
")
    [ -z "$nums" ] && continue

    i=0
    pids=()
    for n in $nums; do
        i=$((i+1))
        outf="$RESULTS/${digits}_${i}"
        (
            s=$(date +%s%N)
            if [ "$TOOL" = "msieve" ]; then
                result=$(timeout ${DEADLINE}s ./msieve -q -s "$RESULTS/s_${digits}_${i}.dat" "$n" 2>/dev/null)
            else
                result=$(timeout ${DEADLINE}s ./factor "$n" "$DEADLINE" 2>/dev/null)
            fi
            e=$(date +%s%N)
            t=$(echo "scale=3; ($e - $s) / 1000000000" | bc)
            # Check if factor was found
            if [ "$TOOL" = "msieve" ]; then
                # msieve outputs prime factors
                if echo "$result" | grep -q "^p[0-9]"; then
                    echo "OK $t" > "$outf"
                else
                    echo "FAIL $t" > "$outf"
                fi
            else
                if echo "$result" | grep -q "^[0-9]"; then
                    echo "OK $t" > "$outf"
                else
                    echo "FAIL $t" > "$outf"
                fi
            fi
        ) &
        pids+=($!)
    done
    for pid in "${pids[@]}"; do wait $pid; done

    max=0; ok=1
    for j in $(seq 1 5); do
        r=$(cat "$RESULTS/${digits}_${j}" 2>/dev/null)
        status=$(echo "$r" | awk '{print $1}')
        t=$(echo "$r" | awk '{print $2}')
        if [ "$status" != "OK" ]; then
            ok=0; break
        fi
        if (( $(echo "$t > $max" | bc -l) )); then max=$t; fi
    done
    if [ "$ok" = "1" ]; then
        printf "%3d digits: %8.3fs\n" $digits $max
    else
        printf "%3d digits: INCOMPLETE\n" $digits
    fi
done

rm -rf "$RESULTS"
