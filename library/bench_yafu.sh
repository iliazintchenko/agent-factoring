#!/bin/bash
# Benchmark YAFU on all semiprimes
# Runs 5 per size in parallel, each with THREADS threads

cd /tmp/agent-factoring-2
START=${1:-30}
END=${2:-100}
DEADLINE=${3:-280}
THREADS=${4:-9}
RESULTS="/tmp/yafu_bench_$$"
mkdir -p "$RESULTS"

echo "=== YAFU benchmark, $START-$END digits, ${THREADS}t per instance, ${DEADLINE}s deadline ==="

for digits in $(seq $START $END); do
    nums=$(python3 -c "
import json
d = json.load(open('semiprimes.json'))
for n in d.get('$digits', []):
    print(n)
")
    [ -z "$nums" ] && continue

    i=0; pids=()
    for n in $nums; do
        i=$((i+1))
        outf="$RESULTS/${digits}_${i}"
        (
            s=$(date +%s%N)
            result=$(echo "factor($n)" | timeout "${DEADLINE}s" ./yafu_bin -threads "$THREADS" 2>/dev/null)
            e=$(date +%s%N)
            t=$(echo "scale=3; ($e - $s) / 1000000000" | bc)
            if echo "$result" | grep -q "^P[0-9]"; then
                echo "OK $t" > "$outf"
            else
                echo "FAIL $t" > "$outf"
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
        if [ "$status" != "OK" ]; then ok=0; break; fi
        if (( $(echo "$t > $max" | bc -l) )); then max=$t; fi
    done
    if [ "$ok" = "1" ]; then
        printf "%3d digits: %8.3fs\n" $digits $max
    else
        printf "%3d digits: INCOMPLETE\n" $digits
    fi
done

rm -rf "$RESULTS"
