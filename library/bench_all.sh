#!/bin/bash
# Run MPQS2 benchmark across all sizes, output JSON-compatible results
cd "$(dirname "$0")"

SEMIPRIMES="../semiprimes.json"
LOG="../experiments.log"
RESULTS=""

for digits in $(seq 30 100); do
    nums=$(python3 -c "import json; d=json.load(open('$SEMIPRIMES')); print('\n'.join(d.get('$digits', [])))" 2>/dev/null)
    if [ -z "$nums" ]; then continue; fi

    worst_time=0
    all_ok=1

    for n in $nums; do
        start_time=$(date +%s.%N)
        result=$(timeout 295 ./mpqs2 "$n" 2>/dev/null)
        end_time=$(date +%s.%N)
        elapsed=$(echo "$end_time - $start_time" | bc)

        if [ -z "$result" ]; then
            echo "$digits: FAIL on $n (timeout)" >&2
            elapsed_str="FAIL"
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] size: $digits | approach: mpqs2 | time: FAIL | notes: timeout on $n" >> "$LOG"
            all_ok=0
            break
        fi

        p1=$(echo "$result" | awk '{print $1}')
        p2=$(echo "$result" | awk '{print $2}')
        check=$(python3 -c "print($p1 * $p2)" 2>/dev/null)
        if [ "$check" != "$n" ]; then
            echo "$digits: WRONG on $n" >&2
            all_ok=0
            break
        fi

        cmp=$(echo "$elapsed > $worst_time" | bc -l)
        if [ "$cmp" -eq 1 ]; then worst_time=$elapsed; fi
    done

    if [ "$all_ok" -eq 1 ]; then
        worst_time=$(printf "%.3f" "$worst_time")
        echo "$digits: ${worst_time}s"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] size: $digits | approach: mpqs2 | time: $worst_time | notes: worst-case of 5 semiprimes" >> "$LOG"
    else
        echo "$digits: STOPPED"
        break
    fi
done
