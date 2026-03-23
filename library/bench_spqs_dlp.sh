#!/bin/bash
# Benchmark SPQS-DLP across all sizes, all 5 semiprimes per size
# Usage: bash library/bench_spqs_dlp.sh [start_digits] [end_digits]
# Output: worst-case time per digit size

START=${1:-30}
END=${2:-80}
BINARY=./spqs_dlp
SEMIPRIMES=semiprimes.json

if [ ! -f "$BINARY" ]; then
    echo "Compiling spqs_dlp..."
    gcc -O3 -march=native -o spqs_dlp library/spqs_dlp.c -lgmp -lm || exit 1
fi

echo "Benchmarking SPQS-DLP from ${START}d to ${END}d"
echo "================================================"

for sz in $(seq $START $END); do
    NUMS=$(python3 -c "import json; d=json.load(open('$SEMIPRIMES')); print('\n'.join(d.get('$sz',[])))")
    if [ -z "$NUMS" ]; then continue; fi

    worst=0
    all_ok=1
    times=""
    idx=0

    while IFS= read -r n; do
        idx=$((idx + 1))
        result=$(timeout 295 $BINARY "$n" 2>/tmp/spqs_dlp_err_${sz}_${idx}.txt)

        if [ $? -ne 0 ] || [ "$result" = "FAIL" ] || [ -z "$result" ]; then
            t="FAIL"
            all_ok=0
        else
            # Extract time from stderr
            t=$(grep -oP 'factored \d+d in \K[0-9.]+' /tmp/spqs_dlp_err_${sz}_${idx}.txt)
            if [ -z "$t" ]; then
                t=$(grep -oP 'in \K[0-9.]+(?=s)' /tmp/spqs_dlp_err_${sz}_${idx}.txt | tail -1)
            fi
            if [ -z "$t" ]; then t="FAIL"; all_ok=0; fi
        fi

        if [ "$t" != "FAIL" ]; then
            cmp=$(python3 -c "print(1 if float('$t') > float('$worst') else 0)")
            if [ "$cmp" = "1" ]; then worst=$t; fi
        fi

        times="$times $t"
    done <<< "$NUMS"

    if [ "$all_ok" = "1" ]; then
        echo "${sz}d: worst=${worst}s  [${times}]"
    else
        echo "${sz}d: PARTIAL worst=${worst}s  [${times}]"
    fi
done
