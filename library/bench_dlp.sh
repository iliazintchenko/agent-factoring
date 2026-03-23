#!/bin/bash
# Benchmark DLP-SIQS across sizes, single-threaded, seed=42
# Usage: bench_dlp.sh <size1> <size2> ...

BINARY="./dlp_siqs"
SEMIPRIMES="semiprimes.json"

for SIZE in $@; do
    NUMS=$(python3 -c "
import json
with open('$SEMIPRIMES') as f:
    d = json.load(f)
for n in d.get('$SIZE', []):
    print(n)
")

    WORST=0
    IDX=0
    ALL_OK=1
    for N in $NUMS; do
        START=$(date +%s.%N)
        RESULT=$(timeout 295 $BINARY $N 2>&1)
        END=$(date +%s.%N)
        ELAPSED=$(python3 -c "print(f'{$END - $START:.3f}')")

        STDOUT=$(echo "$RESULT" | grep -v "^DLP-SIQS\|^FB:\|^M=\|^  p=\|^Sieve\|^After\|^LA:\|^Found\|^TIMEOUT")
        if echo "$STDOUT" | grep -qE "^[0-9]+$"; then
            STATUS="OK"
        else
            STATUS="FAIL"
            ELAPSED="FAIL"
            ALL_OK=0
        fi

        echo "[$SIZE] #$IDX: ${ELAPSED}s $STATUS" >&2

        if [ "$STATUS" = "OK" ]; then
            WORSE=$(python3 -c "print(max($WORST, $ELAPSED))")
            WORST=$WORSE
        fi
        IDX=$((IDX + 1))
    done

    if [ "$ALL_OK" = "1" ]; then
        echo "$SIZE $WORST"
    else
        echo "$SIZE FAIL"
    fi
done
