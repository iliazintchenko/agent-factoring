#!/bin/bash
BINARY="./spqs"
SEMIPRIMES="semiprimes.json"
for SIZE in $@; do
    NUMS=$(python3 -c "import json; [print(n) for n in json.load(open('$SEMIPRIMES')).get('$SIZE', [])]")
    WORST=0; IDX=0; ALL_OK=1
    for N in $NUMS; do
        START=$(date +%s.%N)
        RESULT=$(timeout 295 $BINARY $N 2>&1)
        END=$(date +%s.%N)
        ELAPSED=$(python3 -c "print(f'{$END - $START:.3f}')")
        STDOUT=$(echo "$RESULT" | grep -v "^SPQS:\|^Sieve\|^LA:\|^TIMEOUT\|^  p=")
        if echo "$STDOUT" | grep -qE "^[0-9]+$"; then STATUS="OK"; else STATUS="FAIL"; ELAPSED="FAIL"; ALL_OK=0; fi
        echo "[$SIZE] #$IDX: ${ELAPSED}s $STATUS" >&2
        [ "$STATUS" = "OK" ] && WORSE=$(python3 -c "print(max($WORST, $ELAPSED))") && WORST=$WORSE
        IDX=$((IDX + 1))
    done
    [ "$ALL_OK" = "1" ] && echo "$SIZE $WORST" || echo "$SIZE FAIL"
done
