#!/bin/bash
# Benchmark hyper_siqs across all digit sizes
# Usage: bash library/bench_hyper.sh

BINARY="./hyper_siqs"
SEMIPRIMES="semiprimes.json"

if [ ! -f "$BINARY" ]; then
    echo "Compiling hyper_siqs..."
    gcc -O3 -march=native -o hyper_siqs library/hyper_siqs.c -lgmp -lm
fi

for SIZE in 30 35 40 45 50 55 60 65 70 75 80; do
    NUMS=$(python3 -c "import json; sp=json.load(open('$SEMIPRIMES')); print(' '.join(sp.get('$SIZE',[])))")
    if [ -z "$NUMS" ]; then continue; fi

    echo "=== ${SIZE} digits ==="
    WORST=0
    ALL_OK=1
    IDX=0
    for N in $NUMS; do
        IDX=$((IDX+1))
        RESULT=$(timeout 295 $BINARY $N 2>/tmp/hyper_bench_err_${SIZE}_${IDX}.txt)
        TIME=$(grep -oP 'factored.*in \K[0-9.]+' /tmp/hyper_bench_err_${SIZE}_${IDX}.txt 2>/dev/null || echo "FAIL")
        if [ "$RESULT" = "FAIL" ] || [ -z "$RESULT" ] || [ "$TIME" = "FAIL" ]; then
            echo "  [$IDX] FAIL on $N"
            ALL_OK=0
            WORST=999
        else
            echo "  [$IDX] ${TIME}s → factor=$RESULT"
            if [ $(echo "$TIME > $WORST" | bc -l) -eq 1 ]; then
                WORST=$TIME
            fi
        fi
    done
    if [ "$ALL_OK" = "1" ]; then
        echo "  WORST: ${WORST}s"
    fi
    echo ""
done
