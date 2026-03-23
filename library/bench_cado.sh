#!/bin/bash
# Benchmark CADO-NFS across sizes, single-core via taskset
# Usage: bench_cado.sh <size1> <size2> ...
# NOTE: CADO internally uses threads for some steps, so we use taskset -c 0
# and report the total CPU time (not wallclock)

CADO="/tmp/agent-factoring-9/cado-nfs/cado-nfs.py"
SEMIPRIMES="/tmp/agent-factoring-9/semiprimes.json"
cd /tmp/agent-factoring-9/cado-nfs

for SIZE in $@; do
    NUMS=$(python3 -c "import json; [print(n) for n in json.load(open('$SEMIPRIMES')).get('$SIZE', [])]")
    WORST=0; IDX=0; ALL_OK=1
    for N in $NUMS; do
        START=$(date +%s.%N)
        RESULT=$(taskset -c 0 timeout 295 python3 $CADO $N --server-threads 1 --client-threads 1 2>&1)
        END=$(date +%s.%N)
        ELAPSED=$(python3 -c "print(f'{$END - $START:.3f}')")

        # Check for factor in output
        if echo "$RESULT" | grep -qE "^[0-9]+ [0-9]+$"; then
            STATUS="OK"
            # Extract CPU time from CADO output if available
            CPU_TIME=$(echo "$RESULT" | grep "Total cpu/elapsed time" | grep -oP "[\d.]+" | head -1)
            echo "[$SIZE] #$IDX: wall=${ELAPSED}s cpu=${CPU_TIME}s $STATUS" >&2
        else
            STATUS="FAIL"
            ELAPSED="FAIL"
            ALL_OK=0
            echo "[$SIZE] #$IDX: $STATUS" >&2
        fi

        if [ "$STATUS" = "OK" ]; then
            WORSE=$(python3 -c "print(max($WORST, $ELAPSED))")
            WORST=$WORSE
        fi
        IDX=$((IDX + 1))
    done
    [ "$ALL_OK" = "1" ] && echo "$SIZE $WORST" || echo "$SIZE FAIL"
done
