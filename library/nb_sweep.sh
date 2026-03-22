#!/bin/bash
# NB sweep on a single size
# Usage: nb_sweep.sh <size> <NB_values...>
# Example: nb_sweep.sh 74 "" "-siqsNB 8" "-siqsNB 10"

SIZE=$1
shift
YAFU="/tmp/agent-factoring-9/yafu/yafu"

NUMS=$(python3 -c "import json; d=json.load(open('/tmp/agent-factoring-9/semiprimes.json')); print(' '.join(d['$SIZE']))")

for PARAMS in "$@"; do
    echo "=== ${SIZE}d NB=[$PARAMS] ==="
    WORST=0
    for i in 0 1 2 3 4; do
        N=$(echo $NUMS | cut -d' ' -f$((i+1)))
        WORKDIR=$(mktemp -d /tmp/yafu_XXXXXX)
        START=$(date +%s%N)
        echo "siqs($N)" | timeout 295 $YAFU -threads 1 -seed 42 $PARAMS 2>/dev/null > "$WORKDIR/out" &
        BGPID=$!
        cd /tmp
        wait $BGPID
        END=$(date +%s%N)
        ELAPSED=$(echo "scale=1; ($END - $START) / 1000000000" | bc)
        FACTOR=$(grep -oP 'P\d+ = \K\d+' "$WORKDIR/out" 2>/dev/null | head -1)
        rm -rf "$WORKDIR"
        if [ -z "$FACTOR" ]; then
            echo "  [$i] TIMEOUT/FAIL"
            WORST=999
        else
            echo "  [$i] ${ELAPSED}s"
            if (( $(echo "$ELAPSED > $WORST" | bc -l) )); then
                WORST=$ELAPSED
            fi
        fi
    done
    echo "  WORST: ${WORST}s"
done
