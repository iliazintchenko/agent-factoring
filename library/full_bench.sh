#!/bin/bash
# Full benchmark: test all sizes 30-89 with best known params
# Usage: ./full_bench.sh [start_size] [end_size]
START=${1:-30}
END=${2:-89}
YAFU="/tmp/agent-factoring-4/yafu/yafu"
SP="/tmp/agent-factoring-4/semiprimes.json"

echo "Benchmarking sizes $START to $END"
echo "============================================"

for SIZE in $(seq $START $END); do
    CMD="siqs"
    EXTRA=""
    
    # Use smallmpqs for 30-44d
    if [ $SIZE -le 44 ]; then CMD="smallmpqs"; fi
    
    # NB tuning based on size
    if [ $SIZE -ge 73 ] && [ $SIZE -le 79 ]; then
        if [ $SIZE -le 73 ]; then EXTRA="-siqsNB 10"
        else EXTRA="-siqsNB 11"; fi
    elif [ $SIZE -ge 80 ] && [ $SIZE -le 81 ]; then
        EXTRA="-siqsNB 12"
    elif [ $SIZE -ge 82 ] && [ $SIZE -le 84 ]; then
        EXTRA="-siqsNB 14"
    elif [ $SIZE -ge 85 ]; then
        EXTRA="-siqsNB 18"
        case $SIZE in
            85) EXTRA="$EXTRA -siqsB 70000" ;;
            86) EXTRA="$EXTRA -siqsB 80000" ;;
            87) EXTRA="$EXTRA -siqsB 90000" ;;
            88) EXTRA="$EXTRA -siqsB 80000" ;;
            89) EXTRA="$EXTRA -siqsB 100000" ;;
        esac
    fi
    
    WORST=0
    ALL_OK=true
    for i in 0 1 2 3 4; do
        N=$(python3 -c "import json; print(json.load(open('$SP'))['$SIZE'][$i])")
        WORKDIR=$(mktemp -d /tmp/yafu_fb_XXXXXX)
        cd $WORKDIR
        START_T=$(date +%s.%N)
        echo "$CMD($N)" | timeout 295 $YAFU -threads 1 -seed 42 $EXTRA 2>&1 | grep -q "^P" || ALL_OK=false
        END_T=$(date +%s.%N)
        ELAPSED=$(echo "$END_T - $START_T" | bc)
        if (( $(echo "$ELAPSED > $WORST" | bc -l) )); then WORST=$ELAPSED; fi
        cd /tmp && rm -rf $WORKDIR
    done
    
    echo "${SIZE}d: worst=${WORST}s all_ok=${ALL_OK} params='$CMD $EXTRA'"
done
