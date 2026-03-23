#!/bin/bash
BINARY=./turbo_siqs2
SEMIPRIMES=semiprimes.json
OUTDIR=/tmp/turbo2_results
mkdir -p "$OUTDIR"
run_one() {
    local size=$1 idx=$2 N=$3
    local start_time=$(date +%s%N)
    result=$(timeout 295 $BINARY "$N" 2>"$OUTDIR/${size}d_${idx}.err")
    local exit_code=$?
    local end_time=$(date +%s%N)
    local elapsed=$(echo "scale=3; ($end_time - $start_time) / 1000000000" | bc)
    if [ $exit_code -eq 0 ] && [ -n "$result" ]; then echo "$elapsed" > "$OUTDIR/${size}d_${idx}.txt"
    elif [ $exit_code -eq 124 ]; then echo "TIMEOUT" > "$OUTDIR/${size}d_${idx}.txt"
    else echo "FAIL" > "$OUTDIR/${size}d_${idx}.txt"; fi
}
for size in "$@"; do
    echo "=== ${size}d ==="
    NUMS=$(python3 -c "import json; d=json.load(open('$SEMIPRIMES')); [print(n) for n in d.get('$size', [])]")
    idx=0; pids=()
    for N in $NUMS; do run_one $size $idx "$N" & pids+=($!); idx=$((idx+1)); done
    for pid in "${pids[@]}"; do wait $pid; done
    worst=0; all_ok=1
    for i in 0 1 2 3 4; do
        f="$OUTDIR/${size}d_${i}.txt"; [ -f "$f" ] || continue
        content=$(cat "$f")
        if [[ "$content" == TIMEOUT* ]] || [[ "$content" == FAIL* ]]; then echo "  [$i]: $content"; all_ok=0
        else echo "  [$i]: ${content}s"; (( $(echo "$content > $worst" | bc -l) )) && worst=$content; fi
    done
    [ "$all_ok" = "1" ] && echo "  WORST: ${worst}s"
    echo ""
done
