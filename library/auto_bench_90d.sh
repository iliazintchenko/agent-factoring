#!/bin/bash
# auto_bench_90d.sh - Automatically benchmark 90d when load drops below threshold
# Runs all 5 90d numbers sequentially with GNFS (pre-computed poly)
# Usage: ./auto_bench_90d.sh [max_load] [output_file]

MAX_LOAD="${1:-5}"
OUTFILE="${2:-/tmp/90d_benchmark_results.txt}"
YAFU="/tmp/agent-factoring-4/yafu/yafu"
SIEVER_DIR="/tmp/agent-factoring-1/yafu_mod/factor/lasieve5_64"
POLY_DIR="/tmp/agent-factoring-1/library/gnfs_polys"

NUMBERS=(
  "315547728763353682629201910311839572981703304736061985377687053517447821670148261179477299"
  "144009255506966132983765710705715101627946917437406468019915946276667323993959592666335633"
  "321739159513091262170255410982298427140058934463728910158638687584271260451816679219749747"
  "149689871933523898487085003856984468870482390037601728900188208806500682729361812646669171"
  "226492946512013642353358853937990356716879227059976248121251466901503301655504057622409363"
)
POLYS=("90d_0.poly" "90d_1.poly" "90d_2.poly" "90d_3.poly" "90d_4.poly")

echo "Waiting for load to drop below $MAX_LOAD..." >> "$OUTFILE"

# Wait for low load
while true; do
    LOAD=$(awk '{print $1}' /proc/loadavg)
    if (( $(echo "$LOAD < $MAX_LOAD" | bc -l) )); then
        echo "Load $LOAD < $MAX_LOAD - starting benchmark at $(date)" >> "$OUTFILE"
        break
    fi
    sleep 30
done

WORST_TIME=0
ALL_PASS=1

for i in 0 1 2 3 4; do
    N="${NUMBERS[$i]}"
    POLY="$POLY_DIR/${POLYS[$i]}"

    echo "=== 90d[$i] ===" >> "$OUTFILE"

    # Check load before each run
    LOAD=$(awk '{print $1}' /proc/loadavg)
    echo "Load: $LOAD" >> "$OUTFILE"

    if (( $(echo "$LOAD > $MAX_LOAD + 3" | bc -l) )); then
        echo "Load too high ($LOAD), waiting..." >> "$OUTFILE"
        sleep 60
    fi

    WORKDIR=$(mktemp -d /tmp/bench90_XXXXXX)
    cd "$WORKDIR"
    for f in "$SIEVER_DIR"/gnfs-lasieve4I*e; do ln -sf "$f" .; done
    echo "ggnfs_dir=$WORKDIR/" > yafu.ini

    if [ -f "$POLY" ]; then
        cp "$POLY" nfs.job
        echo "210000" > "nfs.job.$(hostname).last_spq0"

        START=$(date +%s)
        RESULT=$(echo "nfs($N)" | timeout 295 "$YAFU" -threads 1 -seed 42 -xover 85 -R 2>&1)
        END=$(date +%s)
        WALL=$((END - START))
    else
        # Fallback to SIQS
        START=$(date +%s)
        RESULT=$(echo "siqs($N)" | timeout 295 "$YAFU" -threads 1 -seed 42 -siqsNB 20 -siqsB 120000 2>&1)
        END=$(date +%s)
        WALL=$((END - START))
    fi

    if echo "$RESULT" | grep -q "prp\|^P[0-9]"; then
        echo "PASS: ${WALL}s" >> "$OUTFILE"
        echo "$RESULT" | grep "prp\|Factorization" >> "$OUTFILE"
        if [ $WALL -gt $WORST_TIME ]; then WORST_TIME=$WALL; fi
    else
        echo "FAIL: timeout at ${WALL}s" >> "$OUTFILE"
        ALL_PASS=0
    fi

    rm -rf "$WORKDIR"
done

echo "" >> "$OUTFILE"
echo "=== SUMMARY ===" >> "$OUTFILE"
echo "All pass: $ALL_PASS" >> "$OUTFILE"
echo "Worst time: ${WORST_TIME}s" >> "$OUTFILE"
echo "Finished at $(date)" >> "$OUTFILE"
