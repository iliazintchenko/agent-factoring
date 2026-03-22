#!/bin/bash
# Comprehensive 90d testing script
# Tests SIQS (original YAFU) on all 5 90d semiprimes
# Run when machine load is < 5 for reliable results

YAFU="/tmp/agent-factoring-9/yafu/yafu"
LOGFILE="/tmp/agent-factoring-9/test_90d_results_$(date +%Y%m%d_%H%M%S).txt"

NUMS=(
"315547728763353682629201910311839572981703304736061985377687053517447821670148261179477299"
"144009255506966132983765710705715101627946917437406468019915946276667323993959592666335633"
"321739159513091262170255410982298427140058934463728910158638687584271260451816679219749747"
"149689871933523898487085003856984468870482390037601728900188208806500682729361812646669171"
"226492946512013642353358853937990356716879227059976248121251466901503301655504057622409363"
)

echo "=== 90d SIQS Test Run $(date) ===" | tee $LOGFILE
echo "Load at start: $(uptime)" | tee -a $LOGFILE
echo "" | tee -a $LOGFILE

WORST=0
ALL_PASS=true

for i in "${!NUMS[@]}"; do
    N="${NUMS[$i]}"
    WORKDIR=$(mktemp -d /tmp/yafu_XXXXXX)
    cd $WORKDIR
    echo -n "90d[$i]: " | tee -a $LOGFILE

    result=$(echo "siqs($N)" | LD_LIBRARY_PATH=/usr/local/lib timeout 295 $YAFU -threads 1 -seed 42 -siqsNB 20 -siqsB 120000 2>&1)
    elapsed=$(echo "$result" | grep "elapsed time" | head -1 | grep -oP '[\d.]+')
    factor=$(echo "$result" | grep "^P" | head -1)
    rels_line=$(echo "$result" | grep "rels found" | tail -1)

    if [ -z "$elapsed" ]; then
        echo "TIMEOUT (>295s) | $rels_line" | tee -a $LOGFILE
        ALL_PASS=false
    else
        echo "${elapsed}s | $factor" | tee -a $LOGFILE
        # Track worst case
        COMP=$(python3 -c "print(1 if $elapsed > $WORST else 0)")
        if [ "$COMP" = "1" ]; then
            WORST=$elapsed
        fi
    fi
    rm -rf $WORKDIR
done

echo "" | tee -a $LOGFILE
echo "Worst case: ${WORST}s" | tee -a $LOGFILE
echo "All pass: $ALL_PASS" | tee -a $LOGFILE
echo "Load at end: $(uptime)" | tee -a $LOGFILE
