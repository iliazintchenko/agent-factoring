#!/bin/bash
# Auto-test 90d semiprimes when load is low
# Runs continuously, testing one number at a time when load < threshold
# Usage: bash library/auto_90d.sh

YAFU_MOD=/tmp/agent-factoring-6/yafu_mod/yafu
YAFU_GNFS=/tmp/agent-factoring-4/yafu/yafu
GGNFS_DIR=/tmp/agent-factoring-4/yafu/factor/lasieve5_64/bin/
LOAD_THRESHOLD=8
LOGFILE=/tmp/agent-factoring-6/experiments.log
SEED=42
TIMEOUT=295

declare -a NUMS=(
  "315547728763353682629201910311839572981703304736061985377687053517447821670148261179477299"
  "144009255506966132983765710705715101627946917437406468019915946276667323993959592666335633"
  "321739159513091262170255410982298427140058934463728910158638687584271260451816679219749747"
  "149689871933523898487085003856984468870482390037601728900188208806500682729361812646669171"
  "226492946512013642353358853937990356716879227059976248121251466901503301655504057622409363"
)

wait_for_low_load() {
    while true; do
        LOAD=$(awk '{print $1}' /proc/loadavg)
        if (( $(echo "$LOAD < $LOAD_THRESHOLD" | bc -l) )); then
            return
        fi
        sleep 15
    done
}

test_siqs() {
    local idx=$1
    local N=${NUMS[$idx]}
    local WORKDIR=$(mktemp -d /tmp/auto_siqs_XXXXXX)
    cd $WORKDIR

    local LOAD=$(awk '{print $1}' /proc/loadavg)
    local START=$SECONDS
    timeout $TIMEOUT $YAFU_MOD -threads 1 -seed $SEED -siqsNB 20 -siqsB 120000 <<< "siqs($N)" > output.txt 2>&1
    local EXIT=$?
    local ELAPSED=$((SECONDS - START))

    if [ $EXIT -eq 0 ] && grep -q "factors found" output.txt; then
        local TIME=$(grep "elapsed time" output.txt | grep -oP '[\d.]+')
        echo "SIQS 90d[$idx] PASS: ${TIME}s (load=$LOAD)"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] size: 90 | approach: yafu_mod SIQS NB=20 B=120K 90d[$idx] (load $LOAD) | time: $TIME | notes: AUTO TEST PASS" >> $LOGFILE
        rm -rf $WORKDIR
        return 0
    else
        echo "SIQS 90d[$idx] FAIL: ${ELAPSED}s (load=$LOAD)"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] size: 90 | approach: yafu_mod SIQS NB=20 B=120K 90d[$idx] (load $LOAD) | time: FAIL ($ELAPSED) | notes: AUTO TEST" >> $LOGFILE
        rm -rf $WORKDIR
        return 1
    fi
}

test_gnfs() {
    local idx=$1
    local N=${NUMS[$idx]}
    local WORKDIR=$(mktemp -d /tmp/auto_gnfs_XXXXXX)
    cd $WORKDIR
    cat > yafu.ini << EOF
ggnfs_dir=$GGNFS_DIR
EOF

    local LOAD=$(awk '{print $1}' /proc/loadavg)
    local START=$SECONDS
    timeout $TIMEOUT $YAFU_GNFS -threads 1 -seed $SEED -xover 85 <<< "nfs($N)" > output.txt 2>&1
    local EXIT=$?
    local ELAPSED=$((SECONDS - START))

    if [ $EXIT -eq 0 ] && grep -q "factors found" output.txt; then
        echo "GNFS 90d[$idx] PASS: ${ELAPSED}s (load=$LOAD)"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] size: 90 | approach: YAFU GNFS 90d[$idx] (load $LOAD) | time: $ELAPSED | notes: AUTO TEST PASS" >> $LOGFILE
        rm -rf $WORKDIR
        return 0
    else
        echo "GNFS 90d[$idx] FAIL: ${ELAPSED}s (load=$LOAD)"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] size: 90 | approach: YAFU GNFS 90d[$idx] (load $LOAD) | time: FAIL ($ELAPSED) | notes: AUTO TEST" >> $LOGFILE
        rm -rf $WORKDIR
        return 1
    fi
}

echo "=== Auto 90d Test Starting at $(date) ==="
echo "Load threshold: $LOAD_THRESHOLD"

# Test each number with SIQS first, then GNFS if SIQS fails
declare -a SIQS_RESULTS
declare -a GNFS_RESULTS
MAX_SIQS_TIME=0

for i in 0 1 2 3 4; do
    echo ""
    echo "--- Testing 90d[$i] ---"
    wait_for_low_load
    echo "Load OK: $(awk '{print $1}' /proc/loadavg)"

    if test_siqs $i; then
        SIQS_RESULTS[$i]="PASS"
    else
        SIQS_RESULTS[$i]="FAIL"
        # Try GNFS as fallback
        echo "Trying GNFS fallback..."
        wait_for_low_load
        if test_gnfs $i; then
            GNFS_RESULTS[$i]="PASS"
        else
            GNFS_RESULTS[$i]="FAIL"
        fi
    fi
done

echo ""
echo "=== Results ==="
for i in 0 1 2 3 4; do
    echo "90d[$i]: SIQS=${SIQS_RESULTS[$i]:-N/A} GNFS=${GNFS_RESULTS[$i]:-N/A}"
done
