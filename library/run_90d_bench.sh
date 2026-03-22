#!/bin/bash
# 90d benchmark - test all approaches
# Usage: bash library/run_90d_bench.sh
# Only run when load < 5

set -e
TIMEOUT=295
SEED=42

# Binaries
YAFU_ORIG=/tmp/agent-factoring-4/yafu/yafu
YAFU_MOD=/tmp/agent-factoring-7/yafu_mod/yafu

# 90-digit semiprimes
NUMS=(
  "315547728763353682629201910311839572981703304736061985377687053517447821670148261179477299"
  "144009255506966132983765710705715101627946917437406468019915946276667323993959592666335633"
  "321739159513091262170255410982298427140058934463728910158638687584271260451816679219749747"
  "149689871933523898487085003856984468870482390037601728900188208806500682729361812646669171"
  "226492946512013642353358853937990356716879227059976248121251466901503301655504057622409363"
)

check_load() {
    local load=$(awk '{print int($1)}' /proc/loadavg)
    if [ "$load" -gt 5 ]; then
        echo "SKIP: load is $load (>5), waiting..."
        return 1
    fi
    return 0
}

run_test() {
    local label="$1"
    local binary="$2"
    local cmd="$3"
    local params="$4"
    local num="$5"
    local idx="$6"

    # Check load before each test
    if ! check_load; then
        echo "[$label] 90d[$idx]: SKIPPED (high load)"
        return
    fi

    local WORKDIR=$(mktemp -d /tmp/yafu_XXXXXX)
    cd $WORKDIR

    # For GNFS, need sievers
    if [[ "$label" == *"gnfs"* ]]; then
        for f in /tmp/agent-factoring-1/yafu_mod/factor/lasieve5_64/gnfs-lasieve4I*e; do
            [ -f "$f" ] && ln -sf "$f" .
        done
    fi

    export LD_LIBRARY_PATH=/usr/local/lib
    local start=$(date +%s%3N)
    local result=$(echo "$cmd($num)" | timeout $TIMEOUT $binary -threads 1 -seed $SEED $params 2>&1)
    local exit_code=$?
    local end=$(date +%s%3N)
    local elapsed=$(( (end - start) ))

    cd /tmp/agent-factoring-7
    rm -rf $WORKDIR

    if [ $exit_code -eq 124 ]; then
        echo "[$label] 90d[$idx]: TIMEOUT (${elapsed}ms)"
    elif echo "$result" | grep -q "prp\|factor"; then
        local time_s=$(echo "scale=1; $elapsed / 1000" | bc)
        echo "[$label] 90d[$idx]: ${time_s}s ✓"
    else
        echo "[$label] 90d[$idx]: FAIL (exit=$exit_code, ${elapsed}ms)"
    fi
}

echo "=== 90d Benchmark ==="
echo "Load: $(cat /proc/loadavg | cut -d' ' -f1-3)"
echo "Time: $(date)"
echo ""

# Test 1: Original YAFU, NB=20 B=120K (best known params for 90d)
echo "--- Original YAFU SIQS (NB=20 B=120K) ---"
for i in 0 1 2 3 4; do
    run_test "orig-siqs" "$YAFU_ORIG" "siqs" "-siqsNB 20 -siqsB 120000" "${NUMS[$i]}" "$i"
done

echo ""

# Test 2: Modified YAFU (closnuf), NB=20 B=120K
echo "--- Modified YAFU SIQS (closnuf, NB=20 B=120K) ---"
for i in 0 1 2 3 4; do
    run_test "mod-siqs" "$YAFU_MOD" "siqs" "-siqsNB 20 -siqsB 120000" "${NUMS[$i]}" "$i"
done

echo ""

# Test 3: Original YAFU, NB=18 B=100K (best for 89d)
echo "--- Original YAFU SIQS (NB=18 B=100K) ---"
for i in 0 1 2 3 4; do
    run_test "orig-nb18" "$YAFU_ORIG" "siqs" "-siqsNB 18 -siqsB 100000" "${NUMS[$i]}" "$i"
done

echo ""
echo "=== Done ==="
