#!/bin/bash
# Comprehensive 90d benchmark - run when load is low
# Tests both SIQS (yafu_mod with closnuf) and GNFS approaches
# Usage: bash library/bench_90d.sh [siqs|gnfs|both]

set -e

MODE=${1:-both}
YAFU_SIQS=/tmp/agent-factoring-6/yafu_mod/yafu
YAFU_GNFS=/tmp/agent-factoring-4/yafu/yafu
GGNFS_DIR=/tmp/agent-factoring-4/yafu/factor/lasieve5_64/bin/
SEMIPRIMES_FILE=/tmp/agent-factoring-6/semiprimes.json
LOGFILE=/tmp/agent-factoring-6/experiments.log
TIMEOUT=295
SEED=42

# 90-digit semiprimes
declare -a NUMS=(
  "315547728763353682629201910311839572981703304736061985377687053517447821670148261179477299"
  "144009255506966132983765710705715101627946917437406468019915946276667323993959592666335633"
  "321739159513091262170255410982298427140058934463728910158638687584271260451816679219749747"
  "149689871933523898487085003856984468870482390037601728900188208806500682729361812646669171"
  "226492946512013642353358853937990356716879227059976248121251466901503301655504057622409363"
)

# Wait for load to drop
echo "Waiting for load to drop below 5..."
while true; do
  LOAD=$(awk '{print $1}' /proc/loadavg)
  if (( $(echo "$LOAD < 5" | bc -l) )); then
    echo "Load is $LOAD - starting benchmarks"
    break
  fi
  echo "Load is $LOAD - waiting 30s..."
  sleep 30
done

TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')

run_siqs() {
  local idx=$1
  local N=${NUMS[$idx]}
  local WORKDIR=$(mktemp -d /tmp/bench_siqs_${idx}_XXXXXX)

  echo "=== SIQS 90d[$idx] ==="
  local START=$(date +%s.%N)
  cd $WORKDIR
  local OUTPUT=$(echo "siqs($N)" | timeout $TIMEOUT $YAFU_SIQS -threads 1 -seed $SEED -siqsNB 20 -siqsB 120000 2>&1)
  local EXIT=$?
  local END=$(date +%s.%N)
  local ELAPSED=$(echo "$END - $START" | bc)

  if [ $EXIT -eq 0 ] && echo "$OUTPUT" | grep -q "factors found"; then
    echo "PASS: 90d[$idx] SIQS = ${ELAPSED}s"
    echo "[$TIMESTAMP] size: 90 | approach: yafu_mod SIQS (closnuf+VBITS512) NB=20 B=120K 90d[$idx] | time: $ELAPSED | notes: closnuf=digits_n+1, VBITS=512" >> $LOGFILE
  else
    echo "FAIL: 90d[$idx] SIQS timeout"
    echo "[$TIMESTAMP] size: 90 | approach: yafu_mod SIQS (closnuf+VBITS512) NB=20 B=120K 90d[$idx] | time: FAIL | notes: timeout at ${TIMEOUT}s" >> $LOGFILE
  fi

  rm -rf $WORKDIR
  echo "$ELAPSED"
}

run_gnfs() {
  local idx=$1
  local N=${NUMS[$idx]}
  local WORKDIR=$(mktemp -d /tmp/bench_gnfs_${idx}_XXXXXX)

  echo "=== GNFS 90d[$idx] ==="
  cd $WORKDIR
  cat > yafu.ini << EOF
ggnfs_dir=$GGNFS_DIR
EOF

  local START=$(date +%s.%N)
  local OUTPUT=$(echo "nfs($N)" | timeout $TIMEOUT $YAFU_GNFS -threads 1 -seed $SEED -xover 85 2>&1)
  local EXIT=$?
  local END=$(date +%s.%N)
  local ELAPSED=$(echo "$END - $START" | bc)

  if [ $EXIT -eq 0 ] && echo "$OUTPUT" | grep -q "factors found"; then
    echo "PASS: 90d[$idx] GNFS = ${ELAPSED}s"
    echo "[$TIMESTAMP] size: 90 | approach: YAFU GNFS (GGNFS sievers) 90d[$idx] | time: $ELAPSED | notes: default params, -xover 85" >> $LOGFILE
  else
    echo "FAIL: 90d[$idx] GNFS timeout"
    echo "[$TIMESTAMP] size: 90 | approach: YAFU GNFS (GGNFS sievers) 90d[$idx] | time: FAIL | notes: timeout at ${TIMEOUT}s" >> $LOGFILE
  fi

  rm -rf $WORKDIR
  echo "$ELAPSED"
}

echo "============================================"
echo "90-digit Benchmark - $(date)"
echo "Mode: $MODE"
echo "Load: $(awk '{print $1, $2, $3}' /proc/loadavg)"
echo "============================================"

MAX_TIME=0

if [ "$MODE" = "siqs" ] || [ "$MODE" = "both" ]; then
  echo ""
  echo "--- SIQS Tests ---"
  for i in 0 1 2 3 4; do
    T=$(run_siqs $i)
    if (( $(echo "$T > $MAX_TIME" | bc -l) )); then
      MAX_TIME=$T
    fi
  done
  echo "SIQS worst case: ${MAX_TIME}s"
fi

if [ "$MODE" = "gnfs" ] || [ "$MODE" = "both" ]; then
  echo ""
  echo "--- GNFS Tests ---"
  MAX_TIME=0
  for i in 0 1 2 3 4; do
    T=$(run_gnfs $i)
    if (( $(echo "$T > $MAX_TIME" | bc -l) )); then
      MAX_TIME=$T
    fi
  done
  echo "GNFS worst case: ${MAX_TIME}s"
fi

echo ""
echo "Benchmark complete at $(date)"
