#!/bin/bash
# Benchmark 90d semiprimes with modified YAFU (VBITS=512, closnuf +2)
# Runs sequentially, one at a time, to avoid cache contention
# Usage: bash library/run_90d_bench.sh

YAFU_MOD=/tmp/agent-factoring-4/yafu_mod/yafu
YAFU_STD=/tmp/agent-factoring-4/yafu/yafu
TIMEOUT=295
SEED=42

declare -a NUMS=(
  "315547728763353682629201910311839572981703304736061985377687053517447821670148261179477299"
  "144009255506966132983765710705715101627946917437406468019915946276667323993959592666335633"
  "321739159513091262170255410982298427140058934463728910158638687584271260451816679219749747"
  "149689871933523898487085003856984468870482390037601728900188208806500682729361812646669171"
  "226492946512013642353358853937990356716879227059976248121251466901503301655504057622409363"
)

# NB/B parameter sets to test
declare -a NB_VALS=(18 20 22)
declare -a B_VALS=(100000 110000 120000)

echo "Waiting for load < 8..."
while true; do
  LOAD=$(awk '{printf "%.0f", $1}' /proc/loadavg)
  if [ "$LOAD" -lt 8 ]; then
    echo "Load is $(cat /proc/loadavg | awk '{print $1}') - starting"
    break
  fi
  sleep 15
done

TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')
echo "=== 90d Benchmark at $TIMESTAMP ==="
echo "Load: $(cat /proc/loadavg)"

MAX_TIME=0
PASS=0
FAIL=0

for idx in 0 1 2 3 4; do
  N=${NUMS[$idx]}
  BEST_TIME=999
  BEST_PARAMS=""

  for NB in "${NB_VALS[@]}"; do
    for B in "${B_VALS[@]}"; do
      # Check load before each run
      CUR_LOAD=$(awk '{printf "%.0f", $1}' /proc/loadavg)
      if [ "$CUR_LOAD" -gt 15 ]; then
        echo "Load too high ($CUR_LOAD), waiting..."
        while [ "$(awk '{printf "%.0f", $1}' /proc/loadavg)" -gt 10 ]; do sleep 10; done
      fi

      WORKDIR=$(mktemp -d /tmp/yafu_bench_XXXXXX)
      cd $WORKDIR

      START=$(date +%s.%N)
      OUTPUT=$(timeout $TIMEOUT bash -c "echo 'siqs($N)' | LD_LIBRARY_PATH=/usr/local/lib $YAFU_MOD -threads 1 -seed $SEED -siqsNB $NB -siqsB $B -noopt 2>&1")
      EXIT=$?
      END=$(date +%s.%N)
      ELAPSED=$(echo "$END - $START" | bc)

      rm -rf $WORKDIR

      if [ $EXIT -eq 0 ] && echo "$OUTPUT" | grep -q "factors found"; then
        echo "  90d[$idx] NB=$NB B=$B: ${ELAPSED}s PASS"
        echo "[$TIMESTAMP] size: 90 | approach: yafu_mod SIQS NB=$NB B=$B 90d[$idx] | time: $ELAPSED | notes: VBITS=512 closnuf+2 noopt" >> /tmp/agent-factoring-4/experiments.log
        if (( $(echo "$ELAPSED < $BEST_TIME" | bc -l) )); then
          BEST_TIME=$ELAPSED
          BEST_PARAMS="NB=$NB B=$B"
        fi
        break 2  # Found a working param set, move to next number
      else
        echo "  90d[$idx] NB=$NB B=$B: timeout"
      fi
    done
  done

  if (( $(echo "$BEST_TIME < 999" | bc -l) )); then
    PASS=$((PASS + 1))
    echo "90d[$idx]: PASS ${BEST_TIME}s ($BEST_PARAMS)"
    if (( $(echo "$BEST_TIME > $MAX_TIME" | bc -l) )); then
      MAX_TIME=$BEST_TIME
    fi
  else
    FAIL=$((FAIL + 1))
    echo "90d[$idx]: FAIL (all params timed out)"
    echo "[$TIMESTAMP] size: 90 | approach: yafu_mod SIQS 90d[$idx] all params | time: FAIL | notes: NB 18-22 B 100K-120K all timeout" >> /tmp/agent-factoring-4/experiments.log
  fi
done

echo ""
echo "=== Results: $PASS/5 pass, worst time: ${MAX_TIME}s ==="
echo "[$TIMESTAMP] size: 90 | approach: yafu_mod SIQS summary | time: ${MAX_TIME} | notes: $PASS/5 pass" >> /tmp/agent-factoring-4/experiments.log
