#!/bin/bash
# Wait for low system load then run YAFU SIQS on a specific number
# Usage: ./wait_and_test.sh <N> <max_load> <NB> <B> [yafu_binary]
#
# Example: ./wait_and_test.sh 31554... 10 18 100000

N="$1"
MAX_LOAD="${2:-10}"
NB="${3:-18}"
B="${4:-100000}"
YAFU="${5:-/tmp/agent-factoring-10/yafu_mod/yafu}"

if [ -z "$N" ]; then
    echo "Usage: $0 <N> [max_load] [NB] [B] [yafu_binary]"
    exit 1
fi

echo "Waiting for load < $MAX_LOAD..."

while true; do
    LOAD=$(cat /proc/loadavg | awk '{print int($1)}')
    if [ "$LOAD" -lt "$MAX_LOAD" ]; then
        echo "Load is $LOAD < $MAX_LOAD. Starting test at $(date)."
        break
    fi
    sleep 10
done

WORKDIR=$(mktemp -d /tmp/yafu_wait_XXXXXX)
cd "$WORKDIR"

START=$(date +%s)
timeout 295 "$YAFU" -threads 1 -seed 42 -siqsNB $NB -siqsB $B -noopt <<< "siqs($N)" 2>&1 | tail -20
EXIT=$?
ELAPSED=$(( $(date +%s) - START ))

echo ""
echo "Exit: $EXIT, Elapsed: ${ELAPSED}s"
echo "Load at end: $(cat /proc/loadavg | awk '{print $1}')"

cd /tmp
rm -rf "$WORKDIR"
