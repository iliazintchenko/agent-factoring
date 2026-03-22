#!/bin/bash
# Race multiple YAFU SIQS instances with different seeds, fall back to msieve
# Usage: ./race_factor.sh <N> [timeout]
N="$1"
TIMEOUT="${2:-280}"
DIR="$(cd "$(dirname "$0")/.." && pwd)"
DIGITS=${#N}
TMPDIR=$(mktemp -d)
RESULT_FILE="$TMPDIR/result"
trap "rm -rf $TMPDIR; kill 0 2>/dev/null" EXIT

# Phase 1: Race YAFU instances with different seeds
# YAFU timeout per instance
YAFU_TIMEOUT=$((TIMEOUT * 55 / 100))
if [ "$YAFU_TIMEOUT" -lt 20 ]; then YAFU_TIMEOUT=20; fi
if [ "$YAFU_TIMEOUT" -gt 200 ]; then YAFU_TIMEOUT=200; fi

# Determine number of parallel YAFU instances and threads per instance
# Use 6 instances with 8 threads each = 48 cores
NINST=6
NTHREADS=8

for i in $(seq 1 $NINST); do
    (
        SEED=$((i * 1000 + RANDOM))
        result=$(echo "siqs($N)" | LD_LIBRARY_PATH=/usr/local/lib timeout "${YAFU_TIMEOUT}s" "$DIR/yafu/yafu" -threads "$NTHREADS" -seed "$SEED" -siqsT "$YAFU_TIMEOUT" 2>/dev/null)
        factor=$(echo "$result" | grep "^P[0-9]" | head -1 | sed 's/^P[0-9]* = //')
        if [ -n "$factor" ]; then
            echo "$factor" > "$RESULT_FILE"
        fi
    ) &
done

# Wait for first result or all to finish
while true; do
    if [ -f "$RESULT_FILE" ]; then
        cat "$RESULT_FILE"
        exit 0
    fi
    if ! jobs -r | grep -q .; then break; fi
    sleep 0.2
done

# Phase 2: Fallback to msieve
REMAINING=$((TIMEOUT - YAFU_TIMEOUT - 5))
if [ "$REMAINING" -lt 30 ]; then REMAINING=30; fi
SESSION="/tmp/msieve_race_${N}_$$.dat"
result=$(timeout "${REMAINING}s" "$DIR/msieve" -q -s "$SESSION" "$N" 2>/dev/null)
rm -f "$SESSION"
factor=$(echo "$result" | grep "^p[0-9]" | head -1 | sed 's/^p[0-9]*: //')
if [ -n "$factor" ]; then
    echo "$factor"
    exit 0
fi
exit 1
