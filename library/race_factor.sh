#!/bin/bash
# Fast factoring using YAFU SIQS with tuned parameters
# Usage: ./race_factor.sh <N> [timeout]
N="$1"
TIMEOUT="${2:-280}"
DIR="$(cd "$(dirname "$0")/.." && pwd)"
DIGITS=${#N}

# Select factor base size based on digit count
# Default YAFU overestimates, causing "matrix corrupt" errors
if [ "$DIGITS" -le 60 ]; then
    SIQSB=""  # default is fine for small numbers
elif [ "$DIGITS" -le 70 ]; then
    SIQSB="-siqsB 8000"
elif [ "$DIGITS" -le 80 ]; then
    SIQSB="-siqsB 15000"
elif [ "$DIGITS" -le 85 ]; then
    SIQSB="-siqsB 25000"
elif [ "$DIGITS" -le 90 ]; then
    SIQSB="-siqsB 40000"
elif [ "$DIGITS" -le 95 ]; then
    SIQSB="-siqsB 60000"
else
    SIQSB="-siqsB 90000"
fi

result=$(echo "siqs($N)" | LD_LIBRARY_PATH=/usr/local/lib timeout "${TIMEOUT}s" "$DIR/yafu/yafu" -threads 48 -seed 42 $SIQSB 2>/dev/null)
factor=$(echo "$result" | grep "^P[0-9]" | head -1 | sed 's/^P[0-9]* = //')
if [ -n "$factor" ]; then
    echo "$factor"
    exit 0
fi

# Fallback: msieve
SESSION="/tmp/msieve_race_${N}_$$.dat"
REMAINING=$((TIMEOUT - 10))
if [ "$REMAINING" -lt 30 ]; then REMAINING=30; fi
result=$(timeout "${REMAINING}s" "$DIR/msieve" -q -s "$SESSION" "$N" 2>/dev/null)
rm -f "$SESSION"
factor=$(echo "$result" | grep "^p[0-9]" | head -1 | sed 's/^p[0-9]*: //')
if [ -n "$factor" ]; then
    echo "$factor"
    exit 0
fi
exit 1
