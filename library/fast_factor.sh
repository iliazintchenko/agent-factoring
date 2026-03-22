#!/bin/bash
# Fast factoring: try YAFU SIQS first, fall back to msieve
# Usage: ./fast_factor.sh <N> [timeout]
N="$1"
TIMEOUT="${2:-280}"
DIR="$(cd "$(dirname "$0")/.." && pwd)"
DIGITS=${#N}

# Try YAFU SIQS with a limited timeout (YAFU is fast when it works, hangs when buggy)
YAFU_TIMEOUT=$((TIMEOUT / 3))
if [ "$YAFU_TIMEOUT" -lt 20 ]; then YAFU_TIMEOUT=20; fi
if [ "$YAFU_TIMEOUT" -gt 180 ]; then YAFU_TIMEOUT=180; fi

result=$(echo "siqs($N)" | LD_LIBRARY_PATH=/usr/local/lib timeout "${YAFU_TIMEOUT}s" "$DIR/yafu/yafu" -threads 48 2>/dev/null)
factors=$(echo "$result" | grep "^P[0-9]" | sed 's/^P[0-9]* = //' | sort -n)
smallest=$(echo "$factors" | head -1)
if [ -n "$smallest" ]; then
    echo "$smallest"
    exit 0
fi

# Fallback: msieve
SESSION="/tmp/msieve_fb_${N}_$$.dat"
result=$(timeout "${TIMEOUT}s" "$DIR/msieve" -q -s "$SESSION" "$N" 2>/dev/null)
rm -f "$SESSION"
factor=$(echo "$result" | grep "^p[0-9]" | head -1 | sed 's/^p[0-9]*: //')
if [ -n "$factor" ]; then
    echo "$factor"
    exit 0
fi

exit 1
