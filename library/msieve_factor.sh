#!/bin/bash
# Factor a number using msieve with a unique session file
# Usage: ./msieve_factor.sh <N> [timeout_seconds]
# Prints the smallest prime factor to stdout

N="$1"
TIMEOUT="${2:-280}"
DIR="$(cd "$(dirname "$0")/.." && pwd)"
SESSION="/tmp/msieve_${N}_$$.dat"

result=$(timeout "${TIMEOUT}s" "$DIR/msieve" -q -s "$SESSION" "$N" 2>/dev/null)
rm -f "$SESSION"

# Extract first factor
factor=$(echo "$result" | grep "^p[0-9]" | head -1 | sed 's/^p[0-9]*: //')
if [ -n "$factor" ]; then
    echo "$factor"
    exit 0
fi
exit 1
