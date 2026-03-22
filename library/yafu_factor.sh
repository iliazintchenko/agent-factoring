#!/bin/bash
# Factor using YAFU with multi-threaded SIQS
# Usage: ./yafu_factor.sh <N> [threads] [deadline]
# Prints smallest factor to stdout

N="$1"
THREADS="${2:-8}"
DEADLINE="${3:-280}"
DIR="$(cd "$(dirname "$0")/.." && pwd)"

result=$(echo "factor($N)" | timeout "${DEADLINE}s" "$DIR/yafu_bin" -threads "$THREADS" 2>/dev/null)

# Extract factor from YAFU output
factor=$(echo "$result" | grep "^P[0-9]" | head -1 | sed 's/^P[0-9]* = //')
if [ -n "$factor" ]; then
    echo "$factor"
    exit 0
fi
exit 1
