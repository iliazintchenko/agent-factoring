#!/bin/bash
# Factor using YAFU's SIQS directly (skipping ECM for balanced semiprimes)
# Usage: ./yafu_factor.sh <N> [threads] [deadline]
# Prints the smallest factor to stdout.
N="$1"
THREADS="${2:-48}"
DEADLINE="${3:-280}"
DIR="$(cd "$(dirname "$0")/.." && pwd)"

result=$(echo "siqs($N)" | LD_LIBRARY_PATH=/usr/local/lib timeout "${DEADLINE}s" "$DIR/yafu/yafu" -threads "$THREADS" 2>/dev/null)

# Extract factors (format: P37 = 12345...)
factors=$(echo "$result" | grep "^P[0-9]" | sed 's/^P[0-9]* = //' | sort -n)
smallest=$(echo "$factors" | head -1)
if [ -n "$smallest" ]; then
    echo "$smallest"
    exit 0
fi
exit 1
