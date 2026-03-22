#!/bin/bash
# Factorize a number using the best available method
# Usage: ./factorize.sh <N> [deadline_seconds]
# Prints smallest prime factor

N=$1
DEADLINE=${2:-280}
DIGITS=${#N}
DIR="$(cd "$(dirname "$0")/.." && pwd)"

# For all sizes, msieve SIQS is fastest
result=$("$DIR/msieve" -q "$N" 2>/dev/null)
factor=$(echo "$result" | grep "^p[0-9]" | head -1 | sed 's/^p[0-9]*: //')
if [ -n "$factor" ]; then
    echo "$factor"
    exit 0
fi

# Fallback: ECM
if [ -f "$DIR/factor" ]; then
    result=$("$DIR/factor" "$N" "$DEADLINE" 2>/dev/null)
    if echo "$result" | grep -q "^[0-9]"; then
        echo "$result" | head -1
        exit 0
    fi
fi

echo "FAIL" >&2
exit 1
