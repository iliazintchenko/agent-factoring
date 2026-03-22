#!/bin/bash
# Hybrid factoring: race msieve SIQS against parallel ECM
# Usage: ./hybrid_factor.sh <N> [deadline_seconds]
# For each number, runs:
#   1. msieve SIQS (1 core)
#   2. parallel ECM with remaining cores
# First to find a factor wins.

N="$1"
DEADLINE="${2:-280}"
DIGITS=${#N}
DIR="$(cd "$(dirname "$0")/.." && pwd)"
TMPDIR="/tmp/hybrid_$$_${N:0:10}"
mkdir -p "$TMPDIR"

# FIFO for result communication
FIFO="$TMPDIR/result"
mkfifo "$FIFO"

cleanup() {
    # Kill all child processes
    kill $(jobs -p) 2>/dev/null
    wait 2>/dev/null
    rm -rf "$TMPDIR"
}
trap cleanup EXIT

# Method 1: msieve SIQS (single-threaded, most efficient per-core)
(
    result=$("$DIR/msieve" -q -s "$TMPDIR/siqs.dat" "$N" 2>/dev/null)
    factor=$(echo "$result" | grep "^p[0-9]" | head -1 | sed 's/^p[0-9]*: //')
    [ -n "$factor" ] && echo "$factor" > "$FIFO" 2>/dev/null
) &

# Method 2: parallel ECM (good use of remaining cores)
if [ -f "$DIR/factor" ]; then
    (
        result=$("$DIR/factor" "$N" 2>/dev/null)
        factor=$(echo "$result" | head -1)
        echo "$result" | grep -q "^[0-9]" && echo "$factor" > "$FIFO" 2>/dev/null
    ) &
fi

# Wait for first result
factor=$(timeout "${DEADLINE}s" head -1 "$FIFO" 2>/dev/null)

if [ -n "$factor" ]; then
    echo "$factor"
    exit 0
fi
exit 1
