#!/bin/bash
# Race multiple msieve instances with different random seeds
# Usage: ./race_factor.sh <N> [deadline_seconds] [num_racers]
# Prints the first factor found

N="$1"
DEADLINE="${2:-280}"
RACERS="${3:-5}"
DIR="$(cd "$(dirname "$0")/.." && pwd)"
TMPDIR="/tmp/race_$$"
mkdir -p "$TMPDIR"

# Create a named pipe for results
PIPE="$TMPDIR/result"
mkfifo "$PIPE"

# Launch racers
for i in $(seq 1 $RACERS); do
    (
        SESSION="$TMPDIR/s_${i}.dat"
        result=$("$DIR/msieve" -q -s "$SESSION" "$N" 2>/dev/null)
        factor=$(echo "$result" | grep "^p[0-9]" | head -1 | sed 's/^p[0-9]*: //')
        if [ -n "$factor" ]; then
            echo "$factor" > "$PIPE" 2>/dev/null
        fi
    ) &
done

# Wait for first result with timeout
factor=$(timeout "${DEADLINE}s" head -1 "$PIPE" 2>/dev/null)

# Kill all remaining msieve processes for this race
pkill -P $$ 2>/dev/null
wait 2>/dev/null

# Cleanup
rm -rf "$TMPDIR"

if [ -n "$factor" ]; then
    echo "$factor"
    exit 0
fi
exit 1
