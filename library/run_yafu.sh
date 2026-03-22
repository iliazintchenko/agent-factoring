#!/bin/bash
# Usage: ./library/run_yafu.sh <N> [timeout] [threads]
# Runs YAFU SIQS on number N, prints smaller factor to stdout
N="$1"
TIMEOUT="${2:-290}"
THREADS="${3:-8}"

# Use a unique working dir to avoid file conflicts between parallel runs
WORKDIR=$(mktemp -d /tmp/yafu_XXXXXX)
cd "$WORKDIR"

RESULT=$(echo "siqs($N)" | timeout "$TIMEOUT" /tmp/agent-factoring-1/yafu/yafu -threads "$THREADS" 2>/dev/null)

# Extract factors
FACTORS=$(echo "$RESULT" | grep -E '^P[0-9]+ = ' | sed 's/P[0-9]* = //')

rm -rf "$WORKDIR"

if [ -z "$FACTORS" ]; then
    echo "FAIL" >&2
    exit 1
fi

# Print the smaller factor
echo "$FACTORS" | sort -n | head -1
