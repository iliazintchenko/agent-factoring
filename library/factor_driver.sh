#!/bin/bash
# Optimal factoring driver for balanced semiprimes 30-89 digits.
# Automatically selects the best YAFU command and parameters.
# Usage: ./factor_driver.sh <N> [yafu_binary]
# All processes are single-threaded with seed=42 and timeout 295.

N="$1"
YAFU="${2:-/tmp/agent-factoring-4/yafu/yafu}"

if [ -z "$N" ]; then
    echo "Usage: $0 <N> [yafu_binary]" >&2
    exit 1
fi

# Count digits
DIGITS=$(echo -n "$N" | wc -c)

# Set up working directory
WORKDIR=$(mktemp -d /tmp/yafu_XXXXXX)
cd "$WORKDIR" || exit 1

# Select command and parameters based on digit count
CMD="siqs"
EXTRA=""

if [ "$DIGITS" -le 44 ]; then
    CMD="smallmpqs"
elif [ "$DIGITS" -le 72 ]; then
    CMD="siqs"
elif [ "$DIGITS" -le 73 ]; then
    EXTRA="-siqsNB 10"
elif [ "$DIGITS" -le 76 ]; then
    CMD="siqs"
elif [ "$DIGITS" -le 79 ]; then
    EXTRA="-siqsNB 11"
elif [ "$DIGITS" -le 81 ]; then
    EXTRA="-siqsNB 12"
elif [ "$DIGITS" -le 84 ]; then
    EXTRA="-siqsNB 14"
elif [ "$DIGITS" -eq 85 ]; then
    EXTRA="-siqsNB 18 -siqsB 70000"
elif [ "$DIGITS" -eq 86 ]; then
    EXTRA="-siqsNB 18 -siqsB 80000"
elif [ "$DIGITS" -eq 87 ]; then
    EXTRA="-siqsNB 18 -siqsB 90000"
elif [ "$DIGITS" -eq 88 ]; then
    EXTRA="-siqsNB 18 -siqsB 80000"
elif [ "$DIGITS" -ge 89 ]; then
    EXTRA="-siqsNB 18 -siqsB 100000"
fi

echo "${CMD}($N)" | timeout 295 "$YAFU" -threads 1 -seed 42 $EXTRA 2>&1
EXIT_CODE=$?

cd /
rm -rf "$WORKDIR"
exit $EXIT_CODE
