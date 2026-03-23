#!/bin/bash
# Meta-factorer: picks the best SIQS implementation based on digit count
# This script selects the fastest implementation for each size range
# based on benchmarked worst-case times:
#
# 30-55d: siqs_bucket2 (Gray code + small M, fastest overall)
# 56-70d: turbo_pp (interleaved sieve + DLP graph + structured GE)

N="$1"
if [ -z "$N" ]; then echo "Usage: $0 <N>" >&2; exit 1; fi

DIGITS=$(echo -n "$N" | wc -c)
DIR=$(dirname "$0")

if [ "$DIGITS" -le 55 ]; then
    exec "$DIR/../siqs_bucket2" "$N"
elif [ "$DIGITS" -le 90 ]; then
    exec "$DIR/../turbo_pp" "$N"
else
    echo "FAIL: number too large" >&2
    exit 1
fi
