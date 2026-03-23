#!/bin/bash
# Run benchmark for a single semiprime
# Usage: run_bench.sh <binary> <approach> <digits> <index> <number>
# Outputs: <digits> <index> <time_or_FAIL> <factor>

BINARY="$1"
APPROACH="$2"
DIGITS="$3"
INDEX="$4"
NUM="$5"

start_time=$(date +%s.%N)
result=$(timeout 295 "$BINARY" "$NUM" 2>/dev/null)
exit_code=$?
end_time=$(date +%s.%N)

elapsed=$(echo "$end_time - $start_time" | bc)

if [ $exit_code -eq 0 ] && [ -n "$result" ]; then
    result_trimmed=$(echo "$result" | tr -d '[:space:]')
    is_valid=$(python3 -c "
n = int('$NUM')
f = int('$result_trimmed')
if f > 1 and f < n and n % f == 0:
    print('YES')
else:
    print('NO')
" 2>/dev/null)
    if [ "$is_valid" = "YES" ]; then
        echo "$DIGITS $INDEX $elapsed $result_trimmed"
    else
        echo "$DIGITS $INDEX FAIL invalid_factor"
    fi
elif [ $exit_code -eq 124 ]; then
    echo "$DIGITS $INDEX TIMEOUT 295"
else
    echo "$DIGITS $INDEX FAIL exit_$exit_code"
fi
