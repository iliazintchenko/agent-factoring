#!/bin/bash
# bench_all.sh - Benchmark YAFU on all semiprimes
# Usage: ./bench_all.sh <yafu_binary> <start_digits> <end_digits> [extra_args...]
# Output: timing results to stdout, one line per number

YAFU_BIN="${1:-/tmp/agent-factoring-10/yafu_build/yafu}"
START_D="${2:-30}"
END_D="${3:-89}"
shift 3 2>/dev/null
EXTRA_ARGS="$*"

SEMIPRIMES="/tmp/agent-factoring-10/semiprimes.json"

# Extract semiprimes for a given digit count
get_semiprimes() {
    python3 -c "
import json
with open('$SEMIPRIMES') as f:
    data = json.load(f)
nums = data.get('$1', [])
for n in nums:
    print(n)
"
}

# Get YAFU command and NB/B params for a digit count
get_params() {
    local d=$1
    if [ $d -le 44 ]; then
        echo "smallmpqs" ""
    elif [ $d -le 72 ]; then
        echo "siqs" ""
    elif [ $d -le 76 ]; then
        echo "siqs" ""
    elif [ $d -le 79 ]; then
        echo "siqs" "-siqsNB 11"
    elif [ $d -le 81 ]; then
        echo "siqs" "-siqsNB 12"
    elif [ $d -le 84 ]; then
        echo "siqs" "-siqsNB 14"
    elif [ $d -eq 85 ]; then
        echo "siqs" "-siqsNB 18 -siqsB 70000"
    elif [ $d -eq 86 ]; then
        echo "siqs" "-siqsNB 18 -siqsB 80000"
    elif [ $d -eq 87 ]; then
        echo "siqs" "-siqsNB 18 -siqsB 90000"
    elif [ $d -eq 88 ]; then
        echo "siqs" "-siqsNB 18 -siqsB 80000"
    elif [ $d -eq 89 ]; then
        echo "siqs" "-siqsNB 18 -siqsB 100000"
    else
        echo "siqs" "-siqsNB 18 -siqsB 100000"
    fi
}

for d in $(seq $START_D $END_D); do
    read cmd params <<< "$(get_params $d)"

    max_time=0
    all_times=""
    idx=0

    while IFS= read -r N; do
        WORKDIR=$(mktemp -d /tmp/yafu_bench_XXXXXX)

        start_t=$(date +%s%N)
        cd "$WORKDIR"
        echo "${cmd}(${N})" | timeout 295 "$YAFU_BIN" -threads 1 -seed 42 $params $EXTRA_ARGS > /dev/null 2>&1
        exit_code=$?
        end_t=$(date +%s%N)

        elapsed=$(echo "scale=3; ($end_t - $start_t) / 1000000000" | bc)

        if [ $exit_code -ne 0 ]; then
            elapsed="TIMEOUT"
            all_times="${all_times}TIMEOUT "
        else
            all_times="${all_times}${elapsed} "
            # Track max time
            is_greater=$(echo "$elapsed > $max_time" | bc -l 2>/dev/null)
            if [ "$is_greater" = "1" ]; then
                max_time=$elapsed
            fi
        fi

        rm -rf "$WORKDIR"
        idx=$((idx + 1))
    done < <(get_semiprimes $d)

    echo "SIZE=$d | MAX=$max_time | ALL=[$all_times] | PARAMS=$cmd $params"
done
