#!/bin/bash
# Benchmark a factoring binary against semiprimes.json
# Usage: ./benchmark.sh <binary> <approach_name> [min_digits] [max_digits]
#
# Runs the binary on all semiprimes in the specified digit range,
# records results to experiments.log and updates algo-scaling.json

BINARY="$1"
APPROACH="$2"
MIN_DIGITS="${3:-30}"
MAX_DIGITS="${4:-100}"
SEMIPRIMES="../semiprimes.json"
LOG="../experiments.log"
SCALING="../algo-scaling.json"

if [ -z "$BINARY" ] || [ -z "$APPROACH" ]; then
    echo "Usage: $0 <binary> <approach_name> [min_digits] [max_digits]"
    exit 1
fi

if [ ! -f "$BINARY" ]; then
    echo "Binary not found: $BINARY"
    exit 1
fi

echo "Benchmarking $APPROACH ($BINARY) from $MIN_DIGITS to $MAX_DIGITS digits"

for digits in $(seq $MIN_DIGITS $MAX_DIGITS); do
    # Extract semiprimes for this digit count using python
    nums=$(python3 -c "
import json
with open('$SEMIPRIMES') as f:
    data = json.load(f)
key = str($digits)
if key in data:
    for n in data[key]:
        print(n)
")
    if [ -z "$nums" ]; then
        continue
    fi

    worst_time=0
    all_passed=true

    while IFS= read -r num; do
        start_time=$(date +%s.%N)
        result=$(timeout 295 $BINARY "$num" 2>/tmp/bench_stderr_$$)
        exit_code=$?
        end_time=$(date +%s.%N)

        elapsed=$(echo "$end_time - $start_time" | bc)
        stderr_out=$(cat /tmp/bench_stderr_$$ 2>/dev/null)

        if [ $exit_code -eq 0 ] && [ -n "$result" ]; then
            # Verify the factor
            is_valid=$(python3 -c "
n = int('$num')
f = int('$result'.strip())
if f > 1 and f < n and n % f == 0:
    print('YES')
else:
    print('NO')
")
            if [ "$is_valid" = "YES" ]; then
                timestamp=$(date '+%Y-%m-%d %H:%M:%S')
                echo "[$timestamp] size: $digits | approach: $APPROACH | time: ${elapsed}s | notes: factor=$result" >> "$LOG"

                # Track worst time
                is_worse=$(echo "$elapsed > $worst_time" | bc)
                if [ "$is_worse" -eq 1 ]; then
                    worst_time=$elapsed
                fi
            else
                timestamp=$(date '+%Y-%m-%d %H:%M:%S')
                echo "[$timestamp] size: $digits | approach: $APPROACH | time: FAIL | notes: invalid factor $result for $num" >> "$LOG"
                all_passed=false
            fi
        else
            timestamp=$(date '+%Y-%m-%d %H:%M:%S')
            if [ $exit_code -eq 124 ]; then
                echo "[$timestamp] size: $digits | approach: $APPROACH | time: FAIL | notes: timeout (295s) on $num" >> "$LOG"
            else
                echo "[$timestamp] size: $digits | approach: $APPROACH | time: FAIL | notes: exit=$exit_code on $num" >> "$LOG"
            fi
            all_passed=false
        fi
    done <<< "$nums"

    if [ "$all_passed" = true ]; then
        # Round to 3 decimal places
        worst_time_rounded=$(printf "%.3f" "$worst_time")
        echo "  $digits digits: worst time = ${worst_time_rounded}s"

        # Update algo-scaling.json
        python3 -c "
import json
try:
    with open('$SCALING') as f:
        data = json.load(f)
except:
    data = {}
if '$APPROACH' not in data:
    data['$APPROACH'] = {}
data['$APPROACH'][str($digits)] = float($worst_time_rounded)
with open('$SCALING', 'w') as f:
    json.dump(data, f, indent=2)
"
    else
        echo "  $digits digits: FAILED (some semiprimes not factored)"
        # Stop at this digit size - won't improve
        break
    fi
done

rm -f /tmp/bench_stderr_$$
echo "Benchmarking complete."
