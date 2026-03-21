#!/bin/bash
# Parallel factoring: runs multiple ECM workers in parallel
# Usage: ./par_factor.sh <N> [digits]
# Prints factor to stdout on success

N="$1"
DIGITS="${2:-$(echo -n "$N" | wc -c)}"
NCORES=$(nproc)
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKER="$SCRIPT_DIR/../ecm_worker"
export LD_LIBRARY_PATH=/usr/local/lib

# Trial division first
for p in 2 3 5 7 11 13 17 19 23 29 31 37 41 43 47; do
    if python3 -c "exit(0 if $N % $p == 0 else 1)" 2>/dev/null; then
        echo "$p"
        exit 0
    fi
done

# Determine ECM parameters based on digit count
# For balanced semiprimes, factor has ~digits/2 digits
HALF=$((DIGITS / 2))

# ECM B1 levels and curves per worker
# We'll run NCORES workers in parallel at each level
declare -a B1_LEVELS
declare -a CURVES_PER_WORKER

if [ "$HALF" -le 15 ]; then
    B1_LEVELS=(500 2000 11000)
    CURVES_PER_WORKER=(5 5 10)
elif [ "$HALF" -le 20 ]; then
    B1_LEVELS=(2000 11000 50000)
    CURVES_PER_WORKER=(5 10 20)
elif [ "$HALF" -le 25 ]; then
    B1_LEVELS=(11000 50000 250000)
    CURVES_PER_WORKER=(10 20 30)
elif [ "$HALF" -le 30 ]; then
    B1_LEVELS=(50000 250000 1000000)
    CURVES_PER_WORKER=(15 25 40)
elif [ "$HALF" -le 35 ]; then
    B1_LEVELS=(250000 1000000 3000000)
    CURVES_PER_WORKER=(20 40 60)
elif [ "$HALF" -le 40 ]; then
    B1_LEVELS=(1000000 3000000 11000000)
    CURVES_PER_WORKER=(30 50 80)
elif [ "$HALF" -le 45 ]; then
    B1_LEVELS=(3000000 11000000 43000000)
    CURVES_PER_WORKER=(40 80 100)
else
    B1_LEVELS=(11000000 43000000 110000000 260000000)
    CURVES_PER_WORKER=(50 100 150 200)
fi

TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR; kill 0 2>/dev/null" EXIT

for lvl_idx in "${!B1_LEVELS[@]}"; do
    B1=${B1_LEVELS[$lvl_idx]}
    CPW=${CURVES_PER_WORKER[$lvl_idx]}

    # Launch NCORES parallel workers
    for i in $(seq 1 $NCORES); do
        (
            result=$($WORKER "$N" "$B1" "$CPW" 2>/dev/null)
            if [ $? -eq 0 ] && [ -n "$result" ]; then
                echo "$result" > "$TMPDIR/result"
                # Signal success
                kill -USR1 $$ 2>/dev/null
            fi
        ) &
    done

    # Wait for all workers at this level
    wait

    # Check if any worker found a factor
    if [ -f "$TMPDIR/result" ]; then
        cat "$TMPDIR/result"
        exit 0
    fi
done

echo "FAILED" >&2
exit 1
