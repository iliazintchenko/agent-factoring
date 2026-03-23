#!/bin/bash
# gnfs_90d.sh - Direct GGNFS sieve + YAFU post-processing for 90d semiprimes
#
# Usage: ./gnfs_90d.sh <N> <poly_file> [timeout_seconds]
#
# Uses pre-computed polynomial, direct GGNFS siever (not YAFU wrapper),
# and YAFU for post-processing (filter + BL + sqrt).
#
# Requires:
#   - GGNFS siever at /tmp/agent-factoring-4/yafu/factor/lasieve5_64/bin/
#   - YAFU at /tmp/agent-factoring-4/yafu/yafu
#   - Pre-computed polynomial in YAFU .job format

set -e

N="$1"
POLY_FILE="$2"
TIMEOUT="${3:-295}"

if [ -z "$N" ] || [ -z "$POLY_FILE" ]; then
    echo "Usage: $0 <N> <poly_file> [timeout_seconds]"
    exit 1
fi

YAFU="/tmp/agent-factoring-4/yafu/yafu"
SIEVER="/tmp/agent-factoring-4/yafu/factor/lasieve5_64/bin/gnfs-lasieve4I12e"
WORKDIR=$(mktemp -d /tmp/gnfs_90d_XXXXXX)

echo "GNFS 90d pipeline: N=$N"
echo "Workdir: $WORKDIR"
echo "Timeout: ${TIMEOUT}s"

START=$(date +%s%N)

cd $WORKDIR

# Setup siever links
for f in /tmp/agent-factoring-4/yafu/factor/lasieve5_64/bin/gnfs-lasieve4I*e; do
    ln -sf "$f" .
done

# Setup YAFU config
echo "ggnfs_dir=$WORKDIR/" > yafu.ini

# Copy polynomial
cp "$POLY_FILE" nfs.job

# Phase 1: Direct GGNFS sieve (bypass YAFU wrapper for speed)
TARGET_RELS=1460000  # 1.46M relations needed
START_Q=250000       # Start at higher Q for better yield
QRANGE=5000          # Process in batches

echo "Phase 1: Sieving (target ${TARGET_RELS} rels, Q starting at ${START_Q})"

TOTAL_RELS=0
CURRENT_Q=$START_Q
SIEVE_START=$(date +%s)

while [ $TOTAL_RELS -lt $TARGET_RELS ]; do
    ELAPSED=$(( $(date +%s) - $(date -d @$((${START}000000000 / 1000000000)) +%s 2>/dev/null || echo $SIEVE_START) ))
    # Check time budget: reserve 35s for post-processing
    NOW=$(date +%s%N)
    ELAPSED_TOTAL=$(( (NOW - START) / 1000000000 ))
    if [ $ELAPSED_TOTAL -gt $((TIMEOUT - 35)) ]; then
        echo "Time budget exceeded for sieve (${ELAPSED_TOTAL}s used, need 35s for post)"
        break
    fi

    # Run siever batch
    timeout $((TIMEOUT - ELAPSED_TOTAL - 30)) $SIEVER \
        -f $CURRENT_Q -c $QRANGE -o rels_batch.out -n 0 -a nfs.job 2>/dev/null

    if [ -f rels_batch.out ]; then
        BATCH_RELS=$(wc -l < rels_batch.out)
        cat rels_batch.out >> rels_all.out
        TOTAL_RELS=$((TOTAL_RELS + BATCH_RELS))
        rm rels_batch.out
    fi

    CURRENT_Q=$((CURRENT_Q + QRANGE))
    SIEVE_ELAPSED=$(($(date +%s) - SIEVE_START))
    RATE=$((TOTAL_RELS / (SIEVE_ELAPSED > 0 ? SIEVE_ELAPSED : 1)))

    echo "  Q=${CURRENT_Q}: ${TOTAL_RELS}/${TARGET_RELS} rels (${RATE}/sec, ${SIEVE_ELAPSED}s)"
done

SIEVE_TIME=$(($(date +%s) - SIEVE_START))
echo "Sieve done: ${TOTAL_RELS} rels in ${SIEVE_TIME}s"

if [ $TOTAL_RELS -lt $((TARGET_RELS * 90 / 100)) ]; then
    echo "FAIL: insufficient relations (${TOTAL_RELS} < ${TARGET_RELS})"
    rm -rf $WORKDIR
    exit 1
fi

# Phase 2: YAFU post-processing (filter + BL + sqrt)
echo "Phase 2: Post-processing with YAFU"

# Create combined relation file in YAFU format
# YAFU expects relations in nfs.job.dat
cp rels_all.out nfs.job.dat

# Set last_spq to avoid re-sieving
echo "$CURRENT_Q" > "nfs.job.$(hostname).last_spq0"

# Run YAFU with -R (resume) -nc (no sieve continuation)
REMAINING=$((TIMEOUT - $(( ($(date +%s%N) - START) / 1000000000 ))))
echo "Running YAFU post-processing (${REMAINING}s budget)"

export LD_LIBRARY_PATH=/usr/local/lib
echo "nfs($N)" | timeout $REMAINING $YAFU -threads 1 -seed 42 -xover 85 -R -nc 2>&1 | tr '\r' '\n' | grep -E "(P[0-9]+ =|elapsed|factor)" | head -10

END=$(date +%s%N)
TOTAL_MS=$(( (END - START) / 1000000 ))
echo "Total: ${TOTAL_MS}ms (sieve ${SIEVE_TIME}s + post $((TOTAL_MS/1000 - SIEVE_TIME))s)"

rm -rf $WORKDIR
