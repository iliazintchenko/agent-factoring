#!/bin/bash
# GNFS Pipeline: Pre-computed poly + GGNFS sieve + YAFU post-processing
# Usage: ./gnfs_pipeline.sh <N> <job_file> [timeout_seconds]
#
# This script:
# 1. Copies the pre-computed polynomial job file
# 2. Runs GGNFS direct sieving (I=12 for 90d)
# 3. When enough relations collected, runs YAFU for filter+LA+sqrt
# 4. Outputs the factors

set -e

N="${1:?Usage: $0 <N> <job_file> [timeout]}"
JOB_FILE="${2:?Usage: $0 <N> <job_file> [timeout]}"
TIMEOUT="${3:-295}"

SIEVER="/tmp/agent-factoring-4/yafu/factor/lasieve5_64/bin/gnfs-lasieve4I12e"
YAFU="/tmp/agent-factoring-4/yafu/yafu"

if [ ! -f "$SIEVER" ]; then
    echo "ERROR: GGNFS siever not found: $SIEVER"
    exit 1
fi

if [ ! -x "$SIEVER" ]; then
    chmod +x "$SIEVER"
fi

# Create isolated work directory
WORKDIR=$(mktemp -d /tmp/gnfs_work_XXXXXX)
cd "$WORKDIR"

# Setup YAFU
for f in /tmp/agent-factoring-4/yafu/factor/lasieve5_64/bin/gnfs-lasieve4I*e; do
    ln -sf "$f" .
    chmod +x "$(basename $f)" 2>/dev/null
done
echo "ggnfs_dir=$WORKDIR/" > yafu.ini

# Copy job file
cp "$JOB_FILE" nfs.job

# Target relations (1.46M for 90d with standard params)
TARGET_RELS=1460000
START_Q=250000
QRANGE=200000

START_TIME=$(date +%s)

echo "=== GNFS Pipeline ==="
echo "N: $N"
echo "Target: $TARGET_RELS relations"
echo "Timeout: ${TIMEOUT}s"

# Phase 1: Sieving
echo "Phase 1: Sieving (Q=$START_Q, range=$QRANGE)..."
timeout $((TIMEOUT - 35)) "$SIEVER" \
    -f $START_Q -c $QRANGE -o rels.out -n 0 -a nfs.job 2>sieve.log &
SIEVE_PID=$!

# Monitor sieve progress
while kill -0 $SIEVE_PID 2>/dev/null; do
    sleep 5
    NOW=$(date +%s)
    ELAPSED=$((NOW - START_TIME))
    REMAINING=$((TIMEOUT - ELAPSED))

    if [ -f rels.out ]; then
        RELS=$(wc -l < rels.out)
    else
        RELS=0
    fi

    RATE=0
    if [ $ELAPSED -gt 0 ]; then
        RATE=$((RELS / ELAPSED))
    fi

    echo "  t=${ELAPSED}s: ${RELS}/${TARGET_RELS} rels (${RATE}/sec), ${REMAINING}s remaining"

    # Stop sieving when we have enough relations (leave time for post-processing)
    if [ $RELS -ge $TARGET_RELS ] && [ $REMAINING -gt 30 ]; then
        echo "  Enough relations! Stopping sieve."
        kill $SIEVE_PID 2>/dev/null
        wait $SIEVE_PID 2>/dev/null
        break
    fi

    # If running low on time, stop even if not enough rels
    if [ $REMAINING -lt 35 ]; then
        echo "  Time running low (${REMAINING}s left). Stopping sieve."
        kill $SIEVE_PID 2>/dev/null
        wait $SIEVE_PID 2>/dev/null
        break
    fi
done

wait $SIEVE_PID 2>/dev/null || true

SIEVE_END=$(date +%s)
SIEVE_TIME=$((SIEVE_END - START_TIME))

if [ ! -f rels.out ]; then
    echo "ERROR: No relations file produced"
    rm -rf "$WORKDIR"
    exit 1
fi

TOTAL_RELS=$(wc -l < rels.out)
echo "Phase 1 complete: $TOTAL_RELS rels in ${SIEVE_TIME}s"

if [ $TOTAL_RELS -lt $((TARGET_RELS * 80 / 100)) ]; then
    echo "WARNING: Only ${TOTAL_RELS}/${TARGET_RELS} relations. Post-processing may fail."
fi

# Phase 2: Post-processing with YAFU
# YAFU's -R flag resumes NFS, -nc skips sieving
# We need to set up the nfs.job with relation count info
NOW=$(date +%s)
REMAINING=$((TIMEOUT - (NOW - START_TIME)))
echo "Phase 2: Post-processing (${REMAINING}s budget)..."

# Add relations to nfs.dat (YAFU format)
# YAFU expects relations in its own format, but can read GGNFS output
cp rels.out msieve.dat.r 2>/dev/null || true

# Try YAFU post-processing
RESULT=$(echo "nfs($N)" | timeout $((REMAINING - 2)) LD_LIBRARY_PATH=/usr/local/lib "$YAFU" \
    -threads 1 -seed 42 -xover 85 -R -nc 2>&1) || true

END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))

echo "=== Result ==="
echo "$RESULT" | grep -i "factor\|prp\|p[0-9]"
echo "Total time: ${TOTAL_TIME}s"
echo "Workdir: $WORKDIR"

# Don't clean up on success so we can inspect
# rm -rf "$WORKDIR"
