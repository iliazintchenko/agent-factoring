#!/bin/bash
# YAFU SIQS with resume support for hard numbers
# Usage: bash yafu_resume.sh <N> [phase1_timeout] [total_timeout]
# Runs SIQS with save/resume to handle numbers that need >300s in a single run
# Single-threaded, seed 42

N=$1
P1_TIMEOUT=${2:-200}
TOTAL_TIMEOUT=${3:-295}
YAFU="/tmp/agent-factoring-1/yafu/yafu"

WORKDIR=$(mktemp -d /tmp/yafu_resume_XXXXXX)
cd "$WORKDIR"

TOTAL_START=$(date +%s%N)

# Phase 1: initial sieving
echo "siqs($N)" | timeout $P1_TIMEOUT "$YAFU" -threads 1 -seed 42 > /dev/null 2>&1
RC=$?

ELAPSED_NS=$(($(date +%s%N) - TOTAL_START))
ELAPSED_S=$(echo "scale=3; $ELAPSED_NS / 1000000000" | bc)

if [ $RC -eq 0 ]; then
    # Factored in phase 1
    # Extract factors from factor.log
    FACTORS=$(grep -oP 'P\d+ = \d+' factor.log 2>/dev/null | tail -2)
    echo "DONE phase=1 time=$ELAPSED_S factors=$FACTORS"
    rm -rf "$WORKDIR"
    exit 0
fi

# Phase 2+: resume until done or total timeout exceeded
PHASE=2
while true; do
    REMAINING_NS=$((TOTAL_TIMEOUT * 1000000000 - ($(date +%s%N) - TOTAL_START)))
    REMAINING_S=$(echo "scale=0; $REMAINING_NS / 1000000000" | bc)

    if [ "$REMAINING_S" -le 5 ]; then
        echo "FAIL time=$ELAPSED_S (no time left for phase $PHASE)"
        rm -rf "$WORKDIR"
        exit 1
    fi

    echo "siqs($N)" | timeout $REMAINING_S "$YAFU" -threads 1 -seed 42 > /dev/null 2>&1
    RC=$?

    ELAPSED_NS=$(($(date +%s%N) - TOTAL_START))
    ELAPSED_S=$(echo "scale=3; $ELAPSED_NS / 1000000000" | bc)

    if [ $RC -eq 0 ]; then
        FACTORS=$(grep -oP 'P\d+ = \d+' factor.log 2>/dev/null | tail -2)
        echo "DONE phase=$PHASE time=$ELAPSED_S factors=$FACTORS"
        rm -rf "$WORKDIR"
        exit 0
    fi

    PHASE=$((PHASE + 1))
    if [ $PHASE -gt 10 ]; then
        echo "FAIL time=$ELAPSED_S (too many phases)"
        rm -rf "$WORKDIR"
        exit 1
    fi
done
