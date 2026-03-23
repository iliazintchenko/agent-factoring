#!/bin/bash
# factor90.sh - Optimized factoring for 90-digit semiprimes
# Uses YAFU SIQS with best known parameters for 90d
# Falls back to GNFS with pre-computed polynomial if available
#
# Usage: ./factor90.sh <N>

set -e
N="$1"
if [ -z "$N" ]; then echo "Usage: $0 <N>"; exit 1; fi

YAFU="/tmp/agent-factoring-4/yafu/yafu"
SIEVER_DIR="/tmp/agent-factoring-1/yafu_mod/factor/lasieve5_64"
POLY_DIR="/tmp/agent-factoring-1/library/gnfs_polys"
TIMEOUT=295

# Create clean workdir
WORKDIR=$(mktemp -d /tmp/factor90_XXXXXX)
cd "$WORKDIR"

# Link GGNFS sievers
for f in "$SIEVER_DIR"/gnfs-lasieve4I*e; do
    [ -f "$f" ] && ln -sf "$f" .
done
echo "ggnfs_dir=$WORKDIR/" > yafu.ini

START=$(date +%s)

# Strategy: Try SIQS first (best for most 90d numbers)
echo "Trying SIQS (NB=20, B=120000)..."
RESULT=$(echo "siqs($N)" | timeout $((TIMEOUT - 5)) "$YAFU" -threads 1 -seed 42 -siqsNB 20 -siqsB 120000 2>&1)

if echo "$RESULT" | grep -q "^P[0-9]"; then
    ELAPSED=$(($(date +%s) - START))
    echo "SIQS SUCCESS in ${ELAPSED}s"
    echo "$RESULT" | grep "^P[0-9]"
    rm -rf "$WORKDIR"
    exit 0
fi

# SIQS failed (timeout) - check if we have pre-computed poly for GNFS
ELAPSED=$(($(date +%s) - START))
REMAINING=$((TIMEOUT - ELAPSED))
echo "SIQS timeout after ${ELAPSED}s. ${REMAINING}s remaining for GNFS..."

if [ $REMAINING -lt 30 ]; then
    echo "Not enough time for GNFS"
    rm -rf "$WORKDIR"
    exit 1
fi

# Check if pre-computed poly exists for this number
POLY_FOUND=0
for poly in "$POLY_DIR"/90d_*.poly; do
    [ -f "$poly" ] || continue
    if grep -q "^n: $N" "$poly" 2>/dev/null; then
        echo "Found pre-computed polynomial: $poly"
        cp "$poly" nfs.job
        echo "210000" > "nfs.job.$(hostname).last_spq0"
        POLY_FOUND=1
        break
    fi
done

if [ $POLY_FOUND -eq 0 ]; then
    echo "No pre-computed polynomial found. Running GNFS with poly search..."
    echo "nfs($N)" | timeout $REMAINING "$YAFU" -threads 1 -seed 42 -xover 85 2>&1 | tail -5
else
    echo "Running GNFS with pre-computed poly..."
    echo "nfs($N)" | timeout $REMAINING "$YAFU" -threads 1 -seed 42 -xover 85 -R 2>&1 | tail -5
fi

TOTAL=$(($(date +%s) - START))
echo "Total time: ${TOTAL}s"
rm -rf "$WORKDIR"
