#!/bin/bash
# factor_smart.sh - Smart factoring for 90d+ semiprimes
# Tries GNFS with pre-computed poly if available, falls back to SIQS
# Single-core, seed 42, timeout 295s
#
# Usage: ./factor_smart.sh <N>

N="$1"
if [ -z "$N" ]; then echo "Usage: $0 <N>"; exit 1; fi

YAFU="/tmp/agent-factoring-4/yafu/yafu"
SIEVER_DIR="/tmp/agent-factoring-1/yafu_mod/factor/lasieve5_64"
POLY_DIR="/tmp/agent-factoring-1/library/gnfs_polys"

WORKDIR=$(mktemp -d /tmp/factorsmart_XXXXXX)
cd "$WORKDIR"

# Link GGNFS sievers
for f in "$SIEVER_DIR"/gnfs-lasieve4I*e; do
    [ -f "$f" ] && ln -sf "$f" .
done
echo "ggnfs_dir=$WORKDIR/" > yafu.ini

# Check if pre-computed GNFS poly exists for this number
POLY_FOUND=""
for poly in "$POLY_DIR"/90d_*.poly "$POLY_DIR"/90d_*.job; do
    [ -f "$poly" ] || continue
    if grep -q "^n: $N" "$poly" 2>/dev/null; then
        POLY_FOUND="$poly"
        break
    fi
done

if [ -n "$POLY_FOUND" ]; then
    # Use GNFS with pre-computed poly (best for hard 90d numbers)
    cp "$POLY_FOUND" nfs.job
    echo "210000" > "nfs.job.$(hostname).last_spq0"
    echo "nfs($N)" | timeout 295 "$YAFU" -threads 1 -seed 42 -xover 85 -R 2>&1
else
    # Use SIQS with optimal 90d parameters
    echo "siqs($N)" | timeout 295 "$YAFU" -threads 1 -seed 42 -siqsNB 20 -siqsB 120000 2>&1
fi

rm -rf "$WORKDIR"
