#!/bin/bash
# try_90d_gnfs.sh - Attempt GNFS factoring of 90d semiprimes
# Only runs when load is low enough
# Usage: ./try_90d_gnfs.sh <index 0-4>

IDX=${1:-1}
YAFU="/tmp/agent-factoring-4/yafu/yafu"
SIEVER_DIR="/tmp/agent-factoring-1/yafu_mod/factor/lasieve5_64"
POLY_DIR="/tmp/agent-factoring-8/library/gnfs_polys"

# Get the number
N=$(python3 -c "import json; d=json.load(open('/tmp/agent-factoring-8/semiprimes.json')); print(d['90'][$IDX])")
echo "Attempting GNFS on 90d[$IDX]: $N"

# Check load
LOAD=$(cat /proc/loadavg | awk '{print $1}')
echo "Current load: $LOAD"

# Create workdir
WORKDIR=$(mktemp -d /tmp/gnfs_90d_XXXXXX)
cd "$WORKDIR"

# Setup sievers
for f in "$SIEVER_DIR"/gnfs-lasieve4I*e; do
    [ -f "$f" ] && ln -sf "$f" .
done
echo "ggnfs_dir=$WORKDIR/" > yafu.ini

# Copy polynomial
cp "$POLY_DIR/90d_${IDX}.job" nfs.job 2>/dev/null || true

# Run YAFU NFS
export LD_LIBRARY_PATH=/usr/local/lib
echo "Starting GNFS..."
START=$(date +%s)
echo "nfs($N)" | timeout 295 "$YAFU" -threads 1 -seed 42 -xover 85 2>&1 | tee /tmp/gnfs_90d_${IDX}.log | tail -20
END=$(date +%s)
ELAPSED=$((END - START))
echo "Elapsed: ${ELAPSED}s"

# Cleanup
cd /tmp
rm -rf "$WORKDIR"
