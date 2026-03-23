#!/bin/bash
# Standalone GNFS: GGNFS direct sieve + YAFU post-processing
# Bypasses YAFU's batch overhead for maximum sieve throughput
#
# Usage: timeout 295 bash library/gnfs_standalone.sh <N> <poly_idx>

set -e
N="$1"
IDX="${2:-0}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
POLY="$SCRIPT_DIR/library/gnfs_polys/90d_${IDX}.job"
SIEVER="/tmp/agent-factoring-1/yafu_mod/factor/lasieve5_64/gnfs-lasieve4I12e"
YAFU="/tmp/agent-factoring-9/yafu/yafu"

[ ! -f "$POLY" ] && echo "Error: $POLY not found" && exit 1

WORKDIR=$(mktemp -d /tmp/gnfs_sa_XXXXXX)
cd "$WORKDIR"

echo "[$(date +%H:%M:%S)] Standalone GNFS 90d[$IDX]"
START=$SECONDS

# ============ SIEVE (one big batch) ============
cp "$POLY" job.poly

# Run GGNFS in 100K-Q mega-batches for minimum overhead
Q=210000
TOTAL=0
TARGET=1460000  # Standard YAFU target
SIEVE_BUDGET=$((295 - 30))  # Leave 30s for post

while [ $TOTAL -lt $TARGET ] && [ $((SECONDS - START)) -lt $SIEVE_BUDGET ]; do
    REMAINING=$((SIEVE_BUDGET - SECONDS + START))
    BATCH=$((REMAINING > 15 ? 100000 : 10000))
    timeout $((REMAINING + 5)) $SIEVER -f $Q -c $BATCH -o rels_${Q}.out -a job.poly -n 0 2>/dev/null
    NEW=$(wc -l < rels_${Q}.out 2>/dev/null || echo 0)
    TOTAL=$((TOTAL + NEW))
    Q=$((Q + BATCH))
    ELAPSED=$((SECONDS - START))
    RATE=$((TOTAL / (ELAPSED + 1)))
    echo "  q=$Q rels=$TOTAL/${TARGET} ${RATE}/sec ${ELAPSED}s"
done

SIEVE_TIME=$((SECONDS - START))
echo "[$(date +%H:%M:%S)] Sieve: ${SIEVE_TIME}s, ${TOTAL} rels"

[ $TOTAL -lt 500000 ] && echo "FAIL" && rm -rf "$WORKDIR" && exit 1

# ============ POST-PROCESSING with YAFU ============
# Combine relations into msieve.dat for YAFU
cat rels_*.out > msieve.dat

# Setup YAFU environment
for f in /tmp/agent-factoring-1/yafu_mod/factor/lasieve5_64/gnfs-lasieve4I*e; do
    ln -sf "$f" . 2>/dev/null
done

# Create .job file for YAFU
cp job.poly nfs.job

cat > yafu.ini << EOF
ggnfs_dir=$(pwd)/
EOF

# YAFU expects relations in 'nfs.dat' and polynomial in 'nfs.job'
mv msieve.dat nfs.dat
mv job.poly nfs.job

# Tell YAFU to skip sieving and just do filtering + LA + sqrt
echo "[$(date +%H:%M:%S)] Post-processing with YAFU (-R -nc)..."
REMAINING=$((295 - SECONDS + START))
echo "nfs($N)" | LD_LIBRARY_PATH=/usr/local/lib timeout $REMAINING $YAFU -threads 1 -seed 42 -xover 85 -R -nc 2>&1 | tail -15

TOTAL_TIME=$((SECONDS - START))
echo "[$(date +%H:%M:%S)] Total: ${TOTAL_TIME}s"
rm -rf "$WORKDIR"
