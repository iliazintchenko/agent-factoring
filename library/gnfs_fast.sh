#!/bin/bash
# Fast GNFS factoring for 90-digit semiprimes using pre-computed polynomials
# Uses GGNFS siever + YAFU post-processing (filter + BL + sqrt)
#
# Key optimization: skip polynomial selection entirely (saves 50s!)
# Pre-computed polynomials in library/gnfs_polys/90d_{0-4}.job
#
# Usage: timeout 295 bash library/gnfs_fast.sh <N> [poly_file]
# Example: timeout 295 bash library/gnfs_fast.sh 1440... library/gnfs_polys/90d_1.job

set -e

N="$1"
POLY_FILE="$2"

if [ -z "$N" ] || [ -z "$POLY_FILE" ]; then
    echo "Usage: $0 <number> <poly_file>"
    echo "Example: $0 1440... library/gnfs_polys/90d_1.job"
    exit 1
fi

if [ ! -f "$POLY_FILE" ]; then
    echo "Error: polynomial file $POLY_FILE not found"
    exit 1
fi

# Configuration
SIEVER="/tmp/agent-factoring-1/yafu_mod/factor/lasieve5_64/gnfs-lasieve4I12e"
YAFU="/tmp/agent-factoring-9/yafu/yafu"
MSIEVE="/tmp/agent-factoring-1/library/msieve"
WORKDIR=$(mktemp -d /tmp/gnfs_fast_XXXXXX)
cd "$WORKDIR"

echo "[$(date +%H:%M:%S)] GNFS fast pipeline for $(echo $N | head -c 20)..."
echo "Load: $(uptime | awk '{print $NF}')"
START=$SECONDS

# ============ STEP 1: SIEVING ============
echo "[$(date +%H:%M:%S)] Sieving with GGNFS..."

# Copy polynomial file
cp "$POLY_FILE" job.poly

# Get start_q from polynomial file
START_Q=$(grep "^alim:" job.poly | awk '{print int($2/4)}')
if [ -z "$START_Q" ] || [ "$START_Q" -lt 100000 ]; then
    START_Q=210000  # Default for 90d Gimarel params
fi

Q_RANGE=2000
Q_CURRENT=$START_Q
TOTAL_RELS=0
TARGET_RELS=1460000  # Standard GGNFS target for 90d
SIEVE_BUDGET=$((295 - 35))  # Leave 35s for post-processing

while [ $((SECONDS - START)) -lt $SIEVE_BUDGET ]; do
    Q_END=$((Q_CURRENT + Q_RANGE))

    $SIEVER -f $Q_CURRENT -c $Q_RANGE -o rels_${Q_CURRENT}.out -a job.poly -n 0 2>&1 | tail -1

    NEW_RELS=$(wc -l < rels_${Q_CURRENT}.out 2>/dev/null || echo 0)
    TOTAL_RELS=$((TOTAL_RELS + NEW_RELS))
    Q_CURRENT=$Q_END

    ELAPSED=$((SECONDS - START))
    RATE=$((TOTAL_RELS / (ELAPSED + 1)))
    ETA=$(( (TARGET_RELS - TOTAL_RELS) / (RATE + 1) ))

    echo "  q=$Q_CURRENT, rels=$TOTAL_RELS/$TARGET_RELS (${RATE}/sec), elapsed=${ELAPSED}s, ETA=${ETA}s"

    if [ $TOTAL_RELS -ge $TARGET_RELS ]; then
        echo "  Sufficient relations collected!"
        break
    fi
done

SIEVE_TIME=$((SECONDS - START))
echo "[$(date +%H:%M:%S)] Sieve: ${SIEVE_TIME}s, ${TOTAL_RELS} relations"

if [ $TOTAL_RELS -lt 500000 ]; then
    echo "FAIL: insufficient relations ($TOTAL_RELS)"
    rm -rf "$WORKDIR"
    exit 1
fi

# Combine relation files
cat rels_*.out > all_rels.out
echo "Combined $(wc -l < all_rels.out) relation lines"

# ============ STEP 2: POST-PROCESSING ============
echo "[$(date +%H:%M:%S)] Post-processing..."

# Use msieve for filtering, LA, sqrt
# Convert relations to msieve format and process
# TODO: msieve integration

echo ""
echo "[$(date +%H:%M:%S)] Total elapsed: $((SECONDS - START))s"
echo "WORKDIR: $WORKDIR (relations available for manual post-processing)"

# For now, report sieve rate and leave WORKDIR for manual completion
echo ""
echo "To complete factoring manually:"
echo "  1. cd $WORKDIR"
echo "  2. Create msieve.dat from all_rels.out"
echo "  3. Run: $MSIEVE -s msieve.dat -l msieve.log -t 1 -nf job.poly $N"
