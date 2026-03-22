#!/bin/bash
# Production GNFS factoring for 90-digit semiprimes
# Uses pre-computed polynomials + GGNFS I12 siever + msieve post-processing
# Total target: <295s on idle machine, <295s under moderate load with 1.2M rels
#
# Usage: timeout 295 bash library/gnfs_90d.sh <N> <poly_index>
# Example: timeout 295 bash library/gnfs_90d.sh <90d_number> 1

set -e

N="$1"
IDX="${2:-0}"

if [ -z "$N" ]; then
    echo "Usage: $0 <number> <poly_index>"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
POLY_FILE="$SCRIPT_DIR/library/gnfs_polys/90d_${IDX}.job"
SIEVER="/tmp/agent-factoring-1/yafu_mod/factor/lasieve5_64/gnfs-lasieve4I12e"
MSIEVE="/tmp/agent-factoring-1/library/msieve"

if [ ! -f "$POLY_FILE" ]; then
    echo "Error: polynomial file $POLY_FILE not found"
    exit 1
fi

WORKDIR=$(mktemp -d /tmp/gnfs_90d_XXXXXX)
cd "$WORKDIR"

echo "[$(date +%H:%M:%S)] GNFS 90d[$IDX] - $(echo $N | head -c 20)..."
START=$SECONDS

# ============ SIEVING ============
cp "$POLY_FILE" job.poly

# Get alim from poly for start_q
START_Q=$(grep "^alim:" job.poly | awk '{print int($2/4)}')
[ -z "$START_Q" ] && START_Q=210000

Q_RANGE=2000
Q_CURRENT=$START_Q
TOTAL_RELS=0
TARGET_RELS=1200000  # Reduced target (agent-7 finding)

# Leave 35s for post-processing
SIEVE_BUDGET=$((295 - 35))

while [ $((SECONDS - START)) -lt $SIEVE_BUDGET ]; do
    timeout 15 $SIEVER -f $Q_CURRENT -c $Q_RANGE -o rels_${Q_CURRENT}.out -a job.poly -n 0 2>/dev/null
    NEW_RELS=$(wc -l < rels_${Q_CURRENT}.out 2>/dev/null || echo 0)
    TOTAL_RELS=$((TOTAL_RELS + NEW_RELS))
    Q_CURRENT=$((Q_CURRENT + Q_RANGE))

    ELAPSED=$((SECONDS - START))
    RATE=$((TOTAL_RELS / (ELAPSED + 1)))

    # Print progress every ~30s
    if [ $((Q_CURRENT % 10000)) -lt $Q_RANGE ]; then
        ETA=$(( (TARGET_RELS - TOTAL_RELS) / (RATE + 1) ))
        echo "  q=$Q_CURRENT rels=$TOTAL_RELS/${TARGET_RELS} ${RATE}/sec ${ELAPSED}s ETA:${ETA}s"
    fi

    [ $TOTAL_RELS -ge $TARGET_RELS ] && break
done

SIEVE_TIME=$((SECONDS - START))
echo "[$(date +%H:%M:%S)] Sieve done: ${SIEVE_TIME}s, ${TOTAL_RELS} rels"

if [ $TOTAL_RELS -lt 500000 ]; then
    echo "FAIL: only $TOTAL_RELS rels (need 500K+)"
    rm -rf "$WORKDIR"
    exit 1
fi

# Combine relations
cat rels_*.out > msieve.dat
echo "Combined $(wc -l < msieve.dat) relation lines"

# ============ POST-PROCESSING (filter + LA + sqrt) ============
echo "[$(date +%H:%M:%S)] Post-processing with msieve..."

# Create msieve polynomial file
cat > msieve.fb << EOF
N $N
SKEW $(grep "^skew:" job.poly | awk '{print $2}')
R0 $(grep "^Y0:" job.poly | awk '{print $2}')
R1 $(grep "^Y1:" job.poly | awk '{print $2}')
A0 $(grep "^c0:" job.poly | awk '{print $2}')
A1 $(grep "^c1:" job.poly | awk '{print $2}')
A2 $(grep "^c2:" job.poly | awk '{print $2}')
A3 $(grep "^c3:" job.poly | awk '{print $2}')
A4 $(grep "^c4:" job.poly | awk '{print $2}')
EOF

echo "$N" > worktodo.ini

REMAINING=$((295 - SECONDS + START))
echo "Budget remaining: ${REMAINING}s"

timeout $REMAINING $MSIEVE -s msieve.dat -l msieve.log -nf msieve.fb -t 1 -v -nc1 -nc2 -ncr $N 2>&1 | tail -20

echo "---"
echo "[$(date +%H:%M:%S)] Total: $((SECONDS - START))s"

# Check if factoring succeeded
if grep -q "p[0-9]*:" msieve.log 2>/dev/null; then
    echo "*** FACTORS FOUND ***"
    grep "p[0-9]*:" msieve.log
fi

rm -rf "$WORKDIR"
