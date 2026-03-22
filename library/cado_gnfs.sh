#!/bin/bash
# Custom GNFS pipeline using CADO-NFS tools
# Single-threaded, optimized for 90-digit semiprimes
#
# Pipeline: polyselect → las sieving → filter → linalg → sqrt
#
# Usage: timeout 295 bash library/cado_gnfs.sh <N>

set -e

N="$1"
if [ -z "$N" ]; then
    echo "Usage: $0 <number>"
    exit 1
fi

CADO="/tmp/agent-factoring-2/cado-nfs/build/build"
WORKDIR=$(mktemp -d /tmp/cado_gnfs_XXXXXX)
cd "$WORKDIR"

echo "[$(date +%H:%M:%S)] CADO GNFS for $N"
START=$SECONDS

# ============ POLYNOMIAL SELECTION ============
echo "[$(date +%H:%M:%S)] Step 1: Polynomial selection..."

$CADO/polyselect/polyselect \
    -N $N -degree 4 -P 10000 \
    -admin 0 -admax 50000 -incr 60 \
    -nq 256 -t 1 -keep 5 \
    > poly_raw.txt 2>&1

# Extract best polynomial (last one in output that has n: line)
python3 << 'PYEOF' > best_poly.txt
import re
with open('poly_raw.txt') as f:
    content = f.read()

# Find all size-optimized polynomials
blocks = content.split('# Size-optimized polynomial:')
best_e = -1
best_poly = None
for block in blocks[1:]:  # skip first (before any optimized poly)
    lines = block.strip().split('\n')
    poly_lines = []
    exp_e = None
    for line in lines:
        if line.startswith('n:') or line.startswith('Y0:') or line.startswith('Y1:') or line.startswith('c'):
            poly_lines.append(line)
        m = re.search(r'exp_E\s+([\d.]+)', line)
        if m:
            exp_e = float(m.group(1))
    if exp_e and exp_e > best_e and poly_lines:
        best_e = exp_e
        best_poly = poly_lines

if best_poly:
    print(f'# Best polynomial: exp_E = {best_e}')
    for line in best_poly:
        print(line)
else:
    print('# ERROR: no polynomial found')
    exit(1)
PYEOF

cat best_poly.txt
POLY_TIME=$((SECONDS - START))
echo "[$(date +%H:%M:%S)] Poly select: ${POLY_TIME}s"

# Create CADO-NFS format poly file (just extract the polynomial lines directly)
grep -v '^#' best_poly.txt > gnfs.poly
# Add skewness if missing
if ! grep -q 'skew:' gnfs.poly; then
    echo "skew: 1.0" >> gnfs.poly
fi

echo ""
echo "=== Polynomial file ==="
cat gnfs.poly

# ============ FACTOR BASE GENERATION ============
echo ""
echo "[$(date +%H:%M:%S)] Step 2: Factor base generation..."

# Parameters for 90d (GGNFS-style, larger bounds for faster single-core)
LIM0=1200000
LIM1=800000
LPB0=25
LPB1=26
MFB0=50
MFB1=52
NCURVES0=12
NCURVES1=14
I=11

# Generate factor base for algebraic side (side 1)
$CADO/sieve/makefb -poly gnfs.poly -lim $LIM1 -side 1 -out fb1.gz -t 1 2>&1 | tail -3
echo "[$(date +%H:%M:%S)] Factor base generated"

# ============ SIEVING ============
echo ""
echo "[$(date +%H:%M:%S)] Step 3: Sieving..."

SIEVE_START=$SECONDS
Q_START=$((LIM1 * 3 / 4))
Q_RANGE=5000
Q_CURRENT=$Q_START
TOTAL_RELS=0
TARGET_RELS=100000  # Rough target

while [ $TOTAL_RELS -lt $TARGET_RELS ] && [ $((SECONDS - START)) -lt 260 ]; do
    Q_END=$((Q_CURRENT + Q_RANGE))

    $CADO/sieve/las \
        -poly gnfs.poly \
        -fb1 fb1.gz \
        -lim0 $LIM0 -lim1 $LIM1 \
        -lpb0 $LPB0 -lpb1 $LPB1 \
        -mfb0 $MFB0 -mfb1 $MFB1 \
        -ncurves0 $NCURVES0 -ncurves1 $NCURVES1 \
        -I $I \
        -q0 $Q_CURRENT -q1 $Q_END \
        -sqside 1 \
        -out rels_${Q_CURRENT}.gz \
        -t 1 \
        2>&1 | tail -1

    # Count relations
    NEW_RELS=$(zcat rels_${Q_CURRENT}.gz 2>/dev/null | grep -c "^[0-9a-f]" || echo 0)
    TOTAL_RELS=$((TOTAL_RELS + NEW_RELS))
    Q_CURRENT=$Q_END

    echo "  q=$Q_CURRENT, rels=$TOTAL_RELS (target $TARGET_RELS), elapsed=$((SECONDS - START))s"
done

SIEVE_TIME=$((SECONDS - SIEVE_START))
echo "[$(date +%H:%M:%S)] Sieve: ${SIEVE_TIME}s, ${TOTAL_RELS} relations"

if [ $TOTAL_RELS -lt 50000 ]; then
    echo "FAIL: insufficient relations ($TOTAL_RELS < 50000)"
    rm -rf "$WORKDIR"
    exit 1
fi

echo ""
echo "[$(date +%H:%M:%S)] Total elapsed: $((SECONDS - START))s"
echo "WORKDIR: $WORKDIR (not cleaned for debugging)"
