#!/bin/bash
# GNFS factoring driver for 90d semiprimes
# Uses pre-computed polynomial, GGNFS lattice siever, msieve post-processing
# Usage: timeout 295 ./gnfs_run.sh <N>
# Single-threaded.

set -e
N=$1
if [ -z "$N" ]; then echo "Usage: $0 <N>"; exit 1; fi

MSIEVE="/tmp/agent-factoring-3/msieve"
SIEVER="/tmp/ggnfs_build/gnfs-lasieve4I12e"
POLY_DIR="$(dirname "$0")/gnfs_polys"
WORKDIR=$(mktemp -d /tmp/gnfs_run_XXXXXX)
trap "rm -rf $WORKDIR" EXIT

cd "$WORKDIR"

START_NS=$(date +%s%N)
elapsed() { echo $(( ($(date +%s%N) - START_NS) / 1000000000 )); }

echo "=== GNFS: $N ===" >&2

# Step 1: Find or generate polynomial
POLY_FILE=""
for f in "$POLY_DIR"/*.job; do
    if grep -q "n: $N" "$f" 2>/dev/null; then
        POLY_FILE="$f"
        break
    fi
done

if [ -n "$POLY_FILE" ]; then
    echo "[$(elapsed)s] Using pre-computed polynomial from $POLY_FILE" >&2
    cp "$POLY_FILE" gnfs.job
else
    echo "[$(elapsed)s] No pre-computed poly, running msieve poly select..." >&2
    timeout 25 $MSIEVE -s msieve.dat -l msieve.log -t 1 -np \
        "polydegree=4" "poly_deadline=20" "$N" 2>/dev/null

    if [ ! -f msieve.dat.fb ]; then
        echo "ERROR: Polynomial selection failed" >&2
        exit 1
    fi

    # Convert msieve .fb to GGNFS job format
    python3 -c "
import re
with open('msieve.dat.fb') as f: content = f.read()
job = 'n: $N\n'
m = re.search(r'SKEW\s+([\d.e+-]+)', content)
if m: job += f'skew: {m.group(1)}\n'
for i in range(6):
    m = re.search(rf'A{i}\s+(-?\d+)', content)
    if m: job += f'c{i}: {m.group(1)}\n'
for i in range(2):
    m = re.search(rf'R{i}\s+(-?\d+)', content)
    if m: job += f'Y{i}: {m.group(1)}\n'
job += '''rlim: 350000
alim: 840000
lpbr: 25
lpba: 25
mfbr: 50
mfba: 50
rlambda: 2.400
alambda: 2.400
'''
with open('gnfs.job', 'w') as f: f.write(job)
"
    echo "[$(elapsed)s] Polynomial ready" >&2
fi

# Also create msieve fb file for post-processing
python3 -c "
with open('gnfs.job') as f: lines = f.readlines()
fb = ''
for line in lines:
    line = line.strip()
    if line.startswith('n:'): fb += f'N {line[3:]}\n'
    elif line.startswith('skew:'): fb += f'SKEW {line[6:]}\n'
    elif line.startswith('c') and ':' not in line:
        idx = line[1]
        val = line.split(':')[0][3:].strip() if ':' in line else line.split()[1]
        fb += f'A{idx} {val}\n'
    elif line.startswith('Y'):
        idx = line[1]
        val = line.split()[1]
        fb += f'R{idx} {val}\n'
# Handle c0-c4 and Y0-Y1 properly
import re
with open('gnfs.job') as f: content = f.read()
fb = ''
m = re.search(r'n:\s*(\d+)', content)
if m: fb += f'N {m.group(1)}\n'
m = re.search(r'skew:\s*([\d.]+)', content)
if m: fb += f'SKEW {m.group(1)}\n'
for i in range(6):
    m = re.search(rf'c{i}:\s*(-?\d+)', content)
    if m: fb += f'A{i} {m.group(1)}\n'
for i in range(2):
    m = re.search(rf'Y{i}:\s*(-?\d+)', content)
    if m: fb += f'R{i} {m.group(1)}\n'
with open('msieve.dat.fb', 'w') as f: f.write(fb)
"

# Step 2: GGNFS lattice sieving
echo "[$(elapsed)s] Starting GGNFS sieving..." >&2
STARTQ=210000
QRANGE=3000
BATCH=0
TOTAL_RELS=0
TARGET=1460000

touch msieve.dat

while [ $(elapsed) -lt 270 ] && [ $TOTAL_RELS -lt $TARGET ]; do
    CUR_Q=$((STARTQ + BATCH * QRANGE))
    REMAINING=$((275 - $(elapsed)))
    [ $REMAINING -lt 3 ] && break

    timeout $REMAINING $SIEVER -f $CUR_Q -c $QRANGE -o rels_${BATCH}.dat -a gnfs.job 2>/dev/null

    if [ -f rels_${BATCH}.dat ]; then
        NEW=$(wc -l < rels_${BATCH}.dat)
        cat rels_${BATCH}.dat >> msieve.dat
        rm rels_${BATCH}.dat
        TOTAL_RELS=$((TOTAL_RELS + NEW))
        echo "[$(elapsed)s] q=$CUR_Q: +$NEW rels (total: $TOTAL_RELS/$TARGET)" >&2
    fi

    BATCH=$((BATCH + 1))
done

echo "[$(elapsed)s] Sieving done: $TOTAL_RELS relations" >&2

if [ $TOTAL_RELS -lt 100000 ]; then
    echo "ERROR: Too few relations ($TOTAL_RELS)" >&2
    exit 1
fi

# Step 3: msieve post-processing (filtering + linear algebra + square root)
REMAINING=$((293 - $(elapsed)))
if [ $REMAINING -lt 10 ]; then
    echo "ERROR: Not enough time for post-processing ($REMAINING s)" >&2
    exit 1
fi

echo "[$(elapsed)s] Post-processing..." >&2
timeout $REMAINING $MSIEVE -s msieve.dat -l msieve.log -t 1 -nc "$N" 2>/dev/null
echo "[$(elapsed)s] Post-processing done" >&2

# Extract factors
FACTORS=$(grep "prp" msieve.log 2>/dev/null | grep -o '[0-9]\{10,\}' | sort -u)
if [ -n "$FACTORS" ]; then
    echo "$FACTORS"
    echo "[$(elapsed)s] SUCCESS" >&2
else
    echo "No factors found" >&2
    cat msieve.log >&2
    exit 1
fi
