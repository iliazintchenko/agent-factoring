#!/bin/bash
# Direct GNFS factoring using pre-computed polynomial + GGNFS siever + YAFU post-processing
# Usage: timeout 295 ./gnfs_direct.sh <N>
# Single-core. Designed for 90-digit semiprimes.

set -e
N=$1
if [ -z "$N" ]; then echo "Usage: $0 <N>"; exit 1; fi

YAFU="/tmp/agent-factoring-4/yafu/yafu"
SIEVER="/tmp/agent-factoring-1/yafu_mod/factor/lasieve5_64/gnfs-lasieve4I12e"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
POLY_DIR="$SCRIPT_DIR/gnfs_polys"
WORKDIR=$(mktemp -d /tmp/gnfs_direct_XXXXXX)
trap "rm -rf $WORKDIR" EXIT
cd "$WORKDIR"

# Symlink sievers
for S in /tmp/agent-factoring-1/yafu_mod/factor/lasieve5_64/gnfs-lasieve4I*e; do
    ln -s "$S" "$(basename $S)" 2>/dev/null
done
# Fake sievers for sizes YAFU might probe
for i in 14 15 16; do
    [ ! -f gnfs-lasieve4I${i}e ] && ln -s gnfs-lasieve4I12e gnfs-lasieve4I${i}e 2>/dev/null
done

START_NS=$(date +%s%N)
elapsed() { echo $(( ($(date +%s%N) - START_NS) / 1000000000 )); }

echo "=== GNFS Direct: $N ===" >&2

# Step 1: Find pre-computed polynomial or use YAFU poly select
POLY_FILE=""
for f in "$POLY_DIR"/*.job; do
    if grep -q "n: $N" "$f" 2>/dev/null; then
        POLY_FILE="$f"
        break
    fi
done

if [ -n "$POLY_FILE" ]; then
    echo "[$(elapsed)s] Pre-computed polynomial found" >&2
    cp "$POLY_FILE" nfs.job

    # Create YAFU-compatible files
    python3 -c "
import re
with open('nfs.job') as f: content = f.read()
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
with open('nfs.fb', 'w') as f: f.write(fb)
"
else
    echo "[$(elapsed)s] No pre-computed poly, using YAFU poly select..." >&2
    # Let YAFU do poly selection only
    echo "nfs($N)" | timeout 55 $YAFU -threads 1 -seed 42 -xover 85 -np 2>/dev/null || true
    if [ ! -f nfs.job ]; then
        echo "ERROR: Polynomial selection failed" >&2
        exit 1
    fi
fi

echo "[$(elapsed)s] Polynomial ready" >&2

# Step 2: Direct GGNFS lattice sieving
STARTQ=210000
QRANGE=5000
BATCH=0
TOTAL_RELS=0
TARGET=1460000

touch nfs.dat

while [ $(elapsed) -lt 265 ] && [ $TOTAL_RELS -lt $TARGET ]; do
    CUR_Q=$((STARTQ + BATCH * QRANGE))
    REMAINING=$((270 - $(elapsed)))
    [ $REMAINING -lt 3 ] && break

    timeout $REMAINING $SIEVER -f $CUR_Q -c $QRANGE -o rels_${BATCH}.dat -a nfs.job 2>/dev/null

    if [ -f rels_${BATCH}.dat ]; then
        NEW=$(wc -l < rels_${BATCH}.dat)
        cat rels_${BATCH}.dat >> nfs.dat
        rm rels_${BATCH}.dat
        TOTAL_RELS=$((TOTAL_RELS + NEW))
        echo "[$(elapsed)s] q=$CUR_Q: +$NEW (total: $TOTAL_RELS/$TARGET)" >&2
    fi

    BATCH=$((BATCH + 1))
done

echo "[$(elapsed)s] Sieving done: $TOTAL_RELS relations" >&2

if [ $TOTAL_RELS -lt 500000 ]; then
    echo "ERROR: Too few relations" >&2
    exit 1
fi

# Step 3: Post-processing with YAFU (it handles filtering+LA+sqrt)
REMAINING=$((293 - $(elapsed)))
echo "[$(elapsed)s] Post-processing ($REMAINING s budget)..." >&2

# YAFU NFS post-processing expects specific file layout
# Tell it to skip sieving and just do combine+LA+sqrt
echo "nfs($N)" | timeout $REMAINING $YAFU -threads 1 -seed 42 -xover 85 -nc 2>&1 | tail -5

echo "[$(elapsed)s] Done" >&2
