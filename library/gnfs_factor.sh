#!/bin/bash
# gnfs_factor.sh - Complete GNFS factoring pipeline
# Uses: YAFU poly select + GGNFS sieve + YAFU post-processing
# Designed for 85-100 digit balanced semiprimes, single-core
#
# Usage: ./gnfs_factor.sh <N> [timeout]
# Output: prints "p1\np2" on success

set -euo pipefail

N="${1:?Usage: $0 <N> [timeout]}"
TIMEOUT="${2:-295}"

YAFU="/tmp/agent-factoring-4/yafu/yafu"
SIEVER_DIR="/tmp/agent-factoring-4/yafu/factor/lasieve5_64/bin"
SIEVER="$SIEVER_DIR/gnfs-lasieve4I12e"

# Ensure sievers are executable
chmod +x "$SIEVER_DIR"/gnfs-lasieve4I*e 2>/dev/null || true

WORKDIR=$(mktemp -d /tmp/gnfs_XXXXXX)
trap "rm -rf $WORKDIR" EXIT
cd "$WORKDIR"

# Setup YAFU environment
for f in "$SIEVER_DIR"/gnfs-lasieve4I*e; do ln -sf "$f" .; done
echo "ggnfs_dir=$WORKDIR/" > yafu.ini

START=$(date +%s%N)
elapsed_ms() { echo $(( ($(date +%s%N) - START) / 1000000 )); }
elapsed_s() { echo $(( ($(date +%s%N) - START) / 1000000000 )); }

DIGITS=${#N}
echo "GNFS factoring ${DIGITS}d number" >&2

# Phase 1: Polynomial selection (budget: 20% of timeout, max 50s)
POLY_BUDGET=$(( TIMEOUT * 20 / 100 ))
[ $POLY_BUDGET -gt 50 ] && POLY_BUDGET=50

echo "Phase 1: Polynomial selection (${POLY_BUDGET}s budget)" >&2

# Run YAFU just for poly selection
echo "nfs($N)" | timeout $POLY_BUDGET "$YAFU" -threads 1 -seed 42 -xover 85 2>/dev/null || true

# Check if we got a polynomial
if [ ! -f nfs.job ]; then
    # YAFU stores polys in nfs.dat.p, need to check
    if [ -f nfs.dat.p ]; then
        # Convert msieve format to GGNFS format
        # Take the last (best) poly from nfs.dat.p
        echo "n: $N" > nfs.job
        # Extract best poly (last complete block)
        python3 -c "
import re, sys
with open('nfs.dat.p') as f:
    content = f.read()
# Find all poly blocks with e values
blocks = re.findall(r'(# norm.*?e\s+([\d.e+-]+).*?(?=# norm|\Z))', content, re.DOTALL)
if blocks:
    # Sort by e value (descending)
    best = max(blocks, key=lambda b: float(b[1]))
    block = best[0]
    for line in block.split('\n'):
        line = line.strip()
        if line.startswith(('skew:', 'c0:', 'c1:', 'c2:', 'c3:', 'c4:', 'c5:', 'Y0:', 'Y1:')):
            print(line)
" >> nfs.job 2>/dev/null || true

        # Add sieving parameters for 90d
        cat >> nfs.job <<'PARAMS'
rlim: 1200000
alim: 1200000
lpbr: 25
lpba: 25
mfbr: 50
mfba: 50
rlambda: 2.5
alambda: 2.5
PARAMS
    fi
fi

if [ ! -f nfs.job ] || ! grep -q "^c4:" nfs.job; then
    echo "ERROR: Polynomial selection failed" >&2
    exit 1
fi

POLY_TIME=$(elapsed_s)
echo "Poly select done in ${POLY_TIME}s" >&2

# Phase 2: Sieving with GGNFS
SIEVE_BUDGET=$((TIMEOUT - POLY_TIME - 35))  # Reserve 35s for post-processing
[ $SIEVE_BUDGET -lt 60 ] && { echo "ERROR: Not enough time for sieving" >&2; exit 1; }

# Target relations for 90d
TARGET_RELS=1460000
START_Q=250000

echo "Phase 2: Sieving (${SIEVE_BUDGET}s budget, target ${TARGET_RELS} rels)" >&2

timeout $SIEVE_BUDGET "$SIEVER" -f $START_Q -c 200000 -o rels.out -n 0 -a nfs.job 2>/dev/null &
SIEVE_PID=$!

# Monitor and stop when enough relations collected
while kill -0 $SIEVE_PID 2>/dev/null; do
    sleep 3
    REMAINING=$((TIMEOUT - $(elapsed_s)))
    if [ -f rels.out ]; then
        RELS=$(wc -l < rels.out)
        if [ $RELS -ge $TARGET_RELS ] && [ $REMAINING -gt 30 ]; then
            kill $SIEVE_PID 2>/dev/null || true
            break
        fi
    fi
    if [ $REMAINING -lt 35 ]; then
        kill $SIEVE_PID 2>/dev/null || true
        break
    fi
done
wait $SIEVE_PID 2>/dev/null || true

if [ ! -f rels.out ]; then
    echo "ERROR: No relations produced" >&2
    exit 1
fi

TOTAL_RELS=$(wc -l < rels.out)
SIEVE_END=$(elapsed_s)
echo "Sieve done: ${TOTAL_RELS} rels in $((SIEVE_END - POLY_TIME))s" >&2

# Phase 3: Post-processing with YAFU
REMAINING=$((TIMEOUT - $(elapsed_s)))
echo "Phase 3: Post-processing (${REMAINING}s budget)" >&2

# YAFU needs the relation data in its expected format
# Create a fresh nfs.job for YAFU resumption
echo "nfs($N)" | timeout $((REMAINING - 2)) "$YAFU" \
    -threads 1 -seed 42 -xover 85 -R -nc 2>yafu_post.log | \
    grep -E "^P[0-9]|factor" || true

# Check for factors in YAFU output
if grep -q "^P" yafu_post.log 2>/dev/null; then
    grep "^P" yafu_post.log | sed 's/^P[0-9]* = //'
elif grep -q "factor" factor.log 2>/dev/null; then
    grep "factor" factor.log
else
    echo "FAIL: Post-processing did not find factors" >&2
    exit 1
fi

TOTAL_TIME=$(elapsed_s)
echo "Total: ${TOTAL_TIME}s (poly=${POLY_TIME}s sieve=$((SIEVE_END-POLY_TIME))s post=$((TOTAL_TIME-SIEVE_END))s)" >&2
