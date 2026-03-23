#!/bin/bash
# Direct GNFS factoring pipeline for 90-digit semiprimes
# Uses precomputed polynomials + GGNFS sievers + YAFU post-processing
# Usage: bash library/gnfs_factor_direct.sh <N> [timeout_secs] [poly_file]
#
# Designed for single-core execution within 300s budget.
# Pre-computed polys skip ~50s of poly selection.
# Budget: 0s poly + ~260s sieve + ~28s post = ~288s total

set -e

N="$1"
TIMEOUT="${2:-295}"
POLY_FILE="$3"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"

# GGNFS siever locations (try multiple known paths)
SIEVER=""
for path in \
    "$REPO_DIR/yafu/factor/lasieve5_64/bin/gnfs-lasieve4I12e" \
    "/tmp/agent-factoring-1/yafu_mod/factor/lasieve5_64/gnfs-lasieve4I12e" \
    "/tmp/agent-factoring-4/yafu/factor/lasieve5_64/bin/gnfs-lasieve4I12e" \
    "/tmp/agent-factoring-8/yafu/factor/lasieve5_64/bin/gnfs-lasieve4I12e"; do
    if [ -x "$path" ]; then
        SIEVER="$path"
        break
    fi
done

if [ -z "$SIEVER" ]; then
    echo "ERROR: No GGNFS siever found"
    exit 1
fi

# YAFU location
YAFU=""
for path in \
    "$REPO_DIR/yafu_mod/yafu" \
    "$REPO_DIR/yafu/yafu" \
    "/tmp/agent-factoring-4/yafu/yafu"; do
    if [ -x "$path" ]; then
        YAFU="$path"
        break
    fi
done

if [ -z "$YAFU" ]; then
    echo "ERROR: No YAFU binary found"
    exit 1
fi

echo "Using siever: $SIEVER"
echo "Using YAFU: $YAFU"
echo "N = $N"
echo "Timeout = $TIMEOUT"

# Create working directory
WORKDIR=$(mktemp -d /tmp/gnfs_XXXXXX)
cd "$WORKDIR"

# Set up YAFU ini pointing to GGNFS sievers
SIEVER_DIR=$(dirname "$SIEVER")
cat > yafu.ini <<EOF
ggnfs_dir=$SIEVER_DIR/
EOF

# Determine or use poly file
if [ -n "$POLY_FILE" ] && [ -f "$POLY_FILE" ]; then
    echo "Using precomputed polynomial from $POLY_FILE"
    cp "$POLY_FILE" nfs.job
else
    # Try to find precomputed poly
    DIGITS=${#N}
    for i in 0 1 2 3 4; do
        PFILE="$SCRIPT_DIR/gnfs_polys/${DIGITS}d_${i}.job"
        if [ -f "$PFILE" ]; then
            # Check if N matches
            POLY_N=$(grep "^n:" "$PFILE" | awk '{print $2}')
            if [ "$POLY_N" = "$N" ]; then
                echo "Found matching precomputed poly: $PFILE"
                cp "$PFILE" nfs.job
                break
            fi
        fi
    done
fi

if [ ! -f nfs.job ]; then
    echo "No precomputed polynomial found. Running YAFU poly select..."
    # Quick poly select with YAFU
    START_TIME=$SECONDS
    echo "nfs($N)" | timeout 60 $YAFU -threads 1 -seed 42 -xover 85 -psearch min -plan none -noecm -v 2>&1 | tee yafu_poly.log
    POLY_TIME=$((SECONDS - START_TIME))
    echo "Poly select took ${POLY_TIME}s"

    # Check if poly was generated
    if [ ! -f nfs.job ]; then
        echo "ERROR: Poly select failed"
        rm -rf "$WORKDIR"
        exit 1
    fi
fi

# Phase 1: Sieving with GGNFS
echo ""
echo "=== Phase 1: GGNFS Sieving ==="
START_SIEVE=$SECONDS
RELS_TARGET=1500000  # Slightly more than YAFU's 1.46M to be safe
Q_START=100000
Q_BATCH=20000
TOTAL_RELS=0
BATCH_NUM=0

# Create combined output file
> all_rels.dat

SIEVE_TIMEOUT=$((TIMEOUT - 35))  # Reserve 35s for post-processing
echo "Sieve timeout: ${SIEVE_TIMEOUT}s (reserving 35s for post-processing)"

while [ $((SECONDS - START_SIEVE)) -lt "$SIEVE_TIMEOUT" ] && [ "$TOTAL_RELS" -lt "$RELS_TARGET" ]; do
    Q_END=$((Q_START + Q_BATCH))
    BATCH_NUM=$((BATCH_NUM + 1))

    # Run siever for this q range
    BATCH_TIMEOUT=$((SIEVE_TIMEOUT - (SECONDS - START_SIEVE)))
    if [ "$BATCH_TIMEOUT" -lt 5 ]; then
        break
    fi

    timeout "$BATCH_TIMEOUT" "$SIEVER" -f "$Q_START" -c "$Q_BATCH" -o "batch_${BATCH_NUM}.dat" -n 0 -a nfs.job 2>/dev/null || true

    if [ -f "batch_${BATCH_NUM}.dat" ]; then
        NEW_RELS=$(wc -l < "batch_${BATCH_NUM}.dat" 2>/dev/null || echo 0)
        cat "batch_${BATCH_NUM}.dat" >> all_rels.dat
        TOTAL_RELS=$((TOTAL_RELS + NEW_RELS))
        ELAPSED=$((SECONDS - START_SIEVE))
        RATE=$((TOTAL_RELS * 1000 / (ELAPSED + 1)))
        echo "Batch $BATCH_NUM: q=$Q_START-$Q_END, +${NEW_RELS} rels, total=${TOTAL_RELS}, ${ELAPSED}s elapsed, ${RATE}/1000 rels/sec"
    fi

    Q_START=$Q_END
done

SIEVE_TIME=$((SECONDS - START_SIEVE))
echo ""
echo "Sieving complete: ${TOTAL_RELS} relations in ${SIEVE_TIME}s"

if [ "$TOTAL_RELS" -lt 100000 ]; then
    echo "ERROR: Not enough relations (${TOTAL_RELS} < 100000)"
    rm -rf "$WORKDIR"
    exit 1
fi

# Phase 2: Post-processing with YAFU
echo ""
echo "=== Phase 2: YAFU Post-Processing ==="

# Copy relations to format YAFU expects
# YAFU looks for spairs.out or the file specified in the job
cp all_rels.dat spairs.out

# Add YAFU-specific fields to job file if missing
if ! grep -q "^type:" nfs.job; then
    echo "type: gnfs" >> nfs.job
fi

REMAINING=$((TIMEOUT - (SECONDS - START_SIEVE)))
echo "Post-processing with ${REMAINING}s remaining"

# Try YAFU NFS post-processing (-nc = no sieve, -R = resume)
echo "nfs($N)" | timeout "$REMAINING" "$YAFU" -threads 1 -seed 42 -xover 85 -R -nc 2>&1 | tee yafu_post.log

# Check for factors
if grep -q "^P[0-9]" yafu_post.log; then
    echo ""
    echo "=== FACTORS FOUND ==="
    grep "^P[0-9]" yafu_post.log
    TOTAL_TIME=$((SECONDS))
    echo "Total time: ~${TOTAL_TIME}s"
else
    echo "Post-processing did not find factors"
fi

# Cleanup
cd /
rm -rf "$WORKDIR"
