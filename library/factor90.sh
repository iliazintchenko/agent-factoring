#!/bin/bash
# Smart 90-digit factoring script
# Tries YAFU SIQS first, falls back to GNFS if available
# Usage: bash library/factor90.sh <N> [method: siqs|gnfs|auto]
#
# Strategy:
# - SIQS: Uses yafu_mod with VBITS=512, closnuf+1, NB=20, B=120K
# - GNFS: Uses precomputed polys + GGNFS sievers + yafu post-processing
# - auto: Checks load and picks the best approach

set -e

N="$1"
METHOD="${2:-auto}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"

# Find best YAFU binary
YAFU_MOD="$REPO_DIR/yafu_mod/yafu"
YAFU_ORIG="$REPO_DIR/yafu/yafu"
if [ -x "$YAFU_MOD" ]; then
    YAFU="$YAFU_MOD"
    echo "Using yafu_mod (VBITS=512, closnuf+1)"
elif [ -x "$YAFU_ORIG" ]; then
    YAFU="$YAFU_ORIG"
    echo "Using original yafu"
else
    echo "ERROR: No YAFU binary found"
    exit 1
fi

# Determine method
if [ "$METHOD" = "auto" ]; then
    LOAD=$(cat /proc/loadavg | awk '{print $1}')
    echo "Current 1-min load: $LOAD"
    # GNFS needs low load for good sieve rates
    if (( $(echo "$LOAD < 8" | bc -l) )); then
        METHOD="gnfs"
        echo "Auto-selected: GNFS (low load)"
    else
        METHOD="siqs"
        echo "Auto-selected: SIQS (moderate/high load)"
    fi
fi

echo ""
echo "=== Factoring $N ==="
echo "Method: $METHOD"
echo ""

case "$METHOD" in
    siqs)
        WORKDIR=$(mktemp -d /tmp/yafu_XXXXXX) && cd "$WORKDIR"
        START=$SECONDS
        echo "siqs($N)" | timeout 295 "$YAFU" -threads 1 -seed 42 -siqsNB 20 -siqsB 120000 -noopt 2>&1
        ELAPSED=$((SECONDS - START))
        echo "SIQS elapsed: ${ELAPSED}s"
        rm -rf "$WORKDIR"
        ;;
    gnfs)
        WORKDIR=$(mktemp -d /tmp/yafu_XXXXXX) && cd "$WORKDIR"

        # Setup GGNFS sievers
        SIEVER_DIR="$REPO_DIR/yafu/factor/lasieve5_64/bin"
        chmod +x $SIEVER_DIR/gnfs-lasieve4I* 2>/dev/null
        cat > yafu.ini <<EOF
ggnfs_dir=$SIEVER_DIR/
EOF

        # Find precomputed poly
        DIGITS=${#N}
        POLY_FOUND=0
        for i in 0 1 2 3 4; do
            PFILE="$SCRIPT_DIR/gnfs_polys/${DIGITS}d_${i}.job"
            if [ -f "$PFILE" ]; then
                POLY_N=$(grep "^n:" "$PFILE" | awk '{print $2}')
                if [ "$POLY_N" = "$N" ]; then
                    cp "$PFILE" nfs.job
                    POLY_FOUND=1
                    echo "Using precomputed poly from $PFILE"
                    break
                fi
            fi
        done

        if [ "$POLY_FOUND" = "0" ]; then
            echo "No precomputed poly found, running with poly select..."
        fi

        START=$SECONDS
        echo "nfs($N)" | timeout 295 "$YAFU" -threads 1 -seed 42 -xover 85 -R -v 2>&1
        ELAPSED=$((SECONDS - START))
        echo "GNFS elapsed: ${ELAPSED}s"
        rm -rf "$WORKDIR"
        ;;
    *)
        echo "Unknown method: $METHOD (use siqs, gnfs, or auto)"
        exit 1
        ;;
esac
