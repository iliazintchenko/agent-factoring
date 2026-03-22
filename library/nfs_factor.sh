#!/bin/bash
# NFS factoring orchestrator for 85-95 digit semiprimes
# msieve poly select + GGNFS sieve + msieve filter/LA/sqrt
# Single-threaded, target sub-300s wallclock
#
# Usage: ./nfs_factor.sh <N> [timeout_secs]

N="$1"
TIMEOUT="${2:-290}"

if [ -z "$N" ]; then
    echo "Usage: $0 <number> [timeout]"
    exit 1
fi

DIGITS=${#N}
WORKDIR=$(mktemp -d /tmp/nfs_XXXXXX)
GGNFS_DIR="/tmp/agent-factoring-3/yafu/factor/lasieve5_64/bin"
MSIEVE="/tmp/agent-factoring-3/msieve"
SIEVER="$GGNFS_DIR/gnfs-lasieve4I11e"
chmod +x "$SIEVER" 2>/dev/null

START=$(date +%s)
elapsed() { echo $(( $(date +%s) - START )); }
remaining() { echo $(( TIMEOUT - $(elapsed) )); }

cd "$WORKDIR"
echo "NFS: ${DIGITS}d, timeout ${TIMEOUT}s, workdir=$WORKDIR"

###############################################
# Phase 1: Polynomial selection via msieve
###############################################
echo "--- Phase 1: Poly select ---"

# msieve poly selection: -np1 (stage 1) + -nps (size opt) + -npr (root opt)
# For 90d degree 4: ~30-60s total
POLY_BUDGET=55

timeout $POLY_BUDGET $MSIEVE -v -np -t 1 \
    -s "$WORKDIR/msieve.dat" \
    -l "$WORKDIR/msieve.log" \
    -nf "$WORKDIR/msieve.fb" \
    "$N" "polydegree=4" 2>&1 | grep -E "coeff|score|found|poly" | tail -5
echo "Poly select: $(elapsed)s"

if [ ! -f "$WORKDIR/msieve.fb" ]; then
    echo "FAIL: No polynomial generated"
    rm -rf "$WORKDIR"
    exit 1
fi

# Display polynomial
echo "Polynomial (msieve format):"
cat "$WORKDIR/msieve.fb"

# Convert msieve.fb to GGNFS .job format for the siever
python3 << 'PYEOF'
import re, sys

with open("msieve.fb") as f:
    lines = f.readlines()

n_val = None
skew = None
coeffs = {}
rat_coeffs = {}

for line in lines:
    line = line.strip()
    if line.startswith("N "):
        n_val = line.split()[1]
    elif line.startswith("SKEW "):
        skew = line.split()[1]
    elif line.startswith("A"):
        m = re.match(r'A(\d+)\s+([-\d]+)', line)
        if m:
            coeffs[int(m.group(1))] = m.group(2)
    elif line.startswith("R"):
        m = re.match(r'R(\d+)\s+([-\d]+)', line)
        if m:
            rat_coeffs[int(m.group(1))] = m.group(2)

if not n_val or not coeffs:
    print("ERROR: couldn't parse msieve.fb", file=sys.stderr)
    sys.exit(1)

degree = max(coeffs.keys())
digits = len(n_val)

# Write GGNFS job file
with open("nfs.job", "w") as f:
    f.write(f"n: {n_val}\n")
    if skew:
        f.write(f"skew: {skew}\n")
    for i in range(degree, -1, -1):
        if i in coeffs:
            f.write(f"c{i}: {coeffs[i]}\n")
    for i in sorted(rat_coeffs.keys(), reverse=True):
        f.write(f"Y{i}: {rat_coeffs[i]}\n")

    # NFS parameters tuned by digit size
    if digits <= 85:
        f.write("rlim: 300000\nalim: 600000\nlpbr: 24\nlpba: 25\n")
        f.write("mfbr: 48\nmfba: 50\nrlambda: 2.3\nalambda: 2.3\n")
    elif digits <= 90:
        f.write("rlim: 400000\nalim: 800000\nlpbr: 25\nlpba: 26\n")
        f.write("mfbr: 50\nmfba: 52\nrlambda: 2.4\nalambda: 2.4\n")
    elif digits <= 95:
        f.write("rlim: 800000\nalim: 1200000\nlpbr: 26\nlpba: 26\n")
        f.write("mfbr: 52\nmfba: 52\nrlambda: 2.5\nalambda: 2.5\n")
    else:
        f.write("rlim: 1200000\nalim: 2000000\nlpbr: 26\nlpba: 27\n")
        f.write("mfbr: 52\nmfba: 54\nrlambda: 2.6\nalambda: 2.6\n")

print(f"GGNFS job file written ({digits}d)")
PYEOF

echo "GGNFS job file:"
cat nfs.job

###############################################
# Phase 2: Lattice sieving with GGNFS
###############################################
echo "--- Phase 2: Sieve ---"

ALIM=$(grep "^alim:" nfs.job | awk '{print $2}')
STARTQ=$(( ALIM * 2 / 5 ))  # Start at 40% of alim
QRANGE=20000

# Target relations based on digit size
if [ "$DIGITS" -le 85 ]; then
    MIN_RELS=400000
elif [ "$DIGITS" -le 90 ]; then
    MIN_RELS=700000
else
    MIN_RELS=1000000
fi

TOTAL_RELS=0
CURRENT_Q=$STARTQ

while true; do
    REM=$(remaining)
    if [ "$REM" -lt 40 ]; then
        echo "Time limit for sieving (${REM}s left)"
        break
    fi

    BATCH_TIME=$(( REM - 35 ))
    if [ "$BATCH_TIME" -gt 60 ]; then BATCH_TIME=60; fi
    if [ "$BATCH_TIME" -lt 5 ]; then break; fi

    timeout $BATCH_TIME "$SIEVER" -f $CURRENT_Q -c $QRANGE \
        -o spairs.out -n 0 -a nfs.job 2>&1 | tail -1

    if [ -f spairs.out ]; then
        NEW=$(wc -l < spairs.out)
        TOTAL_RELS=$((TOTAL_RELS + NEW))
        # Append to msieve data file (GGNFS format is compatible)
        cat spairs.out >> msieve.dat
        rm -f spairs.out
        echo "Rels: $TOTAL_RELS / $MIN_RELS ($(remaining)s left, q=$CURRENT_Q)"

        if [ $TOTAL_RELS -ge $MIN_RELS ]; then
            echo "Target reached"
            break
        fi
    fi

    CURRENT_Q=$((CURRENT_Q + QRANGE))
done

echo "Sieve done: $TOTAL_RELS rels in $(elapsed)s"

if [ $TOTAL_RELS -lt 50000 ]; then
    echo "FAIL: Not enough relations ($TOTAL_RELS)"
    rm -rf "$WORKDIR"
    exit 1
fi

###############################################
# Phase 3: Filter + LA + Sqrt via msieve
###############################################
echo "--- Phase 3: Filter + LA + Sqrt ---"

REM=$(remaining)
echo "Budget: ${REM}s"

# msieve -nc uses:
#   msieve.fb  - polynomial (factor base)
#   msieve.dat - relations
# Both should already exist

# Run filtering
echo "Running filter..."
timeout $((REM - 10)) $MSIEVE -v -nc1 -t 1 \
    -s "$WORKDIR/msieve.dat" \
    -l "$WORKDIR/msieve.log" \
    -nf "$WORKDIR/msieve.fb" \
    "$N" 2>&1 | tail -10

# Run linear algebra
REM=$(remaining)
if [ "$REM" -gt 10 ]; then
    echo "Running LA..."
    timeout $((REM - 5)) $MSIEVE -v -nc2 -t 1 \
        -s "$WORKDIR/msieve.dat" \
        -l "$WORKDIR/msieve.log" \
        -nf "$WORKDIR/msieve.fb" \
        "$N" 2>&1 | tail -5
fi

# Run square root
REM=$(remaining)
if [ "$REM" -gt 3 ]; then
    echo "Running sqrt..."
    timeout $((REM)) $MSIEVE -v -nc3 -t 1 \
        -s "$WORKDIR/msieve.dat" \
        -l "$WORKDIR/msieve.log" \
        -nf "$WORKDIR/msieve.fb" \
        "$N" 2>&1 | tail -5
fi

echo ""
echo "=== Total: $(elapsed)s ==="

# Check for factors
if grep -q "prp" "$WORKDIR/msieve.log" 2>/dev/null; then
    echo "FACTORS FOUND:"
    grep "prp" "$WORKDIR/msieve.log"
else
    echo "No factors found"
    echo "Last log entries:"
    tail -10 "$WORKDIR/msieve.log" 2>/dev/null
fi

# Cleanup
cd /tmp
rm -rf "$WORKDIR"
