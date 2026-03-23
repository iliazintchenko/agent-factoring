#!/bin/bash
# scan_pm1.sh - Scan semiprimes with Pollard p-1 and Williams p+1
# Usage: ./scan_pm1.sh [start_size] [end_size]
# Checks if any semiprime has a factor with smooth p-1 or p+1

START=${1:-80}
END=${2:-100}

echo "Scanning ${START}d to ${END}d semiprimes with p-1/p+1..."

for size in $(seq $START $END); do
    for i in 0 1 2 3 4; do
        N=$(python3 -c "import json; d=json.load(open('semiprimes.json')); print(d['${size}'][$i])")
        # p-1 with B1=1e6 (finds factors with B1-smooth p-1)
        result=$(echo "$N" | timeout 10 ecm -pm1 1e6 2>&1)
        if echo "$result" | grep -q "Factor found"; then
            factor=$(echo "$result" | grep "Factor found" | head -1)
            echo "P-1 HIT: ${size}d[$i] - $factor"
        fi
        # p+1 with B1=1e6
        result=$(echo "$N" | timeout 10 ecm -pp1 1e6 2>&1)
        if echo "$result" | grep -q "Factor found"; then
            factor=$(echo "$result" | grep "Factor found" | head -1)
            echo "P+1 HIT: ${size}d[$i] - $factor"
        fi
    done
    echo "  ${size}d done"
done
echo "Scan complete"
