#!/bin/bash
# best_factor.sh — Use the best available factoring method
# For small N: ECM. For medium N: try LGSH then ECM fallback.
N="$1"
DIGITS=${#N}

# Try ECM first (best for all sizes up to ~65 digits)
timeout 290 ./library/ecm_factor "$N" 2>/dev/null
if [ $? -eq 0 ]; then exit 0; fi

# Try LGSH (MPQS) — good for 30-60 digits
timeout 290 ./library/lgsh "$N" 2>/dev/null
if [ $? -eq 0 ]; then exit 0; fi

# Try HSD — good for 30-55 digits
timeout 290 ./library/hsd "$N" 2>/dev/null
if [ $? -eq 0 ]; then exit 0; fi

exit 1
