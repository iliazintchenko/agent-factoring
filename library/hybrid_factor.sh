#!/bin/bash
# Hybrid ECM+SIQS factoring script
# Usage: ./hybrid_factor.sh <N>
# Strategy: Try ECM briefly, then fall back to SIQS
N=$1
YAFU="/tmp/agent-factoring-2/yafu/yafu"
TIMEOUT=295

# Phase 1: ECM (5 seconds, ~1 curve with B1=100M for ~35d factors, fast)
echo "Phase 1: Quick ECM probe..." >&2
ECM_RESULT=$(echo "$N" | timeout 5 ecm -sigma 3:42 100000000 2>&1)
ECM_EXIT=$?
if echo "$ECM_RESULT" | grep -q "Factor found"; then
  FACTOR=$(echo "$ECM_RESULT" | grep "Factor found" | grep -oP '\d+$')
  echo "ECM found factor: $FACTOR" >&2
  echo "$FACTOR"
  exit 0
fi

# Phase 2: SIQS with remaining budget
REMAINING=$((TIMEOUT - 6))
echo "Phase 2: SIQS with ${REMAINING}s budget..." >&2
WORKDIR=$(mktemp -d /tmp/yafu_hybrid_XXXXXX)
cd $WORKDIR
RESULT=$(echo "siqs($N)" | timeout $REMAINING $YAFU -threads 1 -seed 42 -siqsNB 20 -siqsB 120000 2>&1)
echo "$RESULT" | tail -5
rm -rf $WORKDIR
