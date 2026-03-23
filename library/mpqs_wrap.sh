#!/bin/bash
# Wrapper to make mpqs output match expected format (single factor on stdout)
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
result=$("$SCRIPT_DIR/mpqs" "$1" 2>/dev/null)
echo "$result" | awk '{print $1}'
