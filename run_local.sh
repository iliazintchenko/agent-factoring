#!/bin/bash
# Run PI agent locally
# Detach:   Ctrl-b d
# Reattach: tmux attach -t factoring
# Kill:     tmux kill-session -t factoring; pkill -f 'claude -p'
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
LOG="$SCRIPT_DIR/pi.log"

# Clean up previous investigator dirs
rm -rf /tmp/inv-*

# If already running, just reattach
if tmux has-session -t factoring 2>/dev/null; then
  exec tmux attach -t factoring
fi

tmux new-session -s factoring -d \
  "cd $SCRIPT_DIR && claude -p 'Read program.md and go.' --dangerously-skip-permissions --verbose --output-format stream-json 2>&1 | tee $LOG | jq --unbuffered -r 'select(.type==\"assistant\" and .message.content) | .message.content[] | select(.type==\"text\") | .text'"
tmux attach -t factoring
