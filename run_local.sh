#!/bin/bash
# Run a single agent locally with tmux showing human-readable messages
set -e

LOG="agent.log"

claude -p 'Read program.md and go.' --effort max --dangerously-skip-permissions --verbose --output-format stream-json > "$LOG" 2>&1 &

tmux new-session -s factoring -d \
  "tail -f $LOG | jq --unbuffered -r 'select(.type==\"assistant\" and .message.content) | .message.content[] | select(.type==\"text\") | .text' 2>/dev/null"
tmux attach -t factoring
