#!/bin/bash
# Run the agent locally with tmux: log top, costs bottom-left, agent steps bottom-right
set -e

LOG="agent.log"

tmux new-session -s factoring -d \
  "claude -p 'Read program.md and go.' --effort max --dangerously-skip-permissions --verbose --output-format stream-json 2>&1 | tee $LOG"
tmux split-window -v -t factoring \
  "tail -f $LOG | jq -r 'select(.type==\"assistant\" and .message.usage) | .message.usage | \"in: \\(.input_tokens) cache: \\(.cache_read_input_tokens) out: \\(.output_tokens)\"'"
tmux split-window -h -t factoring:0.1 \
  "tail -f $LOG | jq -r 'select(.type==\"assistant\" and .message.content) | .message.content[] | select(.type==\"text\") | .text' 2>/dev/null"
tmux select-pane -t factoring:0.0
tmux attach -t factoring
