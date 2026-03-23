#!/bin/bash
# EC2 instance: c8a.12xlarge (48 vCPUs, 96GB RAM, AMD EPYC 5th gen)
#   Instance ID: i-097c43774a5e86e69
#   IP: 44.200.192.220
#
# Launch (or reattach):  ./run_ec2.sh --host ec2-user@44.200.192.220 --agents 3
# Detach:                Ctrl-b d
# Reattach:              ssh -t ec2-user@44.200.192.220 'tmux attach -t factoring'
# Switch agent windows:  Ctrl-b n (next) / Ctrl-b p (prev) / Ctrl-b <number>
# Kill:                  ssh ec2-user@44.200.192.220 'pkill -9 -u ec2-user'
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
HOST=""
NUM_AGENTS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --host) HOST="$2"; shift 2 ;;
    --agents) NUM_AGENTS="$2"; shift 2 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

if [ -z "$HOST" ] || [ -z "$NUM_AGENTS" ]; then
  echo "Usage: $0 --host <user@host> --agents <n>"
  exit 1
fi

# Source .env locally so GITHUB_ACCESS_TOKEN is available for URL rewriting
if [ -f "$SCRIPT_DIR/.env" ]; then
  set -a; source "$SCRIPT_DIR/.env"; set +a
fi

REPO_URL="${REPO_URL:-$(git -C "$SCRIPT_DIR" remote get-url origin)}"
GIT_USER_NAME="${GIT_USER_NAME:-$(git config user.name 2>/dev/null)}"
GIT_USER_EMAIL="${GIT_USER_EMAIL:-$(git config user.email 2>/dev/null)}"
if [ -z "$GIT_USER_NAME" ] || [ -z "$GIT_USER_EMAIL" ]; then
  echo "Error: git user.name and user.email must be set."
  echo "  git config --global user.name \"Your Name\""
  echo "  git config --global user.email \"you@example.com\""
  exit 1
fi

# Convert SSH URL to HTTPS and inject token for authenticated access on remote
if [[ "$REPO_URL" == git@github.com:* ]]; then
  REPO_URL="https://github.com/${REPO_URL#git@github.com:}"
fi
if [[ "$REPO_URL" == https://github.com/* && -n "${GITHUB_ACCESS_TOKEN:-}" && "$REPO_URL" != *@* ]]; then
  REPO_URL="${REPO_URL/https:\/\/github.com/https://${GITHUB_ACCESS_TOKEN}@github.com}"
fi

scp "$SCRIPT_DIR/.env" "$HOST":~/

# If already running, just reattach
if ssh "$HOST" "tmux has-session -t factoring 2>/dev/null"; then
  exec ssh -t "$HOST" 'tmux attach -t factoring'
fi

# Provision and launch agents on remote
ssh "$HOST" "bash -s $(printf '%q %q %q %q' "$NUM_AGENTS" "$REPO_URL" "$GIT_USER_NAME" "$GIT_USER_EMAIL")" <<'REMOTE'
set -e
NUM_AGENTS="$1"; REPO_URL="$2"; GIT_USER_NAME="$3"; GIT_USER_EMAIL="$4"

source ~/.env
export CLAUDE_CODE_OAUTH_TOKEN="$CLAUDE_CODE_OAUTH_TOKEN"
export PATH="$HOME/.local/bin:$PATH"

# Install system dependencies (dnf is idempotent, always run to ensure nothing is missing)
sudo dnf install -y gcc gcc-c++ gmp-devel cmake make git hwloc-devel python3-flask python3-requests autoconf automake libtool tmux jq

# Build gmp-ecm from source if not installed (not in Amazon Linux default repos)
if ! pkg-config --exists ecm 2>/dev/null && [ ! -f /usr/local/lib/libecm.a ]; then
  cd /tmp
  git clone --depth 1 https://gitlab.inria.fr/zimmerma/ecm.git gmp-ecm-build
  cd gmp-ecm-build
  autoreconf -i
  ./configure --disable-shared --enable-static
  make -j"$(nproc)"
  sudo make install
  sudo ldconfig
  cd /tmp && rm -rf gmp-ecm-build
fi

# Install Claude Code if not present
if ! command -v claude &> /dev/null; then
  curl -fsSL https://claude.ai/install.sh | bash
fi

# Claude Code settings (clean slate — remove memories from previous runs)
mkdir -p ~/.claude
rm -rf ~/.claude/projects
printf '%s\n' '{"permissions":{"defaultMode":"bypassPermissions"},"model":"opus[1m]","effortLevel":"max","skipDangerousModePermissionPrompt":true}' > ~/.claude/settings.json

# Ensure CLAUDE_CODE_OAUTH_TOKEN is set for all future shells (including tmux panes)
grep -q 'CLAUDE_CODE_OAUTH_TOKEN' ~/.bashrc 2>/dev/null || cat >> ~/.bashrc <<'BASHRC'
source ~/.env 2>/dev/null
export CLAUDE_CODE_OAUTH_TOKEN="$CLAUDE_CODE_OAUTH_TOKEN"
export PATH="$HOME/.local/bin:$PATH"
BASHRC

# Clone repos — each agent gets its own directory
for i in $(seq 1 "$NUM_AGENTS"); do
  REPO_DIR="/tmp/agent-factoring-$i"
  if [ -d "$REPO_DIR/.git" ]; then
    git -C "$REPO_DIR" pull
  else
    git clone "$REPO_URL" "$REPO_DIR"
  fi
  [ -n "$GIT_USER_NAME" ] && git -C "$REPO_DIR" config user.name "$GIT_USER_NAME"
  [ -n "$GIT_USER_EMAIL" ] && git -C "$REPO_DIR" config user.email "$GIT_USER_EMAIL"
done

# Launch agents in background, fully detached from this SSH session
for i in $(seq 1 "$NUM_AGENTS"); do
  REPO_DIR="/tmp/agent-factoring-$i"
  LOG="$REPO_DIR/agent.log"
  nohup bash -c "cd $REPO_DIR && claude -p 'Read program.md and go.' --dangerously-skip-permissions --verbose --output-format stream-json" > "$LOG" 2>&1 &
done

# Write a monitor script that tails all agent logs with prefixed output
cat > /tmp/monitor-agents.sh <<'MONITOR'
#!/bin/bash
for log in /tmp/agent-factoring-*/agent.log; do
  agent=$(echo "$log" | grep -o 'agent-factoring-[0-9]*' | grep -o '[0-9]*')
  (tail -f "$log" | jq --unbuffered -r '
    select(.type=="assistant" and .message.content)
    | .message.content[]
    | select(.type=="text")
    | "[agent-'"$agent"'] " + .text
  ' 2>/dev/null) &
done
wait
MONITOR
chmod +x /tmp/monitor-agents.sh
tmux new-session -s factoring -d /tmp/monitor-agents.sh
REMOTE

ssh -t "$HOST" 'tmux attach -t factoring'
