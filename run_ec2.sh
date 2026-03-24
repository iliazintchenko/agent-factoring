#!/bin/bash
# EC2 instance: c8a.12xlarge (48 vCPUs, 96GB RAM, AMD EPYC 5th gen)
#   Instance ID: i-097c43774a5e86e69
#   IP: 44.200.192.220
#
# Launch (or reattach):  ./run_ec2.sh --host ec2-user@44.200.192.220
# Detach:                Ctrl-b d
# Reattach:              ssh -t ec2-user@44.200.192.220 'tmux attach -t factoring'
# Kill:                  ssh ec2-user@44.200.192.220 'pkill -9 -u ec2-user'
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
HOST=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --host) HOST="$2"; shift 2 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

if [ -z "$HOST" ]; then
  echo "Usage: $0 --host <user@host>"
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

# Provision and launch on remote
ssh "$HOST" "bash -s $(printf '%q %q %q' "$REPO_URL" "$GIT_USER_NAME" "$GIT_USER_EMAIL")" <<'REMOTE'
set -e
REPO_URL="$1"; GIT_USER_NAME="$2"; GIT_USER_EMAIL="$3"

source ~/.env
export CLAUDE_CODE_OAUTH_TOKEN="$CLAUDE_CODE_OAUTH_TOKEN"
export PATH="$HOME/.local/bin:$PATH"

git config --global user.name "$GIT_USER_NAME"
git config --global user.email "$GIT_USER_EMAIL"

# Install system dependencies (dnf is idempotent, always run to ensure nothing is missing)
sudo dnf install -y gcc gcc-c++ gmp-devel cmake make git hwloc-devel autoconf automake libtool tmux jq

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

# Clean up any previous runs
rm -rf /tmp/agent-factoring /tmp/inv-*

# Clone repo
REPO_DIR="/tmp/agent-factoring"
git clone "$REPO_URL" "$REPO_DIR"
git -C "$REPO_DIR" config user.name "$GIT_USER_NAME"
git -C "$REPO_DIR" config user.email "$GIT_USER_EMAIL"

# Launch single PI agent in tmux
LOG="$REPO_DIR/pi.log"
tmux new-session -s factoring -d \
  "cd $REPO_DIR && claude -p 'Read program.md and go.' --dangerously-skip-permissions --verbose --output-format stream-json 2>&1 | tee $LOG"
REMOTE

ssh -t "$HOST" 'tmux attach -t factoring'
