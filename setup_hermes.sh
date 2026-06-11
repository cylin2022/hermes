#!/bin/bash
# Hermes Setup Script
# Ubuntu 24.04 LTS — 160 cores, 2.2 TiB RAM
set -e

HERMES_DIR="$HOME/hermes"
CLAUDE_CONFIG="$HOME/.claude"

echo "================================================"
echo " Hermes Bioinformatics Agent — Installation"
echo "================================================"

# ── 1. Node.js 22 LTS ────────────────────────────────────────────────────────
echo ""
echo "[1/5] Installing Node.js 22 LTS..."
if ! command -v node &>/dev/null; then
    curl -fsSL https://deb.nodesource.com/setup_22.x | sudo -E bash -
    sudo apt-get install -y nodejs
    echo "      Node.js $(node --version) installed"
else
    echo "      Node.js $(node --version) already installed"
fi

# ── 2. Claude Code CLI ───────────────────────────────────────────────────────
echo ""
echo "[2/5] Installing Claude Code..."
if ! command -v claude &>/dev/null; then
    npm install -g @anthropic-ai/claude-code
    echo "      Claude Code installed"
else
    echo "      Claude Code $(claude --version) already installed"
fi

# ── 3. Python dependencies ───────────────────────────────────────────────────
echo ""
echo "[3/5] Installing Python dependencies..."
pip3 install --quiet --upgrade \
    fastmcp \
    snakemake \
    snakemake-executor-plugin-cluster-generic \
    pulp \
    pyyaml
echo "      fastmcp + snakemake installed"

# ── 4. Configure Claude Code MCP server ─────────────────────────────────────
echo ""
echo "[4/5] Configuring Claude Code MCP integration..."
mkdir -p "$CLAUDE_CONFIG"

# Add Hermes MCP server to Claude Code settings
SETTINGS="$CLAUDE_CONFIG/settings.json"
if [ -f "$SETTINGS" ]; then
    # Merge into existing settings with Python
    python3 - <<PYEOF
import json, sys
with open("$SETTINGS") as f:
    cfg = json.load(f)
cfg.setdefault("mcpServers", {})["hermes"] = {
    "command": "python3",
    "args": ["$HERMES_DIR/mcp_server.py"],
    "type": "stdio"
}
with open("$SETTINGS", "w") as f:
    json.dump(cfg, f, indent=2)
print("      MCP server added to existing settings.json")
PYEOF
else
    cat > "$SETTINGS" <<JSON
{
  "mcpServers": {
    "hermes": {
      "command": "python3",
      "args": ["$HERMES_DIR/mcp_server.py"],
      "type": "stdio"
    }
  }
}
JSON
    echo "      settings.json created"
fi

# ── 5. Systemd service for MCP server (optional background mode) ─────────────
echo ""
echo "[5/5] Installing systemd service (hermes-mcp)..."
sudo tee /etc/systemd/system/hermes-mcp.service > /dev/null <<SERVICE
[Unit]
Description=Hermes Bioinformatics MCP Server
After=network.target

[Service]
Type=simple
User=$USER
WorkingDirectory=$HERMES_DIR
ExecStart=/usr/bin/python3 $HERMES_DIR/mcp_server.py --http 8765
Restart=on-failure
RestartSec=5

[Install]
WantedBy=multi-user.target
SERVICE

sudo systemctl daemon-reload
echo "      Service installed (use: sudo systemctl enable --now hermes-mcp)"

# ── Done ─────────────────────────────────────────────────────────────────────
echo ""
echo "================================================"
echo " Installation complete!"
echo ""
echo " Next steps:"
echo "   1. Set API key:  export ANTHROPIC_API_KEY=sk-ant-..."
echo "   2. Start Claude: claude"
echo "   3. Try:          /hermes-status"
echo "================================================"
