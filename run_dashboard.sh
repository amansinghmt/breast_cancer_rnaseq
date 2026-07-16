#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

PYTHON_BIN=".venv/bin/python"
ONCORNA_PORT="${ONCORNA_PORT:-8502}"

if [[ ! -x "$PYTHON_BIN" ]]; then
  printf 'OncoRNA virtual environment not found at %s. Create it and install dashboard/requirements.txt.\n' "$PYTHON_BIN" >&2
  exit 1
fi

if [[ ! "$ONCORNA_PORT" =~ ^[0-9]+$ ]] || ((ONCORNA_PORT < 1 || ONCORNA_PORT > 65535)); then
  printf 'ONCORNA_PORT must be an integer from 1 to 65535; received: %s\n' "$ONCORNA_PORT" >&2
  exit 1
fi

"$PYTHON_BIN" -c 'import streamlit' 2>/dev/null || {
  printf 'OncoRNA dashboard dependencies are missing. Run: %s -m pip install -r dashboard/requirements.txt\n' "$PYTHON_BIN" >&2
  exit 1
}

printf 'Starting OncoRNA at http://127.0.0.1:%s\n' "$ONCORNA_PORT"
exec "$PYTHON_BIN" -m streamlit run dashboard/app.py \
  --server.address=127.0.0.1 \
  --server.port="$ONCORNA_PORT"
