#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

if [[ -x ".venv/bin/python" ]]; then
  PYTHON_BIN=".venv/bin/python"
else
  PYTHON_BIN="python3"
fi

"$PYTHON_BIN" -c 'import streamlit' 2>/dev/null || {
  printf 'Streamlit is not installed. Run: %s -m pip install -r dashboard/requirements.txt\n' "$PYTHON_BIN" >&2
  exit 1
}

exec "$PYTHON_BIN" -m streamlit run dashboard/app.py --server.address=127.0.0.1
