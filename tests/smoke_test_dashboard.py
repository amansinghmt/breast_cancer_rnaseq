#!/usr/bin/env python3

import os
import hashlib
import subprocess
import socket
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CANONICAL_FILES = sorted(
    [
        *list((ROOT / "results_v2/deseq2").glob("*")),
        *list((ROOT / "results_v2/enrichment").glob("*")),
        *list((ROOT / "results_v2/robustness").glob("*")),
        *list((ROOT / "results_v2/metadata").glob("*")),
        *list((ROOT / "figures_v2/final").glob("*")),
        *list((ROOT / "figures_v2/vector").glob("*")),
    ]
)
DOCUMENT_FILES = [
    ROOT / "docs/ONCORNA_FINAL_SCIENTIFIC_REPORT.md",
    ROOT / "docs/ONCORNA_MSC_PORTFOLIO_SUMMARY.md",
    ROOT / "docs/ONCORNA_VIVA_SHEET.md",
    ROOT / "docs/ONCORNA_FUTURE_STUDY_GUIDE.md",
]


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def available_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
        sock.bind(("127.0.0.1", 0))
        return int(sock.getsockname()[1])


def main() -> int:
    port = available_port()
    health_url = f"http://127.0.0.1:{port}/_stcore/health"
    app_source = (ROOT / "dashboard/app.py").read_text()
    launcher_source = (ROOT / "run_dashboard.sh").read_text()
    config_source = (ROOT / ".streamlit/config.toml").read_text()
    if "use_container_width" in app_source:
        raise RuntimeError("Deprecated use_container_width call remains in dashboard/app.py.")
    if 'ONCORNA_PORT="${ONCORNA_PORT:-8502}"' not in launcher_source:
        raise RuntimeError("run_dashboard.sh does not default to port 8502.")
    for setting in ('toolbarMode = "minimal"', 'showErrorDetails = "none"', "gatherUsageStats = false"):
        if setting not in config_source:
            raise RuntimeError(f"Missing Streamlit config setting: {setting}")
    for path in DOCUMENT_FILES:
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing dashboard document download: {path.name}")

    before = {path: md5(path) for path in CANONICAL_FILES}
    bare_run = subprocess.run(
        [sys.executable, "dashboard/app.py"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    if bare_run.returncode != 0:
        raise RuntimeError(f"Dashboard script execution failed:\n{bare_run.stdout}\n{bare_run.stderr}")
    env = os.environ.copy()
    env.setdefault("STREAMLIT_BROWSER_GATHER_USAGE_STATS", "false")
    process = subprocess.Popen(
        [
            sys.executable,
            "-m",
            "streamlit",
            "run",
            "dashboard/app.py",
            "--server.headless=true",
            f"--server.port={port}",
            "--server.address=127.0.0.1",
        ],
        cwd=ROOT,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    try:
        for _ in range(40):
            if process.poll() is not None:
                output = process.stdout.read() if process.stdout else ""
                raise RuntimeError(f"Dashboard exited before health check:\n{output}")
            try:
                with urllib.request.urlopen(health_url, timeout=1) as response:
                    body = response.read().decode("utf-8").strip()
                if response.status == 200 and body == "ok":
                    after = {path: md5(path) for path in CANONICAL_FILES}
                    if before != after:
                        raise RuntimeError("Dashboard modified one or more canonical outputs.")
                    print("DASHBOARD SMOKE TEST PASSED (canonical outputs unchanged)")
                    return 0
            except (OSError, urllib.error.URLError):
                time.sleep(0.25)
        raise RuntimeError("Dashboard did not become healthy within 10 seconds.")
    finally:
        process.terminate()
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait(timeout=5)


if __name__ == "__main__":
    raise SystemExit(main())
