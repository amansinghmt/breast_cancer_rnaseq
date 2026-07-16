#!/usr/bin/env python3

import os
import hashlib
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PORT = 8502
HEALTH_URL = f"http://127.0.0.1:{PORT}/_stcore/health"
CANONICAL_FILES = [
    ROOT / "results_v2/deseq2/deseq2_paired_v2_results.tsv",
    ROOT / "results_v2/enrichment/hallmark_gsea_paired_v2.tsv",
    ROOT / "results_v2/enrichment/go_bp_ora_paired_v2.tsv",
    ROOT / "results_v2/enrichment/go_bp_ora_tumor_higher_paired_v2.tsv",
    ROOT / "results_v2/enrichment/go_bp_ora_normal_higher_paired_v2.tsv",
    *[ROOT / f"figures_v2/final/F{i:02d}.png" for i in range(1, 8)],
]


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
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
            f"--server.port={PORT}",
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
                with urllib.request.urlopen(HEALTH_URL, timeout=1) as response:
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
