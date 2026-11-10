"""Shared pytest helpers: run the bin/ scripts on the synthetic fixtures."""
import json
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
BIN = REPO / "bin"
FIXTURES = REPO / "tests" / "fixtures"
GOLDEN = REPO / "tests" / "golden"
FIXED_TIMESTAMP = "2026-01-15T08:00:00+00:00"

sys.path.insert(0, str(BIN))
sys.path.insert(0, str(REPO / "tests"))

MANIFEST = json.loads((FIXTURES / "samples.json").read_text())
SAMPLES = sorted(s for s in MANIFEST if s != "EMPTY01")


def script_env():
    env = dict(os.environ)
    env["KPFHIR_FIXED_TIMESTAMP"] = FIXED_TIMESTAMP
    env["PYTHONPATH"] = str(BIN)
    return env


def run_script(name, *args, cwd=None, check=True):
    proc = subprocess.run([sys.executable, str(BIN / name), *map(str, args)],
                          cwd=cwd, env=script_env(), capture_output=True, text=True)
    if check and proc.returncode != 0:
        raise AssertionError(f"{name} exited {proc.returncode}\nstdout:\n{proc.stdout}\n"
                             f"stderr:\n{proc.stderr}")
    return proc


def load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))
