
import json
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "bin" / "assembly_qc.py"


def run_qc(tmp_path, coverage, polishing):
    out = tmp_path / "qc.json"
    subprocess.run([sys.executable, SCRIPT, "--sample_id", "S1", "--platform", "ont", "--coverage", coverage,
                    "--min_coverage", "30", "--polishing", polishing, "--output", out], check=True)
    return json.loads(out.read_text())


@pytest.mark.parametrize("coverage,polishing,status", [
    ("64.2", "polished", "PASS"),
    ("11.4", "polished", "WARNING"),
    ("64.2", "failed", "FAIL"),
    ("11.4", "failed", "FAIL"),
])
def test_status(tmp_path, coverage, polishing, status):
    assert run_qc(tmp_path, coverage, polishing)["status"] == status
