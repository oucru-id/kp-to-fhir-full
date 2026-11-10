"""parse_plasmids.py: unchanged v1 output, plus the confidence label of fasta input."""
import json

import pytest

from chain import run_plasmid, run_typing
from conftest import FIXTURES, GOLDEN, MANIFEST, run_script

PLASMID_SAMPLES = sorted(s for s, info in MANIFEST.items() if info["plasmid"])


@pytest.mark.parametrize("sample", PLASMID_SAMPLES)
def test_plasmid_json_matches_v1(sample, tmp_path):
    produced = run_plasmid(sample, run_typing(sample, tmp_path), tmp_path)
    assert produced.read_bytes() == (GOLDEN / "v1" / f"{sample}.plasmid.json").read_bytes()


def test_submitted_assembly_confidence(tmp_path):
    typing = run_typing("CRE01", tmp_path)
    recon = FIXTURES / "plasmid" / "CRE01"
    out = tmp_path / "fasta.json"
    run_script("parse_plasmids.py", "--sample_id", "CRE01", "--typing_json", typing,
               "--assembly", FIXTURES / "assembly" / "CRE01.fasta", "--assembly_type", "fasta",
               "--contig_report", recon / "contig_report.txt", "--mobtyper", recon / "mobtyper_results.txt",
               "--paf", recon / "amr_hits.paf", "--recon_status", recon / "recon_status.txt",
               "--output", out, cwd=tmp_path)
    data = json.loads(out.read_text())
    assert data["assembly_type"] == "fasta"
    assert data["confidence"] == "predicted (submitted assembly)"
    assert data["gene_locations"]["KPC-2"]["location"] == "plasmid"
