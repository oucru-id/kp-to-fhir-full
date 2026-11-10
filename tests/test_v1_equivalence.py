"""P2 gate: --fhir_profile satusehat reproduces v1.4.2 except for the allowlist.

The v1 goldens (tests/golden/v1/) were produced by tests/compat/run_v1.py
from the baseline commit. KP_001 / KP_002 are excluded: v1 merged both into
sample "KP" (see tests/test_ids.py).
"""
import pytest

from chain import read, run_all
from compat.allowlist import expected_v2, without_new_resources
from conftest import GOLDEN, SAMPLES
from normalize import dumps, normalize_bundle

BUGGY_IN_V1 = {"KP_001", "KP_002"}
WITH_CLINICAL = {"CRE01", "HVKP01", "KPC258", "NOVEL01"}


@pytest.fixture(scope="module")
def satusehat_run(tmp_path_factory):
    return run_all(SAMPLES, tmp_path_factory.mktemp("satusehat"), profile="satusehat")


def observations(bundle):
    return {"entry": [e for e in bundle["entry"] if e["resource"]["resourceType"] == "Observation"]}


@pytest.mark.parametrize("sample", sorted(WITH_CLINICAL))
def test_linked_bundle_matches_v1_within_allowlist(sample, satusehat_run):
    actual = normalize_bundle(read(satusehat_run / "bundle" / f"{sample}.bundle.fhir.json"))
    expected = expected_v2(read(GOLDEN / "v1" / f"{sample}.merged.json"))
    assert dumps(without_new_resources(actual)) == dumps(expected)


@pytest.mark.parametrize("sample", sorted(set(SAMPLES) - BUGGY_IN_V1))
def test_genomic_observations_match_v1(sample, satusehat_run):
    genomic = normalize_bundle(read(satusehat_run / "genomic" / f"{sample}.genomic.fhir.json"))
    v1 = expected_v2(read(GOLDEN / "v1" / f"{sample}.genomic.json"))
    assert dumps(observations(genomic)) == dumps(observations(v1))


@pytest.mark.parametrize("sample", sorted(set(SAMPLES) - BUGGY_IN_V1))
def test_report_matches_v1_report(sample, satusehat_run):
    genomic = normalize_bundle(read(satusehat_run / "genomic" / f"{sample}.genomic.fhir.json"))
    v1 = expected_v2(read(GOLDEN / "v1" / f"{sample}.merged.json"))
    report = next(e["resource"] for e in genomic["entry"] if e["resource"]["resourceType"] == "DiagnosticReport")
    v1_report = next(e["resource"] for e in v1["entry"] if e["resource"]["resourceType"] == "DiagnosticReport")
    for key in ("conclusion", "conclusionCode", "presentedForm", "code", "category", "result"):
        assert dumps(report[key]) == dumps(v1_report[key]), key


def test_samples_without_clinical_rows_are_not_linked(satusehat_run):
    summary = (satusehat_run / "bundle" / "link_summary.tsv").read_text().splitlines()
    status = dict(line.split("\t")[:2] for line in summary[1:])
    assert status["SUS01"] == "missing_clinical"
    assert not (satusehat_run / "bundle" / "SUS01.bundle.fhir.json").exists()
