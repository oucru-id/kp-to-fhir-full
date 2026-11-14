
import json

import pytest

from chain import read, run_clinical, run_genomic, run_link
from conftest import SAMPLES

WITH_CLINICAL = ["CRE01", "HVKP01", "KPC258", "NOVEL01", "KP_001", "KP_002"]


@pytest.fixture(scope="module")
def parts(tmp_path_factory):
    root = tmp_path_factory.mktemp("link")
    genomic = [run_genomic(s, root / "genomic", "satusehat") for s in SAMPLES]
    run_clinical(SAMPLES, root / "clinical", "satusehat")
    return root, genomic


def all_references(node, out):
    if isinstance(node, dict):
        for key, value in node.items():
            if key == "reference":
                out.append(value)
            else:
                all_references(value, out)
    elif isinstance(node, list):
        for item in node:
            all_references(item, out)
    return out


def test_every_reference_resolves_inside_the_bundle(parts, tmp_path):
    root, genomic = parts
    run_link(genomic, root / "clinical", tmp_path, "satusehat")
    bundles = list(tmp_path.glob("*.bundle.fhir.json"))
    assert len(bundles) == len(WITH_CLINICAL)
    for path in bundles:
        bundle = read(path)
        full_urls = {e["fullUrl"] for e in bundle["entry"]}
        assert len(full_urls) == len(bundle["entry"])
        refs = all_references([e["resource"] for e in bundle["entry"]], [])
        dangling = [r for r in refs if r not in full_urls]
        assert not dangling, f"{path.name}: {dangling}"
        assert all(e["request"]["method"] == "PUT" for e in bundle["entry"])


def test_report_links_to_order_encounter_and_practitioner(parts, tmp_path):
    root, genomic = parts
    run_link(genomic, root / "clinical", tmp_path, "satusehat")
    bundle = read(tmp_path / "CRE01.bundle.fhir.json")
    by_url = {e["fullUrl"]: e["resource"] for e in bundle["entry"]}
    report = next(r for r in by_url.values() if r["resourceType"] == "DiagnosticReport")
    assert by_url[report["basedOn"][0]["reference"]]["resourceType"] == "ServiceRequest"
    assert by_url[report["encounter"]["reference"]]["resourceType"] == "Encounter"
    assert [by_url[p["reference"]]["resourceType"] for p in report["performer"]] == ["Organization", "Practitioner"]
    plasmid = next(r for r in by_url.values() if r.get("id") == "CRE01-plasmid-AA001")
    members = [by_url[m["reference"]]["id"] for m in plasmid["hasMember"]]
    assert members and all(m.startswith("CRE01-resistance-") for m in members)


def test_missing_clinical_warns_by_default_and_fails_when_required(parts, tmp_path):
    root, genomic = parts
    proc = run_link(genomic, root / "clinical", tmp_path / "warn", "satusehat")
    assert "SUS01" in proc.stdout
    proc = run_link(genomic, root / "clinical", tmp_path / "strict", "satusehat",
                    extra=["--require_clinical"], check=False)
    assert proc.returncode != 0
    assert "SUS01" in proc.stderr


def test_post_mode(parts, tmp_path):
    root, genomic = parts
    run_link(genomic, root / "clinical", tmp_path, "satusehat", extra=["--bundle_request_method", "POST"])
    bundle = read(tmp_path / "CRE01.bundle.fhir.json")
    for entry in bundle["entry"]:
        assert entry["request"]["method"] == "POST"
        assert "id" not in entry["resource"]
    patient = next(e for e in bundle["entry"] if e["resource"]["resourceType"] == "Patient")
    assert patient["request"]["ifNoneExist"] == "identifier=https://fhir.kemkes.go.id/id/nik|9901010000000001"


def test_external_patient_reference(parts, tmp_path):
    root, genomic = parts
    refs = tmp_path / "external.json"
    refs.write_text(json.dumps({"CRE01": {"Patient": "P0123456789"}}))
    run_link(genomic, root / "clinical", tmp_path / "out", "satusehat", extra=["--external_refs", refs])
    bundle = read(tmp_path / "out" / "CRE01.bundle.fhir.json")
    assert not any(e["resource"]["resourceType"] == "Patient" for e in bundle["entry"])
    refs_found = all_references([e["resource"] for e in bundle["entry"]], [])
    assert "Patient/P0123456789" in refs_found
    assert not any("CRE01-patient" in r for r in refs_found)
