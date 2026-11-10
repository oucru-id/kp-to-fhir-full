"""Profile switch and the PHI-free guarantee of the genomic bundles."""
import base64
import csv
import json

import pytest

from chain import read, run_all, run_genomic
from conftest import FIXTURES, SAMPLES

PHI_COLUMNS = ["family_name", "given_name", "birth_date", "nik", "contact_number", "address",
               "latitude", "longitude", "village_code"]


@pytest.fixture(scope="module")
def generic_run(tmp_path_factory):
    return run_all(SAMPLES, tmp_path_factory.mktemp("generic"), profile="generic")


@pytest.fixture(scope="module")
def satusehat_run(tmp_path_factory):
    return run_all(SAMPLES, tmp_path_factory.mktemp("satusehat"), profile="satusehat")


def all_text(path):
    text = path.read_text(encoding="utf-8")
    bundle = json.loads(text)
    for entry in bundle.get("entry", []):
        for form in entry["resource"].get("presentedForm", []):
            text += base64.b64decode(form["data"]).decode("utf-8")
    return text


def synthetic_phi_values():
    values = set()
    for name, columns in (("patient_clinical_metadata", PHI_COLUMNS),
                          ("practitioner_metadata", ("nik", "name", "phone", "str_kki_number"))):
        with open(FIXTURES / "metadata" / f"{name}.csv", encoding="utf-8") as handle:
            lines = [ln for ln in handle if not ln.startswith("#")]
        for row in csv.DictReader(lines):
            values.update(row[c] for c in columns if row.get(c))
    return {v for v in values if len(v) > 3}


@pytest.mark.parametrize("profile", ["generic", "satusehat"])
def test_genomic_bundles_carry_no_phi(profile, generic_run, satusehat_run):
    root = generic_run if profile == "generic" else satusehat_run
    phi = synthetic_phi_values()
    assert phi, "fixture PHI list is empty"
    paths = list((root / "genomic").glob("*.genomic.fhir.json"))
    assert len(paths) == len(SAMPLES)
    for path in paths:
        text = all_text(path)
        leaked = sorted(v for v in phi if v in text)
        assert not leaked, f"{path.name} contains clinical values: {leaked}"
        types = {e["resource"]["resourceType"] for e in json.loads(path.read_text())["entry"]}
        assert types <= {"DiagnosticReport", "Observation", "Device", "Provenance"}


def test_generic_output_has_no_kemkes_content(generic_run):
    for folder in ("genomic", "clinical", "bundle"):
        for path in (generic_run / folder).glob("*.json"):
            text = path.read_text(encoding="utf-8").lower()
            assert "kemkes" not in text, f"{folder}/{path.name} mentions kemkes"
            assert "sys-ids" not in text
            assert "genenames.org" not in text, f"{folder}/{path.name} codes bacterial genes as HGNC"
            assert "kaptive.holtlab.net" not in text


def test_satusehat_output_keeps_kemkes_content(satusehat_run):
    bundle = read(satusehat_run / "bundle" / "CRE01.bundle.fhir.json")
    patient = next(e["resource"] for e in bundle["entry"] if e["resource"]["resourceType"] == "Patient")
    assert patient["meta"]["profile"] == ["https://fhir.kemkes.go.id/r4/StructureDefinition/Patient"]
    assert any(i["system"] == "https://fhir.kemkes.go.id/id/nik" for i in patient["identifier"])
    assert bundle["meta"]["profile"] == ["https://fhir.kemkes.go.id/r4/StructureDefinition/Bundle"]
    capsule = next(e["resource"] for e in bundle["entry"] if e["resource"]["id"] == "CRE01-capsule")
    assert capsule["code"]["coding"][0]["code"] == "SP000678"


def report_of(bundle):
    return next(e["resource"] for e in bundle["entry"] if e["resource"]["resourceType"] == "DiagnosticReport")


def test_generic_report_uses_ig3_profile_and_code(generic_run):
    report = report_of(read(generic_run / "genomic" / "CRE01.genomic.fhir.json"))
    assert report["meta"]["profile"] == [
        "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/genomic-report"]
    assert [c["code"] for c in report["code"]["coding"]] == ["51969-4", "81247-9"]
    assert "tag" not in report.get("meta", {})


def test_satusehat_report_keeps_v1_profile_and_code(satusehat_run):
    report = report_of(read(satusehat_run / "genomic" / "CRE01.genomic.fhir.json"))
    assert report["meta"]["profile"] == [
        "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/genomics-report"]
    assert report["code"]["coding"][0]["code"] == "81247-9"


@pytest.mark.parametrize("sample,code", [("CRE01", "1098201000112108"), ("KPC258", None),
                                         ("NOVEL01", None), ("HVKP01", None)])
def test_generic_conclusion_codes_are_snomed_or_text_only(sample, code, generic_run):
    first = report_of(read(generic_run / "genomic" / f"{sample}.genomic.fhir.json"))["conclusionCode"][0]
    if code is None:
        assert "coding" not in first
    else:
        assert first["coding"][0] == {"system": "http://snomed.info/sct", "code": code,
                                      "display": "Carbapenemase-producing Klebsiella pneumoniae (organism)"}


def test_generic_project_codes_use_published_codesystems(generic_run):
    from conftest import REPO
    published = {json.loads(p.read_text())["url"]: {c["code"] for c in json.loads(p.read_text()).get("concept", [])}
                 for p in (REPO / "assets" / "fhir").glob("CodeSystem-*.json")}
    used = {}
    for path in (generic_run / "genomic").glob("*.json"):
        def walk(node):
            if isinstance(node, dict):
                if node.get("system", "").startswith("http://terminology.spheres.id/"):
                    used.setdefault(node["system"], set()).add(node.get("code"))
                for value in node.values():
                    walk(value)
            elif isinstance(node, list):
                for value in node:
                    walk(value)
        walk(read(path))
    assert used
    for system, codes in used.items():
        assert system in published, f"{system} has no CodeSystem in assets/fhir/"
        if published[system]:
            assert codes <= published[system], (system, codes - published[system])


def test_genomic_bundle_without_performer_has_no_organization_reference(tmp_path):
    path = run_genomic("CRE01", tmp_path, profile="generic", performer=False)
    assert "Organization/" not in path.read_text()


def test_output_is_deterministic(tmp_path):
    first = run_genomic("CRE01", tmp_path / "a", profile="generic").read_text()
    second = run_genomic("CRE01", tmp_path / "b", profile="generic").read_text()
    assert first == second
