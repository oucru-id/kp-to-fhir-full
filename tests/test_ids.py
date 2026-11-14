
import pytest

from chain import read, run_all, run_typing
from conftest import FIXTURES, run_script
from kpfhir_common import assembly_sample_id, fhir_id, lookup_sample


def test_lookup_is_exact():
    records = {"KP_001": "a", "KP_002": "b", "X": "x"}
    assert lookup_sample(records, "KP_001") == "a"
    assert lookup_sample(records, "KP_002") == "b"
    assert lookup_sample(records, "KP") is None
    assert lookup_sample(records, "KP_003") is None


def test_lookup_strips_platform_suffix_and_sanitises():
    assert lookup_sample({"S1": 1}, "S1_ont") == 1
    assert lookup_sample({"S1": 1}, "S1_illumina") == 1
    assert lookup_sample({"S1": 1}, "S1_ont.fastq") == 1
    assert lookup_sample({"KP-001": 1}, "KP_001") == 1


def test_lookup_rejects_ambiguous_sanitised_ids():
    with pytest.raises(ValueError):
        lookup_sample({"KP_001": 1, "KP.001": 2}, "KP-001")


def test_fhir_id_limit():
    assert fhir_id("KP_001", "patient") == "KP-001-patient"
    with pytest.raises(ValueError):
        fhir_id("x" * 70)


@pytest.mark.parametrize("name,sample", [
    ("KP_001_contigs", "KP_001"), ("KP_001_assembly", "KP_001"), ("KP_001_contigs.fasta", "KP_001"),
    ("ERR10367345_assembly", "ERR10367345"), ("KP_0010_contigs", "KP_0010"), ("contigs", "contigs"),
])
def test_cgmlst_row_names(name, sample):
    assert assembly_sample_id(name) == sample


def test_cgmlst_row_is_matched_exactly(tmp_path):
    table = tmp_path / "results_alleles.tsv"
    table.write_text("FILE\tL1\tL2\nKP_0010_contigs\t1\t2\nKP_001_assembly\t3\tLNF\n")
    typing = run_typing("KP_001", tmp_path / "typing")
    run_script("kpfhir_genomic.py", "--typing_json", typing, "--sample_id", "KP_001", "--cgmlst_file", table,
               "--output", tmp_path / "g.json", cwd=tmp_path)
    cgmlst = next(e["resource"] for e in read(tmp_path / "g.json")["entry"]
                  if e["resource"]["id"] == "KP-001-cgmlst")
    alleles = {c["code"]["coding"][0]["code"]: c.get("valueInteger", c.get("valueString"))
               for c in cgmlst["component"][3:]}
    assert alleles == {"L1": 3, "L2": "LNF"}


def report_of(bundle):
    return next(e["resource"] for e in bundle["entry"] if e["resource"]["resourceType"] == "DiagnosticReport")


@pytest.fixture(scope="module")
def kp_run(tmp_path_factory):
    return run_all(["KP_001", "KP_002"], tmp_path_factory.mktemp("kp"), profile="satusehat")


def test_underscore_samples_keep_their_own_results(kp_run):
    r1 = report_of(read(kp_run / "genomic" / "KP_001.genomic.fhir.json"))
    r2 = report_of(read(kp_run / "genomic" / "KP_002.genomic.fhir.json"))
    assert r1["id"] == "KP-001-genomic-report" and r2["id"] == "KP-002-genomic-report"
    assert r1["conclusionCode"][0]["text"] == "Extended-Spectrum Beta-Lactamase (ESBL)"
    assert r2["conclusionCode"][0]["text"] == "Susceptible"
    assert "MLST ST15" in r1["conclusion"] and "MLST ST45" in r2["conclusion"]


def test_underscore_samples_get_their_own_patient(kp_run):
    ids1 = {e["resource"]["id"] for e in read(kp_run / "bundle" / "KP-001.bundle.fhir.json")["entry"]}
    ids2 = {e["resource"]["id"] for e in read(kp_run / "bundle" / "KP-002.bundle.fhir.json")["entry"]}
    assert {"KP-001-patient", "KP-001-genomic-report"} <= ids1
    assert {"KP-002-patient", "KP-002-genomic-report"} <= ids2
    assert not any(i.startswith("KP-002") for i in ids1)
    assert not any(i.startswith("KP-001") for i in ids2)


def test_typing_file_for_another_sample_is_rejected(tmp_path):
    typing = run_typing("KP_002", tmp_path)
    proc = run_script("kpfhir_genomic.py", "--typing_json", typing, "--sample_id", "KP_001",
                      "--output", tmp_path / "x.json", cwd=tmp_path, check=False)
    assert proc.returncode != 0
    assert "not 'KP-001'" in proc.stderr


def test_qc_file_for_another_sample_is_rejected(tmp_path):
    typing = run_typing("NOVEL01", tmp_path)
    proc = run_script("kpfhir_genomic.py", "--typing_json", typing, "--sample_id", "NOVEL01",
                      "--qc_json", FIXTURES / "qc" / "HVKP01_assembly_qc.json",
                      "--output", tmp_path / "x.json", cwd=tmp_path, check=False)
    assert proc.returncode != 0
    assert "HVKP01" in proc.stderr
