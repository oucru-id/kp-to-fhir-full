
import pytest

from conftest import FIXTURES, GOLDEN, SAMPLES, run_script


@pytest.mark.parametrize("sample", SAMPLES)
def test_typing_and_lineage_match_v1_byte_for_byte(sample, tmp_path):
    run_script("kp_parse_kleborate.py", "--kleborate", FIXTURES / "kleborate" / f"{sample}_kleborate_full.txt",
               "--sample_id", sample, cwd=tmp_path)
    for kind in ("typing", "lineage"):
        produced = (tmp_path / f"{sample}_{kind}.json").read_bytes()
        assert produced == (GOLDEN / "v1" / f"{sample}.{kind}.json").read_bytes(), kind


def test_header_only_table_is_an_error_not_a_susceptible_result(tmp_path):
    proc = run_script("kp_parse_kleborate.py", "--kleborate", FIXTURES / "kleborate" / "EMPTY01_kleborate_full.txt",
                      "--sample_id", "EMPTY01", cwd=tmp_path, check=False)
    assert proc.returncode == 2
    assert "no result row" in proc.stderr
    assert not (tmp_path / "EMPTY01_typing.json").exists()


def test_table_without_strain_column_is_an_error(tmp_path):
    table = tmp_path / "x.tsv"
    table.write_text("a\tb\n1\t2\n")
    proc = run_script("kp_parse_kleborate.py", "--kleborate", table, "--sample_id", "X", cwd=tmp_path, check=False)
    assert proc.returncode == 2
    assert "not a Kleborate results table" in proc.stderr


def test_amr_only_table_has_no_point_mutation_crash(tmp_path):
    table = tmp_path / "amr.tsv"
    table.write_text("Input_file_name\tGene_symbol\tDrug_class\nX\tKPC-2\tbeta-lactam\n")
    run_script("kp_parse_kleborate.py", "--kleborate", table, "--sample_id", "X", cwd=tmp_path)
    typing = (tmp_path / "X_typing.json").read_text()
    assert '"point_mutations": []' in typing
    assert '"carbapenemase": "KPC-2"' in typing
