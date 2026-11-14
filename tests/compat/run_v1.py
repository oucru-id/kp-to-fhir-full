#!/usr/bin/env python3

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tests"))
from normalize import dumps, normalize_bundle 

FIXTURES = REPO / "tests" / "fixtures"
SCRIPTS = ["annotated_to_fhir.py", "merge_clinical_fhir.py", "clinical_metadata_parser.py", "kp_codes.py",
           "parse_plasmids.py"]
PLASMID_ARGS = ["--min_identity", "90", "--min_coverage", "80", "--chrom_min_contig", "20000"]


def baseline_commit():
    return Path(__file__).with_name("BASELINE").read_text().split()[0].strip()


def git_show(commit, path):
    return subprocess.run(["git", "-C", str(REPO), "show", f"{commit}:{path}"],
                          check=True, capture_output=True).stdout


def parse_results_source(commit):
    nf = git_show(commit, "workflows/typing.nf").decode("utf-8")
    block = nf.split("process parse_results", 1)[1]
    block = block.split("#!/usr/bin/env python3", 1)[1].split('\n    """', 1)[0]
    return textwrap.dedent(block).replace("\\\\", "\\")


def strip_marker(src, dst):
    """v1 reads CSVs with pandas, which cannot skip the '# SYNTHETIC' line."""
    lines = Path(src).read_text(encoding="utf-8").splitlines(keepends=True)
    Path(dst).write_text("".join(ln for ln in lines if not ln.startswith("#")), encoding="utf-8")


def run(cmd, cwd):
    return subprocess.run([str(c) for c in cmd], cwd=cwd, capture_output=True, text=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--python", default=sys.executable,
                        help="interpreter with pandas for the baseline scripts")
    parser.add_argument("--commit", default=None, help="baseline commit (default: tests/compat/BASELINE)")
    parser.add_argument("--out", default=str(REPO / "tests" / "golden" / "v1"))
    args = parser.parse_args()

    commit = args.commit or baseline_commit()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((FIXTURES / "samples.json").read_text())

    with tempfile.TemporaryDirectory(prefix="kpfhir-v1-") as tmp:
        tmp = Path(tmp)
        scripts = tmp / "bin"
        scripts.mkdir()
        for name in SCRIPTS:
            (scripts / name).write_bytes(git_show(commit, f"bin/{name}"))
        parse_source = parse_results_source(commit)

        meta = tmp / "metadata"
        meta.mkdir()
        for name in ("patient_clinical_metadata", "organization_metadata", "practitioner_metadata"):
            strip_marker(FIXTURES / "metadata" / f"{name}.csv", meta / f"{name}.csv")

        summary = {}
        for sample in sorted(manifest):
            info = manifest[sample]
            work = tmp / "work" / sample
            work.mkdir(parents=True)
            status = summary.setdefault(sample, {})

            kleborate = work / f"{sample}_kleborate_full.txt"
            shutil.copy(FIXTURES / "kleborate" / kleborate.name, kleborate)
            script = work / "parse_results.py"
            script.write_text(parse_source.replace("${kleborate_txt}", kleborate.name)
                              .replace("${sample_id}", sample))
            proc = run([args.python, script], work)
            status["parse_exit"] = proc.returncode
            if proc.returncode != 0:
                continue
            typing = work / f"{sample}_typing.json"
            outputs = [(typing, f"{sample}.typing.json"),
                       (work / f"{sample}_lineage.json", f"{sample}.lineage.json")]

            plasmid = None
            if info["plasmid"]:
                recon = FIXTURES / "plasmid" / sample
                plasmid = work / f"{sample}_plasmid.json"
                proc = run([args.python, scripts / "parse_plasmids.py", "--sample_id", sample,
                            "--typing_json", typing.name, "--assembly", FIXTURES / "assembly" / f"{sample}.fasta",
                            "--assembly_type", info["platform"],
                            "--contig_report", recon / "contig_report.txt",
                            "--mobtyper", recon / "mobtyper_results.txt", "--paf", recon / "amr_hits.paf",
                            "--recon_status", recon / "recon_status.txt", *PLASMID_ARGS,
                            "--output", plasmid.name], work)
                status["plasmid_exit"] = proc.returncode
                if proc.returncode != 0:
                    print(proc.stdout, proc.stderr, file=sys.stderr)
                    sys.exit(f"baseline parse_plasmids.py failed for {sample}")
                outputs.append((plasmid, f"{sample}.plasmid.json"))

            genomic = work / f"{sample}.fhir.json"
            cmd = [args.python, scripts / "annotated_to_fhir.py", "--input", typing.name, "--sample_id", sample,
                   "--output", genomic.name, "--organization_metadata", meta / "organization_metadata.csv",
                   "--cgmlst_file", FIXTURES / "cgmlst" / "results_alleles.tsv"]
            if plasmid:
                cmd += ["--plasmid_json", plasmid.name]
            if info["qc"]:
                cmd += ["--qc_json", FIXTURES / "qc" / f"{sample}_assembly_qc.json"]
            proc = run(cmd, work)
            status["genomic_exit"] = proc.returncode
            if proc.returncode != 0:
                print(proc.stdout, proc.stderr, file=sys.stderr)
                sys.exit(f"baseline annotated_to_fhir.py failed for {sample}")

            merged = work / f"{sample}.merged.fhir.json"
            proc = run([args.python, scripts / "merge_clinical_fhir.py", "--input", genomic.name,
                        "--output", merged.name,
                        "--patient_metadata", meta / "patient_clinical_metadata.csv",
                        "--organization_metadata", meta / "organization_metadata.csv",
                        "--practitioner_metadata", meta / "practitioner_metadata.csv"], work)
            status["merge_exit"] = proc.returncode
            if proc.returncode != 0:
                print(proc.stdout, proc.stderr, file=sys.stderr)
                sys.exit(f"baseline merge_clinical_fhir.py failed for {sample}")

            for src, name in outputs:
                shutil.copy(src, out_dir / name)
            for src, suffix in ((genomic, "genomic"), (merged, "merged")):
                bundle = normalize_bundle(json.loads(src.read_text()))
                (out_dir / f"{sample}.{suffix}.json").write_text(dumps(bundle))

        (out_dir / "summary.json").write_text(
            json.dumps({"baseline_commit": commit, "samples": summary}, indent=2, sort_keys=True) + "\n")
    print(f"v1 goldens written to {out_dir} from {commit}")


if __name__ == "__main__":
    main()
