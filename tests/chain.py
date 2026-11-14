
import json
from pathlib import Path

from conftest import FIXTURES, MANIFEST, run_script

ORG_ID = "99000001"
PIPELINE_VERSION = "2.0.0"
PLASMID_ARGS = ["--min_identity", 90, "--min_coverage", 80, "--chrom_min_contig", 20000]


def run_typing(sample, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    run_script("kp_parse_kleborate.py", "--kleborate", FIXTURES / "kleborate" / f"{sample}_kleborate_full.txt",
               "--sample_id", sample, cwd=outdir)
    return outdir / f"{sample}_typing.json"


def run_plasmid(sample, typing, outdir):
    recon = FIXTURES / "plasmid" / sample
    out = Path(outdir) / f"{sample}_plasmid.json"
    run_script("parse_plasmids.py", "--sample_id", sample, "--typing_json", typing,
               "--assembly", FIXTURES / "assembly" / f"{sample}.fasta",
               "--assembly_type", MANIFEST[sample]["platform"], "--contig_report", recon / "contig_report.txt",
               "--mobtyper", recon / "mobtyper_results.txt", "--paf", recon / "amr_hits.paf",
               "--recon_status", recon / "recon_status.txt", *PLASMID_ARGS, "--output", out, cwd=outdir)
    return out


def genomic_inputs(sample, workdir):
    typing = run_typing(sample, workdir)
    args = ["--typing_json", typing, "--cgmlst_file", FIXTURES / "cgmlst" / "results_alleles.tsv"]
    if MANIFEST[sample]["plasmid"]:
        args += ["--plasmid_json", run_plasmid(sample, typing, workdir)]
    if MANIFEST[sample]["qc"]:
        args += ["--qc_json", FIXTURES / "qc" / f"{sample}_assembly_qc.json"]
    return args


def run_genomic(sample, outdir, profile="generic", performer=True, extra=()):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    args = [*genomic_inputs(sample, outdir / "work"), "--sample_id", sample,
            "--output", outdir / f"{sample}.genomic.fhir.json", "--fhir_profile", profile,
            "--pipeline_version", PIPELINE_VERSION, *extra]
    if performer:
        args += ["--performer_org_id", ORG_ID]
    run_script("kpfhir_genomic.py", *args, cwd=outdir)
    return outdir / f"{sample}.genomic.fhir.json"


def run_clinical(samples, outdir, profile="generic", organization=True, practitioner=True, check=True):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    ids = outdir / "sample_ids.txt"
    ids.write_text("\n".join(samples) + "\n")
    args = ["--patient_metadata", FIXTURES / "metadata" / "patient_clinical_metadata.csv",
            "--sample_ids", ids, "--fhir_profile", profile, "--outdir", outdir]
    if organization:
        args += ["--organization_metadata", FIXTURES / "metadata" / "organization_metadata.csv"]
    if practitioner:
        args += ["--practitioner_metadata", FIXTURES / "metadata" / "practitioner_metadata.csv"]
    return run_script("clinical_to_fhir.py", *args, cwd=outdir, check=check)


def run_link(genomic_files, clinical_dir, outdir, profile="generic", extra=(), check=True):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    return run_script("link_fhir_bundle.py", "--genomic", *genomic_files, "--clinical_dir", clinical_dir,
                      "--fhir_profile", profile, "--outdir", outdir, *extra, cwd=outdir, check=check)


def run_all(samples, root, profile="generic"):
    root = Path(root)
    genomic = [run_genomic(s, root / "genomic", profile) for s in samples]
    run_clinical(samples, root / "clinical", profile)
    run_link(genomic, root / "clinical", root / "bundle", profile)
    return root


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))
