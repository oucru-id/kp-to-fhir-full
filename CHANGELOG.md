# Changelog

## [2.0.0] — Unreleased

2.0 restructures KPtoFHIR along the lines of TBtoFHIR 2.0. Kleborate, MOB-suite and chewBBACA
are called exactly as in 1.4.2, and `--fhir_profile satusehat` reproduces the 1.4.2 FHIR output
apart from the fixes listed below (`tests/compat/allowlist.py`). Illumina reads are now
assembled with SPAdes instead of MEGAHIT, so results from Illumina reads differ from 1.4.2.

### Breaking changes

- **Inputs.** `--input samplesheet.csv` (`sample,platform,fastq_1,fastq_2,fasta`) is the main
  input. `--reads_dir` keeps the 1.4 file naming. The default `data/NGS` folder is gone.
- **Outputs.** `--results_dir` is now `--outdir`. The layout follows TBtoFHIR:
  - `assembly/`, `qc/`, `typing/`, `plasmid/`, `cgmlst/` and `reports/`;
  - `fhir/{genomic,clinical,bundle,validation,upload}`;
  - `pipeline_info/`.
- **FHIR structure.** One merged bundle per sample is replaced by three layers:
  - a PHI-free genomic `collection` Bundle (always);
  - clinical resources (with `--patient_metadata`);
  - a linked `transaction` Bundle (with `--patient_metadata`).
  The default profile is `generic`; use `--fhir_profile satusehat` for the 1.4 Kemkes content.
- **No clinical row, no placeholder.** A sample without a patient row gets the genomic bundle
  only (`missing_clinical`), not a "Patient KP-<id>" with the hardcoded organisation 100007732.
- **Credentials.** Upload credentials are Nextflow secrets (`KPFHIR_CLIENT_ID`,
  `KPFHIR_CLIENT_SECRET`, `KPFHIR_API_KEY`). `data/access_token.json` and `data/input_sso.json`
  are no longer read, and `--api_key` is gone. Upload runs only when `--fhir_server_url` is set.
- **Parameters.**
  - Replaced: `--skip_validation` by `--skip_fhir_validation`; `--skip_upload` by "upload only with
    `--fhir_server_url`"; `genome_size_bp` and the hardcoded Flye `5.5m` by `--genome_size`.
  - Removed: `--nanoplot_env` and `--mobsuite_env`, because tools come from containers or conda.
- **cgMLST** runs only with `--cgmlst_schema`. The schema is no longer expected in the
  repository; `tools/cgmlst/fetch_schema.sh` builds it.

### Fixed

- **Sample identity.**
  - v1 cut sample ids at the first underscore, so `KP_001` and `KP_002` were both reported as
    sample `KP`, with one patient and one report id.
  - Sample ids are now passed explicitly and matched exactly, in the typing, plasmid, QC and
    clinical inputs and in the cgMLST rows (v1 matched those by substring).
  - An input file for another sample is an error.
- **Encounter.** v1's DiagnosticReport and ServiceRequest referenced an Encounter that was never
  created. The clinical step now creates it, or references an existing one through
  `encounter_id`.
- **Placeholders.** No `Unknown` names, birth dates or NIKs, no run time as a specimen date, and
  no `unknown-org` / `unknown-practitioner` references.
- **NIK leading zeros.** v1 read the clinical CSV with pandas, which turned NIKs and
  administrative codes into numbers. All values are now read as text.
- **Silent failures.**
  - An unreadable metadata file, a missing `sample_id` column, duplicate or colliding ids, or
    unreadable typing / plasmid / QC JSON stop the run with a message. v1 returned `{}` or exited
    silently.
  - The merge step no longer swallows exceptions.
- **Kleborate parsing.**
  - The inline `parse_results` code is now `bin/kp_parse_kleborate.py`, with byte-identical output.
  - `point_mutations` raised a NameError for AMR-only tables and tables without rows.
  - A header-only or non-Kleborate table is now an error, instead of an all-"Unknown" result that
    reads as susceptible.
- **Determinism.**
  - uuid5 `fullUrl`s and a fixed timestamp via `KPFHIR_FIXED_TIMESTAMP` or `SOURCE_DATE_EPOCH`.
  - Gene lists in conclusions keep their order; v1 joined Python sets, whose order changed from
    run to run.
- **Codes.**
  - The pHGVS system is `http://varnomen.hgvs.org`.
  - The placeholder plasmid and QC CodeSystems on `kp-pipeline-docs.readthedocs.io` are replaced
    by published CodeSystems in `assets/fhir/`.
- **Patient names** appear only in the Patient resource; they are no longer copied into displays,
  the Specimen note or the HTML report.
- **Command-line parameters on Nextflow 26.04.** 26.04 passes every command-line value as text:
  - Every number given on the command line failed schema validation.
  - `false` counted as true. This affected `--plasmid_analysis false`, `--skip_fhir_validation
    false`, `--require_clinical false`, `--publish_clinical false`, `--fail_on_validation_error
    false` and `--fhir_upload_api_key false`.

  Values are now converted to their schema type before validation
  (`subworkflows/local/utils_params.nf`, shared with TBtoFHIR) and where the pipeline reads them.
  This works the same on Nextflow 25.04, 25.10 and 26.04. Resuming on 26.04.0 or 26.04.1 warns
  that those versions leave parameters out of the task hash.
- **Validation results of all but one batch were lost.** Files are validated in batches of
  `--fhir_validation_batch_size` (50). Every batch published `validation.json` and
  `validation_summary.tsv` to the same path, so a run with more than 50 FHIR files kept only one
  batch's results. The summaries are now merged into one `validation_summary.tsv`, and each
  batch's validator output is kept as `validation_<n>.json`.
- **Optional inputs.** An absent plasmid, QC, cgMLST or metadata file was passed to the scripts
  as if it were a real file, because the placeholder was staged in a sub-folder.

### Added

- `fasta` input: assemblies, including gzipped files and URLs, go straight to typing, plasmid
  analysis and cgMLST. Their plasmid confidence is `predicted (submitted assembly)`.
- `--fhir_profile generic`:
  - HL7 Genomics Reporting IG 3.0.0 `genomic-report` profile, with LOINC 51969-4 and 81247-9;
  - SNOMED CT conclusion codes where one exists, otherwise text;
  - text-only bacterial genes; v1 coded them as HGNC, which is for human genes.
- A Device and a Provenance in every genomic bundle, and a narrative on the DiagnosticReport.
- Project CodeSystems in `assets/fhir/`: `kp-plasmid`, `kpfhir-observation-code`,
  `pipeline-property`, `klebsiella-mlst-st`, `kaptive-k-locus` and `kaptive-o-locus`.
- HL7 FHIR validation, offline by default (`fhir/validation/validation_summary.tsv`). Errors stop
  the run for `generic`.
- `--bundle_request_method POST` (`ifNoneExist` on NIK), `--external_refs`, `--require_clinical`
  and `--identifier_system_base`.
- Containers and conda for every process, pinned to the 1.4.2 tool versions. Profiles: `docker`,
  `podman`, `singularity`, `apptainer`, `conda`, `mamba` and `slurm`.
- `--mobsuite_db`, or a cached `mob_init` download into `--mobsuite_db_cache`.
- Test profiles:
  - `test_stub`: wiring only.
  - `test`: the RefSeq HS11286 assembly, downloaded at run time.
- pytest suite:
  - synthetic fixtures;
  - v1 equivalence goldens generated from the baseline commit;
  - tests of profiles and PHI-free output, clinical resources, linking, ids and the
    classification tables.
- GitHub Actions CI:
  - the PHI guard, gitleaks, ruff and `nextflow lint`;
  - pytest on Python 3.10 and 3.12;
  - the stub run on Nextflow 25.04.6, 25.10.2 and the latest release, with command-line
    parameters of every type and rejected values;
  - the test profile under Docker and conda.
- `tools/dev/check_no_phi.py` as a pre-commit hook, `.gitleaks.toml` and Dependabot.

### Changed

- **Illumina assembly uses SPAdes 4.3.0 (`--isolate`) instead of MEGAHIT 1.2.9.**
  - SPAdes is built for single isolates, MEGAHIT for metagenomes.
  - fastp trimming and the bwa / samtools / Pilon polishing are unchanged. Contigs under 500 bp
    are still dropped.
  - SPAdes' working folder, including its corrected reads, is deleted once the contigs are
    copied out.
- Nextflow layout:
  - a thin `main.nf` and `workflows/kpfhir.nf`;
  - six subworkflows (Illumina assembly, Nanopore assembly, typing, plasmid, cgMLST, FHIR);
  - one process per module, each with a stub and `versions.yml`, collated into
    `pipeline_info/software_versions.yml`.
- Pilon runs through the bioconda `pilon` wrapper instead of `/opt/pilon/pilon-1.24.jar`. The
  bwa / samtools / Pilon commands are unchanged but now run as three processes, because no
  container holds bwa 0.7.17 with samtools 1.21.
- The plasmid step runs as two processes: MOB-suite, and minimap2 against Kleborate's CARD FASTA
  (staged from the Kleborate container). The commands are as in 1.4.2.
- `get_access_token.py` moved to `tools/satusehat/` and keeps its token in `~/.kpfhir`, mode 600.
- Requires Nextflow ≥ 25.04, and nf-schema 2.5.1.

### Removed

- `annotated_to_fhir.py`, `merge_clinical_fhir.py` and `kp_codes.py`. They are replaced by
  `kpfhir_genomic.py`, `clinical_to_fhir.py`, `link_fhir_bundle.py` and the `kpfhir_*` modules.
- The unused `create_cgmlst_report`, the commented-out SPAdes step and `get_versions.py`.
- `data/access_token.json`, `data/input_sso.json` and the `test.md` stubs. The header-only
  metadata CSVs moved to `assets/metadata_templates/`.

### Kept on purpose

- **Rifampicin mappings.** The handover listed the rifampicin drug-class mappings as TB
  leftovers. They stay, because Kleborate reports a `Rif_acquired` column (arr genes) for
  *K. pneumoniae*, and the mappings feed the drug-efficacy component and the MDR class count.
- **Unchanged from 1.4.2** (see the README limitations):
  - the clonal-complex table, `high_risk_clone` and the subspecies ST ranges;
  - the medaka model default;
  - Kleborate as the only AMR caller;
  - no species identification.

## [1.4.2] — 2026-11-02 (local, not published to GitHub)

Development notes: [docs/v1/report_pipeline_142_20261031.md](docs/v1/report_pipeline_142_20261031.md)
and [docs/v1/plan_plasmid_20261030.md](docs/v1/plan_plasmid_20261030.md).

- Plasmid analysis with MOB-suite 3.1.9:
  - each Kleborate resistance gene is placed on a plasmid, the chromosome, or reported as
    unresolved (minimap2 against Kleborate's CARD alleles);
  - plasmid Observations and a plasmid-analysis summary Observation in FHIR;
  - a cross-sample `plasmid_sharing.tsv`.
- NanoPlot replaces FastQC for Nanopore reads, with a read-coverage warning carried into a
  sequencing-QC Observation and the DiagnosticReport conclusion.
- medaka:
  - polishing works: a working samtools, and `--medaka_model` (consensus model
    `r941_min_sup_g507`);
  - `--medaka_batch_size`, and a CPU retry when GPU inference fails;
  - failed polishing is reported.
- Susceptibility panel:
  - gentamicin and amikacin rules now match Kleborate's aminoglycoside names;
  - Kleborate 3 point-mutation columns (fluoroquinolone, porin, colistin) are read.
- Typing:
  - `Bla_chr` genes are marked intrinsic;
  - Kleborate's own virulence and resistance scores are reported;
  - capsule-null loci and undetermined STs are reported;
  - clonal complexes come from a fixed table, without invented ones;
  - the unused lineage input to FHIR was removed.
- Layout: `scripts/` moved to `bin/`, the cgMLST schema moved to `db/`, trimmed reads are no
  longer published, and `--skip_validation` was added.

## [1.4.1]

Last version published on GitHub (`v1.4.1` tag).
