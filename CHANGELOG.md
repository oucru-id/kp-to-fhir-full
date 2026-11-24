# Changelog

## [2.0.1]

Panel rules revised against two benchmark datasets: Long et al. 2017 (486 Houston isolates passing
assembly QC, MICs read with current CLSI breakpoints) and Tamma et al. 2019 (40 isolates, broth
microdilution).

- **Amikacin.** `aac(6')-Ib-cr` no longer predicts amikacin resistance (8 of 82 carriers resistant).
  `aac(6')-Ib'` is now matched explicitly and still predicts amikacin resistance (120 of 127
  resistant); before, it only matched because the prime was stripped. Houston specificity 75% to
  97%, sensitivity 91% to 86%.
- **Ciprofloxacin.** `aac(6')-Ib-cr` alone no longer predicts resistance, including when Kleborate
  lists it as a fluoroquinolone gene (1 of 6 resistant). `qnr` genes and QRDR mutations still do.
  Houston specificity 87% to 98%.
- **Meropenem.** The porin rule (OmpK35 and OmpK36 changes with an ESBL) now needs OmpK36 to be
  truncated; the OmpK36 loop-3 insertion (`p.134_135insGlyAsp`, `GD`, `TD`) no longer counts
  (7 of 54 resistant). Houston specificity 77% to 92%, sensitivity 95% to 91%.
- **Trimethoprim-sulfamethoxazole** needs a `dfr` gene; a `sul` gene alone (9 of 24 resistant) no
  longer predicts resistance. Houston specificity 70% to 83%; Tamma 79% to 93%.
- **Ceftazidime** uses LOINC 18893-8 (Ceftazidime [Susceptibility]); it shared 18879-7 with
  cefepime.
- **Unpolished Nanopore assemblies.** When medaka fails, the sequencing QC status is `FAIL`, and
  frameshift, insertion and deletion mutation calls are withheld from the panel and the variant
  Observations, with a QC warning giving their number. Base-call errors in unpolished assemblies
  produced false porin, MgrB and PmrB truncations.

## [2.0.0]

- **Inputs.** `--input samplesheet.csv` (`sample,platform,fastq_1,fastq_2,fasta`) is the main
  input.
- **FHIR structure.** One merged bundle per sample is replaced by three layers:
  - a PHI-free genomic `collection` Bundle (always);
  - clinical resources (with `--patient_metadata`);
  - a linked `transaction` Bundle (with `--patient_metadata`).
  The default profile is `generic`; use `--fhir_profile satusehat` for the Kemkes content.
- **Credentials.** Upload credentials are Nextflow secrets (`KPFHIR_CLIENT_ID`,
  `KPFHIR_CLIENT_SECRET`, `KPFHIR_API_KEY`). `data/access_token.json` and `data/input_sso.json`
  are no longer read, and `--api_key` is gone. Upload runs only when `--fhir_server_url` is set.
- **Parameters.**
  - Replaced: `--skip_validation` by `--skip_fhir_validation`; `--skip_upload` by "upload only with
    `--fhir_server_url`"; `genome_size_bp` and the hardcoded Flye `5.5m` by `--genome_size`.
  - Removed: `--nanoplot_env` and `--mobsuite_env`, because tools come from containers or conda.
- **cgMLST** runs only with `--cgmlst_schema`. The schema is no longer expected in the
  repository; `tools/cgmlst/fetch_schema.sh` builds it.
- `fasta` input: assemblies, including gzipped files and URLs, go straight to typing, plasmid
  analysis and cgMLST. Their plasmid confidence is `predicted (submitted assembly)`.
- `--fhir_profile generic`:
  - HL7 Genomics Reporting IG 3.0.0 `genomic-report` profile, with LOINC 51969-4 and 81247-9;
  - SNOMED CT conclusion codes where one exists, otherwise text;
  - text-only bacterial genes; v1 coded them as HGNC, which is for human genes.
- A Device and a Provenance in every genomic bundle, and a narrative on the DiagnosticReport.
- Containers and conda for every process, pinned to the 1.4.2 tool versions. Profiles: `docker`,
  `podman`, `singularity`, `apptainer`, `conda`, `mamba` and `slurm`.
- `--mobsuite_db`, or a cached `mob_init` download into `--mobsuite_db_cache`.

## [1.4.2]

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
