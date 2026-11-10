# KP Genomics to FHIR Pipeline (KPtoFHIR)

A Nextflow pipeline for *Klebsiella pneumoniae* genomic analysis from raw reads or assemblies. It
assembles the genome, then reports:
- [Kleborate](https://github.com/klebgenomics/Kleborate) typing: MLST, K and O loci, virulence
  and resistance genes;
- the plasmid or chromosomal location of the resistance genes
  ([MOB-suite](https://github.com/phac-nml/mob-suite));
- [chewBBACA](https://github.com/B-UMMI/chewBBACA) cgMLST.

The results are HL7 FHIR R4 resources conforming to the
[Genomics Reporting IG 3.0.0](https://hl7.org/fhir/uv/genomics-reporting/STU3/).

**Version 2.0.0 (unreleased)**. See [CHANGELOG.md](CHANGELOG.md). 2.0 changes the inputs, the
output layout and the FHIR structure. Read the
[breaking changes](CHANGELOG.md#breaking-changes) before upgrading from 1.4.

## Key features

- **Three kinds of input**:
  - Illumina paired-end reads, assembled with fastp, MEGAHIT and Pilon.
  - Oxford Nanopore reads, assembled with chopper, NanoPlot, Flye and medaka.
  - Finished assemblies (`fasta`), which go straight to typing.
- **Kleborate 3 typing**:
  - MLST, capsule (K) and O loci, virulence loci and Kleborate's virulence and resistance scores.
  - Acquired resistance genes and point mutations, and a predicted 12-drug susceptibility panel.
  - An overall resistance classification: CRE, XDR, ESBL, AmpC, MDR and so on.
- **Plasmid analysis**:
  - MOB-suite reconstructs the plasmids, giving replicon, relaxase, mobility and cluster for each.
  - Every resistance gene is placed on a plasmid, the chromosome, or reported as unresolved.
  - A cross-sample plasmid-sharing table.
- **cgMLST** with chewBBACA and the cgMLST.org *K. pneumoniae* sensu lato scheme (2,358 loci).
- **QC**:
  - FastQC and fastp (Illumina) or NanoPlot (ONT), collected in MultiQC.
  - For ONT, a read-coverage warning that is carried into the FHIR report.
- **PHI-free genomic FHIR by default**:
  - The genomic report refers to the patient by id only.
  - Clinical resources are built in a separate, optional step and joined to the genomic report in
    a third step.
- **Two FHIR profiles**:
  - `generic` (default): plain HL7 Genomics Reporting.
  - `satusehat`: reproduces the v1.4 output for the Indonesian Ministry of Health SATUSEHAT
    platform (Kemkes profiles, NIK identifiers, SP and KP-SO/KP-RL codes).

## Requirements

You need [Nextflow](https://www.nextflow.io/) ≥ 25.04 and **one** of the following:

| Profile | Provides the tools via |
|---|---|
| `-profile docker` (or `podman`) | pinned [BioContainers](https://biocontainers.pro/) images |
| `-profile singularity` (or `apptainer`) | the same images, e.g. on HPC |
| `-profile conda` (or `mamba`) | pinned bioconda / conda-forge packages |
| none | tools already on `PATH` |

Every process uses the tool versions that produced the v1.4.2 results:

| Step | Tools |
|---|---|
| Illumina | FastQC 0.12.1, fastp 0.23.4, MEGAHIT 1.2.9, bwa 0.7.17, samtools 1.21, Pilon 1.24 |
| ONT | chopper 0.12.0, NanoPlot 1.48.0, Flye 2.9.6, medaka 2.1.0 |
| typing, plasmid, cgMLST | Kleborate 3.2.4, MOB-suite 3.1.9, minimap2 2.26, chewBBACA 3.4.2 |
| other | MultiQC 1.30, Python 3.12 (standard library only), Java 21 |

Images are pulled only for the processes that run; for example, Illumina-only runs never
download medaka. The largest images are MOB-suite (about 1.9 GB), Kleborate (0.7 GB) and
medaka (0.6 GB).

Notes for container runs:
- The Python image has no `openpyxl`, so give `.xlsx` clinical metadata as CSV, or use
  `-profile conda`.
- medaka runs on the CPU in containers. It retries on the CPU if GPU inference fails.

## Quick start

```bash
git clone https://github.com/oucru-id/kp-to-fhir-full.git
cd kp-to-fhir-full

nextflow run . -profile test_stub -stub-run   # wiring only, no tools needed
nextflow run . -profile test,docker           # HS11286 (ST11, KPC-2) RefSeq assembly, downloaded at run time
```

GitHub Actions (`.github/workflows/ci.yml`) runs on every push and pull request:
- the PHI guard, gitleaks, ruff and `nextflow lint`
- pytest
- the stub run on Nextflow 25.04 and the latest release
- the test profile under Docker and under conda

## Input

### Samplesheet (recommended)

```csv
sample,platform,fastq_1,fastq_2,fasta
KP0001,illumina,/data/KP0001_R1.fastq.gz,/data/KP0001_R2.fastq.gz,
KP0002,ont,/data/KP0002.fastq.gz,,
KP0003,fasta,,,/data/KP0003.fasta
```

```bash
nextflow run oucru-id/kp-to-fhir-full -r v2.0.0 -profile docker --input samplesheet.csv --outdir results
```

Sample ids:
- may use letters, digits, `.`, `_` and `-`, up to 40 characters;
- should be laboratory accessions, not national ID numbers, because they appear in the PHI-free
  genomic output.

A `fasta` sample may be gzipped and may be a URL. It skips assembly and read QC.

### Folder (v1.4 naming)

```bash
nextflow run . -profile docker --reads_dir /data/NGS
```

The folder files must be named `<id>_1_illumina.fastq.gz` / `<id>_2_illumina.fastq.gz` or
`<id>_ont.fastq.gz`.

### Databases

- **MOB-suite.** The database is about 3 GB. Give an initialised one with `--mobsuite_db <dir>`.
  Otherwise `mob_init` downloads it once into `--mobsuite_db_cache` (default `db/mob_suite`).
  `--plasmid_analysis false` skips the plasmid step.
- **cgMLST.** cgMLST runs only with `--cgmlst_schema <dir>`. To build the schema:

  ```bash
  tools/cgmlst/fetch_schema.sh db/cgmlst_schema    # cgMLST.org scheme 2187931 via chewBBACA PrepExternalSchema
  ```

  chewBBACA adds inferred alleles to the schema in place, as in v1, so keep one copy per
  analysis series.

## FHIR output

The output has three layers:

| Step | Output | Contains |
|---|---|---|
| genomic (always) | `fhir/genomic/<id>.genomic.fhir.json` | DiagnosticReport; Observations for resistance genes and mutations, susceptibility panel, MLST, virulence score, capsule, O antigen, cgMLST, plasmids and sequencing QC; Device; Provenance. **No patient data**: the subject is the reference `Patient/<id>-patient`. |
| clinical (with `--patient_metadata`) | `fhir/clinical/` | Patient, Specimen, Encounter and ServiceRequest per sample; Organization, Practitioner and PractitionerRole once. **Personal data.** |
| link (with `--patient_metadata`) | `fhir/bundle/<id>.bundle.fhir.json` | One FHIR `transaction` per sample, with every reference resolved inside the bundle. **Personal data.** |

```bash
nextflow run . -profile docker --input samplesheet.csv \
  --patient_metadata site/patients.csv \
  --organization_metadata site/organization.csv \
  --practitioner_metadata site/practitioner.csv \
  --performer_org_id 100000001
```

The CSV columns are described in [docs/clinical_metadata.md](docs/clinical_metadata.md). Keep
these files, and the `fhir/clinical` and `fhir/bundle` outputs, out of version control.

A sample with no row in the patient table gets a genomic bundle only:
- It is listed as `missing_clinical` in `fhir/bundle/link_summary.tsv`.
- `--require_clinical` turns that into an error.
- v1.4 invented a placeholder Patient instead.

### Profiles

- `--fhir_profile generic` (default) emits only HL7 Genomics Reporting IG 3.0.0 content:
  - the `genomic-report` profile, with LOINC 51969-4 and 81247-9;
  - SNOMED CT conclusion codes where one exists, otherwise text;
  - bacterial genes as text, not HGNC codes;
  - the project CodeSystems in `assets/fhir/` for capsule, O antigen, ST, plasmid and QC codes.
- `--fhir_profile satusehat` keeps v1.4's content:
  - Kemkes StructureDefinitions, `genomics` tags and `sys-ids.kemkes.go.id` identifiers;
  - NIK and the administrative-code and citizenship extensions;
  - SP000678/680/682 and SP000681/683/687, and KP-SO / KP-RL;
  - the v1.4 specimen defaults.

  Its output matches v1.4.2 apart from the documented fixes in `tests/compat/allowlist.py`.

### Validation

Every bundle is checked with the HL7 FHIR validator against `hl7.fhir.uv.genomics-reporting#3.0.0`
and the project CodeSystems in `assets/fhir/`.
- Results: `fhir/validation/validation_summary.tsv`.
- Validation runs offline by default (`--fhir_tx_server n/a`), so no codes are sent to a public
  terminology server. Code checks that need a server are reported as not performed.
- A validation error stops the run for `generic` and only warns for `satusehat`, because the
  validator cannot resolve the Kemkes profiles. Override with `--fail_on_validation_error true|false`.
- `--skip_fhir_validation` skips the step.

### Upload

Upload runs only when `--fhir_server_url` is given. It sends only linked bundles that passed
validation. Credentials are Nextflow secrets, never parameters or files in the repository:

```bash
nextflow secrets set KPFHIR_CLIENT_ID     '<client id>'
nextflow secrets set KPFHIR_CLIENT_SECRET '<client secret>'
nextflow run . -profile docker --input samplesheet.csv --patient_metadata site/patients.csv \
  --fhir_profile satusehat \
  --fhir_server_url https://<server>/fhir --fhir_auth_url https://<sso>/realms/<realm>
```

Further options:
- `X-Api-Key` header: set the `KPFHIR_API_KEY` secret and add `--fhir_upload_api_key`.
- `--bundle_request_method POST`: server-assigned-id transactions, with `ifNoneExist` on the NIK.
- `--external_refs`: maps samples to Patients that already exist on the server.
- `tools/satusehat/get_access_token.py`: fetches a token for manual tests. It keeps the token in
  `~/.kpfhir`, not in the repository.

For SATUSEHAT, test against the sandbox first.

## Configuration

All parameters are described in [nextflow_schema.json](nextflow_schema.json) and listed by
`nextflow run . --help`. The most used:

| Parameter | Default | Purpose |
|---|---|---|
| `--input` | | samplesheet |
| `--outdir` | `results` | output folder |
| `--fhir_profile` | `generic` | `generic` or `satusehat` |
| `--patient_metadata` | | enables the clinical and link steps |
| `--performer_org_id` | | performing laboratory, as the genomic report's performer |
| `--genome_size` | `5.5m` | Flye genome size and Nanopore coverage estimate |
| `--medaka_model` | `r941_min_sup_g507` | must match the flow cell and basecaller |
| `--min_nanopore_coverage` | `30` | below this, the report carries a quality warning |
| `--plasmid_analysis` | `true` | MOB-suite and gene location |
| `--cgmlst_schema` | | enables cgMLST |
| `--max_cpus` / `--max_memory` / `--max_time` | `8` / `16 GB` / `24h` | cap on any single task |

## Output structure

```
results/
├── pipeline_info/      execution report, timeline, trace, DAG, software_versions.yml
├── assembly/           <id>_contigs.fasta (Illumina), <id>_assembly.fasta (ONT, fasta input)
├── qc/                 fastqc/  fastp/  nanoplot/  assembly/ (<id>_assembly_qc.json)  multiqc/
├── typing/             <id>_kleborate_full.txt, <id>_typing.json, <id>_lineage.json
├── plasmid/            <id>/ (MOB-suite reports, amr_hits.paf), <id>_plasmid.json, plasmid_sharing.tsv
├── cgmlst/             results_alleles.tsv, cgmlst_statistics.txt
├── reports/            <id>.summary_report.txt
└── fhir/
    ├── genomic/        <id>.genomic.fhir.json                      (PHI-free)
    ├── clinical/       clinical resources, clinical_summary.tsv   (personal data)
    ├── bundle/         <id>.bundle.fhir.json, link_summary.tsv    (personal data)
    ├── validation/     validation.json, validation_summary.tsv
    └── upload/         <id>.bundle.upload.json
```

## Repository layout

```
main.nf  workflows/kpfhir.nf  subworkflows/local/  modules/local/   Nextflow code
bin/                  Python steps (Kleborate parsing, plasmids, QC, FHIR genomic / clinical / link, upload)
conf/                 resources, output layout, test profiles
assets/               samplesheet schema, FHIR CodeSystems, metadata templates
tools/cgmlst/         build the chewBBACA schema
tools/satusehat/      SATUSEHAT token helper
tools/dev/            PHI / credential guard (also a pre-commit hook)
tests/                pytest suite, synthetic fixtures, v1 equivalence goldens
docs/                 clinical metadata guide; docs/v1/ the v1.4.2 development notes
```

Developers can run `pytest`. `tests/compat/` checks that `--fhir_profile satusehat` reproduces
v1.4.2 output, apart from the documented allowlist in `tests/compat/allowlist.py`.

## Known limitations

These are kept from v1.4 on purpose; 2.0 changes the structure, not the science.
- **No species identification.** Every sample is assumed to be *K. pneumoniae*. Check
  Kleborate's species call for other *Klebsiella* or contaminants.
- **Kleborate is the only resistance caller.** The susceptibility panel is a genotypic prediction
  from a fixed rule set (`bin/kpfhir_resistance.py`), not a phenotype.
- **Hand-written tables**:
  - The clonal-complex table covers about 20 STs; other STs report `Unknown`.
  - The `high_risk_clone` list (ST11, ST258, ST147, ST15) and the subspecies ST ranges in
    `lineage.json` are incomplete. They are used only in the text report.
- **medaka model.** The default `r941_min_sup_g507` suits R9.4.1 flow cells. Set
  `--medaka_model` for R10 data.
- **Plasmid placement from short reads is a prediction.** MEGAHIT breaks plasmids at repeats,
  so many genes are `unresolved`. Nanopore assemblies place genes more reliably.
- **Point-mutation HGVS.** Kleborate's porin and truncation notations (e.g. `OmpK36GD`,
  `MgrB-62%`) are not HGVS. The protein-change component is filled only when the notation looks
  like a substitution, and v1.4's parsing is kept.
- **Every resistance Observation repeats the sample's virulence loci** as components (v1.4
  layout).
- **SATUSEHAT upload has not been verified** in 2.0 against the sandbox or production platform.

## Support

[GitHub Issues](https://github.com/oucru-id/kp-to-fhir-full/issues)
