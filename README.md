# KP Genomics to FHIR Pipeline (KPtoFHIR)

A platform-agnostic Nextflow pipeline for *Klebsiella pneumoniae* genomic analysis from raw sequencing data, producing HL7 FHIR R4 genomics bundles (IG v3.0.0). [Full documentation](https://kp-pipeline-docs.readthedocs.io/)

## Key Features

- **Multi-platform:** Illumina paired-end short reads and Oxford Nanopore (ONT) long reads.
- **Comprehensive Typing:** MLST, Capsule type (K locus), O-antigen (O locus), virulence factor scoring, and resistance gene detection using [Kleborate](https://github.com/klebgenomics/Kleborate).
- **cgMLST:** Core-genome MLST using [chewBBACA](https://github.com/B-UMMI/chewBBACA) with a [Ridom](https://www.cgmlst.org/ncs) schema.
- **Quality Control:** Per-sample FastQC reports aggregated into MultiQC.
- **FHIR Compliance:** HL7 FHIR R4 bundles with Susceptibility Panel, MLST, Capsule Type, O-Antigen, Virulence Score, cgMLST, and DiagnosticReport resources.
- **Clinical Integration:** Merges genomic results with patient, organization, and practitioner metadata.

## Installation

### Setup

```bash
git clone https://github.com/oucru-id/kp-to-fhir-full.git
cd kp-to-fhir
```

### Dependencies

```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash

# Verify
nextflow -v
```

Required tools:

| Tool | Purpose |
|------|---------|
| `kleborate` | MLST, capsule, virulence, resistance typing |
| `chewBBACA.py` | Core-genome MLST |
| `fastp` | Illumina adapter trimming |
| `megahit`/`SPAdes` | Illumina assembly |
| `pilon` | Illumina assembly polishing |
| `chopper` | Nanopore quality filtering |
| `flye` | Nanopore assembly |
| `medaka` | Nanopore assembly polishing |
| `fastqc` | Per-sample QC |
| `multiqc` | Aggregated QC report |
| Java + `fhir-validator.jar` | FHIR bundle validation |

## Directory Structure

```
kp-to-fhir
├── main.nf                             # Main workflow
├── nextflow.config                     # Configuration and parameters
├── workflows/
│   ├── illumina.nf                     # Illumina sub-workflow
│   ├── nanopore.nf                     # Nanopore sub-workflow
│   ├── typing.nf                       # Kleborate typing sub-workflow
│   ├── cgmlst.nf                       # cgMLST sub-workflow
│   ├── fhir.nf                         # FHIR bundle generation
│   ├── validate_fhir.nf                # FHIR validation
│   ├── merge_clinical_data.nf          # Clinical metadata merge
│   ├── upload_fhir.nf                  # FHIR server upload
│   ├── report.nf                       # QC and sample report generation
│   └── utils.nf                        # Utility functions
├── scripts/
│   ├── annotated_to_fhir.py            # Kleborate JSON-to-FHIR converter
│   ├── clinical_metadata_parser.py     # clinical metadata parser
│   ├── generate_sample_report.py       # Per-sample text report
│   ├── merge_clinical_fhir.py          # FHIR genomics + clinical data merger
│   ├── upload_fhir.py                  # FHIR uploader (OAuth 2.0)
│   ├── get_access_token.py             # Standalone token fetcher
│   └── get_versions.py                 # Software version collector
├── data/
│   ├── NGS/                            # Input FASTQ files
│   ├── cgmlst_schema/                  # chewBBACA cgMLST schema
│   ├── patient_clinical_metadata.csv   # Patient metadata
│   ├── organization_metadata.csv       # Organization metadata
│   └── practitioner_metadata.csv       # Practitioner metadata
└── tools/
    └── fhir-validator.jar              # HL7 FHIR validator
```

## Input Data

### Illumina Reads

Place paired-end FASTQ files in `data/NGS/`:

```
data/NGS/SAMPLE_1_illumina.fastq.gz
data/NGS/SAMPLE_2_illumina.fastq.gz
```

### Nanopore Reads

Place single-end FASTQ files in `data/NGS/`:

```
data/NGS/SAMPLE_ont.fastq.gz
```

### cgMLST Schema

Place Klebsiella cgMLST allele FASTA files in `data/cgmlst_schema/`:

```
data/cgmlst_schema/KP1_RS00005.fasta
data/cgmlst_schema/KP1_RS00010.fasta
...
```

## Usage

### Get Access Token (FHIR Upload)

```bash
python scripts/get_access_token.py
```

### Basic Run

```bash
nextflow run main.nf
```

### Run with FHIR Upload

> Get the access token first before running with upload.

```bash
nextflow run main.nf \
  --fhir_server_url "https://<BASE_URL>/fhir"
```

## Output Structure

```
results/
├── qc/
│   └── multiqc_report.html             # Aggregated QC report
├── assemblies/
│   └── *.fasta                         # Polished genome assemblies
├── typing/
│   └── *.typing.json                   # Kleborate typing results per sample
├── cgmlst/
│   └── results_alleles.tsv             # cgMLST allele calls
├── fhir/
│   └── *.fhir.json                     # FHIR genomics bundles
├── fhir_merged/
│   └── *.merged.fhir.json              # FHIR bundles with clinical data
├── fhir_validated/
│   └── *.validation.txt                # FHIR validation results
├── fhir_upload/
│   └── *.upload.json                   # FHIR upload results
├── reports/
│   └── *.summary_report.txt            # Per-sample summary reports
├── runningstat/
│   ├── execution.html                  # Nextflow execution report
│   ├── timeline.html                   # Timeline report
│   └── dag.html                        # Workflow DAG
└── software_versions.yml               # Software version manifest
```

## Support

[GitHub Issues](https://github.com/oucru-id/kp-to-fhir-full/issues)
