# KPtoFHIR pipeline v1.4.2 — technical report

**Date:** 31 October 2026
**Pipeline:** `~/.nextflow/assets/KP_Mut_Test` (Nextflow 26.04.1, local executor, WSL2 Ubuntu)
**Covers:** work done on 30–31 October 2026, plus the 2 November update (section 0)

## 0. Update — 2 November 2026

The high-priority issues H1, H2, H5 and the silent-fallback part of H3 (section 6) are fixed. H4 (FHIR validation) and the chemistry question in H3 are still open.

| Issue | Fix | Where |
|---|---|---|
| H2: `^` treated as chromosomal | `parse_results` sets `"intrinsic": true` only for genes in Kleborate's `Bla_chr` column, and leaves them out of `acquired_resistance`. The FHIR narrative says "(intrinsic chromosomal)" only for those genes, so the merge step now classifies `^` genes (aac(3)-IIa^, qnrB1^, aadA2^, aph3-Ia^, …) as acquired. The sample report uses the same flag. Its "clinically significant" check now also understands Kleborate's class codes (`agly`, `flq`, `tet`). | `typing.nf`, `annotated_to_fhir.py`, `generate_sample_report.py` |
| H3: silent medaka fallback | `medaka_polish` reports `POLISH_STATUS`, and the run log shows `WARN:` when polishing fails. | `nanopore.nf` |
| H3: medaka still failed in the re-run | **CUDA out of memory.** Medaka's default batch size of 100 needs about 5.8 GB of GPU memory; this machine has a 4 GB RTX 2050. Fix: new `medaka_batch_size = 20` (about 40 s per sample on the GPU), `maxForks = 1` for `medaka_polish` because samples share one GPU, and a CPU retry if the GPU run fails. | `nanopore.nf`, `nextflow.config` |
| H5: low-coverage samples | New `assembly_qc` step and `bin/assembly_qc.py`. They write `results/qc/assembly/<sample>_assembly_qc.json` (coverage, coverage status, polishing status, warnings) for each Nanopore sample. This record feeds: a "Sequencing quality control" Observation in FHIR (placeholder code `sequencing-qc`, CodeSystem `…/CodeSystem/kp-qc` in `kp_codes.py`, interpretation `A` when there are warnings); a "Quality warning: …; interpret these results with caution" sentence and a "Sequencing quality warning" `conclusionCode` in the DiagnosticReport, plus an HTML banner; and a banner and SEQUENCING QC section in the sample report. | `nanopore.nf`, `fhir.nf`, `report.nf`, `main.nf`, `annotated_to_fhir.py`, `merge_clinical_fhir.py`, `generate_sample_report.py` |
| Minor | Removed the double full stop in DiagnosticReport conclusions. | `merge_clinical_fhir.py` |
| H1: full re-run | Pre-fix outputs moved to `results_pre_fix_20261031/`. Full run 2 Nov 09:24–10:00 (`lethal_joliot`, medaka out of memory), then `-resume` 10:03–10:10 with the GPU fix. Exit code 0 both times; upload and validation skipped. | — |

**Results of the 2 November run** (panel order: AMP AMC TZP CTX CAZ FEP MEM CIP GEN AMK SXT COL):

| Sample | QC | Panel (old → new) | Notes |
|---|---|---|---|
| ERR10367345 | PASS: 39.6×, polished | `RSSSSSSSSSSS` (no change) | cgMLST call rate 94.1% → 94.8%. Conclusion: no acquired genes; intrinsic ampicillin resistance only. Capsule is still "null" (KL28). |
| SRR25759560 | Illumina, no QC record | `RSRRRRSRSSRS` (no change) | The conclusion now lists aadA2 and aph3-Ia, which were dropped before. |
| SRR25759561 | Illumina | `RSRRRRSRSRRS` → `RSRRRRSRRRRS` | Gentamicin R (aac(3)-IIa). |
| SRR25759568 | Illumina | `RRRRRRRRSRRS` → `RRRRRRRRRRRS` | Gentamicin R (aac(3)-IIa + rmtG). GyrA S83I + ParC S80I now reported. |
| SRR32156747 | **WARNING: 8.8×**, polished | `RSSSSSSSSSSS` → `RSSSSSSSSSSR` | cgMLST 8.3% → 19.0%. CTX-M-15, OXA-1, aac(6')-Ib-cr and others are still truncated. Colistin R comes from a pmrB frameshift that is probably an assembly error. |
| SRR32156754 | **WARNING: 13.6×**, polished | `RSSSSSSSSSRS` → `RSSSSSSSSSRR` | cgMLST 16.1% → 32.0%. ST405-4LV → ST405-3LV. Truncated genes and colistin R as for SRR32156747. |

Polishing can't make up for 9–14× coverage. The resistance profiles of SRR32156747 and SRR32156754 should not be used. The quality warning now appears in their FHIR conclusions and sample reports, but the pipeline still reports their genotype-based predictions.

**Medium-priority fixes M1–M3 (2 November):**

| Issue | Fix | Where |
|---|---|---|
| M1: virulence score | `typing.json` now stores Kleborate's own 0–5 score (`klebsiella_pneumo_complex__virulence_score__virulence_score`); it is `Unknown` if the column is missing. A score of 0 now produces an Observation. For these six samples the values don't change (1, 0, 1, 0, 1, 1), but SRR25759560 and SRR25759568 now have a virulence-score Observation. | `typing.nf`, `create_virulence_score_observation` |
| M2: Capsule null | When Kaptive reports "Capsule null", the capsule Observation now carries the locus, e.g. `KL28 (capsule null)`, code `KL28`. It is left out only when both locus and type are unknown. | `create_capsule_observation` |
| M2: unknown ST | The MLST Observation is kept and given only a text value: `ST not determined (inexact or missing alleles)`, or `Novel ST (new allele combination)` when all 7 alleles are exact. Every MLST Observation now carries 7 "MLST allele <locus>" components, which the merge step already reads. | `create_mlst_observation` |
| M3: clonal complex | The `CC` + ST-number fallback is removed. Only STs in a fixed table get a clonal complex; locus variants (`ST405-3LV`), novel STs and STs not in the table stay `Unknown`. The table now puts ST11, ST512, ST340 and ST437 in CC258, which is standard *K. pneumoniae* nomenclature (CG258; ST11 was `CC11` before). Effect here: ST1774, ST45, ST629 and ST405-3LV go from CC1774, CC45, CC629 and CC4053 to `Unknown`. The clonal complex only appears in `typing.json` and `lineage.json`; FHIR never carried it. | `typing.nf` |

Two new issues found while fixing M1–M3:

- **Resistance score was overwritten (fixed 2 November).** `parse_results` read Kleborate's resistance score, then replaced it with its own formula (carbapenemase 3, ESBL 2, more than 2 classes 1). It now keeps Kleborate's score, which is `Unknown` if the column is missing. Kleborate's scale is: 0 no ESBL or carbapenemase, 1 ESBL, 2 carbapenemase, 3 carbapenemase + colistin resistance. The sample report now prints the meaning next to the score, e.g. `2/3 (Kleborate: carbapenemase, no colistin resistance)`, and prints `Not determined` when the score is missing. Effect: SRR25759560 and SRR25759561 go from 2 to 1, SRR25759568 from 3 to 2, and the other samples stay at 0. The score is not in the FHIR bundles. (`typing.nf`, `generate_sample_report.py`)
- **Lineage data was never loaded into FHIR (input removed 2 November).** `load_lineage_data` in `annotated_to_fhir.py` looked for `*.lineage.json`, but the files are named `<sample>_lineage.json`, so `sample_lineage_info` was always empty. Little was lost, because `lineage.json` is built from the same Kleborate row as `typing.json`, and its `lineage_designation` is just the ST string. If the file had been found, each resistance Observation would have gained:
  - `KP lineage designation` = the ST. This duplicates the MLST Observation and the existing "KP MLST: …" text.
  - `High risk clone` = Yes only for ST11, ST258, ST147 and ST15. That list is too short: it would mark ST307 and ST101 as "No".
  - The narrative suffix " - Lineage: ST629", which also duplicates "KP MLST: ST629".

  The lineage input is now removed from the FHIR step:
  - `annotated_to_fhir.py`: removed `--lineage_dir`, `load_lineage_data` and the lineage branches in `add_standard_kp_components`, `create_enhanced_annotation` and the no-resistance narrative.
  - `fhir.nf`: `CREATE_FHIR` and the `FHIR` workflow no longer take lineage files.
  - `main.nf`: the call no longer passes `typing_results.lineage_report`.

  The rebuilt bundles are identical to the previous ones for all six samples, ignoring UUIDs and timestamps. `CREATE_FHIR` no longer receives every sample's lineage file, so changing one sample no longer reruns FHIR creation for all samples. `lineage.json` is still produced and still used by the sample report (risk level). If high-risk clones should be reported, add them later as one sample-level Observation with a full clone list.

New follow-ups from the re-run:

- Consider not reporting resistance predictions, or marking them preliminary, when coverage is LOW.
- Report Kleborate's `truncated_resistance_hits` (M5): they would have shown the CTX-M-15 and OXA-1 genes broken by assembly errors.
- DiagnosticReport text: point mutations appear under "Resistance genes detected" (e.g. "OMPK35 P.SER110FS"), and one conclusion still has a double space before "Predicted".

## 1. Summary

v1.4.2 adds plasmid analysis to the pipeline. Plasmid results are reported inside the existing genomic FHIR DiagnosticReport and in the sample text report.

A six-sample test run (run `trusting_ramanujan`, 31 Oct 13:42–14:16) found these problems:

- **Medaka never polished any Nanopore assembly.** The pipeline carried on with unpolished Flye assemblies and reported nothing.
- **The FHIR susceptibility panel always reported gentamicin as susceptible.** Aminoglycoside genes never reached the gentamicin rules.
- **Resistance point mutations were never applied to the panel.** Kleborate 3's fluoroquinolone and porin mutations were not read at all.

All three are fixed, and each fix was tested on its own (section 7). The pipeline has **not yet been re-run** with the fixes. Everything currently in `results/` comes from the run before the fixes and must not be used for reporting or upload.

Other changes:

- NanoPlot replaces FastQC for Nanopore reads and warns about low coverage.
- Scripts moved from `scripts/` to `bin/`, and the cgMLST schema moved to `db/`.
- Large or unneeded files are no longer copied to `results/`.
- FHIR validation can now be skipped with `--skip_validation`.

## 2. Changes

### 2.1 Plasmid analysis (new)

The design is in `doc/plan_plasmid_20261030.md`.

| Component | File | What it does |
|---|---|---|
| `PLASMID_RECON` | `workflows/plasmid.nf` | Runs `mob_recon` (MOB-suite 3.1.9) on each assembly. Also aligns Kleborate's own resistance gene reference (`CARD_v3.2.9.fasta`) to the assembly with `minimap2 --end-bonus=10 --eqx -c -x map-ont`, the same settings as Kleborate's resistance module. If MOB-suite fails, the step writes empty placeholder files and `status=failed` in `recon_status.txt` instead of stopping the run. |
| `PARSE_PLASMIDS` | `bin/parse_plasmids.py` | Joins the MOB-suite contig report, the MOB-typer results and the gene alignments with `typing.json`, and writes `<sample>_plasmid.json`. |
| `PLASMID_SHARING` | `bin/plasmid_sharing.py` | Writes `results/plasmid/plasmid_sharing.tsv`: MOB-suite cluster, replicon types, mobility, number of samples, sample IDs and resistance genes. |
| Code constants | `bin/kp_codes.py` | Placeholder FHIR CodeSystem `https://kp-pipeline-docs.readthedocs.io/fhir/CodeSystem/kp-plasmid` with codes `plasmid-analysis`, `plasmid` and `genomic-location`. Sequence Ontology codes for the location: `SO:0000155` (plasmid) and `SO:0000340` (chromosome). |

**How each resistance gene is located** (`parse_plasmids.py`):

- A gene alignment counts if it has at least 90% identity and 80% coverage (`amr_locate_min_identity`, `amr_locate_min_coverage`). Overlapping alignments are reduced to the best one, as Kleborate does.
- Each Kleborate gene is matched first to its exact allele, then to another allele in the same CARD cluster.
- A gene is `plasmid` if any copy sits on a plasmid contig.
- It is `chromosome` only if the contig is at least 20 kb long (`plasmid_chrom_min_contig`). Otherwise it is `unresolved`.
- Point mutations are always `chromosome`.
- Results from Illumina assemblies are labelled `predicted (short-read)`; Nanopore results are labelled `long-read`.

**FHIR output** (`bin/annotated_to_fhir.py`):

- Each resistance-gene Observation gets a "Genomic location" component, coded with Sequence Ontology or given as text "Unresolved". Plasmid-borne genes also get a "Plasmid" component, e.g. `AA276 IncFIB, IncFII`.
- A plasmid-analysis summary Observation is always present. Its value is the plasmid count, or a `dataAbsentReason` if the analysis failed.
- There is one Observation per plasmid, with components for: cluster, size, contigs, relaxase, mating-pair formation type, oriT, mobility, nearest reference plasmid, host range, carried genes and confidence.
- `hasMember` links each plasmid Observation to the gene Observations for the genes it carries.

**DiagnosticReport** (`bin/merge_clinical_fhir.py`):

- Plasmid Observations are left out of the resistance classification.
- The conclusion gets an extra sentence when a carbapenemase or ESBL gene is on a plasmid, and a text `conclusionCode` is added ("Plasmid-borne carbapenemase" or "Plasmid-borne ESBL").
- The HTML report gets a "Plasmid Content" table, and each gene in the gene list is tagged with its location.

**Sample text report** (`bin/generate_sample_report.py`): a new PLASMID ANALYSIS section.

**Switch:** `--plasmid_analysis false` turns the feature off; the FHIR bundles then contain no plasmid content.

### 2.2 Medaka polishing fix

**Cause 1: `samtools`.**

- `~/medaka-env/bin/samtools` was a link to Ubuntu's `/usr/bin/samtools` 1.19.2.
- That build's `--version` text contains bytes that aren't valid UTF-8 (`\xabBUILDPATH\xbb`).
- `medaka_version_report` decodes that text, fails, and reports samtools as "Not found". `medaka_consensus` then exits within seconds.
- The `set +e` fallback in `medaka_polish` copied the unpolished Flye assembly into place, so all three Nanopore samples were affected and nothing was reported.

**Fix 1:**

- bioconda `samtools` 1.21 is installed in `~/miniforge3/envs/samtools`.
- `~/medaka-env/bin/samtools` now points to it, and `medaka_version_report` passes.

**Cause 2: model.** The model was hard-coded as `r941_e81_sup_variant_g514`. That is a variant-calling model; `medaka_consensus` needs a consensus model.

**Fix 2:** a new setting, `params.medaka_model`, with default `r941_min_sup_g507`.

- ENA records only "MinION" for all three runs and gives no chemistry.
- Median read quality (Q13.7 for ERR10367345, Q10.9 for SRR32156747 and SRR32156754) fits R9.4.1, so the default is a best guess.
- Override it per run with `--medaka_model`. These consensus models are installed locally: `r941_min_{sup,hac,fast}_g507` and `r941_e81_sup_g514`.
- The model is recorded in `software_versions.yml`.

### 2.3 Susceptibility panel and point mutations

| Problem | Where | Fix |
|---|---|---|
| Kleborate writes drug classes as short codes (`agly`, `flq`, `sul`, `tmt`, `col`). The panel looked for the word `aminoglycoside`, so no aminoglycoside gene reached the gentamicin or amikacin rules. | `create_susceptibility_panel` | Added `DRUG_CLASS_ALIASES`. |
| Gene names come in two styles (`aac(3)-IIa.v1^` and `aph3-Ia.v1^`). The old first-match prefix rule depended on the order of the rules. | same | `normalise_aminoglycoside_name()` removes `( ) ' ^ * ?`, then the longest matching rule wins. This also stops aph(3')-Ia being called gentamicin-resistant by the catch-all rule. |
| `parse_results` looked for `GyrA_mutations`, `ParC_mutations`, `OmpK35_mutations` and `OmpK36_mutations`. Kleborate 3.2.4 doesn't output those columns, so fluoroquinolone and porin mutations were lost. | `workflows/typing.nf` | It now reads `Flq_mutations`, `Omp_mutations` and `Col_mutations`, and splits each entry into protein and change, e.g. `{"gene": "GyrA", "mutation": "p.Ser83Ile"}`. |
| Colistin mutations were stored as `"gene": "Colistin"`, but the panel checked for `PmrB`, `MgrB` and similar names. | `annotated_to_fhir.py` | Added a shared `COLISTIN_MUTATION_GENES` list, used by both the panel and the Observation builder. It also accepts the old `Colistin` label. |

Effect when re-tested on the current `typing.json` files (tested in isolation, not by a full pipeline run):

| Sample | Change | Reason |
|---|---|---|
| SRR25759561 | Gentamicin S → R | aac(3)-IIa |
| SRR25759568 | Gentamicin S → R | aac(3)-IIa + rmtG. GyrA S83I and ParC S80I are now reported as well. |
| SRR32156747, SRR32156754 | Colistin S → R | pmrB frameshifts, probably Nanopore assembly errors (see section 5). Recheck after the re-run. |

### 2.4 Nanopore read QC: NanoPlot and coverage warning

- A new `nanoplot` step replaces the Nanopore `fastqc` step. It runs on the reads after `chopper`, which are the reads Flye and medaka use. Illumina still uses FastQC.
- Coverage = total bases ÷ `genome_size_bp` (default 5,500,000). Below `min_nanopore_coverage` (default 30×), the run log shows a `WARN:` line. The sample is still processed.
- `<sample>_coverage_mqc.tsv` adds a "Nanopore coverage" table to the MultiQC report, marking each sample PASS, LOW or UNKNOWN.
- Outputs go to `results/qc/nanoplot/`: `<sample>_NanoPlot-report.html`, `<sample>_NanoStats.txt` and `<sample>_coverage_mqc.tsv`. MultiQC 1.30 reads `NanoStats.txt` through its `nanostat` module.
- NanoPlot 1.48.0 is installed in `~/miniforge3/envs/nanoplot` (897 MB), kept separate so it can't change medaka-env's numpy or pandas.

### 2.5 Layout and output changes

| Change | Detail |
|---|---|
| `scripts/` → `bin/` | All Python scripts are called by name (Nextflow adds `bin/` to each step's PATH). Windows line endings were converted to Unix, a shebang was added to `get_versions.py`, and the scripts were made executable (`kp_codes.py` stays non-executable because it is only imported). |
| cgMLST schema → `db/` | `db/cgmlst_schema/` (2,361 files, 3.4 GB) plus `cgmlst_schema_invalid_alleles.txt` and `cgmlst_schema_summary_stats.tsv`. `params.cgmlst_schema` was updated. |
| `publishDir` fix | `PLASMID_RECON` builds a per-sample output folder from `sample_id`. Nextflow 26.04 rejected the plain string ("No such variable: sample_id"), so the path is now a closure: `path: { ... }`. |
| Trimmed reads not published | `chopper` no longer copies `*_trimmed.fastq.gz` (about 200 MB per sample) to `results/trimming/`. |
| `versions.yml` removed | Removed from `CREATE_FHIR`, `MERGE_CLINICAL_FHIR` and `kleborate`; nothing used them. Tool versions are recorded in `results/software_versions.yml`. |
| `--skip_validation` | When set, `VALIDATE` is skipped and the merged bundles are passed straight to upload. |
| MultiQC input folder | `MULTIQC` now reads from `qc_results/` (was `fastqc_results/`) and also accepts `_NanoStats.txt` and `_mqc.tsv` files. |

## 3. Configuration reference (new or changed settings)

| Setting | Default | Purpose |
|---|---|---|
| `medaka_model` | `r941_min_sup_g507` | Medaka consensus model; must match the flow cell and basecaller. |
| `medaka_batch_size` | `20` | Medaka batch size; the default of 100 needs more than 5 GB of GPU memory (added 2 Nov). |
| `nanoplot_env` | `null` (= `$HOME/miniforge3/envs/nanoplot`) | NanoPlot environment. |
| `genome_size_bp` | `5500000` | Genome size used for the coverage estimate. Flye still uses its own hard-coded `5.5m`. |
| `min_nanopore_coverage` | `30` | Coverage warning threshold. |
| `plasmid_analysis` | `true` | Turns plasmid analysis on or off. |
| `mobsuite_env` | `null` (= `$HOME/miniforge3/envs/mobsuite`) | MOB-suite environment. |
| `mobsuite_db` | `null` (= database installed by `mob_init`) | MOB-suite database. |
| `plasmid_chrom_min_contig` | `20000` | Minimum contig length for a `chromosome` call. |
| `amr_locate_min_identity` / `amr_locate_min_coverage` | `90` / `80` | Thresholds for gene alignments. |
| `skip_upload` | `true` | Upload is off by default. |
| `skip_validation` | `true` | FHIR validation is off by default. |
| `cgmlst_schema` | `$baseDir/db/cgmlst_schema` | Moved from `data/`. |
| `version` | `1.4.2` | Pipeline version. |

New resource settings: `PLASMID_RECON` (4 CPUs, 8 GB) and `nanoplot` (4 CPUs, 4 GB).

## 4. External environments

These live outside the pipeline folder and must exist on any machine that runs it.

| Environment | Path | Used by | Notes |
|---|---|---|---|
| medaka-env (venv) | `~/medaka-env` | most steps; activate before running | `bin/samtools` → `~/miniforge3/envs/samtools/bin/samtools` |
| samtools | `~/miniforge3/envs/samtools` | medaka, `pilon_polish` | samtools 1.21 (bioconda) |
| mobsuite | `~/miniforge3/envs/mobsuite` | `PLASMID_RECON` | MOB-suite 3.1.9, Python 3.11. Database from Zenodo record 10304948, about 2.9 GB installed. |
| nanoplot | `~/miniforge3/envs/nanoplot` | `nanoplot` | NanoPlot 1.48.0 |

To undo the samtools change: `ln -sfn /usr/bin/samtools ~/medaka-env/bin/samtools`. That would break medaka again.

## 5. Findings from the six-sample test run (before the fixes)

All six isolates are *K. pneumoniae* (Kleborate species check run by hand, strong match).

| Sample | Data | ST | Key resistance (location) | Reliability |
|---|---|---|---|---|
| SRR25759568 | Illumina, 640 contigs | ST629, K10 | KPC-2 + CTX-M-12 together on a transferable IncL/M plasmid (73 kb); rmtG; GyrA S83I + ParC S80I | Gentamicin and amikacin should be R (fixed) |
| SRR25759561 | Illumina, 87 contigs | ST45, K62 | CTX-M-15, OXA-1, aac(6')-Ib-cr, qnrB1, aac(3)-IIa, tet(A), dfrA14, sul2 on one transferable 187 kb IncFIB/IncFII plasmid | Gentamicin should be R (fixed) |
| SRR25759560 | Illumina, 714 contigs | ST1774, KL107 | CTX-M-15 on a small untyped plasmid | Plasmid placement is a short-read prediction |
| ERR10367345 | Nanopore, 1 contig, 44.6× (raw) | ST20, KL28 | No acquired genes; chromosomal SHV-187 | Unpolished; 94.1% of cgMLST genes called |
| SRR32156747 | Nanopore, 7 contigs, 11.8× (raw) | not typable | Reported as none, but Kleborate lists CTX-M-15, OXA-1, aac(6')-Ib-cr, qnrB1, dfrA14, aac(3)-IIa as **truncated** | Unpolished, low coverage, median read Q10.9; 8.3% of cgMLST genes called |
| SRR32156754 | Nanopore, 4 contigs, 17.9× (raw) | ST405-4LV | Reported as dfrA14 only; truncated CTX-M-15, OXA-1, aac(6')-Ib-cr, qnrB1, sul2, tet(A), strAB | Unpolished, low coverage; 16.1% of cgMLST genes called |

What these results show:

- **No clonal transmission.** Six different STs. The four genomes with usable cgMLST differ at more than 1,600 genes.
- **A shared resistance gene set.** SRR25759561 (ST45) and SRR25759568 (ST629) carry nearly the same OXA-1 / aac(6')-Ib-cr / qnrB1 / tet(A) / dfrA14 gene set. The two SRR32156xxx isolates probably carry it too. Long-read data would show whether it is the same plasmid.
- **Plasmid cluster AA277 is a broad family.** `plasmid_sharing.tsv` lists it in three samples, but the gene content differs in each, so the shared cluster is not evidence of transmission. Only SRR32156747 and SRR32156754 share the closer sub-cluster AI089.
- **No hypervirulence markers.** No aerobactin, salmochelin, colibactin or rmpA/rmpA2 in any sample; three carry yersiniabactin only.
- **The SRR32156xxx results are mostly assembly errors.** Their "Capsule null" calls, pmrB and OmpK35 frameshifts, inexact MLST alleles and truncated genes are expected from unpolished, low-coverage Nanopore assemblies. Polishing will improve them, but at 11–18× coverage it may not fully fix them.

## 6. Remaining issues

### High priority

Status on 2 November: H1, H2 and H5 are fixed, and the silent fallback in H3 is fixed (section 0). H4 and the H3 chemistry question are still open.

| # | Issue | Where | Impact | Suggested action |
|---|---|---|---|---|
| H1 | The full pipeline hasn't been re-run with the fixes, and `work/` was deleted. | — | `results/` holds pre-fix outputs, including wrong gentamicin calls and unpolished Nanopore assemblies. | Run the full pipeline from scratch, then recheck sections 2.3 and 5. |
| H2 | A `^` marker is treated as "chromosomal". In Kleborate, `^` only means the protein matches exactly but the DNA doesn't. | `create_enhanced_annotation` (`annotated_to_fhir.py`), `classify_resistance_profile` (`merge_clinical_fhir.py`), `is_intrinsic_gene` (`generate_sample_report.py`) | Genes such as aac(3)-IIa.v1^, qnrB1.v2^, aadA2^ and aph3-Ia.v1^ are dropped from the DiagnosticReport resistance classification. Narratives can read "(chromosomal) … located on plasmid". | Decide intrinsic status from the actual location (`parse_plasmids.py`) or a list of intrinsic genes, not from `^`. |
| H3 | The medaka model's chemistry is unconfirmed, and a medaka failure still silently falls back to the unpolished assembly. | `medaka_polish` | Assemblies may be unpolished, or polished with the wrong model, without any record in the outputs. | Confirm chemistry with the data providers. Add a per-sample model option if runs differ. Write a polishing status into the outputs and the run log. |
| H4 | The new FHIR resources have not been checked with the FHIR validator. | `VALIDATE` | Possible profile errors in the plasmid Observations, location components and placeholder codes. | Run with `--skip_validation false` and compare error counts with earlier bundles. |
| H5 | Low-coverage Nanopore samples (SRR32156747 about 11.8×, SRR32156754 about 17.9×). | data | Results stay unreliable even after polishing. | Resequence, or flag or exclude them; the new coverage warning reports them. |

### Medium priority

Status on 2 November: M1, M2 and M3 are fixed (section 0).

| # | Issue | Where |
|---|---|---|
| M1 | `virulence_score` in `typing.json` is a count of virulence loci, not Kleborate's 0–5 score. The code reads Kleborate's score, then overwrites it. A score of 0 produces no Observation. | `parse_results` (`typing.nf`), `create_virulence_score_observation` |
| M2 | A "Capsule null" result drops the capsule Observation, so the KL locus is lost; an unknown ST drops the MLST Observation. | `create_capsule_observation`, `create_mlst_observation` |
| M3 | The clonal complex is just "CC" plus the ST number, and `ST405-4LV` becomes `CC4054`. | `parse_results` (`typing.nf`) |
| M4 | The pipeline never checks species; FHIR assumes *K. pneumoniae*. | Add `enterobacterales__species` to the Kleborate modules. |
| M5 | Kleborate columns not read: `Bla_inhR_acquired`, `klebsiella_pneumo_complex__amr___acquired` (TEM-1D appeared here), `truncated_resistance_hits` (would flag broken genes) and `SHV_mutations`. | `parse_results` (`typing.nf`) |
| M6 | The plasmid codes are placeholders; the Sequence Ontology system URI still needs confirming. | `bin/kp_codes.py` |
| M7 | Some panel rules overcall resistance: ciprofloxacin R from a qnr gene or aac(6')-Ib-cr alone, and piperacillin-tazobactam R from any ESBL. | `create_susceptibility_panel` |

### Low priority

| # | Issue | Where |
|---|---|---|
| L1 | No amino-acid-change (pHGVS) component for Kleborate 3 notation such as `p.Ser83Ile`; the regex expects one-letter codes. | Observation builder (`annotated_to_fhir.py`) |
| L2 | Every resistance-gene Observation also carries a yersiniabactin component. | `create_enhanced_annotation` |
| L3 | Flye's genome size (`5.5m`) is hard-coded separately from `genome_size_bp`. | `flye` (`nanopore.nf`) |
| L4 | Machine-specific paths (`$HOME/miniforge3/...`, `/opt/pilon/pilon-1.24.jar`). No containers. The folder isn't a git repository. | various |
| L5 | `PLASMID_RECON` still writes a `versions.yml` in `work/` (it isn't published). | `plasmid.nf` |
| L6 | 13 lint warnings in 6 files, all from before this work. | `nextflow lint` |
| L7 | Editing `bin/` scripts from Windows removes the executable flag and can add Windows line endings. | After editing, run `chmod +x bin/*.py` (except `kp_codes.py`) and check line endings. |

## 7. Verification performed

| Check | Result |
|---|---|
| `nextflow lint main.nf workflows/` | No errors; 13 warnings, all from before this work. |
| `nextflow run main.nf -preview` | Pass. |
| `publishDir` closure fix | Test pipeline published per-sample folders. |
| `medaka_version_report` after the samtools change | Pass (samtools 1.21). |
| `parse_results` on six Kleborate outputs | Only `point_mutations` changed; everything else in `typing.json` is identical. |
| Susceptibility panel, six samples plus 14 aminoglycoside name cases | Results as in section 2.3. |
| `annotated_to_fhir.py` end-to-end on SRR25759568 and SRR32156747 | Runs. GyrA/ParC → fluoroquinolone, OmpK35 → carbapenem, PmrB → colistin. |
| `nanoplot` and `MULTIQC` in a test pipeline (raw reads of two samples) | Low-coverage warning raised for SRR32156747 (11.8×). MultiQC shows the "Nanopore coverage" and NanoStat sections. |
| Not done | Full pipeline run, FHIR validator, upload. |

## 8. How to run

```bash
source ~/medaka-env/bin/activate
cd ~/.nextflow/assets/KP_Mut_Test

nextflow run main.nf                                    # validation and upload are off by default
nextflow run main.nf --skip_validation false            # also validate the FHIR bundles
nextflow run main.nf --medaka_model r941_min_hac_g507   # different medaka consensus model
nextflow run main.nf --plasmid_analysis false           # leave out plasmid analysis
```

`nextflow run main.nf -preview` starts a new, empty run session. A plain `-resume` afterwards resumes that empty session and reruns everything. Either resume a named run (`nextflow log` lists them; then `-resume <run name>`), or don't run `-preview` just before `-resume`.

Nextflow only adds or overwrites files in `results/`; it never deletes old ones. To get a `results/` folder with only the current outputs, move or delete it before running.

## Appendix: files changed

| File | Change |
|---|---|
| `main.nf` | Plasmid sub-workflow wiring; `skip_validation`; MultiQC input filter |
| `nextflow.config` | New settings (section 3); schema path; resource settings |
| `workflows/plasmid.nf` | New |
| `workflows/nanopore.nf` | `fastqc` replaced by `nanoplot` + coverage warning; `medaka_model`; trimmed reads not published |
| `workflows/typing.nf` | Kleborate 3 mutation columns; `versions.yml` removed |
| `workflows/fhir.nf` | Plasmid JSON input; `versions.yml` removed; scripts called by name |
| `workflows/merge_clinical_data.nf` | `versions.yml` removed; script called by name |
| `workflows/report.nf` | Plasmid input to sample report; MultiQC input folder |
| `workflows/utils.nf` | Plasmid, NanoPlot and medaka model in `software_versions.yml` |
| `workflows/upload_fhir.nf` | Script called by name |
| `bin/kp_codes.py`, `bin/parse_plasmids.py`, `bin/plasmid_sharing.py` | New |
| `bin/assembly_qc.py` | New (2 Nov): per-sample assembly QC record |
| `workflows/fhir.nf`, `workflows/report.nf`, `main.nf` | 2 Nov: assembly QC input (Nanopore samples only, optional) |
| `bin/annotated_to_fhir.py` | Plasmid Observations and location components; panel fixes; shared colistin gene list |
| `bin/merge_clinical_fhir.py` | Plasmid skip list; conclusion, `conclusionCode` and HTML plasmid table |
| `bin/generate_sample_report.py` | PLASMID ANALYSIS section |
| `bin/get_versions.py` | MOB-suite and NanoPlot versions |
| `doc/plan_plasmid_20261030.md` | Paths updated to `bin/` |
