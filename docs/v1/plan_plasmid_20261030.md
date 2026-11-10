# Plan: add plasmid analysis to the KP pipeline (FHIR + text report)

**Status, 2026-10-30:** approved, not yet implemented. The MOB-suite install was started and then stopped, so nothing was installed. No pipeline files have changed.

**Updated after approval**, based on what the install attempt showed:
- MOB-suite goes in its own conda environment, not `~/medaka-env`. A pip dry run showed that installing it there would downgrade numpy (1.24.4 → 1.23.4) and pandas (2.0.3 → 1.5.3), which medaka and Kleborate depend on.
- The database download is about 0.55 GB, not "a few GB".
- Free space is limited by the Windows C: drive (20.5 GB free on 2026-10-30), not the 756 GB WSL reports, because WSL's disk is a single file on C:.
- Kleborate's exact minimap2 settings and gene-name markers are now filled in (Changes, steps 2 and 3).
- An existing bug found along the way is listed under "Known issue".

## Context

The pipeline (`KP_Mut_Test`, v1.4.1) reports genomic results only. Kleborate finds resistance genes, but nothing records whether a gene sits on a plasmid or the chromosome. The only location hint is the "(chromosomal)" label at `annotated_to_fhir.py:462`, and that is based on a misreading of Kleborate's `^` marker (see Known issue). The FHIR bundles in `results/fhir_validated/` never mention plasmids.

**Is it feasible?** Yes. Almost everything needed is already installed:
- `/usr/bin` already has BLAST+ 2.12, Mash 2.3 and minimap2 2.26.
- Kleborate 3.2.4 ships the CARD reference alleles it uses (`kleborate/modules/klebsiella_pneumo_complex__amr/data/CARD_v3.2.9.fasta`). The FASTA headers have the form `>cluster__geneclass__ALLELE__seqID`, and the allele field matches Kleborate's output names exactly (`CTX-M-15`, `aph3-Ia.v1`).

The only new pieces are **MOB-suite** (`mob_recon`, from bioconda) and its database.

**Does it add real value?** Yes:
- It shows where each resistance gene is (plasmid, chromosome or unresolved). Carbapenemase or ESBL genes on a conjugative plasmid are the main way resistance spreads between strains and species, and that matters for infection control.
- It reports each plasmid's replicon type (IncFII(K), IncX3, IncL…), relaxase/MPF type and predicted mobility (conjugative, mobilizable or non-mobilizable).
- MOB-suite cluster IDs let you compare plasmids across samples. This catches plasmid-driven outbreaks that cgMLST misses, because cgMLST only tracks clonal spread.

**Limitation to state in every report:** Illumina samples are assembled with MEGAHIT, which breaks plasmids at repeated IS elements, so plasmid assignment there is a *prediction*. Genes on short contigs come out as "unresolved". Nanopore (Flye) assemblies give much more reliable placement. Each report labels its confidence by assembly type.

**Decisions made:**
- Plasmid results go **inside the existing genomic DiagnosticReport**, not a separate report.
- New codes are **placeholders** held in one shared module until kemkes SP codes are issued.
- Include the **cross-sample plasmid-sharing table**.
- Flye circularity parsing is out of scope.

## Data flow

```
assemblies ──┬─ kleborate ─ parse_results ─ typing.json ──────────────┐
(illumina|ont)│                                                        ▼
             └─ PLASMID_RECON (mob_recon + minimap2 vs Kleborate CARD) ─ PARSE_PLASMIDS ─ <sample>_plasmid.json
                                                                        │
              ┌──────────────────────┬─────────────────────────┬──────┴──────────────┐
              ▼                      ▼                         ▼                     ▼
        CREATE_FHIR           GENERATE_SAMPLE_REPORT     PLASMID_SHARING      (merge → validate → upload,
   (+location, plasmid obs)    (+PLASMID section)      plasmid_sharing.tsv     unchanged wiring)
```

## Step 0: install MOB-suite in its own conda environment

Don't install it into `~/medaka-env`. A pip dry run there would downgrade numpy 1.24.4 → 1.23.4 and pandas 2.0.3 → 1.5.3.

Check the free space on C: first. The download is about 0.55 GB in total:
- the MOB-suite database `data.tar.gz`: 472,851,578 bytes (v3.1.8, Zenodo record 10304948);
- the NCBI taxonomy `taxdump.tar.gz`: 79,954,233 bytes (size on 2026-10-30).

Both take more space once unpacked and indexed.

```bash
~/miniforge3/bin/mamba env remove -n mobsuite    # removes the empty stub left by the stopped install
~/miniforge3/bin/mamba create -y -n mobsuite -c conda-forge -c bioconda --strict-channel-priority mob_suite=3.1.9
~/miniforge3/bin/conda run -n mobsuite mob_init  # downloads the database
~/miniforge3/bin/conda run -n mobsuite mob_recon --version
```

## Changes

### 1. `nextflow.config`
- New params:
  - `plasmid_analysis = true`
  - `mobsuite_env = "${System.getenv('HOME')}/miniforge3/envs/mobsuite"`
  - `mobsuite_db = null` (null means the `mob_init` default)
  - `plasmid_chrom_min_contig = 20000`
  - `amr_locate_min_identity = 90`, `amr_locate_min_coverage = 80` (Kleborate's AMR defaults)
  - `skip_upload = false`
- `withName: PLASMID_RECON { cpus = 4; memory = '8 GB' }`.

### 2. `workflows/plasmid.nf` (new)
- **`PLASMID_RECON`**: input `tuple val(sample_id), path(assembly), val(asm_type)`. It does two things:
  - Runs `mob_recon` from the mobsuite environment: `PATH="${params.mobsuite_env}/bin:$PATH" mob_recon --infile <asm> --outdir mob_out -n <cpus> [-d <db>]`. Only this one command gets the environment's PATH, so `python3` still resolves to `medaka-env`, where Kleborate is installed. The call is wrapped in `set +e`, using the same failure pattern as `medaka_polish` in `workflows/nanopore.nf`. A failure must not kill the run.
  - Aligns Kleborate's CARD alleles to the assembly with the same settings Kleborate uses (`kleborate/shared/alignment.py:131`): `minimap2 --end-bonus=10 --eqx -c -x map-ont <assembly> <CARD_v*.fasta> > amr_hits.paf`. The FASTA path is resolved at runtime through `python3 -c "import kleborate…"`, so Kleborate upgrades are followed automatically.
  - Emits `contig_report.txt`, `mobtyper_results.txt` and `amr_hits.paf`. If any is missing it emits an empty placeholder plus a status file. Results are published to `results/plasmid/<sample>/`.
- **`PARSE_PLASMIDS`**: joins the outputs above with `typing.json` on `sample_id` and runs `bin/parse_plasmids.py`. It always emits `<sample>_plasmid.json`, even on failure (`status: "failed"`), because Nextflow `.join` silently drops samples that have no match.
- **`PLASMID_SHARING`**: collects every `*_plasmid.json` and writes `results/plasmid/plasmid_sharing.tsv` with these columns: cluster_id, replicon_types, mobility, n_samples, samples, amr_genes.
- **`workflow PLASMID`**: takes `(sample_id, assembly, asm_type)` and emits `plasmid_report` and `versions`.

### 3. `bin/parse_plasmids.py` (new). It does all the matching, so every consumer sees the same answer.
- It parses MOB-suite TSVs **by column header**, not by position.
- **Joining contigs.** It matches the PAF target name against `contig_report.contig_id`. If an ID doesn't match, it falls back to the `md5` of the contig sequence. MEGAHIT and Pilon names such as `k141_12_pilon` depend on taking only the first token of the header.
- **Kleborate's gene-name markers.** Kleborate appends these in a fixed order (`kleborate/shared/resMinimap.py:170-180`). Strip all of them to get the allele name:
  1. `^`: the protein matches exactly but the DNA doesn't. Kleborate then renames the hit to the allele with the exact protein match.
  2. Otherwise `*` (identity below 100%) and `?` (partial coverage).
  3. A truncation suffix such as `-56%`.
  4. ` +<mutation>` for class-changing SHV mutations.
- **Locating each gene.** For each `detailed_resistance[].gene` in `typing.json`:
  1. Strip the markers above to get the allele name.
  2. Look for PAF hits on that allele that pass the identity and coverage thresholds. Identity is PAF column 10 / column 11; coverage is the aligned part of the allele divided by its length.
  3. If none pass, take the best hit in the same CARD cluster at a locus that no other gene has already claimed. This is the usual path for `^` and `*` genes.
  4. Map each hit's contig to its molecule type and cluster ID.
- **Location rules:**
  - Any hit on a plasmid contig makes the gene `plasmid`.
  - `chromosome` is used only when the contig is at least `plasmid_chrom_min_contig` long.
  - Everything else is `unresolved`.
  - All copies are kept.
  - Point mutations (GyrA, ParC, OmpK35/36, MgrB/PmrB) are always `chromosome`.
- **Output schema.** The `gene_locations` keys are the **raw Kleborate gene strings**, exactly as they appear in `typing.json`, so no re-normalising happens downstream.
  ```json
  {"sample_id","status":"ok|no_plasmids|failed","assembly_type":"illumina|ont",
   "confidence":"predicted (short-read)|long-read","tools":{"mob_suite","card_db"},
   "plasmids":[{"cluster_id","secondary_cluster_id","replicon_types":[],"relaxase_types":[],
     "mpf_type","orit_types":[],"predicted_mobility","size_bp","gc","num_contigs","contigs":[],
     "mash_nearest_neighbor","mash_distance","predicted_host_range","amr_genes":[]}],
   "gene_locations":{"<kleborate gene>":{"location":"plasmid|chromosome|unresolved",
     "copies":[{"molecule","cluster_id","contig","start","end","identity","coverage"}]}}}
  ```

### 4. `bin/kp_codes.py` (new, shared constants)
- `PLASMID_CODE_SYSTEM` is a placeholder URI. It has codes `plasmid-analysis`, `plasmid` and `genomic-location`. Swap in kemkes SP codes here once they're issued.
- Sequence Ontology codes for the location field: `SO:0000155` plasmid, `SO:0000340` chromosome. Confirm the canonical SO system URI against the validator output.
- Both FHIR scripts import this module. That works because Python puts the running script's own folder (`bin/`) on `sys.path`.

### 5. `main.nf`
- Build a tagged channel:
  ```
  illumina_out.assembly.map{ id,f -> tuple(id,f,'illumina') }
      .mix(nanopore_out.assembly.map{ id,f -> tuple(id,f,'ont') })
  ```
- When `params.plasmid_analysis` is on, call `PLASMID(tagged, typing_results.typing_report)`.
- When it's off, map typing reports to `tuple(id, [])`. This is the optional-path idiom; the bundle then stays identical to today's.
- Pass the plasmid channel to `FHIR(...)` and `GENERATE_REPORT(...)`.
- Gate `UPLOAD_FHIR` on `!params.skip_upload`.

### 6. `workflows/fhir.nf` and `bin/annotated_to_fhir.py`
- `CREATE_FHIR` input changes to `tuple val(sample_id), path(typing_json), path(plasmid_json)`, built with `typing_ch.join(plasmid_ch)`. Add `--plasmid_json` only when the path is non-empty.
- **Variant Observations** (the loop around lines 827–866). Add a `Genomic location` component (SO coding, or text `Unresolved`). For plasmid genes, also add a `Plasmid` component (`AA123 IncFII(K)`, valueString). Append "- located on plasmid AA123 (IncFII(K), conjugative)" to the narrative.
  **Never write "chromosomal" or "intrinsic" in this text.** `merge_clinical_fhir.py:710` matches those words and would drop chromosomally integrated acquired genes from the resistance classification. Keep the existing `^` logic in `create_enhanced_annotation` unchanged.
- **New `create_plasmid_observations()`**, written in the same style as `create_mlst_observation()`:
  - **One Observation per plasmid.** Its value is the replicon type(s). Components: cluster ID, size (`valueQuantity`, `bp`), number of contigs, relaxase, MPF, oriT, predicted mobility, nearest reference plasmid and distance, host range, and carried AMR genes (as one valueString). `hasMember` links to the variant Observations for the genes it carries.
  - **One plasmid-analysis summary Observation**, always present. Its value is the number of plasmids as `valueInteger`, or `dataAbsentReason` when status is failed. Components: status, assembly type/confidence, MOB-suite version. `hasMember` links to the plasmid Observations.
  - **Do not** use the LOINC 48018-6 "Gene studied" component on plasmid Observations. The merge step counts that component on every Observation, so the genes would be counted twice.
  - IDs are `{sample}-plasmid-{cluster}` and `{sample}-plasmid-summary`, created through the existing `sanitize_id()`.

### 7. `bin/merge_clinical_fhir.py`
- **Line 715.** Add the plasmid-analysis and plasmid codes to the skip list next to `'612-2', 'SP000678', 'SP000680', 'SP000682'`. This keeps the resistance classification unchanged.
- **`create_diagnostic_report()`** (around lines 1117–1157):
  - Add a "Plasmid Content" HTML table: cluster, replicons, size, mobility, genes, plus the confidence note.
  - Tag items in the "Detected Resistance Genes" list with their location.
  - Add to the conclusion only when a carbapenemase or ESBL gene is on a plasmid, e.g. "KPC-2 on conjugative IncFII(K) plasmid". Also add a text `conclusionCode` entry such as "Plasmid-borne carbapenemase".
- `result[]` already includes every Observation that has an `id` (line 1232), so nothing more is needed there.

### 8. `workflows/report.nf` and `bin/generate_sample_report.py`
- `GENERATE_REPORT` takes the plasmid channel. Use `typing.join(lineage).join(plasmid)`.
- Add an optional `--plasmid-json` argument. It produces a "PLASMID ANALYSIS" section listing each plasmid, the location of each resistance gene and the confidence note. When the JSON is absent, write "not run".

### 9. Versions
- `bin/get_versions.py`: add MOB-suite to the `tools` dict (around line 76), calling `mob_recon --version` from `params.mobsuite_env`.
- `workflows/utils.nf`: add `plasmid: tool 'mob_suite (mob_recon)'` and `amr_location: 'minimap2 vs Kleborate CARD'`.

## Known issue (existing code, out of scope here)

In Kleborate, `^` after a gene name means the protein matches exactly but the DNA doesn't (`kleborate/shared/resMinimap.py:170-171`). It says nothing about location. Even so, `annotated_to_fhir.py:462` labels every `^` gene "(chromosomal)". `merge_clinical_fhir.py:710` then leaves those genes out of the resistance classification unless they are carbapenemase, ESBL, AmpC or colistin genes. A `sul1^`, for example, would be dropped.

Fix this separately. This plan keeps that logic unchanged, so the classification doesn't shift while plasmid reporting is added.

## Verification

1. **Install:** `conda run -n mobsuite mob_recon --version` works, and `~/medaka-env` is unchanged.
2. **Scripts on their own:** run `PLASMID_RECON`'s commands and `parse_plasmids.py` by hand on two samples.
   - `results/assemblies/SRR25759560_contigs.fasta` (Illumina; CTX-M-15, qnrS1, sul1… which are usually plasmid-borne).
   - `ERR10367345_assembly.fasta` (ONT; SHV-187, where chromosome is expected).
   - Check: every Kleborate gene gets a location, contig IDs join, and the JSON matches the schema.
3. **Pipeline:** from WSL, run `source ~/medaka-env/bin/activate && nextflow run main.nf -resume --skip_upload`.
   - Assemblies and Kleborate come from the cache.
   - Expect `results/plasmid/<sample>/`, `results/plasmid/plasmid_sharing.tsv` and new Observations in `results/fhir_merged/*.json`.
4. **Nothing regressed:**
   - The DiagnosticReport resistance classification and conclusion match today's outputs, apart from the added plasmid sentence.
   - Run with `--plasmid_analysis false`; the bundle should match today's apart from UUIDs and timestamps.
5. **FHIR validation:** `results/fhir_validated/*.validation.txt` should have no **new** `Error` lines compared with the current files (count them before and after). Warnings about the placeholder code system are expected.
6. **Upload:** only after review, re-run without `--skip_upload` against the dev gateway.
