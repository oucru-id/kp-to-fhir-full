# Test data

## `-profile test`: HS11286 (downloaded at run time)

| | |
|---|---|
| Organism | *Klebsiella pneumoniae* subsp. *pneumoniae* HS11286 |
| Assembly | RefSeq `GCF_000240185.1` (ASM24018v2), complete genome: chromosome and 6 plasmids |
| URL | https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/185/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.fna.gz |
| Size / md5 | 1,679,403 bytes, `d5ca66f37a3e464d1fa922748cfb87ff` (NCBI `md5checksums.txt`, checked 2026-11-10) |
| Why | ST11 with *bla*<sub>KPC-2</sub> on a plasmid: CI asserts ST11, the KPC-2 call, the classification and 0 FHIR validation errors |
| Expected classification | **XDR**, not CRE. Kleborate 3.2.4 also reports a PmrB `p.Gly20fs` frameshift as a colistin-resistance mutation and gives the genome resistance score 3 (carbapenemase with colistin resistance). KPtoFHIR follows Kleborate, as v1.4 did (see the README limitations). |

The file is not stored in this repository. NCBI data are in the public domain
(NCBI website and data usage policies). The test profile skips plasmid analysis
and cgMLST, so CI does not download the MOB-suite database or a cgMLST schema.

## `-profile test_stub`

`tests/data/stub/` holds empty placeholder reads, a 4-base FASTA and an empty
cgMLST schema directory. Stub blocks never read them; they only let
`-stub-run` exercise the wiring of every process. They are written by
`tests/fixtures/make_fixtures.py`.

## Synthetic fixtures

`tests/fixtures/` (Kleborate tables, MOB-suite reports, contigs, chewBBACA calls,
assembly QC and clinical CSVs) is invented data for the pytest suite; see
`tests/fixtures/make_fixtures.py`. Nothing in it comes from a patient.
