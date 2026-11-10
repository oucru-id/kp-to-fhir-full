#!/usr/bin/env python3
"""Regenerate the synthetic test fixtures.

Every value here is invented. Kleborate rows use the real Kleborate 3.2.4
column layout with made-up calls; contig sequences are pseudo-random DNA;
names, NIKs (prefix 99, which is not an Indonesian province code, plus one
leading-zero NIK), addresses and dates are synthetic.

Run from the repository root:  python3 tests/fixtures/make_fixtures.py
"""
import csv
import hashlib
import json
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
TEST_DATA = REPO / "tests" / "data"

KLEBORATE_COLUMNS = """strain
klebsiella_pneumo_complex__mlst__ST klebsiella_pneumo_complex__mlst__gapA klebsiella_pneumo_complex__mlst__infB
klebsiella_pneumo_complex__mlst__mdh klebsiella_pneumo_complex__mlst__pgi klebsiella_pneumo_complex__mlst__phoE
klebsiella_pneumo_complex__mlst__rpoB klebsiella_pneumo_complex__mlst__tonB
klebsiella_pneumo_complex__amr__AGly_acquired klebsiella_pneumo_complex__amr__Col_acquired
klebsiella_pneumo_complex__amr__Fcyn_acquired klebsiella_pneumo_complex__amr__Flq_acquired
klebsiella_pneumo_complex__amr__Gly_acquired klebsiella_pneumo_complex__amr__MLS_acquired
klebsiella_pneumo_complex__amr__Phe_acquired klebsiella_pneumo_complex__amr__Rif_acquired
klebsiella_pneumo_complex__amr__Sul_acquired klebsiella_pneumo_complex__amr__Tet_acquired
klebsiella_pneumo_complex__amr__Tgc_acquired klebsiella_pneumo_complex__amr__Tmt_acquired
klebsiella_pneumo_complex__amr__Bla_acquired klebsiella_pneumo_complex__amr__Bla_inhR_acquired
klebsiella_pneumo_complex__amr__Bla_ESBL_acquired klebsiella_pneumo_complex__amr__Bla_ESBL_inhR_acquired
klebsiella_pneumo_complex__amr__Bla_Carb_acquired klebsiella_pneumo_complex__amr___acquired
klebsiella_pneumo_complex__amr__Bla_chr klebsiella_pneumo_complex__amr__SHV_mutations
klebsiella_pneumo_complex__amr__Omp_mutations klebsiella_pneumo_complex__amr__Col_mutations
klebsiella_pneumo_complex__amr__Flq_mutations klebsiella_pneumo_complex__amr__truncated_resistance_hits
klebsiella_pneumo_complex__amr__spurious_resistance_hits
klebsiella_pneumo_complex__resistance_score__resistance_score
klebsiella_pneumo_complex__resistance_gene_count__num_resistance_genes
klebsiella_pneumo_complex__resistance_class_count__num_resistance_classes
klebsiella__ybst__YbST klebsiella__ybst__Yersiniabactin klebsiella__ybst__ybtS klebsiella__ybst__ybtX
klebsiella__ybst__ybtQ klebsiella__ybst__ybtP klebsiella__ybst__ybtA klebsiella__ybst__irp2 klebsiella__ybst__irp1
klebsiella__ybst__ybtU klebsiella__ybst__ybtT klebsiella__ybst__ybtE klebsiella__ybst__fyuA
klebsiella__ybst__spurious_ybt_hits
klebsiella__cbst__CbST klebsiella__cbst__Colibactin klebsiella__cbst__clbA klebsiella__cbst__clbB
klebsiella__cbst__clbC klebsiella__cbst__clbD klebsiella__cbst__clbE klebsiella__cbst__clbF klebsiella__cbst__clbG
klebsiella__cbst__clbH klebsiella__cbst__clbI klebsiella__cbst__clbL klebsiella__cbst__clbM klebsiella__cbst__clbN
klebsiella__cbst__clbO klebsiella__cbst__clbP klebsiella__cbst__clbQ klebsiella__cbst__spurious_clb_hits
klebsiella__abst__AbST klebsiella__abst__Aerobactin klebsiella__abst__iucA klebsiella__abst__iucB
klebsiella__abst__iucC klebsiella__abst__iucD klebsiella__abst__iutA klebsiella__abst__spurious_abst_hits
klebsiella__smst__SmST klebsiella__smst__Salmochelin klebsiella__smst__iroB klebsiella__smst__iroC
klebsiella__smst__iroD klebsiella__smst__iroN klebsiella__smst__spurious_smst_hits
klebsiella__rmst__RmST klebsiella__rmst__RmpADC klebsiella__rmst__rmpA klebsiella__rmst__rmpD klebsiella__rmst__rmpC
klebsiella__rmst__spurious_rmst_hits klebsiella__rmpa2__rmpA2
klebsiella_pneumo_complex__virulence_score__virulence_score
klebsiella_pneumo_complex__virulence_score__spurious_virulence_hits
klebsiella_pneumo_complex__kaptive__K_locus klebsiella_pneumo_complex__kaptive__K_type
klebsiella_pneumo_complex__kaptive__K_locus_confidence klebsiella_pneumo_complex__kaptive__K_locus_problems
klebsiella_pneumo_complex__kaptive__K_locus_identity klebsiella_pneumo_complex__kaptive__K_Missing_expected_genes
klebsiella_pneumo_complex__kaptive__O_locus klebsiella_pneumo_complex__kaptive__O_type
klebsiella_pneumo_complex__kaptive__O_locus_confidence klebsiella_pneumo_complex__kaptive__O_locus_problems
klebsiella_pneumo_complex__kaptive__O_locus_identity klebsiella_pneumo_complex__kaptive__O_Missing_expected_genes
klebsiella_pneumo_complex__wzi__wzi""".split()

ZERO_COLUMNS = ("YbST", "CbST", "AbST", "SmST", "RmST", "virulence_score", "resistance_score",
                "num_resistance_genes", "num_resistance_classes")


def column(key):
    matches = [c for c in KLEBORATE_COLUMNS if c.endswith("__" + key)]
    assert len(matches) == 1, (key, matches)
    return matches[0]


def mlst(st, alleles):
    loci = ("gapA", "infB", "mdh", "pgi", "phoE", "rpoB", "tonB")
    return {"ST": st, **dict(zip(loci, alleles))}


# --------------------------------------------------------------------------
# Samples: platform, Kleborate calls, and which other fixtures exist
# --------------------------------------------------------------------------
SAMPLES = {
    "CRE01": {
        "platform": "illumina", "plasmid": True,
        "calls": {**mlst("ST11", ("3", "3", "1", "1", "1", "1", "4")),
                  "AGly_acquired": "aac(6')-Ib-cr.v2;rmtB.v1",
                  "Flq_acquired": "aac(6')-Ib-cr.v2;qnrS1",
                  "Sul_acquired": "sul1", "Tmt_acquired": "dfrA12",
                  "Bla_acquired": "OXA-1", "Bla_ESBL_acquired": "CTX-M-15",
                  "Bla_Carb_acquired": "KPC-2", "Bla_chr": "SHV-11",
                  "Flq_mutations": "GyrA:p.Ser83Ile;ParC:p.Ser80Ile",
                  "resistance_score": "2", "num_resistance_genes": "9", "num_resistance_classes": "6",
                  "YbST": "485", "Yersiniabactin": "ybt 9; ICEKp3", "virulence_score": "1",
                  "K_locus": "KL64", "K_type": "K64", "K_locus_confidence": "Typeable",
                  "O_locus": "O1/O2v1", "O_type": "O2a", "O_locus_confidence": "Typeable",
                  "wzi": "wzi64"},
    },
    "HVKP01": {
        "platform": "ont", "plasmid": True,
        "qc": {"coverage": 64.2, "polishing": "polished"},
        "calls": {**mlst("ST23", ("2", "1", "1", "1", "9", "4", "12")),
                  "Bla_chr": "SHV-1",
                  "YbST": "326", "Yersiniabactin": "ybt 1; ICEKp10",
                  "CbST": "29-1LV", "Colibactin": "clb 2",
                  "AbST": "1", "Aerobactin": "iuc 1", "SmST": "18-1LV", "Salmochelin": "iro 1",
                  "RmST": "1", "RmpADC": "rmp 1; KpVP-1", "rmpA2": "rmpA2_3-47%",
                  "virulence_score": "5",
                  "K_locus": "KL1", "K_type": "K1", "K_locus_confidence": "Typeable",
                  "O_locus": "O1/O2v2", "O_type": "O1v2", "O_locus_confidence": "Typeable",
                  "wzi": "wzi1"},
    },
    "KPC258": {
        "platform": "illumina", "plasmid": False,
        "calls": {**mlst("ST258", ("3", "3", "1", "1", "1", "1", "79")),
                  "AGly_acquired": "aac(6')-Ib.v1;aadA2", "Sul_acquired": "sul1", "Tmt_acquired": "dfrA12",
                  "_acquired": "TEM-1D.v1^", "Bla_Carb_acquired": "KPC-3", "Bla_chr": "SHV-11^",
                  "Omp_mutations": "OmpK35-25%;OmpK36GD", "Col_mutations": "MgrB-62%",
                  "resistance_score": "3", "num_resistance_genes": "6", "num_resistance_classes": "5",
                  "K_locus": "KL106", "K_type": "unknown (KL106)", "K_locus_confidence": "Typeable",
                  "O_locus": "O2a", "O_type": "O2a", "O_locus_confidence": "Typeable",
                  "wzi": "wzi29"},
    },
    "SUS01": {
        "platform": "illumina", "plasmid": True,
        "calls": {**mlst("ST45", ("2", "1", "1", "6", "7", "1", "12")),
                  "Bla_chr": "SHV-1",
                  "K_locus": "KL62", "K_type": "K62", "K_locus_confidence": "Typeable",
                  "O_locus": "O1/O2v1", "O_type": "O1v1", "O_locus_confidence": "Typeable",
                  "wzi": "wzi101"},
    },
    "NOVEL01": {
        "platform": "ont", "plasmid": True, "plasmid_failed": True,
        "qc": {"coverage": 11.4, "polishing": "polished"},
        "calls": {"ST": "NA", "gapA": "2", "infB": "1", "mdh": "2*", "pgi": "1", "phoE": "9", "rpoB": "4",
                  "tonB": "-",
                  "AGly_acquired": "aph(3')-Ia.v1;strA.v1;strB.v1", "Sul_acquired": "sul2",
                  "Tet_acquired": "tet(A).v1", "Tmt_acquired": "dfrA14.v2*",
                  "Bla_ESBL_acquired": "CTX-M-15", "Bla_chr": "SHV-28",
                  "truncated_resistance_hits": "OXA-1-56%",
                  "resistance_score": "1", "num_resistance_genes": "7", "num_resistance_classes": "5",
                  "K_locus": "KL102", "K_type": "Capsule null", "K_locus_confidence": "Typeable",
                  "O_locus": "O3b", "O_type": "O3b", "O_locus_confidence": "Typeable"},
    },
    "KP_001": {
        "platform": "illumina", "plasmid": False,
        "calls": {**mlst("ST15", ("1", "1", "1", "1", "1", "1", "1")),
                  "Bla_ESBL_acquired": "CTX-M-15", "Bla_chr": "SHV-28",
                  "resistance_score": "1", "num_resistance_genes": "2", "num_resistance_classes": "1",
                  "K_locus": "KL24", "K_type": "K24", "K_locus_confidence": "Typeable",
                  "O_locus": "O1/O2v1", "O_type": "O1v1", "O_locus_confidence": "Typeable",
                  "wzi": "wzi24"},
    },
    "KP_002": {
        "platform": "illumina", "plasmid": False,
        "calls": {**mlst("ST45", ("2", "1", "1", "6", "7", "1", "12")),
                  "Bla_chr": "SHV-1",
                  "K_locus": "KL62", "K_type": "K62", "K_locus_confidence": "Typeable",
                  "O_locus": "O1/O2v1", "O_type": "O1v1", "O_locus_confidence": "Typeable",
                  "wzi": "wzi101"},
    },
}

EMPTY = {"EMPTY01": {"platform": "illumina"}}


def kleborate_table(sample_id, info):
    suffix = "assembly" if info["platform"] == "ont" else "contigs"
    header = "\t".join(KLEBORATE_COLUMNS) + "\n"
    if "calls" not in info:
        return header
    row = {c: "-" for c in KLEBORATE_COLUMNS}
    row["strain"] = f"{sample_id}_{suffix}"
    for key in ZERO_COLUMNS:
        row[column(key)] = "0"
    for key, value in info["calls"].items():
        row[column(key)] = value
    return header + "\t".join(row[c] for c in KLEBORATE_COLUMNS) + "\n"


# --------------------------------------------------------------------------
# MOB-suite reconstruction: contigs, MOB-suite reports and minimap2 hits of
# Kleborate's CARD alleles, laid out so every gene has a known location.
# --------------------------------------------------------------------------
CARD_IDS = {
    "KPC-2": ("117", "KPC_Bla", "1201"), "CTX-M-15": ("31", "CTXM1_Bla", "331"),
    "OXA-1": ("162", "OXA1_Bla", "1508"), "SHV-11": ("205", "SHV_Bla", "1794"),
    "SHV-1": ("205", "SHV_Bla", "1788"), "aac(6')-Ib-cr.v2": ("4", "AAC(6')-Ib_AGly_Flq", "41"),
    "rmtB.v1": ("230", "RmtB_AGly", "2046"), "qnrS1": ("221", "QnrS_Flq", "1965"),
    "sul1": ("246", "Sul1_Sul", "2161"), "dfrA12": ("57", "DfrA12_Tmt", "512"),
}

PLASMID_LAYOUT = {
    "CRE01": {
        "contig_prefix": "k141_{}_pilon",
        "contigs": [("chromosome", None, 22000), ("chromosome", None, 21000), ("plasmid", "AA001", 6000),
                    ("plasmid", "AB002", 5000), ("chromosome", None, 2500), ("plasmid", "AA001", 3000)],
        "plasmids": {
            "AA001": {"rep": "IncFII(K),IncFIB(K)", "relaxase": "MOBF", "mpf": "MPF_F", "orit": "MOBF",
                      "mobility": "conjugative", "neighbor": "CP011985", "distance": "0.0021",
                      "host": "Enterobacteriaceae", "secondary": "AA001_A"},
            "AB002": {"rep": "IncR", "relaxase": "MOBP", "mpf": "-", "orit": "MOBP",
                      "mobility": "mobilizable", "neighbor": "MN310379", "distance": "0.0105",
                      "host": "Enterobacteriaceae", "secondary": "-"},
        },
        "hits": [("KPC-2", 3, 1000), ("sul1", 3, 3000), ("dfrA12", 3, 4500), ("CTX-M-15", 4, 1200),
                 ("qnrS1", 4, 3200), ("OXA-1", 5, 500), ("SHV-11", 1, 8000),
                 ("aac(6')-Ib-cr.v2", 6, 200), ("rmtB.v1", 6, 1500)],
    },
    "HVKP01": {
        "contig_prefix": "contig_{}",
        "contigs": [("chromosome", None, 30000), ("plasmid", "AA552", 9000)],
        "plasmids": {
            "AA552": {"rep": "IncHI1B(pNDM-MAR),IncFIB(K)", "relaxase": "-", "mpf": "-", "orit": "-",
                      "mobility": "non-mobilizable", "neighbor": "AP006726", "distance": "0.0008",
                      "host": "Klebsiella pneumoniae", "secondary": "AA552_B"},
        },
        "hits": [("SHV-1", 1, 15000)],
    },
    "SUS01": {
        "contig_prefix": "k141_{}_pilon",
        "contigs": [("chromosome", None, 21000), ("chromosome", None, 1500)],
        "plasmids": {},
        "hits": [("SHV-1", 2, 300)],
    },
    "NOVEL01": {
        "contig_prefix": "contig_{}",
        "contigs": [("chromosome", None, 2000)],
        "plasmids": {},
        "hits": [],
    },
}

CONTIG_REPORT_COLUMNS = ["sample_id", "molecule_type", "primary_cluster_id", "secondary_cluster_id", "contig_id",
                         "size", "gc", "md5", "circularity_status", "rep_type(s)", "rep_type_accession(s)",
                         "relaxase_type(s)", "relaxase_type_accession(s)", "mpf_type", "mpf_type_accession(s)",
                         "orit_type(s)", "orit_accession(s)", "predicted_mobility", "mash_nearest_neighbor",
                         "mash_neighbor_distance", "mash_neighbor_identification", "repetitive_dna_id",
                         "repetitive_dna_type", "filtering_reason"]
MOBTYPER_COLUMNS = ["sample_id", "num_contigs", "size", "gc", "md5", "rep_type(s)", "rep_type_accession(s)",
                    "relaxase_type(s)", "relaxase_type_accession(s)", "mpf_type", "mpf_type_accession(s)",
                    "orit_type(s)", "orit_accession(s)", "predicted_mobility", "mash_nearest_neighbor",
                    "mash_neighbor_distance", "mash_neighbor_identification", "primary_cluster_id",
                    "secondary_cluster_id", "predicted_host_range_overall_rank",
                    "predicted_host_range_overall_name", "observed_host_range_ncbi_rank",
                    "observed_host_range_ncbi_name", "reported_host_range_lit_rank",
                    "reported_host_range_lit_name", "associated_pmid(s)"]
GENE_LENGTH = 870


def dna(seed, length):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def gc_percent(seq):
    return f"{100.0 * sum(seq.count(b) for b in 'GC') / len(seq):.2f}"


def tsv(columns, rows):
    return "\t".join(columns) + "\n" + "".join("\t".join(str(r.get(c, "-")) for c in columns) + "\n"
                                                for r in rows)


def write_plasmid_fixture(sample_id, layout, failed):
    names = [layout["contig_prefix"].format(i) for i in range(1, len(layout["contigs"]) + 1)]
    seqs = {name: dna(f"{sample_id}/{name}", length)
            for name, (_kind, _cluster, length) in zip(names, layout["contigs"])}
    asm_dir = HERE / "assembly"
    asm_dir.mkdir(parents=True, exist_ok=True)
    with open(asm_dir / f"{sample_id}.fasta", "w") as handle:
        for name, seq in seqs.items():
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i + 80] + "\n")

    out = HERE / "plasmid" / sample_id
    out.mkdir(parents=True, exist_ok=True)
    if failed:
        (out / "contig_report.txt").write_text("")
        (out / "mobtyper_results.txt").write_text("")
        (out / "amr_hits.paf").write_text("")
        (out / "recon_status.txt").write_text("status=failed\nmob_suite=3.1.9\ncard_fasta=CARD_v3.2.9.fasta\n")
        return

    contig_rows = []
    for name, (kind, cluster, _length) in zip(names, layout["contigs"]):
        plasmid = layout["plasmids"].get(cluster, {})
        contig_rows.append({
            "sample_id": sample_id, "molecule_type": kind, "primary_cluster_id": cluster or "-",
            "secondary_cluster_id": plasmid.get("secondary", "-"), "contig_id": name,
            "size": len(seqs[name]), "gc": gc_percent(seqs[name]),
            "md5": hashlib.md5(seqs[name].encode()).hexdigest(), "circularity_status": "not tested",
            "rep_type(s)": plasmid.get("rep", "-"), "relaxase_type(s)": plasmid.get("relaxase", "-"),
            "mpf_type": plasmid.get("mpf", "-"), "orit_type(s)": plasmid.get("orit", "-"),
            "predicted_mobility": plasmid.get("mobility", "-"),
            "mash_nearest_neighbor": plasmid.get("neighbor", "-"),
            "mash_neighbor_distance": plasmid.get("distance", "-"), "filtering_reason": "-",
        })
    mobtyper_rows = []
    for cluster, plasmid in layout["plasmids"].items():
        members = [n for n, (_k, c, _l) in zip(names, layout["contigs"]) if c == cluster]
        joined = "".join(seqs[m] for m in members)
        mobtyper_rows.append({
            "sample_id": f"{sample_id}:{cluster}", "num_contigs": len(members), "size": len(joined),
            "gc": gc_percent(joined), "md5": hashlib.md5(joined.encode()).hexdigest(),
            "rep_type(s)": plasmid["rep"], "relaxase_type(s)": plasmid["relaxase"], "mpf_type": plasmid["mpf"],
            "orit_type(s)": plasmid["orit"], "predicted_mobility": plasmid["mobility"],
            "mash_nearest_neighbor": plasmid["neighbor"], "mash_neighbor_distance": plasmid["distance"],
            "primary_cluster_id": cluster, "secondary_cluster_id": plasmid["secondary"],
            "predicted_host_range_overall_rank": "family",
            "predicted_host_range_overall_name": plasmid["host"],
        })
    paf = []
    for gene, contig_index, start in layout["hits"]:
        cluster_no, gene_class, seq_id = CARD_IDS[gene]
        name = names[contig_index - 1]
        query = f"{cluster_no}__{gene_class}__{gene}__{seq_id}"
        paf.append("\t".join(map(str, [query, GENE_LENGTH, 0, GENE_LENGTH, "+", name, len(seqs[name]), start,
                                       start + GENE_LENGTH, GENE_LENGTH, GENE_LENGTH, 60, "NM:i:0",
                                       f"ms:i:{2 * GENE_LENGTH}", f"AS:i:{2 * GENE_LENGTH}", "tp:A:P"])))
    (out / "contig_report.txt").write_text(tsv(CONTIG_REPORT_COLUMNS, contig_rows))
    (out / "mobtyper_results.txt").write_text(tsv(MOBTYPER_COLUMNS, mobtyper_rows) if mobtyper_rows else "")
    (out / "amr_hits.paf").write_text("\n".join(paf) + ("\n" if paf else ""))
    (out / "recon_status.txt").write_text("status=ok\nmob_suite=3.1.9\ncard_fasta=CARD_v3.2.9.fasta\n")


# --------------------------------------------------------------------------
# chewBBACA allele calls (FILE column = assembly basename) and assembly QC
# --------------------------------------------------------------------------
CGMLST_LOCI = [f"KP1_RS{n:05d}" for n in range(5, 45, 5)]
CGMLST_CALLS = {
    "CRE01": ["1", "4", "12", "INF-301", "7", "2", "LNF", "9"],
    "HVKP01": ["3", "1", "1", "2", "ASM", "5", "8", "1"],
    "KPC258": ["1", "4", "12", "6", "7", "2", "3", "9"],
    "SUS01": ["2", "2", "5", "NIPH", "1", "1", "4", "3"],
    "NOVEL01": ["LNF", "LNF", "3", "PLOT5", "LNF", "1", "LNF", "2"],
    "KP_001": ["5", "6", "7", "8", "9", "10", "11", "12"],
    "KP_002": ["2", "2", "5", "4", "1", "1", "4", "3"],
}


def write_cgmlst():
    (HERE / "cgmlst").mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(["FILE"] + CGMLST_LOCI)]
    for sample_id, calls in CGMLST_CALLS.items():
        suffix = "assembly" if SAMPLES[sample_id]["platform"] == "ont" else "contigs"
        lines.append("\t".join([f"{sample_id}_{suffix}"] + calls))
    (HERE / "cgmlst" / "results_alleles.tsv").write_text("\n".join(lines) + "\n")


def assembly_qc(sample_id, coverage, polishing, min_coverage=30.0):
    low = coverage < min_coverage
    warnings = [f"low read coverage ({coverage:g}x; minimum {min_coverage:g}x)"] if low else []
    return {"sample_id": sample_id, "platform": "ont", "estimated_coverage": coverage,
            "min_coverage": min_coverage, "coverage_status": "LOW" if low else "PASS",
            "polishing": polishing, "polishing_tool": "medaka", "polishing_model": "r941_min_sup_g507",
            "status": "WARNING" if warnings else "PASS", "warnings": warnings}


# --------------------------------------------------------------------------
# Clinical metadata (synthetic); same 29-column contract as TBtoFHIR
# --------------------------------------------------------------------------
PATIENT_COLUMNS = [
    "sample_id", "family_name", "given_name", "gender", "birth_date", "age", "nik",
    "diagnosis_date", "treatment_history", "contact_number", "address", "city", "state",
    "province_code", "city_code", "district_code", "village_code", "latitude", "longitude",
    "citizenship_status", "specimen_type_code", "specimen_type_display",
    "specimen_collection_method_code", "specimen_collection_method_display",
    "specimen_collection_method_text", "specimen_quantity_value", "specimen_quantity_unit",
    "specimen_collected_date", "specimen_received_date",
]


def patient_row(sample_id, n, **overrides):
    row = {
        "sample_id": sample_id, "family_name": f"Sintetis{n}", "given_name": "Pasien",
        "gender": "laki-laki" if n % 2 else "perempuan", "birth_date": f"19{70 + n}-0{1 + n % 9}-15",
        "age": str(50 - n), "nik": f"99010100000000{n:02d}", "diagnosis_date": "2026-01-05",
        "treatment_history": "new", "contact_number": f"+62-800-0000-{n:04d}",
        "address": f"Jl. Contoh No. {n}, Kec. Uji", "city": "Kota Contoh", "state": "Provinsi Uji",
        "province_code": "99", "city_code": "9901", "district_code": "990101",
        "village_code": "9901011001", "latitude": "-7.0000", "longitude": "110.0000",
        "citizenship_status": "WNI", "specimen_type_code": "119297000",
        "specimen_type_display": "Blood specimen", "specimen_collection_method_code": "VENIP",
        "specimen_collection_method_display": "Venipuncture",
        "specimen_collection_method_text": "Blood culture collection", "specimen_quantity_value": "10",
        "specimen_quantity_unit": "mL", "specimen_collected_date": "2026-01-10",
        "specimen_received_date": "2026-01-11",
    }
    row.update(overrides)
    return row


PATIENT_ROWS = [
    patient_row("CRE01", 1),
    patient_row("HVKP01", 2, family_name="", given_name="", birth_date="",
                specimen_collected_date="", specimen_received_date="", latitude="", longitude=""),
    patient_row("KPC258", 3, nik="0000012345678901"),
    patient_row("NOVEL01", 4),
    patient_row("KP_001", 5),
    patient_row("KP_002", 6),
    patient_row("ORPHAN01", 7),
]

ORG_ROW = {
    "org_id": "99000001", "name": "Laboratorium Sintetis", "alias": "LabSin",
    "type_code": "prov", "type_display": "Healthcare Provider", "type_text": "Laboratory",
    "phone": "+62-800-0000-0000", "email": "lab@example.org",
    "address_line": "Jl. Contoh No. 99", "city": "Kota Contoh", "state": "Provinsi Uji",
    "country": "ID", "province_code": "99", "city_code": "9901", "district_code": "990101",
    "latitude": "-7.0001", "longitude": "110.0001",
}

PRACTITIONER_ROW = {
    "practitioner_id": "N99999999", "nik": "9900000000000001", "name": "dr. Praktisi Sintetis",
    "phone": "+62-800-0000-0001", "gender": "female", "birth_date": "1980-01-01",
    "str_kki_number": "SYN-STR-0001", "qualification_period_start": "2020-01-01",
    "role_id": "SYN-ROLE-001", "role_code": "309343006", "role_display": "Physician",
    "role_text": "Clinical microbiologist",
}


def write_stub_inputs():
    """Placeholder inputs for `-profile test_stub` (stub blocks never read them)."""
    import gzip
    stub = TEST_DATA / "stub"
    (stub / "cgmlst_schema").mkdir(parents=True, exist_ok=True)
    for name in ("STUBILL_1.fastq.gz", "STUBILL_2.fastq.gz", "STUBONT.fastq.gz"):
        with open(stub / name, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as handle:
            handle.write(b"")
    (stub / "STUBFA.fasta").write_text(">stub\nACGT\n")
    (stub / "cgmlst_schema" / "README.md").write_text("Placeholder chewBBACA schema for -profile test_stub.\n")
    with open(TEST_DATA / "samplesheet_stub.csv", "w") as handle:
        handle.write("sample,platform,fastq_1,fastq_2,fasta\n")
        handle.write("STUBILL,illumina,tests/data/stub/STUBILL_1.fastq.gz,tests/data/stub/STUBILL_2.fastq.gz,\n")
        handle.write("STUBONT,ont,tests/data/stub/STUBONT.fastq.gz,,\n")
        handle.write("STUBFA,fasta,,,tests/data/stub/STUBFA.fasta\n")


def write_csv(path, columns, rows):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        handle.write("# SYNTHETIC test fixture - invented values, not patient data\n")
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    for sub in ("kleborate", "qc", "metadata"):
        (HERE / sub).mkdir(parents=True, exist_ok=True)

    manifest = {}
    for sample_id, info in {**SAMPLES, **EMPTY}.items():
        (HERE / "kleborate" / f"{sample_id}_kleborate_full.txt").write_text(kleborate_table(sample_id, info))
        entry = {"platform": info["platform"], "plasmid": bool(info.get("plasmid")), "qc": "qc" in info}
        if info.get("plasmid"):
            write_plasmid_fixture(sample_id, PLASMID_LAYOUT[sample_id], info.get("plasmid_failed", False))
        if "qc" in info:
            qc = assembly_qc(sample_id, info["qc"]["coverage"], info["qc"]["polishing"])
            (HERE / "qc" / f"{sample_id}_assembly_qc.json").write_text(json.dumps(qc, indent=2) + "\n")
        manifest[sample_id] = entry
    (HERE / "samples.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    write_cgmlst()
    write_csv(HERE / "metadata" / "patient_clinical_metadata.csv", PATIENT_COLUMNS, PATIENT_ROWS)
    write_csv(HERE / "metadata" / "organization_metadata.csv", list(ORG_ROW), [ORG_ROW])
    write_csv(HERE / "metadata" / "practitioner_metadata.csv", list(PRACTITIONER_ROW), [PRACTITIONER_ROW])
    write_stub_inputs()
    print(f"fixtures written under {HERE} and {TEST_DATA}")


if __name__ == "__main__":
    main()
