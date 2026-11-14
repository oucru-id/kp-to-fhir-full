#!/usr/bin/env python3

import argparse
import csv
import json
import sys

MLST_LOCI = ['gapA', 'infB', 'mdh', 'pgi', 'phoE', 'rpoB', 'tonB']

ST_COLUMNS = ['klebsiella_pneumo_complex__mlst__ST', 'ST', 'sequence_type', 'mlst_st']
K_LOCUS_COLUMNS = ['klebsiella_pneumo_complex__kaptive__K_locus', 'K_locus', 'Best_match_locus']
K_TYPE_COLUMNS = ['klebsiella_pneumo_complex__kaptive__K_type', 'K_type', 'Best_match_type']
O_LOCUS_COLUMNS = ['klebsiella_pneumo_complex__kaptive__O_locus', 'klebsiella_pneumo_complex__kaptive__O_type',
                   'O_locus', 'O_type']
WZI_COLUMNS = ['klebsiella_pneumo_complex__wzi__wzi', 'wzi']
VIRULENCE_COLUMNS = {
    'yersiniabactin': ['klebsiella__ybst__Yersiniabactin', 'Yersiniabactin', 'ybt', 'YbST'],
    'colibactin': ['klebsiella__cbst__Colibactin', 'Colibactin', 'clb', 'CbST'],
    'aerobactin': ['klebsiella__abst__Aerobactin', 'Aerobactin', 'iuc', 'AbST'],
    'salmochelin': ['klebsiella__smst__Salmochelin', 'Salmochelin', 'iro', 'SmST'],
    'rmpadc': ['klebsiella__rmst__RmpADC', 'RmpADC', 'rmpADC', 'rmpA_rmpC', 'RmpA'],
    'rmpa2': ['klebsiella__rmpa2__rmpA2', 'RmpA2', 'rmpA2'],
}
VIRULENCE_SCORE_COLUMNS = ['klebsiella_pneumo_complex__virulence_score__virulence_score', 'virulence_score']
RESISTANCE_SCORE_COLUMNS = ['klebsiella_pneumo_complex__resistance_score__resistance_score', 'resistance_score']
BLA_COLUMNS = [
    'klebsiella_pneumo_complex__amr__Bla_Carb_acquired',
    'klebsiella_pneumo_complex__amr__Bla_ESBL_acquired',
    'klebsiella_pneumo_complex__amr__Bla_ESBL_inhR_acquired',
    'klebsiella_pneumo_complex__amr__Bla_acquired',
    'klebsiella_pneumo_complex__amr__Bla_chr',
    'Bla_Carb_acquired', 'Bla_ESBL_acquired', 'Bla_ESBL_inhR_acquired', 'Bla_acquired', 'Bla_chr',
]
OTHER_RESISTANCE_COLUMNS = [
    'klebsiella_pneumo_complex__amr__AGly_acquired', 'klebsiella_pneumo_complex__amr__Col_acquired',
    'klebsiella_pneumo_complex__amr__Fcyn_acquired', 'klebsiella_pneumo_complex__amr__Flq_acquired',
    'klebsiella_pneumo_complex__amr__Gly_acquired', 'klebsiella_pneumo_complex__amr__MLS_acquired',
    'klebsiella_pneumo_complex__amr__Phe_acquired', 'klebsiella_pneumo_complex__amr__Rif_acquired',
    'klebsiella_pneumo_complex__amr__Sul_acquired', 'klebsiella_pneumo_complex__amr__Tet_acquired',
    'klebsiella_pneumo_complex__amr__Tmt_acquired', 'klebsiella_pneumo_complex__amr__Tgc_acquired',
    'AGly_acquired', 'Col_acquired', 'Fcyn_acquired', 'Flq_acquired', 'Gly_acquired',
    'MLS_acquired', 'Phe_acquired', 'Rif_acquired', 'Sul_acquired', 'Tet_acquired', 'Tgc_acquired',
]
MUTATION_COLUMNS = {
    'Fluoroquinolone': ['klebsiella_pneumo_complex__amr__Flq_mutations', 'Flq_mutations'],
    'Porin': ['klebsiella_pneumo_complex__amr__Omp_mutations', 'Omp_mutations'],
    'Colistin': ['klebsiella_pneumo_complex__amr__Col_mutations', 'Col_mutations'],
}
PROTEIN_NAMES = {
    'gyra': 'GyrA', 'parc': 'ParC', 'ompk35': 'OmpK35', 'ompk36': 'OmpK36',
    'mgrb': 'MgrB', 'pmra': 'PmrA', 'pmrb': 'PmrB', 'phop': 'PhoP', 'phoq': 'PhoQ',
    'crra': 'CrrA', 'crrb': 'CrrB',
}
ST_TO_CC = {
    '258': 'CC258', '11': 'CC258', '512': 'CC258', '340': 'CC258', '437': 'CC258',
    '147': 'CC147', '15': 'CC15', '16': 'CC16', '23': 'CC23', '37': 'CC37', '101': 'CC101',
    '14': 'CC14', '17': 'CC17', '20': 'CC20', '25': 'CC25', '268': 'CC268', '392': 'CC392',
    '307': 'CC307', '67': 'CC67', '91': 'CC91', '3067': 'CC3067',
}
HIGH_RISK_STS = ['11', '258', '147', '15']
RHINOSCLEROMATIS_STS = ['67', '68', '69']
OZAENAE_STS = ['90', '91', '92', '93', '95', '96', '97']


class KleborateError(Exception):
    pass


def first_value(row, columns, excluded):
    for col in columns:
        if col in row and row[col] not in excluded:
            return row[col]
    return None


def split_genes(value):
    return [g.strip() for g in value.split(';') if g.strip()]


class Result:
    def __init__(self):
        self.sequence_type = "Unknown"
        self.mlst_alleles = {locus: 'Unknown' for locus in MLST_LOCI}
        self.mlst_confidence = "Not determined"
        self.resistance_genes = []
        self.resistance_classes = set()
        self.esbl_genes = []
        self.carbapenemases = []
        self.other_resistance = []
        self.resistance_score = "Unknown"
        self.point_mutations = []
        self.k_locus = "Unknown"
        self.k_type = "Unknown"
        self.o_locus = "Unknown"
        self.o_type = "Unknown"
        self.virulence_genes = {}
        self.virulence_score = "Unknown"
        self.wzi = "Unknown"


def parse_kleborate_row(row, result):
    st = first_value(row, ST_COLUMNS, ['-', '', 'Unknown', '0', 'NF', 'NA'])
    if st is not None:
        result.sequence_type = st
        result.mlst_confidence = "High"

    for locus in MLST_LOCI:
        columns = [f'klebsiella_pneumo_complex__mlst__{locus}', locus, locus.upper(), f'mlst_{locus}']
        allele = first_value(row, columns, ['-', '', 'Unknown', 'NF', '0', 'NA'])
        if allele is not None:
            result.mlst_alleles[locus] = allele

    result.k_locus = first_value(row, K_LOCUS_COLUMNS, ['-', '', 'Unknown']) or result.k_locus
    result.k_type = first_value(row, K_TYPE_COLUMNS, ['-', '', 'Unknown']) or result.k_type
    o_locus = first_value(row, O_LOCUS_COLUMNS, ['-', '', 'Unknown', 'NA'])
    if o_locus is not None:
        result.o_locus = result.o_type = o_locus
    result.wzi = first_value(row, WZI_COLUMNS, ['-', '', 'Unknown']) or result.wzi

    for name, columns in VIRULENCE_COLUMNS.items():
        value = first_value(row, columns, ['-', '', '0', 'NA'])
        if value is not None:
            result.virulence_genes[name] = value
    result.virulence_score = first_value(row, VIRULENCE_SCORE_COLUMNS, ['-', '', 'NA']) or result.virulence_score

    for col in BLA_COLUMNS:
        if col in row and row[col] not in ['-', '', '0', 'NA']:
            intrinsic = col.endswith('Bla_chr')
            for gene in split_genes(row[col]):
                result.resistance_genes.append({'gene': gene, 'drug_class': 'beta-lactam', 'identity': '100%',
                                                'coverage': '100%', 'intrinsic': intrinsic})
                result.resistance_classes.add('beta-lactam')
                if 'ESBL' in col:
                    result.esbl_genes.append(gene)
                elif 'Carb' in col:
                    result.carbapenemases.append(gene)
                elif not intrinsic:
                    result.other_resistance.append(gene)

    for col in OTHER_RESISTANCE_COLUMNS:
        if col in row and row[col] not in ['-', '', '0', 'NA']:
            drug_class = col.split('__')[-1].split('_')[0].lower()
            for gene in split_genes(row[col]):
                result.resistance_genes.append({'gene': gene, 'drug_class': drug_class, 'identity': '100%',
                                                'coverage': '100%'})
                result.resistance_classes.add(drug_class)
                result.other_resistance.append(gene)

    result.resistance_score = first_value(row, RESISTANCE_SCORE_COLUMNS, ['-', '', 'NA']) or result.resistance_score

    for mutation_class, columns in MUTATION_COLUMNS.items():
        value = first_value(row, columns, ['-', '', '0', 'NA'])
        if value is None:
            continue
        for mut in split_genes(value):
            gene_target = mutation_class
            clean_mut = mut
            for prefix, protein in PROTEIN_NAMES.items():
                if mut.lower().startswith(prefix):
                    gene_target = protein
                    clean_mut = mut[len(prefix):].lstrip(':').strip() or mut
                    break
            result.point_mutations.append({'gene': gene_target, 'mutation': clean_mut, 'type': 'point_mutation'})


def parse_amr_only_rows(rows, result):
    for row in rows:
        gene_symbol = row.get('Gene_symbol', '')
        drug_class = row.get('Drug_class', '')
        if not gene_symbol or gene_symbol in ['-', '']:
            continue
        result.resistance_genes.append({
            'gene': gene_symbol,
            'drug_class': drug_class if drug_class not in ['-', ''] else 'unknown',
            'identity': row.get('Sequence_identity', '100%'),
            'coverage': row.get('Coverage', '100%'),
        })
        if drug_class not in ['-', '']:
            result.resistance_classes.add(drug_class)
        if any(esbl in gene_symbol.upper() for esbl in ['ESBL', 'CTX', 'SHV', 'TEM']):
            result.esbl_genes.append(gene_symbol)
        elif any(carb in gene_symbol.upper() for carb in ['KPC', 'NDM', 'VIM', 'IMP', 'OXA']):
            result.carbapenemases.append(gene_symbol)
        else:
            result.other_resistance.append(gene_symbol)


def parse_kleborate(path):
    try:
        with open(path, newline='', encoding='utf-8') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            columns = list(reader.fieldnames or [])
            rows = list(reader)
    except (OSError, UnicodeDecodeError, csv.Error) as exc:
        raise KleborateError(f"{path}: cannot read Kleborate output ({exc})") from exc

    result = Result()
    if 'strain' in columns:
        if not rows:
            raise KleborateError(f"{path}: Kleborate table has a header but no result row")
        parse_kleborate_row(rows[0], result)
    elif 'Input_file_name' in columns:
        parse_amr_only_rows(rows, result)
    else:
        raise KleborateError(f"{path}: not a Kleborate results table (no 'strain' column)")
    return result


def clonal_complex(sequence_type):
    if sequence_type.startswith('ST') and sequence_type[2:].isdigit():
        return ST_TO_CC.get(sequence_type[2:], 'Unknown')
    return 'Unknown'


def typing_json(sample_id, result, cc):
    return {
        'sample_id': sample_id,
        'assembly': sample_id,
        'mlst': {
            'sequence_type': result.sequence_type,
            'clonal_complex': cc,
            'confidence': result.mlst_confidence,
            'alleles': result.mlst_alleles,
        },
        'capsule': {
            'k_locus': result.k_locus,
            'k_type': result.k_type,
            'k_confidence': "High" if result.k_type != "Unknown" else "Low",
            'k_problems': 'None',
            'wzi': result.wzi,
        },
        'virulence': {
            'yersiniabactin': result.virulence_genes.get('yersiniabactin', 'Unknown'),
            'colibactin': result.virulence_genes.get('colibactin', 'Unknown'),
            'aerobactin': result.virulence_genes.get('aerobactin', 'Unknown'),
            'salmochelin': result.virulence_genes.get('salmochelin', 'Unknown'),
            'hypermucoidy': result.virulence_genes.get('rmpadc', 'Unknown'),
            'rmpa2': result.virulence_genes.get('rmpa2', 'Unknown'),
            'virulence_score': result.virulence_score,
        },
        'resistance': {
            'acquired_resistance': ', '.join(result.other_resistance) if result.other_resistance else 'None',
            'esbl': ', '.join(result.esbl_genes) if result.esbl_genes else 'None',
            'carbapenemase': ', '.join(result.carbapenemases) if result.carbapenemases else 'None',
            'other_resistance': ', '.join(r['gene'] for r in result.resistance_genes),
            'resistance_score': str(result.resistance_score),
            'num_resistance_classes': str(len(result.resistance_classes)),
            'num_resistance_genes': str(len(result.resistance_genes)),
            'detailed_resistance': result.resistance_genes,
            'point_mutations': result.point_mutations,
        },
        'o_locus': {
            'o_locus': result.o_locus,
            'o_type': result.o_type,
            'o_confidence': "High" if result.o_locus != "Unknown" else "Low",
        },
        'analysis_type': 'K. pneumoniae complete Kleborate analysis',
    }


def lineage_json(sample_id, result, cc):
    st = result.sequence_type
    st_number = st.replace('ST', '').split('-')[0]
    num_classes = len(result.resistance_classes)
    return {
        'sample_id': sample_id,
        'lineage': {
            'sequence_type': st,
            'clonal_complex': cc,
            'capsule_type': result.k_type,
            'o_type': result.o_type,
            'wzi_type': result.wzi,
            'confidence': f'MLST: {result.mlst_confidence}',
            'lineage_designation': st if st != 'Unknown' else 'Unknown ST',
            'high_risk_clone': ('Yes' if st_number in HIGH_RISK_STS
                                else 'Unknown' if st == 'Unknown' else 'No'),
            'subspecies': ('rhinoscleromatis' if st_number in RHINOSCLEROMATIS_STS
                           else 'ozaenae' if st_number in OZAENAE_STS else 'pneumoniae'),
            'resistance_profile': {
                'esbl_producer': 'Yes' if result.esbl_genes else 'No',
                'carbapenemase_producer': 'Yes' if result.carbapenemases else 'No',
                'multidrug_resistant': 'Yes' if num_classes >= 3 else 'No',
                'risk_level': ('Critical' if result.carbapenemases else 'High' if result.esbl_genes
                               else 'Moderate' if num_classes > 1 else 'Low'),
            },
        },
        'analysis_type': 'K. pneumoniae Kleborate complete analysis',
    }


def main(argv=None):
    p = argparse.ArgumentParser(description='Summarise one Kleborate results table as typing and lineage JSON')
    p.add_argument('--kleborate', required=True, help='Kleborate results table (TSV)')
    p.add_argument('--sample_id', required=True)
    p.add_argument('--typing_json', help='output typing JSON (default <sample_id>_typing.json)')
    p.add_argument('--lineage_json', help='output lineage JSON (default <sample_id>_lineage.json)')
    args = p.parse_args(argv)

    try:
        result = parse_kleborate(args.kleborate)
    except KleborateError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)

    cc = clonal_complex(result.sequence_type)
    outputs = ((args.typing_json or f"{args.sample_id}_typing.json", typing_json(args.sample_id, result, cc)),
               (args.lineage_json or f"{args.sample_id}_lineage.json", lineage_json(args.sample_id, result, cc)))
    for path, data in outputs:
        with open(path, 'w') as out:
            json.dump(data, out, indent=2)

    print(f"{args.sample_id}: {result.sequence_type}, K {result.k_type}, O {result.o_type}, "
          f"{len(result.resistance_genes)} resistance gene(s), {len(result.point_mutations)} mutation(s)")


if __name__ == "__main__":
    main()
