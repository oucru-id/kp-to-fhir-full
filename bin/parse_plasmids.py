#!/usr/bin/env python3

import argparse
import csv
import glob
import hashlib
import json
import os
import re
import sys
import traceback

KLEBORATE_MIN_SPURIOUS_COVERAGE = 40.0
KLEBORATE_MIN_OVERLAP = 50

def debug_print(message):
    print(f"DEBUG: {message}", file=sys.stderr)

def find_kleborate_card():
    import kleborate
    data_dir = os.path.join(os.path.dirname(kleborate.__file__), 'modules', 'klebsiella_pneumo_complex__amr', 'data')
    matches = sorted(glob.glob(os.path.join(data_dir, 'CARD_v*.fasta')))
    if not matches:
        sys.exit(f"No CARD_v*.fasta found in {data_dir}")
    return matches[-1]

def read_fasta(path):
    sequences = {}
    name = None
    chunks = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if name is not None:
                    sequences[name] = ''.join(chunks)
                name = line[1:].split()[0]
                chunks = []
            elif name is not None:
                chunks.append(line)
    if name is not None:
        sequences[name] = ''.join(chunks)
    return sequences

def read_tsv(path):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return []
    with open(path, newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def read_status(path):
    status = {}
    if path and os.path.exists(path):
        with open(path) as f:
            for line in f:
                if '=' in line:
                    key, value = line.rstrip('\n').split('=', 1)
                    status[key.strip()] = value.strip()
    return status

def first_value(row, *columns):
    for column in columns:
        value = (row.get(column) or '').strip()
        if value and value not in ['-', 'NA', 'nan']:
            return value
    return None

def split_field(value):
    if not value:
        return []
    return [v.strip() for v in re.split(r'[,;]', value) if v.strip() and v.strip() not in ['-', 'NA']]

def to_int(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None

def to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None

def read_allele_clusters(card_fasta):
    clusters = {}
    if not card_fasta or not os.path.exists(card_fasta):
        return clusters
    with open(card_fasta) as f:
        for line in f:
            if line.startswith('>'):
                parts = line[1:].strip().split('__')
                if len(parts) >= 4:
                    clusters.setdefault('__'.join(parts[2:-1]), parts[0])
    return clusters

def read_paf(path):
    hits = []
    if not path or not os.path.exists(path):
        return hits
    with open(path) as f:
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 12:
                continue
            header_parts = fields[0].split('__')
            if len(header_parts) < 4:
                continue
            query_length, query_start, query_end = int(fields[1]), int(fields[2]), int(fields[3])
            matching_bases, block_length = int(fields[9]), int(fields[10])
            score = next((int(tag[5:]) for tag in fields[12:] if tag.startswith('AS:i:')), 0)
            hits.append({
                'query': fields[0],
                'cluster': header_parts[0],
                'allele': '__'.join(header_parts[2:-1]),
                'contig': fields[5],
                'strand': fields[4],
                'start': int(fields[7]),
                'end': int(fields[8]),
                'identity': 100.0 * matching_bases / block_length if block_length else 0.0,
                'coverage': 100.0 * (query_end - query_start) / query_length if query_length else 0.0,
                'score': score,
            })
    return hits

def hits_overlap(a, b):
    if a['contig'] != b['contig'] or a['strand'] != b['strand']:
        return False
    return min(a['end'], b['end']) - max(a['start'], b['start']) > KLEBORATE_MIN_OVERLAP

def best_hit_per_locus(hits, min_identity, min_coverage):
    candidates = [h for h in hits if h['identity'] >= min_identity and h['coverage'] >= KLEBORATE_MIN_SPURIOUS_COVERAGE and h['score'] > 0]
    candidates.sort(key=lambda h: (-(h['identity'] * h['score'] * h['coverage']), h['query']))
    loci = []
    for hit in candidates:
        if not any(hits_overlap(hit, kept) for kept in loci):
            loci.append(hit)
    return [locus for locus in loci if locus['coverage'] >= min_coverage]

def base_allele(gene):
    name = gene.split(' +')[0].strip()
    name = re.sub(r'-\d+(\.\d+)?%$', '', name)
    return name.rstrip('*?^')

def is_exact_call(gene):
    return not re.search(r'[\*\^]', gene.split(' +')[0])

def assign_loci(genes, loci, allele_clusters):
    claimed = set()
    assignments = {}

    for gene in genes:
        allele = base_allele(gene)
        exact = is_exact_call(gene)
        matched = [i for i, locus in enumerate(loci)
                   if i not in claimed and locus['allele'] == allele and (locus['identity'] >= 100.0) == exact]
        if matched:
            assignments[gene] = matched
            claimed.update(matched)

    for gene in genes:
        if gene in assignments:
            continue
        cluster = allele_clusters.get(base_allele(gene))
        matched = [i for i, locus in enumerate(loci) if i not in claimed and cluster and locus['cluster'] == cluster]
        if matched:
            best = max(matched, key=lambda i: loci[i]['identity'] * loci[i]['score'] * loci[i]['coverage'])
            assignments[gene] = [best]
            claimed.add(best)

    return assignments

def load_contig_molecules(contig_rows, assembly):
    by_id = {}
    by_md5 = {}
    for row in contig_rows:
        info = {
            'molecule': (row.get('molecule_type') or '').strip().lower(),
            'cluster_id': first_value(row, 'primary_cluster_id'),
        }
        contig_id = (row.get('contig_id') or '').strip()
        if contig_id:
            by_id[contig_id.split()[0]] = info
        md5 = (row.get('md5') or '').strip().lower()
        if md5:
            by_md5[md5] = info

    molecules = {}
    for name, sequence in assembly.items():
        info = by_id.get(name)
        if info is None:
            info = by_md5.get(hashlib.md5(sequence.encode()).hexdigest()) or \
                   by_md5.get(hashlib.md5(sequence.upper().encode()).hexdigest())
        molecules[name] = info
    return molecules

def classify_locus(locus, molecules, contig_lengths, chrom_min_contig):
    info = molecules.get(locus['contig'])
    if info and info['molecule'] == 'plasmid':
        return 'plasmid', info['cluster_id']
    if info and info['molecule'] == 'chromosome' and contig_lengths.get(locus['contig'], 0) >= chrom_min_contig:
        return 'chromosome', None
    return 'unresolved', None

def build_plasmids(contig_rows, mobtyper_rows):
    plasmids = {}
    for row in contig_rows:
        if (row.get('molecule_type') or '').strip().lower() != 'plasmid':
            continue
        cluster_id = first_value(row, 'primary_cluster_id') or 'unclustered'
        plasmid = plasmids.setdefault(cluster_id, {
            'cluster_id': cluster_id,
            'secondary_cluster_id': first_value(row, 'secondary_cluster_id'),
            'replicon_types': [],
            'relaxase_types': [],
            'mpf_type': None,
            'orit_types': [],
            'predicted_mobility': None,
            'size_bp': 0,
            'gc': None,
            'num_contigs': 0,
            'contigs': [],
            'mash_nearest_neighbor': None,
            'mash_distance': None,
            'predicted_host_range': None,
            'amr_genes': [],
            'carbapenemase_genes': [],
            'esbl_genes': [],
        })
        plasmid['contigs'].append((row.get('contig_id') or '').split()[0])
        plasmid['num_contigs'] += 1
        plasmid['size_bp'] += to_int(row.get('size')) or 0
        for rep in split_field(row.get('rep_type(s)')):
            if rep not in plasmid['replicon_types']:
                plasmid['replicon_types'].append(rep)

    for row in mobtyper_rows:
        cluster_id = first_value(row, 'primary_cluster_id') or (row.get('sample_id') or '').split(':')[-1]
        plasmid = plasmids.get(cluster_id)
        if plasmid is None:
            continue
        plasmid['secondary_cluster_id'] = first_value(row, 'secondary_cluster_id') or plasmid['secondary_cluster_id']
        plasmid['replicon_types'] = split_field(row.get('rep_type(s)')) or plasmid['replicon_types']
        plasmid['relaxase_types'] = split_field(row.get('relaxase_type(s)'))
        plasmid['mpf_type'] = first_value(row, 'mpf_type')
        plasmid['orit_types'] = split_field(row.get('orit_type(s)'))
        plasmid['predicted_mobility'] = first_value(row, 'predicted_mobility')
        plasmid['size_bp'] = to_int(row.get('size')) or plasmid['size_bp']
        plasmid['gc'] = to_float(row.get('gc'))
        plasmid['num_contigs'] = to_int(row.get('num_contigs')) or plasmid['num_contigs']
        plasmid['mash_nearest_neighbor'] = first_value(row, 'mash_nearest_neighbor')
        plasmid['mash_distance'] = to_float(row.get('mash_neighbor_distance'))
        plasmid['predicted_host_range'] = first_value(row, 'predicted_host_range_overall_name')

    return plasmids

def split_gene_list(value):
    if not value or value == 'None':
        return set()
    return {g.strip() for g in value.split(',') if g.strip()}

def analyse(args):
    with open(args.typing_json) as f:
        typing_data = json.load(f)

    status_info = read_status(args.recon_status)
    card_fasta = args.card_fasta or status_info.get('card_fasta')
    assembly = read_fasta(args.assembly)
    contig_lengths = {name: len(seq) for name, seq in assembly.items()}
    contig_rows = read_tsv(args.contig_report)
    mobtyper_rows = read_tsv(args.mobtyper)

    recon_ok = status_info.get('status') == 'ok' and bool(contig_rows)
    molecules = load_contig_molecules(contig_rows, assembly) if recon_ok else {}
    plasmids = build_plasmids(contig_rows, mobtyper_rows) if recon_ok else {}

    unmatched = [name for name, info in molecules.items() if info is None]
    if unmatched:
        debug_print(f"{len(unmatched)} contigs not found in the MOB-suite contig report")

    resistance = typing_data.get('resistance', {})
    genes = []
    for gene_info in resistance.get('detailed_resistance', []):
        gene = gene_info.get('gene')
        if gene and gene not in genes:
            genes.append(gene)

    loci = best_hit_per_locus(read_paf(args.paf), args.min_identity, args.min_coverage)
    assignments = assign_loci(genes, loci, read_allele_clusters(card_fasta))
    carbapenemases = split_gene_list(resistance.get('carbapenemase'))
    esbls = split_gene_list(resistance.get('esbl'))

    gene_locations = {}
    for gene in genes:
        copies = []
        for i in assignments.get(gene, []):
            locus = loci[i]
            molecule, cluster_id = classify_locus(locus, molecules, contig_lengths, args.chrom_min_contig)
            copies.append({
                'molecule': molecule,
                'cluster_id': cluster_id,
                'contig': locus['contig'],
                'start': locus['start'],
                'end': locus['end'],
                'identity': round(locus['identity'], 2),
                'coverage': round(locus['coverage'], 2),
            })
            if cluster_id in plasmids and gene not in plasmids[cluster_id]['amr_genes']:
                plasmids[cluster_id]['amr_genes'].append(gene)
                if gene in carbapenemases:
                    plasmids[cluster_id]['carbapenemase_genes'].append(gene)
                if gene in esbls:
                    plasmids[cluster_id]['esbl_genes'].append(gene)

        molecules_found = {c['molecule'] for c in copies}
        if 'plasmid' in molecules_found:
            location = 'plasmid'
        elif 'chromosome' in molecules_found:
            location = 'chromosome'
        else:
            location = 'unresolved'
        gene_locations[gene] = {'location': location, 'copies': copies}

    for mut in resistance.get('point_mutations', []):
        gene_locations[f"{mut.get('gene')} {mut.get('mutation')}"] = {'location': 'chromosome', 'copies': []}

    if not recon_ok:
        status = 'failed'
    elif plasmids:
        status = 'ok'
    else:
        status = 'no_plasmids'

    return {
        'sample_id': args.sample_id,
        'status': status,
        'assembly_type': args.assembly_type,
        'confidence': 'long-read' if args.assembly_type == 'ont' else 'predicted (short-read)',
        'tools': {
            'mob_suite': status_info.get('mob_suite', 'unknown'),
            'card_db': os.path.basename(card_fasta).replace('.fasta', '') if card_fasta else 'unknown',
        },
        'plasmids': list(plasmids.values()),
        'gene_locations': gene_locations,
    }

def main():
    parser = argparse.ArgumentParser(description='Locate Kleborate resistance genes on plasmids or the chromosome using MOB-suite results')
    parser.add_argument('--find_card', action='store_true', help="Print the path of Kleborate's CARD reference FASTA and exit")
    parser.add_argument('--sample_id', help='Sample ID')
    parser.add_argument('--typing_json', help='Path to typing JSON from parse_results')
    parser.add_argument('--assembly', help='Assembly FASTA that MOB-suite and minimap2 were run on')
    parser.add_argument('--assembly_type', default='illumina', choices=['illumina', 'ont'])
    parser.add_argument('--contig_report', help='MOB-suite contig_report.txt')
    parser.add_argument('--mobtyper', help='MOB-suite mobtyper_results.txt')
    parser.add_argument('--paf', help='minimap2 PAF of Kleborate CARD alleles against the assembly')
    parser.add_argument('--recon_status', help='Status file written by PLASMID_RECON')
    parser.add_argument('--card_fasta', help="Kleborate's CARD FASTA (default: path recorded in the status file)")
    parser.add_argument('--min_identity', type=float, default=90.0)
    parser.add_argument('--min_coverage', type=float, default=80.0)
    parser.add_argument('--chrom_min_contig', type=int, default=20000)
    parser.add_argument('--output', help='Path to output plasmid JSON')
    args = parser.parse_args()

    if args.find_card:
        print(find_kleborate_card())
        return

    missing = [name for name in ['sample_id', 'typing_json', 'assembly', 'output'] if not getattr(args, name)]
    if missing:
        parser.error(f"missing required arguments: {', '.join('--' + m for m in missing)}")

    try:
        result = analyse(args)
    except Exception:
        traceback.print_exc()
        result = {
            'sample_id': args.sample_id,
            'status': 'failed',
            'assembly_type': args.assembly_type,
            'confidence': 'long-read' if args.assembly_type == 'ont' else 'predicted (short-read)',
            'tools': {'mob_suite': 'unknown', 'card_db': 'unknown'},
            'plasmids': [],
            'gene_locations': {},
        }

    with open(args.output, 'w') as out:
        json.dump(result, out, indent=2)

if __name__ == "__main__":
    main()
