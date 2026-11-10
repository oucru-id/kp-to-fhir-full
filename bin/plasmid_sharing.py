#!/usr/bin/env python3

import argparse
import csv
import json
import sys

def debug_print(message):
    print(f"DEBUG: {message}", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(description='Summarise MOB-suite plasmid clusters shared across samples')
    parser.add_argument('--inputs', nargs='+', required=True, help='Per-sample *_plasmid.json files')
    parser.add_argument('--output', required=True, help='Path to output TSV')
    args = parser.parse_args()

    clusters = {}
    for path in args.inputs:
        try:
            with open(path) as f:
                data = json.load(f)
        except (OSError, json.JSONDecodeError) as e:
            debug_print(f"Skipping {path}: {e}")
            continue

        for plasmid in data.get('plasmids', []):
            cluster_id = plasmid.get('cluster_id')
            if not cluster_id:
                continue
            entry = clusters.setdefault(cluster_id, {'replicons': [], 'mobility': [], 'samples': [], 'genes': []})
            for key, values in [('replicons', plasmid.get('replicon_types') or []),
                                ('mobility', [plasmid['predicted_mobility']] if plasmid.get('predicted_mobility') else []),
                                ('samples', [data.get('sample_id', path)]),
                                ('genes', plasmid.get('amr_genes') or [])]:
                for value in values:
                    if value not in entry[key]:
                        entry[key].append(value)

    rows = sorted(clusters.items(), key=lambda item: (-len(item[1]['samples']), item[0]))

    with open(args.output, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['cluster_id', 'replicon_types', 'mobility', 'n_samples', 'samples', 'amr_genes'])
        for cluster_id, entry in rows:
            writer.writerow([
                cluster_id,
                ','.join(entry['replicons']) or '-',
                ','.join(entry['mobility']) or '-',
                len(entry['samples']),
                ','.join(entry['samples']),
                ','.join(entry['genes']) or '-',
            ])

if __name__ == "__main__":
    main()
