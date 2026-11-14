#!/usr/bin/env python3

import argparse
import json


def main():
    parser = argparse.ArgumentParser(description='Write the per-sample assembly QC record used by the FHIR and text reports')
    parser.add_argument('--sample_id', required=True)
    parser.add_argument('--platform', required=True, choices=['illumina', 'ont'])
    parser.add_argument('--coverage', required=True, help='Estimated read coverage, or NA')
    parser.add_argument('--min_coverage', type=float, required=True)
    parser.add_argument('--polishing', required=True, help='polished or failed')
    parser.add_argument('--polishing_tool', default='medaka')
    parser.add_argument('--polishing_model', default='')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    try:
        coverage = round(float(args.coverage), 1)
    except ValueError:
        coverage = None

    warnings = []
    if coverage is None:
        coverage_status = 'UNKNOWN'
        warnings.append('read coverage could not be estimated')
    elif coverage < args.min_coverage:
        coverage_status = 'LOW'
        warnings.append(f'low read coverage ({coverage:g}x; minimum {args.min_coverage:g}x)')
    else:
        coverage_status = 'PASS'

    if args.polishing != 'polished':
        warnings.append(f'{args.polishing_tool} polishing failed, so the unpolished assembly was used')

    qc = {
        'sample_id': args.sample_id,
        'platform': args.platform,
        'estimated_coverage': coverage,
        'min_coverage': args.min_coverage,
        'coverage_status': coverage_status,
        'polishing': args.polishing,
        'polishing_tool': args.polishing_tool,
        'polishing_model': args.polishing_model or None,
        'status': 'WARNING' if warnings else 'PASS',
        'warnings': warnings,
    }
    with open(args.output, 'w') as f:
        json.dump(qc, f, indent=2)


if __name__ == '__main__':
    main()
