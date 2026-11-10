#!/usr/bin/env python3
import json
import datetime
import argparse

# Kleborate writes drug classes as short codes (AGly_acquired -> 'agly')
DRUG_CLASS_ALIASES = {'agly': 'aminoglycoside', 'flq': 'fluoroquinolone', 'tet': 'tetracycline'}

def is_clinically_significant_resistance(resistance_data):
    detailed_resistance = resistance_data.get('detailed_resistance', [])

    significant_genes = []
    for gene_info in detailed_resistance:
        gene = gene_info.get('gene', '')
        drug_class = gene_info.get('drug_class', '')
        drug_class = DRUG_CLASS_ALIASES.get(drug_class.lower(), drug_class)

        if is_intrinsic_gene(gene_info):
            continue

        if any(marker in gene.upper() for marker in ['CTX-M', 'TEM', 'KPC', 'NDM', 'VIM', 'IMP', 'OXA-48']):
            significant_genes.append(gene_info)
        elif any(tet_gene in gene.upper() for tet_gene in ['TET(A)', 'TET(B)', 'TET(C)', 'TET(D)', 'TET(E)', 'TET(K)', 'TET(M)', 'TET(O)']):
            significant_genes.append(gene_info)
        elif drug_class.lower() in ['esbl', 'carbapenemase', 'fluoroquinolone', 'aminoglycoside', 'tetracycline']:
            significant_genes.append(gene_info)
    
    return significant_genes

def is_intrinsic_gene(gene_info):
    # Set by parse_results for Kleborate's Bla_chr column (chromosomal SHV/LEN/OKP)
    return bool(gene_info.get('intrinsic', False))

def format_plasmid_section(plasmid_json):
    section = '''

PLASMID ANALYSIS:
============================================================
'''
    if not plasmid_json:
        return section + "Not run\n"

    try:
        with open(plasmid_json, 'r') as f:
            plasmid_data = json.load(f)
    except (OSError, ValueError):
        return section + "Failed - plasmid results could not be read\n"

    status = plasmid_data.get('status', 'failed')
    mob_version = plasmid_data.get('tools', {}).get('mob_suite', 'unknown')
    plasmids = plasmid_data.get('plasmids', [])
    section += f"Status: {status} (MOB-suite {mob_version})\n"
    section += f"Confidence: {plasmid_data.get('confidence', 'unknown')}\n"
    if status == 'failed':
        return section + "Plasmid content is unknown.\n"
    section += f"Plasmids detected: {len(plasmids)}\n"

    for plasmid in plasmids:
        replicons = ', '.join(plasmid.get('replicon_types') or []) or 'replicon not typed'
        size_bp = plasmid.get('size_bp')
        size_text = f"{size_bp / 1000:.1f} kb" if size_bp else 'size unknown'
        genes = ', '.join(plasmid.get('amr_genes') or []) or 'none'
        section += f'''
Plasmid {plasmid.get('cluster_id')}: {replicons} | {size_text} | {plasmid.get('predicted_mobility') or 'mobility unknown'}
  Resistance genes: {genes}'''

    gene_locations = plasmid_data.get('gene_locations', {})
    if gene_locations:
        section += "\n\nResistance gene locations:"
        for gene, info in gene_locations.items():
            location = info.get('location', 'unresolved')
            clusters = [c.get('cluster_id') for c in info.get('copies', []) if c.get('molecule') == 'plasmid' and c.get('cluster_id')]
            if location == 'plasmid' and clusters:
                location = f"plasmid {', '.join(dict.fromkeys(clusters))}"
            section += f"\n  - {gene}: {location}"
    return section + "\n"

def load_qc_data(qc_json):
    if not qc_json:
        return None
    try:
        with open(qc_json, 'r') as f:
            return json.load(f)
    except (OSError, ValueError):
        return {'status': 'UNKNOWN', 'warnings': ['assembly QC results could not be read']}

def format_qc_section(qc_data):
    section = '''

SEQUENCING QC:
============================================================
'''
    if qc_data is None:
        return section + "Not available for this sample\n"

    coverage = qc_data.get('estimated_coverage')
    min_coverage = qc_data.get('min_coverage')
    coverage_text = f"{coverage:g}x" if coverage is not None else "unknown"
    min_text = f"{min_coverage:g}x" if min_coverage is not None else "unknown"
    section += f"Status: {qc_data.get('status', 'UNKNOWN')}\n"
    section += f"Platform: {qc_data.get('platform', 'unknown')}\n"
    section += f"Estimated coverage: {coverage_text} (minimum {min_text})\n"
    if qc_data.get('polishing'):
        section += f"Polishing: {qc_data.get('polishing_tool') or 'polishing'} {qc_data['polishing']}\n"
    for warning in qc_data.get('warnings') or []:
        section += f"WARNING: {warning}\n"
    return section

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample-id', required=True)
    parser.add_argument('--typing-json', required=True)
    parser.add_argument('--lineage-json', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--plasmid-json', help='Plasmid JSON from parse_plasmids.py')
    parser.add_argument('--qc-json', help='Assembly QC JSON from assembly_qc.py')
    args = parser.parse_args()

    qc_data = load_qc_data(args.qc_json)
    qc_warnings = (qc_data or {}).get('warnings') or []
    qc_status = (qc_data or {}).get('status', 'Not available')
    qc_banner = f"\nQUALITY WARNING: {'; '.join(qc_warnings)}. Interpret these results with caution." if qc_warnings else ""
    
    with open(args.typing_json, 'r') as f:
        typing_data = json.load(f)
    
    with open(args.lineage_json, 'r') as f:
        lineage_data = json.load(f)
    
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    sample_id = typing_data.get('sample_id', args.sample_id)
    resistance = typing_data.get('resistance', {})
    mlst = typing_data.get('mlst', {})
    o_locus = typing_data.get('o_locus', {})
    lineage = lineage_data.get('lineage', {})

    sequence_type = mlst.get('sequence_type', 'Unknown')
    o_antigen = o_locus.get('o_type', 'Unknown')
    
    esbl = resistance.get('esbl', 'None')
    carbapenemase = resistance.get('carbapenemase', 'None')
    other_resistance = resistance.get('acquired_resistance', 'None')
    resistance_score = str(resistance.get('resistance_score', 'Unknown'))
    num_resistance_classes = resistance.get('num_resistance_classes', '0')
    num_resistance_genes = resistance.get('num_resistance_genes', '0')
    
    significant_resistance = is_clinically_significant_resistance(resistance)
    
    detailed_resistance = resistance.get('detailed_resistance', [])
    point_mutations = resistance.get('point_mutations', [])
    
    intrinsic_genes = []
    for gene_info in detailed_resistance:
        if is_intrinsic_gene(gene_info):
            intrinsic_genes.append(gene_info)
    
    resistance_score_meaning = {
        '0': 'no ESBL or carbapenemase',
        '1': 'ESBL, no carbapenemase',
        '2': 'carbapenemase, no colistin resistance',
        '3': 'carbapenemase with colistin resistance'
    }.get(resistance_score)
    resistance_score_text = (f"{resistance_score}/3 (Kleborate: {resistance_score_meaning})"
                             if resistance_score_meaning else "Not determined")

    resistance_status = "NOT DETECTED"
    if ((resistance_score.isdigit() and int(resistance_score) > 0) or
        carbapenemase != 'None' or 
        esbl != 'None' or 
        len(significant_resistance) > 0 or
        len(point_mutations) > 0):
        resistance_status = "DETECTED"
    
    clinical_action = "MONITORING"
    risk_level = lineage.get('resistance_profile', {}).get('risk_level', 'Low')
    
    if carbapenemase != 'None':
        clinical_action = "URGENT - Carbapenemase detected"
    elif esbl != 'None':
        clinical_action = "REQUIRED - ESBL detected"
    elif len(significant_resistance) >= 3:
        clinical_action = "REQUIRED - Multidrug resistant"
    elif len(significant_resistance) > 0:
        clinical_action = "ATTENTION - Acquired resistance detected"
    elif len(point_mutations) > 0:
        clinical_action = "ATTENTION - Resistance mutations detected"
    
    significant_count = len(significant_resistance)
    intrinsic_count = len(intrinsic_genes)
    mutation_count = len(point_mutations)
    
    report_content = f'''================================================================================
KLEBSIELLA PNEUMONIAE RESISTANCE ANALYSIS SUMMARY REPORT - {sample_id}
================================================================================
Generated: {current_date}{qc_banner}

RESISTANCE ANALYSIS:
==================================================
Sequence Type (ST): {sequence_type}
O-antigen: {o_antigen}
Resistance Status: {resistance_status}
Resistance Score: {resistance_score_text}
Number of Resistance Classes: {num_resistance_classes}
Number of Resistance Genes: {num_resistance_genes}
Clinically Significant Genes: {significant_count}
Point Mutations: {mutation_count}
Intrinsic Genes: {intrinsic_count}
Risk Level: {risk_level}

CONFIRMED RESISTANCE DETERMINANTS:
============================================================
'''
    
    if carbapenemase != 'None':
        report_content += f'''
CARBAPENEMASE: {carbapenemase}
Action: Isolation precautions + carbapenem restriction
'''
    
    if esbl != 'None':
        report_content += f'''
ESBL: {esbl}
Action: Avoid 3rd generation cephalosporins
'''
    
    if significant_resistance:
        report_content += f'''
ACQUIRED RESISTANCE GENES: {len(significant_resistance)} detected
Clinical Impact: Acquired antimicrobial resistance present
'''
        for gene_info in significant_resistance:
            gene = gene_info.get('gene', 'Unknown')
            drug_class = gene_info.get('drug_class', 'Unknown')
            report_content += f'''
  - {gene} ({drug_class} resistance)'''

    if point_mutations:
        report_content += f'''
RESISTANCE POINT MUTATIONS: {len(point_mutations)} detected
Clinical Impact: Chromosomal mutations conferring resistance
'''
        for mut in point_mutations:
            gene = mut.get('gene', 'Unknown')
            mutation = mut.get('mutation', 'Unknown')
            report_content += f'''
  - {gene}: {mutation}'''
    
    if resistance_status == "NOT DETECTED" and not significant_resistance and not point_mutations:
        report_content += '''

No clinically significant resistance determinants detected.
'''
    
    if intrinsic_genes:
        report_content += f'''

INTRINSIC GENES DETECTED: {intrinsic_count}
============================================================
'''
        for gene_info in intrinsic_genes:
            gene = gene_info.get('gene', 'Unknown')
            drug_class = gene_info.get('drug_class', 'Unknown')
            report_content += f'''
Gene: {gene} (chromosomal/intrinsic)
Drug Class: {drug_class}
Clinical Significance: Low - naturally occurring
----------------------------------------'''
    
    if significant_resistance:
        report_content += f'''

DETAILED RESISTANCE PROFILE:
============================================================
'''
        for gene_info in significant_resistance:
            gene = gene_info.get('gene', 'Unknown')
            drug_class = gene_info.get('drug_class', 'Unknown')
            identity = gene_info.get('identity', 'Unknown')
            coverage = gene_info.get('coverage', 'Unknown')
            report_content += f'''
Gene: {gene}
Drug Class: {drug_class}
Identity: {identity}
Coverage: {coverage}
Clinical Significance: High - acquired resistance
----------------------------------------'''

    report_content += format_qc_section(qc_data)
    report_content += format_plasmid_section(args.plasmid_json)

    report_content += f'''

CLINICAL SUMMARY:
==============================
Sample: {sample_id}
Sequence Type (ST): {sequence_type}
O-antigen: {o_antigen}
Resistance Status: {resistance_status}
Clinical Action: {clinical_action}
Risk Level: {risk_level}
ESBL Producer: {'YES' if esbl != 'None' else 'NO'}
Carbapenemase Producer: {'YES' if carbapenemase != 'None' else 'NO'}
Multidrug Resistant: {'YES' if len(significant_resistance) >= 3 else 'NO'}
Clinically Significant Resistance: {'YES' if len(significant_resistance) > 0 else 'NO'}
Sequencing QC: {qc_status}

================================================================================
END OF RESISTANCE ANALYSIS REPORT
================================================================================
'''
    
    with open(args.output, 'w') as f:
        f.write(report_content)
    
if __name__ == "__main__":
    main()