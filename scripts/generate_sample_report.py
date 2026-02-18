#!/usr/bin/env python3
import json
import datetime
import argparse

def is_clinically_significant_resistance(resistance_data):
    detailed_resistance = resistance_data.get('detailed_resistance', [])
    
    significant_genes = []
    for gene_info in detailed_resistance:
        gene = gene_info.get('gene', '')
        drug_class = gene_info.get('drug_class', '')
        
        if '^' in gene:
            if gene.upper().startswith('SHV-1'):
                continue
            if 'esbl' in drug_class.lower():
                significant_genes.append(gene_info)
                continue
            continue
        
        if any(marker in gene.upper() for marker in ['CTX-M', 'TEM', 'KPC', 'NDM', 'VIM', 'IMP', 'OXA-48']):
            significant_genes.append(gene_info)
        elif any(tet_gene in gene.upper() for tet_gene in ['TET(A)', 'TET(B)', 'TET(C)', 'TET(D)', 'TET(E)', 'TET(K)', 'TET(M)', 'TET(O)']):
            significant_genes.append(gene_info)
        elif drug_class.lower() in ['esbl', 'carbapenemase', 'fluoroquinolone', 'aminoglycoside', 'tetracycline']:
            significant_genes.append(gene_info)
    
    return significant_genes

def is_intrinsic_gene(gene, drug_class):
    if '^' in gene:
        if gene.upper().startswith('SHV-1'):
            return True
        if gene.upper().startswith('SHV') and 'esbl' not in drug_class.lower():
            return True
    return False

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample-id', required=True)
    parser.add_argument('--typing-json', required=True)
    parser.add_argument('--lineage-json', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    
    with open(args.typing_json, 'r') as f:
        typing_data = json.load(f)
    
    with open(args.lineage_json, 'r') as f:
        lineage_data = json.load(f)
    
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    sample_id = typing_data.get('sample_id', args.sample_id)
    resistance = typing_data.get('resistance', {})
    lineage = lineage_data.get('lineage', {})
    
    esbl = resistance.get('esbl', 'None')
    carbapenemase = resistance.get('carbapenemase', 'None')
    other_resistance = resistance.get('acquired_resistance', 'None')
    resistance_score = resistance.get('resistance_score', '0')
    num_resistance_classes = resistance.get('num_resistance_classes', '0')
    num_resistance_genes = resistance.get('num_resistance_genes', '0')
    
    significant_resistance = is_clinically_significant_resistance(resistance)
    
    detailed_resistance = resistance.get('detailed_resistance', [])
    point_mutations = resistance.get('point_mutations', [])
    
    intrinsic_genes = []
    for gene_info in detailed_resistance:
        gene = gene_info.get('gene', '')
        drug_class = gene_info.get('drug_class', '')
        if is_intrinsic_gene(gene, drug_class):
            intrinsic_genes.append(gene_info)
    
    resistance_status = "NOT DETECTED"
    if (int(resistance_score) > 0 or 
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
Generated: {current_date}

RESISTANCE ANALYSIS:
==================================================
Resistance Status: {resistance_status}
Resistance Score: {resistance_score}/3
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
    
    report_content += f'''

CLINICAL SUMMARY:
==============================
Sample: {sample_id}
Resistance Status: {resistance_status}
Clinical Action: {clinical_action}
Risk Level: {risk_level}
ESBL Producer: {'YES' if esbl != 'None' else 'NO'}
Carbapenemase Producer: {'YES' if carbapenemase != 'None' else 'NO'}
Multidrug Resistant: {'YES' if len(significant_resistance) >= 3 else 'NO'}
Clinically Significant Resistance: {'YES' if len(significant_resistance) > 0 else 'NO'}

================================================================================
END OF RESISTANCE ANALYSIS REPORT
================================================================================
'''
    
    with open(args.output, 'w') as f:
        f.write(report_content)
    
if __name__ == "__main__":
    main()