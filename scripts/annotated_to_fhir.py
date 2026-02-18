#!/usr/bin/env python3

import json
import argparse
from datetime import datetime, timezone
import os
import sys
import uuid
from collections import defaultdict
import glob
import re

def debug_print(message):
    print(f"DEBUG: {message}", file=sys.stderr)

def load_lineage_data(lineage_dir):
    lineage_data = {}
    
    if not lineage_dir or not os.path.exists(lineage_dir):
        return lineage_data
    
    lineage_files = glob.glob(os.path.join(lineage_dir, "*.lineage.json"))
    
    for lineage_file in lineage_files:
        try:
            with open(lineage_file, 'r') as f:
                lineage_info = json.load(f)
                
                sample_id = lineage_info.get('sample_id', 
                    os.path.basename(lineage_file).replace('.lineage.json', ''))
                
                possible_ids = [
                    sample_id,
                    sample_id.replace('_ont.fastq', ''),
                    sample_id.replace('_illumina', ''),
                    sample_id.replace('.fastq', ''),
                    sample_id.split('_')[0] if '_' in sample_id else sample_id,
                    os.path.basename(lineage_file).replace('.lineage.json', '')
                ]
                
                for pid in set(possible_ids):
                    if pid:
                        lineage_data[pid] = lineage_info
                
        except Exception as e:
            continue
    
    return lineage_data

def extract_sample_id_from_filename(filename):
    basename = os.path.basename(filename)
    
    basename = basename.replace('.typing.json', '')
    basename = basename.replace('.json', '')
    
    clean_basename = basename.replace('_ont', '').replace('_illumina', '').replace('.fastq', '')
    if '_' in clean_basename:
        clean_basename = clean_basename.split('_')[0]
    
    variations = [
        clean_basename, 
        basename,
        basename.replace('_ont', ''),
        basename.replace('_illumina', ''),
        basename.replace('.fastq', ''),
        basename.split('_')[0] if '_' in basename else basename 
    ]
    
    return variations

def get_drug_class_snomed_mapping(drug_class):
    class_mapping = {
        'beta_lactam': {'code': '765422000', 'display': 'Product containing beta-lactam (product)'},
        'bla': {'code': '765422000', 'display': 'Product containing beta-lactam (product)'},
        
        'aminoglycoside': {'code': '324116004', 'display': 'Product containing aminoglycoside (product)'},
        'agly': {'code': '324116004', 'display': 'Product containing aminoglycoside (product)'},
        
        'fluoroquinolone': {'code': '1010205001', 'display': 'Medicinal product containing fluoroquinolone and acting as antibacterial agent'},
        'flq': {'code': '1010205001', 'display': 'Medicinal product containing fluoroquinolone and acting as antibacterial agent'},
        
        'macrolide': {'code': '763878009', 'display': 'Product containing macrolide (product)'},
        'mls': {'code': '763878009', 'display': 'Product containing macrolide (product)'},
        
        'tetracycline': {'code': '66261008', 'display': 'Product containing tetracycline (medicinal product)'},
        'tet': {'code': '66261008', 'display': 'Product containing tetracycline (medicinal product)'},
        
        'sulfonamide': {'code': '763875007', 'display': 'Product containing sulfonamide (product)'},
        'sul': {'code': '763875007', 'display': 'Product containing sulfonamide (product)'},
        
        'trimethoprim': {'code': '32792001', 'display': 'Product containing trimethoprim (medicinal product)'},
        'tmt': {'code': '32792001', 'display': 'Product containing trimethoprim (medicinal product)'},
        'dfr': {'code': '32792001', 'display': 'Product containing trimethoprim (medicinal product)'},
        
        'chloramphenicol': {'code': '57191001', 'display': 'Product containing chloramphenicol (medicinal product)'},
        'phenicol': {'code': '57191001', 'display': 'Product containing chloramphenicol (medicinal product)'},
        
        'colistin': {'code': '73074003', 'display': 'Product containing colistin (medicinal product)'},
        'col': {'code': '73074003', 'display': 'Product containing colistin (medicinal product)'},
        
        'carbapenem': {'code': '350134005', 'display': 'Product containing carbapenem (product)'},
        
        'fosfomycin': {'code': '387065000', 'display': 'Product containing fosfomycin (medicinal product)'},
        'fos': {'code': '387065000', 'display': 'Product containing fosfomycin (medicinal product)'},
        
        'rifampicin': {'code': '77891000', 'display': 'Product containing rifampicin (medicinal product)'},
        'rif': {'code': '77891000', 'display': 'Product containing rifampicin (medicinal product)'}
    }
    
    normalized_class = drug_class.lower().strip().replace('-', '_').replace(' ', '_')
    return class_mapping.get(normalized_class)

def get_capsule_type_coding(capsule_type):
    return {
        'system': 'http://kaptive.holtlab.net/capsule',
        'code': capsule_type,
        'display': f'Klebsiella capsule type {capsule_type}'
    }

def get_mlst_st_coding(sequence_type):
    return {
        'system': 'http://pubmlst.org/klebsiella',
        'code': sequence_type.replace('ST', ''),
        'display': f'Klebsiella pneumoniae MLST {sequence_type}'
    }

def get_susceptibility_interpretation(status):
    if status == 'Resistant':
        return {
            "coding": [{
                "system": "http://loinc.org",
                "code": "LA6676-6",
                "display": "Resistant"
            }],
            "text": "Resistant"
        }
    else:
        return {
            "coding": [{
                "system": "http://loinc.org",
                "code": "LA24225-7",
                "display": "Susceptible"
            }],
            "text": "Susceptible"
        }

def create_susceptibility_panel(file_sample_id, typing_data):
    
    drug_status = {
        "Ampicillin": "Resistant",               
        "Amoxicillin-Clavulanate": "Susceptible",
        "Piperacillin-Tazobactam": "Susceptible",
        "Cefotaxime": "Susceptible",
        "Ceftazidime": "Susceptible",
        "Cefepime": "Susceptible",
        "Meropenem": "Susceptible",
        "Ciprofloxacin": "Susceptible",
        "Gentamicin": "Susceptible",
        "Amikacin": "Susceptible",
        "Trimethoprim-Sulfamethoxazole": "Susceptible",
        "Colistin": "Susceptible"
    }

    resistance = typing_data.get('resistance', {})
    carbapenemase = resistance.get('carbapenemase', 'None')
    esbl = resistance.get('esbl', 'None')
    detailed_resistance = resistance.get('detailed_resistance', [])
    point_mutations = resistance.get('point_mutations', [])

    has_carb = carbapenemase != 'None'
    has_esbl = esbl != 'None'
    
    aminoglycoside_gene_mapping = {
        'aac(3)-i':    {'Gentamicin': True, 'Amikacin': False},
        'aac(3)-ia':   {'Gentamicin': True, 'Amikacin': False},
        'aac(3)-ib':   {'Gentamicin': True, 'Amikacin': False},
        'aac(3)-ii':   {'Gentamicin': True, 'Amikacin': False},
        'aac(3)-iia':  {'Gentamicin': True, 'Amikacin': False},
        'aac(3)-iib':  {'Gentamicin': True, 'Amikacin': False},
        'aac(3)-iii':  {'Gentamicin': True, 'Amikacin': False},
        'aac(3)-iv':   {'Gentamicin': True, 'Amikacin': True},
        'aac(3)-vi':   {'Gentamicin': True, 'Amikacin': False},
        'aac(6\')-i':  {'Gentamicin': False, 'Amikacin': True},
        'aac(6\')-ia': {'Gentamicin': False, 'Amikacin': True},
        'aac(6\')-ib': {'Gentamicin': False, 'Amikacin': True},
        'aac(6\')-ib-cr': {'Gentamicin': False, 'Amikacin': True},
        'aac(6\')-ii': {'Gentamicin': True, 'Amikacin': True},
        'aac(6\')-30': {'Gentamicin': False, 'Amikacin': True},
        'ant(2\'\')-i':  {'Gentamicin': True, 'Amikacin': False},
        'ant(2\'\')-ia': {'Gentamicin': True, 'Amikacin': False},
        'ant(3\'\')-i':  {'Gentamicin': False, 'Amikacin': False},
        'ant(3\'\')-ia': {'Gentamicin': False, 'Amikacin': False},
        'ant(4\')-i':   {'Gentamicin': False, 'Amikacin': True},
        'ant(4\')-ia':  {'Gentamicin': False, 'Amikacin': True},
        'ant(4\')-ii':  {'Gentamicin': False, 'Amikacin': True},
        'ant(4\')-iia': {'Gentamicin': False, 'Amikacin': True},
        'aph(3\')-i':   {'Gentamicin': False, 'Amikacin': False},
        'aph(3\')-ia':  {'Gentamicin': False, 'Amikacin': False},
        'aph(3\')-ii':  {'Gentamicin': False, 'Amikacin': False},
        'aph(3\')-iia': {'Gentamicin': False, 'Amikacin': False},
        'aph(3\')-vi':  {'Gentamicin': False, 'Amikacin': True},
        'aph(3\')-via': {'Gentamicin': False, 'Amikacin': True},
        'aph(3\'\')-i': {'Gentamicin': False, 'Amikacin': False},
        'aph(3\'\')-ib': {'Gentamicin': False, 'Amikacin': False},
        'aph(6)-i':     {'Gentamicin': False, 'Amikacin': False},
        'aph(6)-ia':    {'Gentamicin': False, 'Amikacin': False},
        'aph(6)-id':    {'Gentamicin': False, 'Amikacin': False},
        'arma':   {'Gentamicin': True, 'Amikacin': True},
        'rmta':   {'Gentamicin': True, 'Amikacin': True},
        'rmtb':   {'Gentamicin': True, 'Amikacin': True},
        'rmtc':   {'Gentamicin': True, 'Amikacin': True},
        'rmtd':   {'Gentamicin': True, 'Amikacin': True},
        'rmte':   {'Gentamicin': True, 'Amikacin': True},
        'rmtf':   {'Gentamicin': True, 'Amikacin': True},
        'rmtg':   {'Gentamicin': True, 'Amikacin': True},
        'rmth':   {'Gentamicin': True, 'Amikacin': True},
        'npma':   {'Gentamicin': True, 'Amikacin': True},
    }
    
    has_ompk35_mutation = False
    has_ompk36_mutation = False
    has_flq_gene = False
    has_flq_mutation = False
    has_col_gene = False
    has_col_mutation = False
    has_trimethoprim = False
    has_sulfonamide = False
    
    for gene in detailed_resistance:
        gene_name = gene.get('gene', '').lower().replace('^', '').strip()
        drug_class = gene.get('drug_class', '').lower()
        
        if drug_class == 'aminoglycoside' or 'aminoglycoside' in drug_class:
            gene_lookup = gene_name.lower().replace('_', '-').replace(' ', '-')
            
            matched = False
            for known_gene, resistance_profile in aminoglycoside_gene_mapping.items():
                if known_gene in gene_lookup or gene_lookup.startswith(known_gene.split('-')[0]):
                    if resistance_profile.get('Gentamicin'):
                        drug_status["Gentamicin"] = "Resistant"
                    if resistance_profile.get('Amikacin'):
                        drug_status["Amikacin"] = "Resistant"
                    matched = True
                    break
            
            if not matched:
                rmt_genes = ['arma', 'rmta', 'rmtb', 'rmtc', 'rmtd', 'rmte', 'rmtf', 'rmtg', 'rmth', 'npma']
                if any(rmt in gene_lookup for rmt in rmt_genes):
                    drug_status["Gentamicin"] = "Resistant"
                    drug_status["Amikacin"] = "Resistant"
                    matched = True
            
            if not matched:
                if 'aac' in gene_lookup or 'ant' in gene_lookup or 'aph' in gene_lookup:
                    drug_status["Gentamicin"] = "Resistant"
                    debug_print(f"Unknown aminoglycoside gene {gene_name} - marking Gentamicin resistant (conservative)")
        
        if drug_class == 'fluoroquinolone' or 'qnr' in gene_name:
            has_flq_gene = True
        
        if 'aac(6\')-ib-cr' in gene_name or 'aac(6)-ib-cr' in gene_name:
            has_flq_gene = True
            drug_status["Amikacin"] = "Resistant"
            
        if 'trimethoprim' in drug_class or gene_name.startswith('dfr'):
            has_trimethoprim = True
            
        if 'sulfonamide' in drug_class or gene_name.startswith('sul'):
            has_sulfonamide = True
            
        if 'colistin' in drug_class or gene_name.startswith('mcr'):
            has_col_gene = True
    
    for mut in point_mutations:
        gene = mut.get('gene', '')
        
        if gene in ['GyrA', 'ParC']:
            has_flq_mutation = True
            
        if gene == 'OmpK35':
            has_ompk35_mutation = True
        if gene == 'OmpK36':
            has_ompk36_mutation = True
            
        if gene in ['MgrB', 'PhoP', 'PhoQ', 'PmrA', 'PmrB', 'CrrA', 'CrrB']:
            has_col_mutation = True

    if has_carb:
        drug_status["Meropenem"] = "Resistant"
        drug_status["Amoxicillin-Clavulanate"] = "Resistant"
        drug_status["Piperacillin-Tazobactam"] = "Resistant"
        drug_status["Cefotaxime"] = "Resistant"
        drug_status["Ceftazidime"] = "Resistant"
        drug_status["Cefepime"] = "Resistant"
    
    if has_esbl:
        drug_status["Cefotaxime"] = "Resistant"
        drug_status["Ceftazidime"] = "Resistant"
        drug_status["Cefepime"] = "Resistant"
        drug_status["Piperacillin-Tazobactam"] = "Resistant"
    
    if has_ompk35_mutation and has_ompk36_mutation:
        if has_esbl:
            drug_status["Meropenem"] = "Resistant"
            drug_status["Piperacillin-Tazobactam"] = "Resistant"
    elif has_ompk36_mutation:
        if has_esbl:
            drug_status["Piperacillin-Tazobactam"] = "Resistant"
    
    if has_flq_gene or has_flq_mutation:
        drug_status["Ciprofloxacin"] = "Resistant"
    
    if has_col_gene or has_col_mutation:
        drug_status["Colistin"] = "Resistant"
    
    if has_trimethoprim or has_sulfonamide:
        drug_status["Trimethoprim-Sulfamethoxazole"] = "Resistant"

    components = []
    
    drug_loinc = {
        "Ampicillin": "18864-9",
        "Amoxicillin-Clavulanate": "18861-5",
        "Piperacillin-Tazobactam": "18970-4",
        "Cefotaxime": "18876-3",
        "Ceftazidime": "18879-7",
        "Cefepime": "18879-7",
        "Meropenem": "18943-1",
        "Ciprofloxacin": "18906-8",
        "Gentamicin": "18928-2",
        "Amikacin": "18860-7",
        "Trimethoprim-Sulfamethoxazole": "18998-5",
        "Colistin": "18912-6"
    }
    
    for drug_name, status in drug_status.items():
        interpretation = get_susceptibility_interpretation(status)
        
        component = {
            "code": {
                "coding": [{
                    "system": "http://loinc.org",
                    "code": drug_loinc.get(drug_name),
                    "display": f"{drug_name} [Susceptibility]"
                }],
                "text": drug_name
            },
            "valueCodeableConcept": interpretation
        }
        components.append(component)

    total_drugs = len(drug_status)
    resistant_count = sum(1 for s in drug_status.values() if s == "Resistant")
    susceptible_count = total_drugs - resistant_count
    acquired_resistant_count = resistant_count - 1
    
    div_text = f"""<div xmlns="http://www.w3.org/1999/xhtml">
    <h3>Klebsiella pneumoniae Predicted Susceptibility Panel</h3>
    <p><strong>Sample:</strong> {file_sample_id}</p>
    <p><strong>Method:</strong> Whole Genome Sequencing - Genotypic Prediction</p>
    <p><strong>Panel:</strong> {total_drugs} antimicrobials</p>
    <p><strong>Summary:</strong> {resistant_count} Resistant (including 1 intrinsic), {susceptible_count} Susceptible</p>
    <table border="1">
        <tr><th>Antimicrobial</th><th>Predicted Result</th></tr>
        {"".join(f'<tr><td>{drug}</td><td><strong>{status}</strong></td></tr>' for drug, status in drug_status.items())}
    </table>
    </div>"""

    return {
        "resourceType": "Observation",
        "id": sanitize_id(f"{file_sample_id}-susceptibility-panel"),
        "meta": {
            "profile": ["http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/therapeutic-implication"],
            "tag": [
                {"system": "http://terminology.kemkes.go.id/sp", "code": "genomics", "display": "Genomics"}
            ]
        },
        "text": {
            "status": "generated",
            "div": div_text
        },
        "status": "final",
        "category": [
            {"coding": [{"system": "http://terminology.hl7.org/CodeSystem/observation-category", "code": "laboratory", "display": "Laboratory"}]},
            {"coding": [{"system": "http://terminology.hl7.org/CodeSystem/v2-0074", "code": "GE", "display": "Genetics"}]}
        ],
        "code": {
            "coding": [{
                "system": "http://loinc.org",
                "code": "29576-6",
                "display": "Bacterial susceptibility panel"
            }],
            "text": "Klebsiella pneumoniae Susceptibility Panel (Genotypic Prediction)"
        },
        "subject": {"reference": f"Patient/{file_sample_id}-patient"},
        "specimen": {"reference": f"Specimen/{file_sample_id}-specimen"},
        "effectiveDateTime": datetime.now(timezone.utc).isoformat(),
        "performer": [{"reference": "Organization/100007732"}],
        "interpretation": [{
            "coding": [{
                "system": "http://terminology.hl7.org/CodeSystem/v3-ObservationInterpretation",
                "code": "N" if acquired_resistant_count == 0 else "A",
                "display": "Normal" if acquired_resistant_count == 0 else "Abnormal"
            }],
            "text": f"{'No acquired resistance detected' if acquired_resistant_count == 0 else f'Acquired resistance detected ({acquired_resistant_count} drugs)'}"
        }],
        "component": components
    }

def sanitize_id(id_string):
    sanitized = re.sub(r'[^A-Za-z0-9\-]', '-', str(id_string))
    sanitized = re.sub(r'-+', '-', sanitized)
    sanitized = sanitized.strip('-')
    return sanitized if sanitized else "unknown"

def add_standard_kp_components(components, typing_data, sample_lineage_info):
    
    virulence_info = typing_data.get('virulence', {})

    for vf_name, vf_value in virulence_info.items():
        if vf_name not in ['virulence_score'] and vf_value not in ['Unknown', 'None', '']:
            components.append({
                "code": {"text": f"Virulence factor {vf_name}"},
                "valueString": vf_value
            })

    if sample_lineage_info:
        lineage = sample_lineage_info.get('lineage', {})
        if isinstance(lineage, dict):
            lineage_designation = lineage.get('lineage_designation', '')
            high_risk_clone = lineage.get('high_risk_clone', '')
            
            if lineage_designation:
                components.append({
                    "code": {"text": "KP lineage designation"},
                    "valueString": lineage_designation
                })
            
            if high_risk_clone:
                components.append({
                    "code": {"text": "High risk clone"},
                    "valueString": high_risk_clone
                })

def create_enhanced_annotation(gene, drug_class, typing_data, sample_lineage_info):
    gene_info = {
        'SHV-1': 'narrow-spectrum β-lactamase',
        'CTX-M': 'ESBL family',
        'TEM': 'β-lactamase family',
        'KPC': 'class A carbapenemase',
        'NDM': 'class B metallo-β-lactamase (carbapenemase)',
        'VIM': 'class B metallo-β-lactamase (carbapenemase)',
        'IMP': 'class B metallo-β-lactamase (carbapenemase)',
        'OXA-48': 'class D carbapenemase'
    }
    
    clean_gene = gene.replace('^', '').strip()
    gene_description = gene_info.get(clean_gene, f'{drug_class} resistance gene')
    
    annotation = f"Resistance gene {gene} - {gene_description}"
    if '^' in gene:
        annotation += " (chromosomal)"
    
    mlst_info = typing_data.get('mlst', {})
    sequence_type = mlst_info.get('sequence_type', 'Unknown')
    capsule_info = typing_data.get('capsule', {})
    k_type = capsule_info.get('k_type', 'Unknown')
    
    if sequence_type not in ['Unknown', 'Undetermined']:
        annotation += f" - KP MLST: {sequence_type}"
    if k_type not in ['Unknown', 'Capsule null']:
        annotation += f" - Capsule: {k_type}"
    
    if sample_lineage_info and isinstance(sample_lineage_info.get('lineage'), dict):
        lineage_designation = sample_lineage_info['lineage'].get('lineage_designation', '')
        if lineage_designation:
            annotation += f" - Lineage: {lineage_designation}"
    
    return annotation

def load_cgmlst_data(cgmlst_file, sample_id):
    if not cgmlst_file or not os.path.exists(cgmlst_file):
        return None
        
    try:
        with open(cgmlst_file, 'r') as f:
            header = f.readline().strip().split('\t')
            
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 2: continue
                
                row_sample = parts[0].replace('.fasta', '').replace('_contigs', '')
                
                if row_sample == sample_id or sample_id in row_sample:
                    alleles = parts[1:]
                    total_loci = len(alleles)
                    called_loci = sum(1 for a in alleles if a.isdigit() and a != '0')
                    
                    allele_map = {}
                    for i, allele in enumerate(alleles):
                        if i+1 < len(header):
                            locus = header[i+1]
                            allele_map[locus] = allele
                    
                    return {
                        'total_loci': total_loci,
                        'called_loci': called_loci,
                        'percentage': (called_loci / total_loci) * 100 if total_loci > 0 else 0,
                        'alleles': allele_map
                    }
    except Exception as e:
        debug_print(f"Error parsing cgMLST: {e}")
    return None

def create_mlst_observation(file_sample_id, typing_data):
    mlst_info = typing_data.get('mlst', {})
    sequence_type = mlst_info.get('sequence_type', 'Unknown')
    
    if sequence_type in ['Unknown', 'Undetermined', 'NA']:
        return None
        
    st_coding = get_mlst_st_coding(sequence_type)
    
    div_text = f"<div xmlns=\"http://www.w3.org/1999/xhtml\">Klebsiella pneumoniae MLST {sequence_type}</div>"
    
    return {
        "resourceType": "Observation",
        "id": sanitize_id(f"{file_sample_id}-mlst"),
        "meta": {
            "profile": ["http://hl7.org/fhir/StructureDefinition/Observation"],
            "tag": [{"system": "http://terminology.kemkes.go.id/sp", "code": "genomics", "display": "Genomics"}]
        },
        "text": {
            "status": "generated", 
            "div": div_text
        },
        "status": "final",
        "category": [
            {"coding": [{"system": "http://terminology.hl7.org/CodeSystem/observation-category", "code": "laboratory", "display": "Laboratory"}]},
            {"coding": [{"system": "http://terminology.hl7.org/CodeSystem/v2-0074", "code": "GE", "display": "Genetics"}]}
        ],
        "code": {
            "coding": [{
                "system": "http://loinc.org",
                "code": "612-2",  
                "display": "Bacterial strain [Type] in Isolate by Bacteria subtyping"
            }]
        },
        "valueCodeableConcept": {
            "coding": [st_coding],
            "text": sequence_type
        },
        "subject": {"reference": f"Patient/{file_sample_id}-patient"},
        "specimen": {"reference": f"Specimen/{file_sample_id}-specimen"},
        "effectiveDateTime": datetime.now(timezone.utc).isoformat(),
        "performer": [{"reference": "Organization/100007732"}]
    }

def create_virulence_score_observation(file_sample_id, typing_data):
    virulence_info = typing_data.get('virulence', {})
    virulence_score = virulence_info.get('virulence_score', '0')
    
    if virulence_score == '0':
        return None
        
    try:
        score_value = int(virulence_score)
        div_text = f"<div xmlns=\"http://www.w3.org/1999/xhtml\">Virulence score: {score_value}</div>"
        
        return {
            "resourceType": "Observation",
            "id": sanitize_id(f"{file_sample_id}-virulence-score"),
            "meta": {
                "profile": ["http://hl7.org/fhir/StructureDefinition/Observation"],
                "tag": [{"system": "http://terminology.kemkes.go.id/sp", "code": "genomics", "display": "Genomics"}]
            },
            "text": {
                "status": "generated",
                "div": div_text
            },
            "status": "final",
            "category": [
                {"coding": [{"system": "http://terminology.hl7.org/CodeSystem/observation-category", "code": "laboratory", "display": "Laboratory"}]},
                {"coding": [{"system": "http://terminology.hl7.org/CodeSystem/v2-0074", "code": "GE", "display": "Genetics"}]}
            ],
            "code": {
                "coding": [{
                    "system": "http://terminology.kemkes.go.id/sp",
                    "code": "SP000680",
                    "display": "Klebsiella pneumoniae virulence score [Numeric] by genomic analysis"
                }],
                "text": "Virulence score"
            },
            "valueQuantity": {
                "value": score_value,
                "unit": "score",
                "system": "http://unitsofmeasure.org",
                "code": "score"
            },
            "subject": {"reference": f"Patient/{file_sample_id}-patient"},
            "specimen": {"reference": f"Specimen/{file_sample_id}-specimen"},
            "effectiveDateTime": datetime.now(timezone.utc).isoformat(),
            "performer": [{"reference": "Organization/100007732"}]
        }
    except ValueError:
        return None

def create_capsule_observation(file_sample_id, typing_data):
    capsule_info = typing_data.get('capsule', {})
    k_type = capsule_info.get('k_type', 'Unknown')
    
    if k_type in ['Unknown', 'Capsule null', 'NA']:
        return None
        
    capsule_coding = get_capsule_type_coding(k_type)
    
    div_text = f"<div xmlns=\"http://www.w3.org/1999/xhtml\">Klebsiella pneumoniae capsule type {k_type}</div>"
    
    return {
        "resourceType": "Observation",
        "id": sanitize_id(f"{file_sample_id}-capsule"),
        "meta": {
            "profile": ["http://hl7.org/fhir/StructureDefinition/Observation"],
            "tag": [{"system": "http://terminology.kemkes.go.id/sp", "code": "genomics", "display": "Genomics"}]
        },
        "text": {
            "status": "generated", 
            "div": div_text
        },
        "status": "final",
        "category": [
            {"coding": [{"system": "http://terminology.hl7.org/CodeSystem/observation-category", "code": "laboratory", "display": "Laboratory"}]},
            {"coding": [{"system": "http://terminology.hl7.org/CodeSystem/v2-0074", "code": "GE", "display": "Genetics"}]}
        ],
        "code": {
            "coding": [{
                "system": "http://terminology.kemkes.go.id/sp", 
                "code": "SP000678",
                "display": "Klebsiella pneumoniae capsular type [Identifier] by genomic analysis"
            }],
            "text": "Capsule Type"
        },
        "valueCodeableConcept": {
            "coding": [capsule_coding],
            "text": k_type
        },
        "subject": {"reference": f"Patient/{file_sample_id}-patient"},
        "specimen": {"reference": f"Specimen/{file_sample_id}-specimen"},
        "effectiveDateTime": datetime.now(timezone.utc).isoformat(),
        "performer": [{"reference": "Organization/100007732"}]
    }

def main():
    parser = argparse.ArgumentParser(description='Convert KP typing data to FHIR format')
    parser.add_argument('--input', required=True, help='Path to typing JSON file')
    parser.add_argument('--sample_id', required=True, help='Sample ID')
    parser.add_argument('--output', required=True, help='Path to output FHIR JSON')
    parser.add_argument('--lineage_dir', help='Directory with lineage files')
    parser.add_argument('--cgmlst_file', help='Path to chewBBACA results_alleles.tsv')
    args = parser.parse_args()
    
    lineage_data = {}
    if args.lineage_dir:
        lineage_data = load_lineage_data(args.lineage_dir)
    
    file_sample_id_variations = extract_sample_id_from_filename(args.input)
    matched_sample_id = None
    sample_lineage_info = None
    
    for variation in file_sample_id_variations:
        if variation in lineage_data:
            matched_sample_id = variation
            sample_lineage_info = lineage_data[variation]
            break
    
    if not sample_lineage_info:
        debug_print("No lineage info found")
    
    try:
        if not os.path.exists(args.input):
            sys.exit(1)
        
        if os.path.getsize(args.input) == 0:
            sys.exit(1)

        with open(args.input, 'r') as f:
            typing_data = json.load(f)
        
        file_sample_id = file_sample_id_variations[0]
        bundles = defaultdict(list)
        observation_count = 0
        successful_annotations = 0
        resistance_data = typing_data.get('resistance', {})
        detailed_resistance = resistance_data.get('detailed_resistance', [])
        point_mutations = resistance_data.get('point_mutations', [])

        all_resistance_features = []
        
        for gene in detailed_resistance:
            gene['source_type'] = 'acquired'
            all_resistance_features.append(gene)
            
        for mut in point_mutations:
            mut['source_type'] = 'mutation'
            all_resistance_features.append(mut)

        if all_resistance_features:
            for idx, feature_info in enumerate(all_resistance_features):
                try:
                    components = []
                    
                    source_type = feature_info.get('source_type', 'acquired')
                    
                    if source_type == 'acquired':
                        gene = feature_info.get('gene', 'unknown')
                        drug_class = feature_info.get('drug_class', 'unknown')
                        if gene == 'unknown' or not gene: continue
                        clean_gene = gene.replace('^', '').strip()
                        display_name = clean_gene
                        code_system = "http://www.genenames.org/geneId"
                        
                    else:
                        gene = feature_info.get('gene', 'unknown')
                        mutation = feature_info.get('mutation', 'unknown')
                        drug_class = 'unknown'
                        
                        if gene in ['GyrA', 'ParC']: drug_class = 'fluoroquinolone'
                        elif gene in ['OmpK35', 'OmpK36']: drug_class = 'carbapenem'
                        elif gene == 'Colistin': drug_class = 'colistin'
                        
                        clean_gene = f"{gene} {mutation}"
                        display_name = clean_gene
                        code_system = "http://loinc.org"
                        
                        aa_map = {
                            'A':'Ala','C':'Cys','D':'Asp','E':'Glu','F':'Phe','G':'Gly','H':'His',
                            'I':'Ile','K':'Lys','L':'Leu','M':'Met','N':'Asn','P':'Pro','Q':'Gln',
                            'R':'Arg','S':'Ser','T':'Thr','V':'Val','W':'Trp','Y':'Tyr','*':'Ter'
                        }
                        
                        match = re.search(r'([A-Z\*])(\d+)([A-Z\*]?)', mutation)
                        if match:
                            ref_aa_char, pos, alt_aa_char = match.groups()
                            
                            if ref_aa_char in aa_map and (not alt_aa_char or alt_aa_char in aa_map):
                                ref_3 = aa_map.get(ref_aa_char)
                                alt_3 = aa_map.get(alt_aa_char, '?')
                                
                                phgvs_code = f"p.({ref_3}{pos}{alt_3})"
                                phgvs_display = f"{ref_3}{pos}{alt_3}"
                                
                                components.append({
                                    "code": {
                                        "coding": [{
                                            "system": "http://loinc.org",
                                            "code": "48005-3",
                                            "display": "Amino acid change (pHGVS)"
                                        }]
                                    },
                                    "valueCodeableConcept": {
                                        "coding": [{
                                            "system": "https://varnomen.hgvs.org",
                                            "code": phgvs_code,
                                            "display": phgvs_display
                                        }]
                                    }
                                })
                        
                    successful_annotations += 1
                    
                    components.append({
                        "code": {
                            "coding": [{
                                "system": "http://loinc.org",
                                "code": "48018-6",
                                "display": "Gene studied [ID]"
                            }]
                        },
                        "valueCodeableConcept": {
                            "coding": [{
                                "system": code_system,
                                "code": clean_gene,
                                "display": display_name
                            }],
                            "text": display_name
                        }
                    })
                    
                    if drug_class != 'unknown' and drug_class:
                        drug_mapping = get_drug_class_snomed_mapping(drug_class)
                        if drug_mapping:
                            components.append({
                                "code": {
                                    "coding": [{
                                        "system": "http://loinc.org",
                                        "code": "51961-1",
                                        "display": "Genetic variation's effect on drug efficacy"
                                    }]
                                },
                                "valueCodeableConcept": {
                                    "coding": [{
                                        "system": "http://snomed.info/sct",
                                        "code": drug_mapping['code'],  
                                        "display": drug_mapping['display']
                                    }],
                                    "text": drug_class
                                }
                            })
                    
                    add_standard_kp_components(components, typing_data, sample_lineage_info)

                    if source_type == 'acquired':
                        annotation_text = create_enhanced_annotation(gene, drug_class, typing_data, sample_lineage_info)
                    else:
                        annotation_text = f"Point mutation {mutation} in {gene} conferring resistance"
                        
                    div_text = f"<div xmlns=\"http://www.w3.org/1999/xhtml\">{annotation_text}</div>"
                
                    observation_id = sanitize_id(f"{file_sample_id}-resistance-{idx+1}")
                    
                    observation = {
                        "resourceType": "Observation",
                        "id": observation_id,
                        "meta": {
                            "profile": ["http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/variant"],
                            "tag": [{"system": "http://terminology.kemkes.go.id/sp", "code": "genomics", "display": "Genomics"}]
                        },
                        "text": {
                            "status": "generated", 
                            "div": div_text
                        },
                        "status": "final",
                        "category": [
                            {"coding": [{"system": "http://terminology.hl7.org/CodeSystem/observation-category", "code": "laboratory", "display": "Laboratory"}]},
                            {"coding": [{"system": "http://terminology.hl7.org/CodeSystem/v2-0074", "code": "GE", "display": "Genetics"}]}
                        ],
                        "code": {
                            "coding": [{
                                "system": "http://loinc.org",
                                "code": "69548-6",
                                "display": "Genetic variant assessment"
                            }],
                            "text": f"Resistance Variant: {clean_gene}"
                        },
                        "valueCodeableConcept": {
                            "coding": [{"system": "http://loinc.org", "code": "LA9633-4", "display": "Present"}],
                            "text": "Present"
                        },
                        "subject": {"reference": f"Patient/{file_sample_id}-patient"},
                        "specimen": {"reference": f"Specimen/{file_sample_id}-specimen"},
                        "effectiveDateTime": datetime.now(timezone.utc).isoformat(),
                        "performer": [{"reference": "Organization/100007732"}],
                        "component": components
                    }

                    bundles[file_sample_id].append({
                        "fullUrl": f"urn:uuid:{str(uuid.uuid4()).lower()}",
                        "resource": observation
                    })

                    observation_count += 1
                    
                except Exception as e:
                    import traceback
                    debug_print(f"Error processing resistance feature {idx}: {e}")
                    debug_print(traceback.format_exc())
                    continue
        else:
            debug_print("No resistance genes found - creating susceptible observation")
            try:
                components = []
                
                
                add_standard_kp_components(components, typing_data, sample_lineage_info)
                
                mlst_info = typing_data.get('mlst', {})
                sequence_type = mlst_info.get('sequence_type', 'Unknown')
                capsule_info = typing_data.get('capsule', {})
                k_type = capsule_info.get('k_type', 'Unknown')
                
                annotation_text = "Klebsiella pneumoniae isolate with no antimicrobial resistance genes detected"
                if sequence_type not in ['Unknown', 'Undetermined']:
                    annotation_text += f" - MLST: {sequence_type}"
                if k_type not in ['Unknown', 'Capsule null']:
                    annotation_text += f" - Capsule: {k_type}"
                
                if sample_lineage_info and isinstance(sample_lineage_info.get('lineage'), dict):
                    lineage_designation = sample_lineage_info['lineage'].get('lineage_designation', '')
                    if lineage_designation:
                        annotation_text += f" - Lineage: {lineage_designation}"
                
                div_text = f"<div xmlns=\"http://www.w3.org/1999/xhtml\">{annotation_text}</div>"
                
                susceptible_observation = {
                    "resourceType": "Observation",
                    "id": f"{file_sample_id}-susceptible",
                    "meta": {
                        "profile": [
                            "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/variant"
                        ],
                        "tag": [
                            {
                                "system": "http://terminology.kemkes.go.id/sp",
                                "code": "genomics",
                                "display": "Genomics"
                            }
                        ]
                    },
                    "text": {
                        "status": "generated",
                        "div": div_text
                    },
                    "status": "final",
                    "category": [
                        {
                            "coding": [{
                                "system": "http://terminology.hl7.org/CodeSystem/observation-category",
                                "code": "laboratory",
                                "display": "Laboratory"
                            }]
                        },
                        {
                            "coding": [{
                                "system": "http://terminology.hl7.org/CodeSystem/v2-0074", 
                                "code": "GE", 
                                "display": "Genetics"
                            }]
                        }
                    ],
                    "code": {
                        "coding": [{
                            "system": "http://loinc.org",
                            "code": "69548-6",
                            "display": "Genetic variant assessment"
                        }],
                        "text": "Antimicrobial Susceptibility Assessment"
                    },
                    "valueCodeableConcept": {
                        "coding": [
                            {
                                "system": "http://loinc.org",
                                "code": "LA9634-2", 
                                "display": "Absent"
                            },
                            {
                                "system": "http://snomed.info/sct",
                                "code": "444018003",
                                "display": "No antimicrobial resistance detected"
                            }
                        ],
                        "text": "No resistance genes detected"
                    },
                    "subject": {"reference": f"Patient/{file_sample_id}-patient"},
                    "specimen": {"reference": f"Specimen/{file_sample_id}-specimen"},
                    "effectiveDateTime": datetime.now(timezone.utc).isoformat(),
                    "performer": [{"reference": "Organization/100007732"}],
                    "component": components
                }
                
                bundles[file_sample_id].append({
                    "fullUrl": f"urn:uuid:{str(uuid.uuid4()).lower()}",
                    "resource": susceptible_observation
                })
                observation_count += 1
                
            except Exception as e:
                debug_print(f"Error creating susceptible observation: {e}")
        
        mlst_obs = create_mlst_observation(file_sample_id, typing_data)
        if mlst_obs:
            bundles[file_sample_id].append({
                "fullUrl": f"urn:uuid:{str(uuid.uuid4()).lower()}",
                "resource": mlst_obs
            })
            observation_count += 1
            
        virulence_obs = create_virulence_score_observation(file_sample_id, typing_data)
        if virulence_obs:
            bundles[file_sample_id].append({
                "fullUrl": f"urn:uuid:{str(uuid.uuid4()).lower()}",
                "resource": virulence_obs
            })
            observation_count += 1

        capsule_obs = create_capsule_observation(file_sample_id, typing_data)
        if capsule_obs:
            bundles[file_sample_id].append({
                "fullUrl": f"urn:uuid:{str(uuid.uuid4()).lower()}",
                "resource": capsule_obs
            })
            observation_count += 1

        susceptible_panel = create_susceptibility_panel(file_sample_id, typing_data)
        if susceptible_panel:
             bundles[file_sample_id].append({
                "fullUrl": f"urn:uuid:{str(uuid.uuid4()).lower()}",
                "resource": susceptible_panel
            })
             observation_count += 1
        
        cgmlst_info = None
        if args.cgmlst_file:
            cgmlst_info = load_cgmlst_data(args.cgmlst_file, matched_sample_id or args.sample_id)

        if cgmlst_info:
            try:
                div_text = f"""<div xmlns="http://www.w3.org/1999/xhtml">
                <strong>Core Genome MLST Analysis</strong><br/>
                Called Loci: {cgmlst_info['called_loci']}/{cgmlst_info['total_loci']} ({cgmlst_info['percentage']:.1f}%)<br/>
                Detailed allele profile included in components.
                </div>"""
                
                cgmlst_obs = {
                    "resourceType": "Observation",
                    "id": sanitize_id(f"{file_sample_id}-cgmlst"),
                    "meta": {
                        "profile": ["http://hl7.org/fhir/StructureDefinition/Observation"],
                        "tag": [{"system": "http://terminology.kemkes.go.id/sp", "code": "genomics", "display": "Genomics"}]
                    },
                    "status": "final",
                    "category": [{
                        "coding": [{
                            "system": "http://terminology.hl7.org/CodeSystem/observation-category",
                            "code": "laboratory",
                            "display": "Laboratory"
                        }]
                    }],
                    "code": {
                        "coding": [{
                            "system": "http://terminology.kemkes.go.id/sp",
                            "code": "SP000682",
                            "display": "Core genome MLST [Type] in Isolate by Sequencing"
                        }],
                        "text": "Core Genome MLST Profile"
                    },
                    "subject": {"reference": f"Patient/{file_sample_id}-patient"},
                    "effectiveDateTime": datetime.now(timezone.utc).isoformat(),
                    "performer": [{"reference": "Organization/100007732"}],
                    "component": [
                        {
                            "code": {
                                "coding": [{
                                    "system": "http://terminology.kemkes.go.id/sp",
                                    "code": "SP000690",
                                    "display": "Loci called"
                                }]
                            },
                            "valueQuantity": {
                                "value": cgmlst_info['called_loci'],
                                "unit": "{count}",
                                "system": "http://unitsofmeasure.org",
                                "code": "1"
                            }
                        },
                        {
                            "code": {
                                "coding": [{
                                    "system": "http://terminology.kemkes.go.id/sp",
                                    "code": "SP000691",
                                    "display": "Total loci in scheme"
                                }]
                            },
                            "valueQuantity": {
                                "value": cgmlst_info['total_loci'],
                                "unit": "{count}",
                                "system": "http://unitsofmeasure.org",
                                "code": "1"
                            }
                        },
                        {
                            "code": {
                                "coding": [{
                                    "system": "http://terminology.kemkes.go.id/sp",
                                    "code": "SP000692",
                                    "display": "Call rate"
                                }]
                            },
                            "valueQuantity": {
                                "value": round(cgmlst_info['percentage'], 1),
                                "unit": "%",
                                "system": "http://unitsofmeasure.org",
                                "code": "%"
                            }
                        }
                    ],
                    "text": {"status": "generated", "div": div_text}
                }
                
                for locus, allele in cgmlst_info['alleles'].items():
                    cgmlst_obs['component'].append({
                        "code": {
                            "coding": [{
                                "system": "https://www.cgmlst.org/ncs",
                                "code": locus,
                                "display": locus
                            }]
                        },
                        "valueInteger": allele
                    })
                
                bundles[file_sample_id].append({
                    "fullUrl": f"urn:uuid:{str(uuid.uuid4()).lower()}",
                    "resource": cgmlst_obs
                })
                observation_count += 1
            except Exception as e:
                debug_print(f"Error creating cgMLST observation: {e}")

        fhir_output = {
            "resourceType": "Bundle",
            "type": "collection",
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "total": observation_count,
            "entry": []
        }

        for sample_id, entries in bundles.items():
            fhir_output['entry'].extend(entries)

        with open(args.output, 'w') as out:
            json.dump(fhir_output, out, indent=2)
        
    except Exception as e:
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()