#!/usr/bin/env python3

import json
import argparse
from datetime import datetime, timezone
import os
import sys
import uuid
from collections import defaultdict
import re
import html
from clinical_metadata_parser import load_organization_metadata
from kp_codes import (PLASMID_CODE_SYSTEM, PLASMID_ANALYSIS_CODE, PLASMID_CODE, GENOMIC_LOCATION_CODE,
                      PLASMID_CODE_DISPLAYS, LOCATION_CODINGS,
                      QC_CODE_SYSTEM, SEQUENCING_QC_CODE, SEQUENCING_QC_DISPLAY)

# 'Colistin' is the class label older typing.json files used in place of the protein name
COLISTIN_MUTATION_GENES = ['Colistin', 'MgrB', 'PhoP', 'PhoQ', 'PmrA', 'PmrB', 'CrrA', 'CrrB']

# Kleborate writes drug classes as short codes (AGly_acquired -> 'agly')
DRUG_CLASS_ALIASES = {
    'agly': 'aminoglycoside', 'flq': 'fluoroquinolone', 'sul': 'sulfonamide',
    'tmt': 'trimethoprim', 'col': 'colistin', 'tet': 'tetracycline'
}

def debug_print(message):
    print(f"DEBUG: {message}", file=sys.stderr)

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
        'system': 'http://kaptive.holtlab.net/mlst',
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

def normalise_aminoglycoside_name(name):
    # Kleborate mixes "aac(6')-Ib-cr.v2" and "aph3-Ia.v1"; drop the punctuation so both styles compare
    for char in "()'^*?":
        name = name.replace(char, '')
    return name.lower().replace('_', '-').replace(' ', '-')

def create_susceptibility_panel(file_sample_id, typing_data, org_id='unknown-org'):
    
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
    
    # Longest key first, so aac3-iia is preferred over aac3-i
    aminoglycoside_profiles = sorted(
        ((normalise_aminoglycoside_name(name), profile) for name, profile in aminoglycoside_gene_mapping.items()),
        key=lambda item: len(item[0]), reverse=True)

    for gene in detailed_resistance:
        gene_name = gene.get('gene', '').lower().replace('^', '').strip()
        drug_class = gene.get('drug_class', '').lower()
        drug_class = DRUG_CLASS_ALIASES.get(drug_class, drug_class)

        if 'aminoglycoside' in drug_class:
            gene_lookup = normalise_aminoglycoside_name(gene_name)
            resistance_profile = next((profile for name, profile in aminoglycoside_profiles
                                       if gene_lookup.startswith(name)), None)

            if resistance_profile is not None:
                if resistance_profile.get('Gentamicin'):
                    drug_status["Gentamicin"] = "Resistant"
                if resistance_profile.get('Amikacin'):
                    drug_status["Amikacin"] = "Resistant"
            elif 'aac' in gene_lookup or 'ant' in gene_lookup or 'aph' in gene_lookup:
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
            
        if gene in COLISTIN_MUTATION_GENES:
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
        "performer": [{"reference": f"Organization/{org_id}"}],
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

def add_standard_kp_components(components, typing_data):

    virulence_info = typing_data.get('virulence', {})

    for vf_name, vf_value in virulence_info.items():
        if vf_name not in ['virulence_score'] and vf_value not in ['Unknown', 'None', '']:
            components.append({
                "code": {"text": f"Virulence factor {vf_name}"},
                "valueString": vf_value
            })

def create_enhanced_annotation(gene, drug_class, typing_data, intrinsic=False):
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
    # merge_clinical_fhir.py leaves genes marked like this out of the resistance classification
    if intrinsic:
        annotation += " (intrinsic chromosomal)"
    
    mlst_info = typing_data.get('mlst', {})
    sequence_type = mlst_info.get('sequence_type', 'Unknown')
    capsule_info = typing_data.get('capsule', {})
    k_type = capsule_info.get('k_type', 'Unknown')
    
    if sequence_type not in ['Unknown', 'Undetermined']:
        annotation += f" - KP MLST: {sequence_type}"
    if k_type not in ['Unknown', 'Capsule null']:
        annotation += f" - Capsule: {k_type}"

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

def create_mlst_observation(file_sample_id, typing_data, org_id='unknown-org'):
    mlst_info = typing_data.get('mlst', {})
    sequence_type = mlst_info.get('sequence_type', 'Unknown')
    alleles = {locus: allele for locus, allele in mlst_info.get('alleles', {}).items()
               if allele not in ['Unknown', '-', '', 'NA']}

    if sequence_type in ['Unknown', 'Undetermined', 'NA']:
        if not alleles:
            return None
        # keep the allele profile: "*" = inexact allele match, which is common in error-prone assemblies
        if len(alleles) < 7 or any('*' in a or '?' in a for a in alleles.values()):
            sequence_type = 'ST not determined (inexact or missing alleles)'
        else:
            sequence_type = 'Novel ST (new allele combination)'
        st_value = {"text": sequence_type}
    else:
        st_value = {"coding": [get_mlst_st_coding(sequence_type)], "text": sequence_type}

    allele_text = ', '.join(f"{locus} {allele}" for locus, allele in alleles.items())
    div_text = f"<div xmlns=\"http://www.w3.org/1999/xhtml\">Klebsiella pneumoniae MLST {sequence_type}"
    div_text += f". Alleles: {allele_text}</div>" if allele_text else "</div>"

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
        "valueCodeableConcept": st_value,
        # merge_clinical_fhir.py reads these back as "MLST allele <locus>"
        **({"component": [{"code": {"text": f"MLST allele {locus}"}, "valueString": allele}
                          for locus, allele in alleles.items()]} if alleles else {}),
        "subject": {"reference": f"Patient/{file_sample_id}-patient"},
        "specimen": {"reference": f"Specimen/{file_sample_id}-specimen"},
        "effectiveDateTime": datetime.now(timezone.utc).isoformat(),
        "performer": [{"reference": f"Organization/{org_id}"}]
    }

def create_virulence_score_observation(file_sample_id, typing_data, org_id='unknown-org'):
    virulence_info = typing_data.get('virulence', {})
    # a score of 0 is a result (no ybt, clb or iuc), so it is reported too
    virulence_score = virulence_info.get('virulence_score', 'Unknown')

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
            "performer": [{"reference": f"Organization/{org_id}"}]
        }
    except ValueError:
        return None

def create_capsule_observation(file_sample_id, typing_data, org_id='unknown-org'):
    capsule_info = typing_data.get('capsule', {})
    k_type = capsule_info.get('k_type', 'Unknown')
    k_locus = capsule_info.get('k_locus', 'Unknown')

    if k_type in ['Unknown', 'NA'] and k_locus in ['Unknown', 'NA']:
        return None

    if k_type in ['Capsule null', 'Unknown', 'NA']:
        # Kaptive found the locus but no capsule type: "Capsule null" means an essential gene is missing or truncated
        phenotype = 'capsule null' if k_type == 'Capsule null' else 'type unknown'
        k_type = f"{k_locus} ({phenotype})"
        capsule_coding = {
            'system': 'http://kaptive.holtlab.net/capsule',
            'code': k_locus,
            'display': f'Klebsiella capsule locus {k_locus} ({phenotype})'
        }
    else:
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
        "performer": [{"reference": f"Organization/{org_id}"}]
    }

def create_o_antigen_observation(file_sample_id, typing_data, org_id='unknown-org'):
    o_locus_info = typing_data.get('o_locus', {})
    o_type = o_locus_info.get('o_type', 'Unknown')

    if o_type in ['Unknown', 'NA', '-', '']:
        return None

    div_text = f"<div xmlns=\"http://www.w3.org/1999/xhtml\">Klebsiella pneumoniae O-antigen type {o_type}</div>"

    return {
        "resourceType": "Observation",
        "id": sanitize_id(f"{file_sample_id}-o-antigen"),
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
            "coding": [{
                "system": "http://kaptive.holtlab.net/o-antigen",
                "code": o_type,
                "display": f"Klebsiella pneumoniae {o_type}"
            }],
            "text": o_type
        },
        "subject": {"reference": f"Patient/{file_sample_id}-patient"},
        "specimen": {"reference": f"Specimen/{file_sample_id}-specimen"},
        "effectiveDateTime": datetime.now(timezone.utc).isoformat(),
        "performer": [{"reference": f"Organization/{org_id}"}]
    }

def load_plasmid_data(plasmid_json):
    if not plasmid_json:
        return None
    try:
        with open(plasmid_json, 'r') as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError) as e:
        debug_print(f"Could not read plasmid data {plasmid_json}: {e}")
        return {"status": "failed", "plasmids": [], "gene_locations": {}}

def plasmid_clusters_for(location_info):
    clusters = []
    for copy in location_info.get('copies', []):
        cluster_id = copy.get('cluster_id')
        if copy.get('molecule') == 'plasmid' and cluster_id and cluster_id not in clusters:
            clusters.append(cluster_id)
    return clusters

def describe_plasmid(plasmid):
    details = [', '.join(plasmid.get('replicon_types') or []) or 'replicon not typed']
    if plasmid.get('predicted_mobility'):
        details.append(plasmid['predicted_mobility'])
    return f"{plasmid.get('cluster_id')} ({', '.join(details)})"

def describe_location(location_info, plasmids_by_cluster):
    # Never use the words "chromosomal" or "intrinsic" here: merge_clinical_fhir.py treats them as intrinsic genes.
    location = location_info.get('location', 'unresolved')
    if location == 'plasmid':
        plasmids = [describe_plasmid(plasmids_by_cluster.get(c, {'cluster_id': c})) for c in plasmid_clusters_for(location_info)]
        return f" - located on plasmid {'; '.join(plasmids)}"
    if location == 'chromosome':
        return " - located on chromosome"
    return " - genomic location unresolved"

def create_location_components(location_info, plasmids_by_cluster):
    location = location_info.get('location', 'unresolved')
    coding = LOCATION_CODINGS.get(location)
    components = [{
        "code": {
            "coding": [{
                "system": PLASMID_CODE_SYSTEM,
                "code": GENOMIC_LOCATION_CODE,
                "display": PLASMID_CODE_DISPLAYS[GENOMIC_LOCATION_CODE]
            }],
            "text": PLASMID_CODE_DISPLAYS[GENOMIC_LOCATION_CODE]
        },
        "valueCodeableConcept": {"coding": [coding], "text": location} if coding else {"text": "Unresolved"}
    }]
    for cluster_id in plasmid_clusters_for(location_info):
        replicons = ', '.join(plasmids_by_cluster.get(cluster_id, {}).get('replicon_types') or [])
        components.append({
            "code": {"text": "Plasmid"},
            "valueString": f"{cluster_id} {replicons}".strip()
        })
    return components

def text_component(label, value):
    if isinstance(value, list):
        value = ', '.join(str(v) for v in value)
    if value in [None, '']:
        return None
    return {"code": {"text": label}, "valueString": str(value)}

def plasmid_observation(obs_id, code, div_text, file_sample_id, org_id):
    return {
        "resourceType": "Observation",
        "id": obs_id,
        "meta": {
            "profile": ["http://hl7.org/fhir/StructureDefinition/Observation"],
            "tag": [{"system": "http://terminology.kemkes.go.id/sp", "code": "genomics", "display": "Genomics"}]
        },
        "text": {
            "status": "generated",
            "div": f"<div xmlns=\"http://www.w3.org/1999/xhtml\">{html.escape(div_text, quote=False)}</div>"
        },
        "status": "final",
        "category": [
            {"coding": [{"system": "http://terminology.hl7.org/CodeSystem/observation-category", "code": "laboratory", "display": "Laboratory"}]},
            {"coding": [{"system": "http://terminology.hl7.org/CodeSystem/v2-0074", "code": "GE", "display": "Genetics"}]}
        ],
        "code": {
            "coding": [{"system": PLASMID_CODE_SYSTEM, "code": code, "display": PLASMID_CODE_DISPLAYS[code]}],
            "text": PLASMID_CODE_DISPLAYS[code]
        },
        "subject": {"reference": f"Patient/{file_sample_id}-patient"},
        "specimen": {"reference": f"Specimen/{file_sample_id}-specimen"},
        "effectiveDateTime": datetime.now(timezone.utc).isoformat(),
        "performer": [{"reference": f"Organization/{org_id}"}]
    }

def load_qc_data(qc_json):
    if not qc_json:
        return None
    try:
        with open(qc_json, 'r') as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError) as e:
        debug_print(f"Could not read assembly QC data {qc_json}: {e}")
        return None

def create_sequencing_qc_observation(file_sample_id, qc_data, org_id='unknown-org'):
    status = qc_data.get('status', 'UNKNOWN')
    warnings = qc_data.get('warnings') or []

    components = [text_component("Sequencing platform", qc_data.get('platform'))]
    for label, key in [("Estimated read coverage", 'estimated_coverage'), ("Minimum read coverage", 'min_coverage')]:
        if qc_data.get(key) is not None:
            components.append({"code": {"text": label}, "valueQuantity": {"value": qc_data[key], "unit": "x"}})
    components.append(text_component("Coverage status", qc_data.get('coverage_status')))
    if qc_data.get('polishing'):
        polishing = f"{qc_data.get('polishing_tool') or 'polishing'} {qc_data['polishing']}"
        if qc_data.get('polishing_model'):
            polishing += f" ({qc_data['polishing_model']})"
        components.append(text_component("Assembly polishing", polishing))
    # one component per warning; merge_clinical_fhir.py copies them into the DiagnosticReport
    components.extend(text_component("Quality warning", warning) for warning in warnings)

    div_text = f"Sequencing QC: {status}"
    div_text += f" - {'; '.join(warnings)}. Results for this sample may be unreliable." if warnings else "."

    observation = {
        "resourceType": "Observation",
        "id": sanitize_id(f"{file_sample_id}-sequencing-qc"),
        "meta": {
            "profile": ["http://hl7.org/fhir/StructureDefinition/Observation"],
            "tag": [{"system": "http://terminology.kemkes.go.id/sp", "code": "genomics", "display": "Genomics"}]
        },
        "text": {
            "status": "generated",
            "div": f"<div xmlns=\"http://www.w3.org/1999/xhtml\">{html.escape(div_text, quote=False)}</div>"
        },
        "status": "final",
        "category": [
            {"coding": [{"system": "http://terminology.hl7.org/CodeSystem/observation-category", "code": "laboratory", "display": "Laboratory"}]}
        ],
        "code": {
            "coding": [{"system": QC_CODE_SYSTEM, "code": SEQUENCING_QC_CODE, "display": SEQUENCING_QC_DISPLAY}],
            "text": SEQUENCING_QC_DISPLAY
        },
        "subject": {"reference": f"Patient/{file_sample_id}-patient"},
        "specimen": {"reference": f"Specimen/{file_sample_id}-specimen"},
        "effectiveDateTime": datetime.now(timezone.utc).isoformat(),
        "performer": [{"reference": f"Organization/{org_id}"}],
        "valueString": status,
        "component": [c for c in components if c]
    }
    if warnings:
        observation["interpretation"] = [{"coding": [{
            "system": "http://terminology.hl7.org/CodeSystem/v3-ObservationInterpretation",
            "code": "A", "display": "Abnormal"}]}]
    return observation

def create_plasmid_observations(file_sample_id, plasmid_data, plasmid_members, org_id='unknown-org'):
    status = plasmid_data.get('status', 'failed')
    plasmids = plasmid_data.get('plasmids', []) if status != 'failed' else []
    confidence = plasmid_data.get('confidence', 'unknown')
    tools = plasmid_data.get('tools', {})
    observations = []

    for plasmid in plasmids:
        cluster_id = plasmid.get('cluster_id') or 'unclustered'
        replicons = ', '.join(plasmid.get('replicon_types') or []) or 'Replicon not typed'
        genes = plasmid.get('amr_genes') or []
        size_bp = plasmid.get('size_bp')
        mobility = plasmid.get('predicted_mobility')

        components = []
        cluster_text = cluster_id
        if plasmid.get('secondary_cluster_id'):
            cluster_text += f" (secondary {plasmid['secondary_cluster_id']})"
        components.append(text_component("Plasmid cluster", cluster_text))
        if size_bp:
            components.append({"code": {"text": "Plasmid size"}, "valueQuantity": {"value": size_bp, "unit": "bp"}})
        if plasmid.get('num_contigs'):
            components.append({"code": {"text": "Number of contigs"}, "valueInteger": plasmid['num_contigs']})
        components.append(text_component("Relaxase type", plasmid.get('relaxase_types')))
        components.append(text_component("MPF type", plasmid.get('mpf_type')))
        components.append(text_component("oriT type", plasmid.get('orit_types')))
        components.append(text_component("Predicted mobility", mobility))
        if plasmid.get('mash_nearest_neighbor'):
            neighbor = plasmid['mash_nearest_neighbor']
            if plasmid.get('mash_distance') is not None:
                neighbor += f" (Mash distance {plasmid['mash_distance']})"
            components.append(text_component("Nearest reference plasmid", neighbor))
        components.append(text_component("Predicted host range", plasmid.get('predicted_host_range')))
        components.append(text_component("Carried resistance genes", genes))
        components.append(text_component("Carried carbapenemase genes", plasmid.get('carbapenemase_genes')))
        components.append(text_component("Carried ESBL genes", plasmid.get('esbl_genes')))
        components.append(text_component("Analysis confidence", confidence))

        size_text = f"{size_bp / 1000:.1f} kb" if size_bp else "size unknown"
        div_text = (f"Plasmid {cluster_id}: {replicons}, {size_text}, {mobility or 'mobility unknown'}; "
                    f"carries {', '.join(genes) if genes else 'no resistance genes'}. Confidence: {confidence}.")

        observation = plasmid_observation(sanitize_id(f"{file_sample_id}-plasmid-{cluster_id}"), PLASMID_CODE,
                                          div_text, file_sample_id, org_id)
        observation["valueCodeableConcept"] = {"text": replicons}
        observation["component"] = [c for c in components if c]
        members = plasmid_members.get(cluster_id, [])
        if members:
            observation["hasMember"] = [{"reference": f"Observation/{obs_id}"} for obs_id in members]
        observations.append(observation)

    if status == 'failed':
        div_text = "Plasmid analysis failed; plasmid content is unknown."
    else:
        div_text = (f"Plasmid analysis (MOB-suite {tools.get('mob_suite', 'unknown')}): "
                    f"{len(plasmids)} plasmid(s) detected. Confidence: {confidence}.")

    summary = plasmid_observation(sanitize_id(f"{file_sample_id}-plasmid-summary"), PLASMID_ANALYSIS_CODE,
                                  div_text, file_sample_id, org_id)
    if status == 'failed':
        summary["dataAbsentReason"] = {
            "coding": [{"system": "http://terminology.hl7.org/CodeSystem/data-absent-reason", "code": "error", "display": "Error"}]
        }
    else:
        summary["valueInteger"] = len(plasmids)
    summary_components = [
        text_component("Analysis status", status),
        text_component("Assembly type", plasmid_data.get('assembly_type')),
        text_component("Analysis confidence", confidence),
        text_component("MOB-suite version", tools.get('mob_suite')),
        text_component("Resistance gene reference", tools.get('card_db')),
    ]
    summary["component"] = [c for c in summary_components if c]
    if observations:
        summary["hasMember"] = [{"reference": f"Observation/{obs['id']}"} for obs in observations]

    return [summary] + observations

def main():
    parser = argparse.ArgumentParser(description='Convert KP typing data to FHIR format')
    parser.add_argument('--input', required=True, help='Path to typing JSON file')
    parser.add_argument('--sample_id', required=True, help='Sample ID')
    parser.add_argument('--output', required=True, help='Path to output FHIR JSON')
    parser.add_argument('--cgmlst_file', help='Path to chewBBACA results_alleles.tsv')
    parser.add_argument('--organization_metadata', help='Path to organization metadata CSV/Excel file')
    parser.add_argument('--plasmid_json', help='Path to plasmid JSON from parse_plasmids.py')
    parser.add_argument('--qc_json', help='Path to assembly QC JSON from assembly_qc.py')
    args = parser.parse_args()

    org_data = {}
    if args.organization_metadata and os.path.exists(args.organization_metadata):
        org_data = load_organization_metadata(args.organization_metadata)
    else:
        debug_print(f"Organization metadata not provided or not found: {args.organization_metadata}")
    org_id = org_data.get('org_id', 'unknown-org')

    file_sample_id_variations = extract_sample_id_from_filename(args.input)

    plasmid_data = load_plasmid_data(args.plasmid_json)
    plasmids_by_cluster = {p.get('cluster_id'): p for p in (plasmid_data or {}).get('plasmids', [])}
    gene_locations = (plasmid_data or {}).get('gene_locations', {})
    plasmid_members = defaultdict(list)

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
                        elif gene in COLISTIN_MUTATION_GENES: drug_class = 'colistin'
                        
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

                    location_info = gene_locations.get(gene if source_type == 'acquired' else clean_gene)
                    if location_info:
                        components.extend(create_location_components(location_info, plasmids_by_cluster))

                    add_standard_kp_components(components, typing_data)

                    if source_type == 'acquired':
                        annotation_text = create_enhanced_annotation(gene, drug_class, typing_data,
                                                                     intrinsic=feature_info.get('intrinsic', False))
                    else:
                        annotation_text = f"Point mutation {mutation} in {gene} conferring resistance"
                    if location_info:
                        annotation_text += describe_location(location_info, plasmids_by_cluster)

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
                        "performer": [{"reference": f"Organization/{org_id}"}],
                        "component": components
                    }

                    bundles[file_sample_id].append({
                        "fullUrl": f"urn:uuid:{str(uuid.uuid4()).lower()}",
                        "resource": observation
                    })

                    if location_info:
                        for cluster_id in plasmid_clusters_for(location_info):
                            plasmid_members[cluster_id].append(observation_id)

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
                
                
                add_standard_kp_components(components, typing_data)

                mlst_info = typing_data.get('mlst', {})
                sequence_type = mlst_info.get('sequence_type', 'Unknown')
                capsule_info = typing_data.get('capsule', {})
                k_type = capsule_info.get('k_type', 'Unknown')
                
                annotation_text = "Klebsiella pneumoniae isolate with no antimicrobial resistance genes detected"
                if sequence_type not in ['Unknown', 'Undetermined']:
                    annotation_text += f" - MLST: {sequence_type}"
                if k_type not in ['Unknown', 'Capsule null']:
                    annotation_text += f" - Capsule: {k_type}"

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
                    "performer": [{"reference": f"Organization/{org_id}"}],
                    "component": components
                }
                
                bundles[file_sample_id].append({
                    "fullUrl": f"urn:uuid:{str(uuid.uuid4()).lower()}",
                    "resource": susceptible_observation
                })
                observation_count += 1
                
            except Exception as e:
                debug_print(f"Error creating susceptible observation: {e}")
        
        mlst_obs = create_mlst_observation(file_sample_id, typing_data, org_id)
        if mlst_obs:
            bundles[file_sample_id].append({
                "fullUrl": f"urn:uuid:{str(uuid.uuid4()).lower()}",
                "resource": mlst_obs
            })
            observation_count += 1
            
        virulence_obs = create_virulence_score_observation(file_sample_id, typing_data, org_id)
        if virulence_obs:
            bundles[file_sample_id].append({
                "fullUrl": f"urn:uuid:{str(uuid.uuid4()).lower()}",
                "resource": virulence_obs
            })
            observation_count += 1

        capsule_obs = create_capsule_observation(file_sample_id, typing_data, org_id)
        if capsule_obs:
            bundles[file_sample_id].append({
                "fullUrl": f"urn:uuid:{str(uuid.uuid4()).lower()}",
                "resource": capsule_obs
            })
            observation_count += 1

        o_antigen_obs = create_o_antigen_observation(file_sample_id, typing_data, org_id)
        if o_antigen_obs:
            bundles[file_sample_id].append({
                "fullUrl": f"urn:uuid:{str(uuid.uuid4()).lower()}",
                "resource": o_antigen_obs
            })
            observation_count += 1

        susceptible_panel = create_susceptibility_panel(file_sample_id, typing_data, org_id)
        if susceptible_panel:
             bundles[file_sample_id].append({
                "fullUrl": f"urn:uuid:{str(uuid.uuid4()).lower()}",
                "resource": susceptible_panel
            })
             observation_count += 1
        
        cgmlst_info = None
        if args.cgmlst_file:
            cgmlst_info = load_cgmlst_data(args.cgmlst_file, args.sample_id)

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
                    "performer": [{"reference": f"Organization/{org_id}"}],
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
                    if allele.isdigit():
                        comp_value = {"valueInteger": int(allele)}
                    else:
                        comp_value = {"valueString": allele}
                    cgmlst_obs['component'].append({
                        "code": {
                            "coding": [{
                                "system": "https://www.cgmlst.org/ncs",
                                "code": locus,
                                "display": locus
                            }]
                        },
                        **comp_value
                    })
                
                bundles[file_sample_id].append({
                    "fullUrl": f"urn:uuid:{str(uuid.uuid4()).lower()}",
                    "resource": cgmlst_obs
                })
                observation_count += 1
            except Exception as e:
                debug_print(f"Error creating cgMLST observation: {e}")

        if plasmid_data is not None:
            try:
                for plasmid_obs in create_plasmid_observations(file_sample_id, plasmid_data, plasmid_members, org_id):
                    bundles[file_sample_id].append({
                        "fullUrl": f"urn:uuid:{str(uuid.uuid4()).lower()}",
                        "resource": plasmid_obs
                    })
                    observation_count += 1
            except Exception as e:
                debug_print(f"Error creating plasmid observations: {e}")

        qc_data = load_qc_data(args.qc_json)
        if qc_data is not None:
            bundles[file_sample_id].append({
                "fullUrl": f"urn:uuid:{str(uuid.uuid4()).lower()}",
                "resource": create_sequencing_qc_observation(file_sample_id, qc_data, org_id)
            })
            observation_count += 1

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