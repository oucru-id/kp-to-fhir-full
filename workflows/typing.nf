#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process kleborate {
    publishDir "${params.results_dir}/typing", mode: 'copy'
    
    input:
    tuple val(sample_id), path(assembly)
    
    output:
    tuple val(sample_id), path("${sample_id}_kleborate_full.txt"), emit: kleborate_report
    path "versions.yml", emit: versions
    
    script:
    """
    mkdir -p kleborate_output
    
    kleborate -a ${assembly} \\
              --modules klebsiella_pneumo_complex__mlst,klebsiella_pneumo_complex__amr,klebsiella_pneumo_complex__resistance_score,klebsiella_pneumo_complex__resistance_gene_count,klebsiella_pneumo_complex__resistance_class_count,klebsiella__ybst,klebsiella__cbst,klebsiella__abst,klebsiella__smst,klebsiella__rmst,klebsiella__rmpa2,klebsiella_pneumo_complex__virulence_score,klebsiella_pneumo_complex__kaptive,klebsiella_pneumo_complex__wzi \\
              -o kleborate_output
    
    if [ -f "kleborate_output/results.txt" ]; then
        cp "kleborate_output/results.txt" ${sample_id}_kleborate_full.txt
    elif [ -f "kleborate_output/kleborate_results.txt" ]; then
        cp "kleborate_output/kleborate_results.txt" ${sample_id}_kleborate_full.txt
    else
        output_file=\$(find kleborate_output -name "*.txt" | head -1)
        if [ -f "\$output_file" ]; then
            cp "\$output_file" ${sample_id}_kleborate_full.txt
        else
            ls -la kleborate_output/
            exit 1
        fi
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kleborate: \$(kleborate --version 2>&1 | sed 's/Kleborate v//')
    END_VERSIONS
    """
}

process parse_results {
    publishDir "${params.results_dir}/typing", mode: 'copy'
    
    input:
    tuple val(sample_id), path(kleborate_txt)
    
    output:
    tuple val(sample_id), path("${sample_id}_typing.json"), emit: typing_report
    tuple val(sample_id), path("${sample_id}_lineage.json"), emit: lineage_report
    
    script:
    """
    #!/usr/bin/env python3
    import json
    import csv
    import os
    from collections import defaultdict
    
    sequence_type = "Unknown"
    clonal_complex = "Unknown"
    mlst_alleles = {
        'gapA': 'Unknown', 'infB': 'Unknown', 'mdh': 'Unknown',
        'pgi': 'Unknown', 'phoE': 'Unknown', 'rpoB': 'Unknown', 'tonB': 'Unknown'
    }
    mlst_confidence = "Not determined"
    
    resistance_genes = []
    resistance_classes = set()
    esbl_genes = []
    carbapenemases = []
    other_resistance = []
    
    k_locus = "Unknown"
    k_type = "Unknown" 
    o_locus = "Unknown"
    o_type = "Unknown"
    virulence_genes = {}
    wzi = "Unknown"
    
    try:
        with open('${kleborate_txt}', 'r') as f:
            content = f.read()
            
        with open('${kleborate_txt}', 'r') as f:
            reader = csv.DictReader(f, delimiter='\\t')
            columns = list(reader.fieldnames) if reader.fieldnames else []
            for i, col in enumerate(columns, 1):
                print(f"  {i:2d}. {col}")
            
            if 'strain' in columns:
                print("standard Kleborate output format")
                
                for row_num, row in enumerate(reader):
                    print(f"\\nProcessing row {row_num + 1}:")
                    
                    st_columns = ['klebsiella_pneumo_complex__mlst__ST', 'ST', 'sequence_type', 'mlst_st']
                    for col in st_columns:
                        if col in row and row[col] not in ['-', '', 'Unknown', '0', 'NF', 'NA']:
                            sequence_type = row[col]
                            mlst_confidence = "High"
                            break
                    
                    found_alleles = 0
                    for allele in mlst_alleles.keys():
                        allele_columns = [f'klebsiella_pneumo_complex__mlst__{allele}', allele, allele.upper(), f'mlst_{allele}']
                        for col in allele_columns:
                            if col in row and row[col] not in ['-', '', 'Unknown', 'NF', '0', 'NA']:
                                mlst_alleles[allele] = row[col]
                                found_alleles += 1
                                break
                    
                    if found_alleles > 0:
                        print(f"Found {found_alleles}/7 MLST alleles")
                    
                    k_locus_columns = ['klebsiella_pneumo_complex__kaptive__K_locus', 'K_locus', 'Best_match_locus']
                    for col in k_locus_columns:
                        if col in row and row[col] not in ['-', '', 'Unknown']:
                            k_locus = row[col]
                            break
                    
                    k_type_columns = ['klebsiella_pneumo_complex__kaptive__K_type', 'K_type', 'Best_match_type']
                    for col in k_type_columns:
                        if col in row and row[col] not in ['-', '', 'Unknown']:
                            k_type = row[col]
                            break
                    
                    o_locus_columns = ['klebsiella_pneumo_complex__kaptive__O_locus', 'klebsiella_pneumo_complex__kaptive__O_type', 'O_locus', 'O_type']
                    for col in o_locus_columns:
                        if col in row and row[col] not in ['-', '', 'Unknown', 'NA']:
                            o_locus = row[col]
                            o_type = row[col]  # For O, locus and type are often the same
                            break
                    
                    wzi_columns = ['klebsiella_pneumo_complex__wzi__wzi', 'wzi']
                    for col in wzi_columns:
                        if col in row and row[col] not in ['-', '', 'Unknown']:
                            wzi = row[col]
                            break
                    
                    virulence_mapping = {
                        'yersiniabactin': ['klebsiella__ybst__Yersiniabactin', 'Yersiniabactin', 'ybt', 'YbST'],
                        'colibactin': ['klebsiella__cbst__Colibactin', 'Colibactin', 'clb', 'CbST'], 
                        'aerobactin': ['klebsiella__abst__Aerobactin', 'Aerobactin', 'iuc', 'AbST'],
                        'salmochelin': ['klebsiella__smst__Salmochelin', 'Salmochelin', 'iro', 'SmST'],
                        'rmpadc': ['klebsiella__rmst__RmpADC', 'RmpADC', 'rmpADC', 'rmpA_rmpC', 'RmpA'],
                        'rmpa2': ['klebsiella__rmpa2__rmpA2', 'RmpA2', 'rmpA2']
                    }
                    
                    for vf_name, possible_cols in virulence_mapping.items():
                        for col in possible_cols:
                            if col in row and row[col] not in ['-', '', '0', 'NA']:
                                virulence_genes[vf_name] = row[col]
                                print(f"Found virulence factor {vf_name}: {row[col]} (from column: {col})")
                                break
                    
                    virulence_score_cols = ['klebsiella_pneumo_complex__virulence_score__virulence_score', 'virulence_score']
                    for col in virulence_score_cols:
                        if col in row and row[col] not in ['-', '', 'NA']:
                            virulence_score = row[col]
                            break
                    
                    bla_columns = [
                        'klebsiella_pneumo_complex__amr__Bla_Carb_acquired',
                        'klebsiella_pneumo_complex__amr__Bla_ESBL_acquired', 
                        'klebsiella_pneumo_complex__amr__Bla_ESBL_inhR_acquired',
                        'klebsiella_pneumo_complex__amr__Bla_acquired', 
                        'klebsiella_pneumo_complex__amr__Bla_chr',
                        'Bla_Carb_acquired', 'Bla_ESBL_acquired', 'Bla_ESBL_inhR_acquired', 'Bla_acquired', 'Bla_chr'
                    ]
                    
                    for col in bla_columns:
                        if col in row and row[col] not in ['-', '', '0', 'NA']:
                            genes = [g.strip() for g in row[col].split(';') if g.strip()]
                            for gene in genes:
                                resistance_genes.append({
                                    'gene': gene,
                                    'drug_class': 'beta-lactam',
                                    'identity': '100%',
                                    'coverage': '100%'
                                })
                                resistance_classes.add('beta-lactam')
                                
                                if 'ESBL' in col:
                                    esbl_genes.append(gene)
                                elif 'Carb' in col:
                                    carbapenemases.append(gene)
                                else:
                                    other_resistance.append(gene)
                                                        
                    other_resistance_cols = [
                        'klebsiella_pneumo_complex__amr__AGly_acquired', 'klebsiella_pneumo_complex__amr__Col_acquired', 
                        'klebsiella_pneumo_complex__amr__Fcyn_acquired', 'klebsiella_pneumo_complex__amr__Flq_acquired', 
                        'klebsiella_pneumo_complex__amr__Gly_acquired', 'klebsiella_pneumo_complex__amr__MLS_acquired', 
                        'klebsiella_pneumo_complex__amr__Phe_acquired', 'klebsiella_pneumo_complex__amr__Rif_acquired', 
                        'klebsiella_pneumo_complex__amr__Sul_acquired', 'klebsiella_pneumo_complex__amr__Tet_acquired', 
                        'klebsiella_pneumo_complex__amr__Tmt_acquired', 'klebsiella_pneumo_complex__amr__Tgc_acquired',
                        'AGly_acquired', 'Col_acquired', 'Fcyn_acquired', 'Flq_acquired', 'Gly_acquired', 
                        'MLS_acquired', 'Phe_acquired', 'Rif_acquired', 'Sul_acquired', 'Tet_acquired', 'Tgc_acquired'
                    ]
                    
                    for col in other_resistance_cols:
                        if col in row and row[col] not in ['-', '', '0', 'NA']:
                            drug_class = col.split('__')[-1].split('_')[0].lower() 
                            genes = [g.strip() for g in row[col].split(';') if g.strip()]
                            for gene in genes:
                                resistance_genes.append({
                                    'gene': gene,
                                    'drug_class': drug_class,
                                    'identity': '100%',
                                    'coverage': '100%'
                                })
                                resistance_classes.add(drug_class)
                                other_resistance.append(gene)
                    
                    resistance_score_cols = ['klebsiella_pneumo_complex__resistance_score__resistance_score', 'resistance_score']
                    for col in resistance_score_cols:
                        if col in row and row[col] not in ['-', '', 'NA']:
                            resistance_score = row[col]
                            break
                    
                    point_mutations = []
                    mutation_columns = {
                        'GyrA': ['klebsiella_pneumo_complex__amr__GyrA_mutations', 'GyrA_mutations'],
                        'ParC': ['klebsiella_pneumo_complex__amr__ParC_mutations', 'ParC_mutations'],
                        'OmpK35': ['klebsiella_pneumo_complex__amr__OmpK35_mutations', 'OmpK35_mutations'],
                        'OmpK36': ['klebsiella_pneumo_complex__amr__OmpK36_mutations', 'OmpK36_mutations'],
                        'Colistin': ['klebsiella_pneumo_complex__amr__Col_mutations', 'Col_mutations']
                    }

                    for gene_target, cols in mutation_columns.items():
                        for col in cols:
                            if col in row and row[col] not in ['-', '', '0', 'NA', '']:
                                # Kleborate separates mutations with ;
                                muts = [m.strip() for m in row[col].split(';') if m.strip()]
                                for mut in muts:
                                    # Clean up mutation string (sometimes has gene name prefix)
                                    clean_mut = mut
                                    
                                    point_mutations.append({
                                        'gene': gene_target,
                                        'mutation': clean_mut,
                                        'type': 'point_mutation'
                                    })
                                    print(f"Found point mutation: {gene_target} - {clean_mut}")
                                break

                    print(f"\\nExtracted summary:")
                    print(f"  ST: {sequence_type}")
                    print(f"  K type: {k_type}")
                    print(f"  K locus: {k_locus}")
                    print(f"  wzi: {wzi}")
                    print(f"  O type: {o_type}")
                    print(f"  Resistance genes: {len(resistance_genes)}")
                    print(f"  ESBL genes: {esbl_genes}")
                    print(f"  Carbapenemases: {carbapenemases}")
                    print(f"  Virulence factors: {list(virulence_genes.keys())}")
                    break
                    
            elif 'Input_file_name' in columns:
                print("MLST module did not run properly")
                
                for row_num, row in enumerate(reader):
                    gene_symbol = row.get('Gene_symbol', '')
                    drug_class = row.get('Drug_class', '')
                    
                    if gene_symbol and gene_symbol not in ['-', '']:
                        resistance_genes.append({
                            'gene': gene_symbol,
                            'drug_class': drug_class if drug_class not in ['-', ''] else 'unknown',
                            'identity': row.get('Sequence_identity', '100%'),
                            'coverage': row.get('Coverage', '100%')
                        })
                        
                        if drug_class not in ['-', '']:
                            resistance_classes.add(drug_class)
                        
                        if any(esbl in gene_symbol.upper() for esbl in ['ESBL', 'CTX', 'SHV', 'TEM']):
                            esbl_genes.append(gene_symbol)
                        elif any(carb in gene_symbol.upper() for carb in ['KPC', 'NDM', 'VIM', 'IMP', 'OXA']):
                            carbapenemases.append(gene_symbol)
                        else:
                            other_resistance.append(gene_symbol)
            
            else:
                print("Unknown Kleborate output")
                
    except Exception as e:
        import traceback
        traceback.print_exc()
    
    if sequence_type != "Unknown" and sequence_type.replace('ST', '').replace('-', '').split('L')[0].isdigit():
        st_number = sequence_type.replace('ST', '').replace('-', '').split('L')[0]  # Handle variants like ST258-1LV
        st_to_cc = {
            '11': 'CC11', '258': 'CC258', '147': 'CC147', '15': 'CC15',
            '16': 'CC16', '23': 'CC23', '37': 'CC37', '101': 'CC101',
            '14': 'CC14', '17': 'CC17', '20': 'CC20', '25': 'CC25',
            '268': 'CC268', '392': 'CC392', '307': 'CC307', '67': 'CC67',
            '91': 'CC91', '3067': 'CC3067'
        }
        clonal_complex = st_to_cc.get(st_number, f"CC{st_number}")
    
    num_resistance_genes = len(resistance_genes)
    num_resistance_classes = len(resistance_classes)
    
    resistance_score = 0
    if carbapenemases:
        resistance_score = 3
    elif esbl_genes:
        resistance_score = 2
    elif num_resistance_classes > 2:
        resistance_score = 1
    
    typing_result = {
        'sample_id': '${sample_id}',
        'assembly': '${sample_id}',
        'mlst': {
            'sequence_type': sequence_type,
            'clonal_complex': clonal_complex,
            'confidence': mlst_confidence,
            'alleles': mlst_alleles
        },
        'capsule': {
            'k_locus': k_locus,
            'k_type': k_type,
            'k_confidence': "High" if k_type != "Unknown" else "Low",
            'k_problems': 'None',
            'wzi': wzi
        },
        'virulence': {
            'yersiniabactin': virulence_genes.get('yersiniabactin', 'Unknown'),
            'colibactin': virulence_genes.get('colibactin', 'Unknown'), 
            'aerobactin': virulence_genes.get('aerobactin', 'Unknown'),
            'salmochelin': virulence_genes.get('salmochelin', 'Unknown'),
            'hypermucoidy': virulence_genes.get('rmpadc', 'Unknown'),
            'rmpa2': virulence_genes.get('rmpa2', 'Unknown'),
            'virulence_score': str(len([v for v in virulence_genes.values() if v != 'Unknown']))
        },
        'resistance': {
            'acquired_resistance': ', '.join(other_resistance) if other_resistance else 'None',
            'esbl': ', '.join(esbl_genes) if esbl_genes else 'None', 
            'carbapenemase': ', '.join(carbapenemases) if carbapenemases else 'None',
            'other_resistance': ', '.join([r['gene'] for r in resistance_genes]),
            'resistance_score': str(resistance_score),
            'num_resistance_classes': str(num_resistance_classes),
            'num_resistance_genes': str(num_resistance_genes),
            'detailed_resistance': resistance_genes,
            'point_mutations': point_mutations  # Add this line
        },
        'o_locus': {
            'o_locus': o_locus,
            'o_type': o_type,
            'o_confidence': "High" if o_locus != "Unknown" else "Low"
        },
        'analysis_type': 'K. pneumoniae complete Kleborate analysis'
    }
    
    lineage_result = {
        'sample_id': '${sample_id}',
        'lineage': {
            'sequence_type': sequence_type,
            'clonal_complex': clonal_complex,
            'capsule_type': k_type,
            'o_type': o_type,
            'wzi_type': wzi,
            'confidence': f'MLST: {mlst_confidence}',
            'lineage_designation': f'{sequence_type}' if sequence_type != 'Unknown' else 'Unknown ST',
            'high_risk_clone': 'Yes' if sequence_type.replace('ST', '').split('-')[0] in ['11', '258', '147', '15'] else 'Unknown' if sequence_type == 'Unknown' else 'No',
            'subspecies': 'rhinoscleromatis' if sequence_type.replace('ST', '').split('-')[0] in ['67', '68', '69'] else 'ozaenae' if sequence_type.replace('ST', '').split('-')[0] in ['90', '91', '92', '93', '95', '96', '97'] else 'pneumoniae',
            'resistance_profile': {
                'esbl_producer': 'Yes' if esbl_genes else 'No',
                'carbapenemase_producer': 'Yes' if carbapenemases else 'No', 
                'multidrug_resistant': 'Yes' if num_resistance_classes >= 3 else 'No',
                'risk_level': 'Critical' if carbapenemases else 'High' if esbl_genes else 'Moderate' if num_resistance_classes > 1 else 'Low'
            }
        },
        'analysis_type': 'K. pneumoniae Kleborate complete analysis'
    }
    
    with open('${sample_id}_typing.json', 'w') as out:
        json.dump(typing_result, out, indent=2)
        
    with open('${sample_id}_lineage.json', 'w') as out:
        json.dump(lineage_result, out, indent=2)

    """
}

workflow TYPING_ANALYSIS {
    take:
    assemblies
    
    main:
    kleborate_out = kleborate(assemblies)
    typing_results = parse_results(kleborate_out.kleborate_report)
    
    emit:
    typing_report = typing_results.typing_report
    lineage_report = typing_results.lineage_report
    versions = kleborate_out.versions
}