#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process kleborate {
    publishDir "${params.results_dir}/typing", mode: 'copy'
    
    input:
    tuple val(sample_id), path(assembly)
    
    output:
    tuple val(sample_id), path("${sample_id}_kleborate_full.txt"), emit: kleborate_report
    
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
    kp_parse_kleborate.py --kleborate ${kleborate_txt} --sample_id ${sample_id}
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
}