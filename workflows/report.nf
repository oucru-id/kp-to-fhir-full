#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process MULTIQC {
    publishDir "${params.results_dir}/qc", mode: 'copy'

    input:
    path(qc_files)

    output:
    path "multiqc_report.html", emit: report

    script:
    """
    mkdir -p fastqc_results
    cp -L ${qc_files} fastqc_results/
    multiqc fastqc_results --force --outdir .
    """
}

process GENERATE_SAMPLE_REPORT {
    publishDir "${params.results_dir}/reports", mode: 'copy'
    
    input:
    tuple val(sample_id), path(typing_json), path(lineage_json)
    
    output:
    tuple val(sample_id), path("${sample_id}.summary_report.txt"), emit: sample_report
    
    script:
    """
    python3 $baseDir/scripts/generate_sample_report.py \\
        --sample-id ${sample_id} \\
        --typing-json ${typing_json} \\
        --lineage-json ${lineage_json} \\
        --output ${sample_id}.summary_report.txt
    """
}

workflow GENERATE_REPORT {
    take:
    qc_files
    typing_reports
    lineage_reports

    main:

    multiqc_out = MULTIQC(qc_files)
    
    combined_reports = typing_reports
        .join(lineage_reports)
    
    sample_reports = GENERATE_SAMPLE_REPORT(combined_reports)

    emit:
    sample_reports = sample_reports.sample_report
}
