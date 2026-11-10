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
    mkdir -p qc_results
    cp -L ${qc_files} qc_results/
    multiqc qc_results --force --outdir .
    """
}

process GENERATE_SAMPLE_REPORT {
    publishDir "${params.results_dir}/reports", mode: 'copy'
    
    input:
    tuple val(sample_id), path(typing_json), path(lineage_json), path(plasmid_json), path(qc_json)

    output:
    tuple val(sample_id), path("${sample_id}.summary_report.txt"), emit: sample_report

    script:
    def plasmid_arg = plasmid_json ? "--plasmid-json ${plasmid_json}" : ""
    def qc_arg = qc_json ? "--qc-json ${qc_json}" : ""
    """
    generate_sample_report.py \\
        --sample-id ${sample_id} \\
        --typing-json ${typing_json} \\
        --lineage-json ${lineage_json} \\
        --output ${sample_id}.summary_report.txt \\
        ${plasmid_arg} \\
        ${qc_arg}
    """
}

workflow GENERATE_REPORT {
    take:
    qc_files
    typing_reports
    lineage_reports
    plasmid_reports
    assembly_qc         // Nanopore samples only

    main:

    multiqc_out = MULTIQC(qc_files)

    combined_reports = typing_reports
        .join(lineage_reports)
        .join(plasmid_reports, remainder: true)
        .join(assembly_qc, remainder: true)
        .filter { row -> row.size() > 2 && row[1] != null && row[2] != null }
        .map { row -> tuple(row[0], row[1], row[2],
                            (row.size() > 3 && row[3]) ? row[3] : [],
                            (row.size() > 4 && row[4]) ? row[4] : []) }
    
    sample_reports = GENERATE_SAMPLE_REPORT(combined_reports)

    emit:
    sample_reports = sample_reports.sample_report
}
