#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process CREATE_FHIR {
    publishDir "${params.results_dir}/fhir", mode: 'copy'

    input:
    tuple val(sample_id), path(typing_json), path(plasmid_json), path(qc_json)
    path(cgmlst_matrix)
    path(org_metadata)

    output:
    tuple val(sample_id), path("${sample_id}.fhir.json"), emit: fhir_output

    script:
    def cgmlst_arg = cgmlst_matrix.name != 'NO_CGMLST' ? "--cgmlst_file ${cgmlst_matrix}" : ""
    def plasmid_arg = plasmid_json ? "--plasmid_json ${plasmid_json}" : ""
    def qc_arg = qc_json ? "--qc_json ${qc_json}" : ""
    """
    annotated_to_fhir.py \\
        --input ${typing_json} \\
        --sample_id ${sample_id} \\
        --output ${sample_id}.fhir.json \\
        --organization_metadata ${org_metadata} \\
        ${cgmlst_arg} \\
        ${plasmid_arg} \\
        ${qc_arg}
    """
}

workflow FHIR {
    take:
    typing_ch
    plasmid_ch
    assembly_qc_ch      // Nanopore samples only
    cgmlst_ch
    org_metadata_ch

    main:

    cgmlst_file = cgmlst_ch.ifEmpty(file("NO_CGMLST"))

    typing_with_extras = typing_ch
        .join(plasmid_ch, remainder: true)
        .join(assembly_qc_ch, remainder: true)
        .filter { row -> row.size() > 1 && row[1] != null }
        .map { row -> tuple(row[0], row[1],
                            (row.size() > 2 && row[2]) ? row[2] : [],
                            (row.size() > 3 && row[3]) ? row[3] : []) }

    CREATE_FHIR(typing_with_extras, cgmlst_file, org_metadata_ch)

    emit:
    fhir_output = CREATE_FHIR.out.fhir_output
}