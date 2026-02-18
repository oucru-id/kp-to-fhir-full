#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process CREATE_FHIR {
    publishDir "${params.results_dir}/fhir", mode: 'copy'

    input:
    tuple val(sample_id), path(typing_json)
    path(lineage_files)
    path(cgmlst_matrix)

    output:
    tuple val(sample_id), path("${sample_id}.fhir.json"), emit: fhir_output
    path "versions.yml", emit: versions

    script:
    def cgmlst_arg = cgmlst_matrix.name != 'NO_CGMLST' ? "--cgmlst_file ${cgmlst_matrix}" : ""
    """
    mkdir -p lineage_data
    
    for file in ${lineage_files}; do
        if [ -f "\$file" ]; then
            cp "\$file" lineage_data/
        fi
    done

    python3 $baseDir/scripts/annotated_to_fhir.py \\
        --input ${typing_json} \\
        --sample_id ${sample_id} \\
        --output ${sample_id}.fhir.json \\
        --lineage_dir lineage_data/ \\
        ${cgmlst_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}

workflow FHIR {
    take:
    typing_ch
    lineage_ch
    cgmlst_ch

    main:

    lineage_files = lineage_ch
        .map { sample_id, lineage_file -> lineage_file }
        .collect()

    cgmlst_file = cgmlst_ch.ifEmpty(file("NO_CGMLST"))

    CREATE_FHIR(typing_ch, lineage_files, cgmlst_file)

    emit:
    fhir_output = CREATE_FHIR.out.fhir_output
    versions    = CREATE_FHIR.out.versions
}