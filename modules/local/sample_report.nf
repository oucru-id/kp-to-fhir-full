process SAMPLE_REPORT {
    tag "${meta.id}"
    label 'process_single'
    conda "conda-forge::python=3.12"
    container "quay.io/biocontainers/python:3.12.12"

    input:
    tuple val(meta), path(typing_json), path(lineage_json), path(plasmid_json, stageAs: 'plasmid/*'), path(qc_json, stageAs: 'qc/*')

    output:
    tuple val(meta), path("${meta.id}.summary_report.txt"), emit: report
    path "versions.yml", emit: versions

    script:
    def args = []
    if (!plasmid_json.toString().endsWith('NO_FILE')) args << "--plasmid-json ${plasmid_json}"
    if (!qc_json.toString().endsWith('NO_FILE')) args << "--qc-json ${qc_json}"
    """
    generate_sample_report.py \\
        --sample-id '${meta.id}' \\
        --typing-json ${typing_json} \\
        --lineage-json ${lineage_json} \\
        --output ${meta.id}.summary_report.txt \\
        ${args.join(' ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.summary_report.txt
    printf '"%s":\\n    python: stub\\n' "${task.process}" > versions.yml
    """
}
