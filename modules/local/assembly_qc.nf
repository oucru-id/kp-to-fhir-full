process ASSEMBLY_QC {
    tag "${meta.id}"
    label 'process_single'
    conda "conda-forge::python=3.12"
    container "quay.io/biocontainers/python:3.12.12"

    input:
    tuple val(meta), val(coverage), val(polish_status)

    output:
    tuple val(meta), path("${meta.id}_assembly_qc.json"), emit: json
    path "versions.yml", emit: versions

    script:
    """
    assembly_qc.py \\
        --sample_id '${meta.id}' \\
        --platform ont \\
        --coverage ${coverage} \\
        --min_coverage ${params.min_nanopore_coverage} \\
        --polishing ${polish_status} \\
        --polishing_tool medaka \\
        --polishing_model ${params.medaka_model} \\
        --output ${meta.id}_assembly_qc.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo '{"sample_id": "${meta.id}", "status": "PASS", "warnings": []}' > ${meta.id}_assembly_qc.json
    printf '"%s":\\n    python: stub\\n' "${task.process}" > versions.yml
    """
}
