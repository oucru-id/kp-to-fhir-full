process PARSE_KLEBORATE {
    tag "${meta.id}"
    label 'process_single'
    conda "conda-forge::python=3.12"
    container "quay.io/biocontainers/python:3.12.12"

    input:
    tuple val(meta), path(kleborate_txt)

    output:
    tuple val(meta), path("${meta.id}_typing.json"), emit: typing
    tuple val(meta), path("${meta.id}_lineage.json"), emit: lineage
    path "versions.yml", emit: versions

    script:
    """
    kp_parse_kleborate.py --kleborate ${kleborate_txt} --sample_id '${meta.id}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo '{"sample_id": "${meta.id}", "mlst": {"sequence_type": "ST11"}, "resistance": {}}' > ${meta.id}_typing.json
    echo '{"sample_id": "${meta.id}", "lineage": {}}' > ${meta.id}_lineage.json
    printf '"%s":\\n    python: stub\\n' "${task.process}" > versions.yml
    """
}
