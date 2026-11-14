process PREPARE_ASSEMBLY {
    tag "${meta.id}"
    label 'process_single'
    conda "conda-forge::python=3.12"
    container "quay.io/biocontainers/python:3.12.12"

    input:
    tuple val(meta), path(fasta, stageAs: 'input/*')

    output:
    tuple val(meta), path("${meta.id}_assembly.fasta"), emit: assembly
    path "versions.yml", emit: versions

    script:
    """
    gzip -cdf ${fasta} > ${meta.id}_assembly.fasta
    if ! head -c 1 ${meta.id}_assembly.fasta | grep -q '>'; then
        echo "ERROR: ${fasta.name} for sample ${meta.id} is not a FASTA file" >&2
        exit 1
    fi
    echo "${meta.id}: \$(grep -c '>' ${meta.id}_assembly.fasta) contig(s)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: \$(gzip --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    printf '>contig_1\\nACGT\\n' > ${meta.id}_assembly.fasta
    printf '"%s":\\n    gzip: stub\\n' "${task.process}" > versions.yml
    """
}
