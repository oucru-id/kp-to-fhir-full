process BWA_MEM {
    tag "${meta.id}"
    label 'process_medium'
    conda "bioconda::bwa=0.7.17"
    container "quay.io/biocontainers/bwa:0.7.17--he4a0461_11"

    input:
    tuple val(meta), path(raw_assembly), path(reads)

    output:
    tuple val(meta), path("${meta.id}.sam.gz"), emit: sam
    path "versions.yml", emit: versions

    script:
    """
    bwa index ${raw_assembly}

    bwa mem -t ${task.cpus} ${raw_assembly} ${reads[0]} ${reads[1]} | gzip -1 > ${meta.id}.sam.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | sed -n 's/^Version: //p')
    END_VERSIONS
    """

    stub:
    """
    echo | gzip > ${meta.id}.sam.gz
    printf '"%s":\\n    bwa: stub\\n' "${task.process}" > versions.yml
    """
}
