process SAMTOOLS_SORT {
    tag "${meta.id}"
    label 'process_medium'
    conda "bioconda::samtools=1.21"
    container "quay.io/biocontainers/samtools:1.21--h96c455f_1"

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("aligned.bam"), path("aligned.bam.bai"), emit: bam
    path "versions.yml", emit: versions

    script:
    """
    samtools sort -@ ${task.cpus} -o aligned.bam ${sam}
    samtools index aligned.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch aligned.bam aligned.bam.bai
    printf '"%s":\\n    samtools: stub\\n' "${task.process}" > versions.yml
    """
}
