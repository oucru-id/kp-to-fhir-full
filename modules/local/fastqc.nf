process FASTQC {
    tag "${meta.id}"
    label 'process_low'
    conda "bioconda::fastqc=0.12.1"
    container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_fastqc.zip"), emit: zip
    tuple val(meta), path("*_fastqc.html"), emit: html
    path "versions.yml", emit: versions

    script:
    """
    fastqc ${reads.join(' ')} --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/FastQC v//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_1_fastqc.zip ${meta.id}_1_fastqc.html
    printf '"%s":\\n    fastqc: stub\\n' "${task.process}" > versions.yml
    """
}
