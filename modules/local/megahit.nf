process MEGAHIT {
    tag "${meta.id}"
    label 'process_high'
    conda "bioconda::megahit=1.2.9"
    container "quay.io/biocontainers/megahit:1.2.9--haf24da9_8"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_raw_contigs.fasta"), emit: contigs
    path "versions.yml", emit: versions

    script:
    """
    megahit -1 ${reads[0]} -2 ${reads[1]} \\
        -o megahit_output \\
        --min-contig-len 500 \\
        -t ${task.cpus} \\
        --memory 0.5

    cp megahit_output/final.contigs.fa ${meta.id}_raw_contigs.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megahit: \$(megahit --version 2>&1 | sed 's/MEGAHIT v//')
    END_VERSIONS
    """

    stub:
    """
    printf '>k141_1\\nACGT\\n' > ${meta.id}_raw_contigs.fasta
    printf '"%s":\\n    megahit: stub\\n' "${task.process}" > versions.yml
    """
}
