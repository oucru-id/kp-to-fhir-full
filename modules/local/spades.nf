process SPADES {
    tag "${meta.id}"
    label 'process_high'
    conda "bioconda::spades=4.3.0"
    container "quay.io/biocontainers/spades:4.3.0--hde4eca7_2"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_raw_contigs.fasta"), emit: contigs
    path "versions.yml", emit: versions

    script:
    """
    spades.py --isolate \\
        -1 ${reads[0]} -2 ${reads[1]} \\
        -o spades_output \\
        -t ${task.cpus} \\
        -m ${task.memory.toGiga()}

    awk '/^>/ { split(\$0, f, "_"); keep = f[4] >= 500 } keep' spades_output/contigs.fasta > ${meta.id}_raw_contigs.fasta
    rm -rf spades_output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$(spades.py --version 2>&1 | sed 's/^.*assembler v//')
    END_VERSIONS
    """

    stub:
    """
    printf '>NODE_1_length_4_cov_1.0\\nACGT\\n' > ${meta.id}_raw_contigs.fasta
    printf '"%s":\\n    spades: stub\\n' "${task.process}" > versions.yml
    """
}
