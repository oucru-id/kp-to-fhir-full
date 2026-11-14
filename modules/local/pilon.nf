process PILON {
    tag "${meta.id}"
    label 'process_medium'
    conda "bioconda::pilon=1.24"
    container "quay.io/biocontainers/pilon:1.24--hdfd78af_0"

    input:
    tuple val(meta), path(raw_assembly), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}_contigs.fasta"), emit: assembly
    path "versions.yml", emit: versions

    script:
    def heap = Math.max(1024, (task.memory.toMega() * 0.8).intValue())
    """
    pilon -Xmx${heap}m --genome ${raw_assembly} --frags ${bam} --output polished --changes

    cp polished.fasta ${meta.id}_contigs.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pilon: \$(pilon --version 2>&1 | grep -oE 'version [0-9.]+' | head -1 | sed 's/version //')
    END_VERSIONS
    """

    stub:
    """
    printf '>NODE_1_length_4_cov_1.0_pilon\\nACGT\\n' > ${meta.id}_contigs.fasta
    printf '"%s":\\n    pilon: stub\\n' "${task.process}" > versions.yml
    """
}
