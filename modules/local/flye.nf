process FLYE {
    tag "${meta.id}"
    label 'process_high'
    conda "bioconda::flye=2.9.6"
    container "quay.io/biocontainers/flye:2.9.6--py312h734f728_1"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_flye.fasta"), emit: assembly
    tuple val(meta), path("${meta.id}_flye_info.txt"), emit: info, optional: true
    path "versions.yml", emit: versions

    script:
    """
    flye --nano-raw ${reads} \\
         --out-dir flye_output \\
         --genome-size ${params.genome_size} \\
         --threads ${task.cpus} \\
         --iterations 2

    if [ ! -f flye_output/assembly.fasta ]; then
        echo "Error: Flye assembly failed." >&2
        exit 1
    fi
    cp flye_output/assembly.fasta ${meta.id}_flye.fasta
    if [ -f flye_output/assembly_info.txt ]; then
        cp flye_output/assembly_info.txt ${meta.id}_flye_info.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$(flye --version)
    END_VERSIONS
    """

    stub:
    """
    printf '>contig_1\\nACGT\\n' > ${meta.id}_flye.fasta
    printf '"%s":\\n    flye: stub\\n' "${task.process}" > versions.yml
    """
}
