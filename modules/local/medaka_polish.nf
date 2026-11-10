process MEDAKA_POLISH {
    tag "${meta.id}"
    label 'process_high'
    conda "bioconda::medaka=2.1.0"
    container "quay.io/biocontainers/medaka:2.1.0--py38ha0c3a46_0"

    input:
    tuple val(meta), path(raw_assembly), path(reads)

    output:
    tuple val(meta), path("${meta.id}_assembly.fasta"), emit: assembly
    tuple val(meta), env('POLISH_STATUS'), emit: polish_status
    path "versions.yml", emit: versions

    script:
    """
    mv ${raw_assembly} raw_input.fasta
    mkdir -p medaka_work

    set +e
    medaka_consensus -i ${reads} -d raw_input.fasta -o medaka_work -t ${task.cpus} -m ${params.medaka_model} -b ${params.medaka_batch_size}
    MEDAKA_EXIT_CODE=\$?
    if [ \$MEDAKA_EXIT_CODE -ne 0 ] || [ ! -s medaka_work/consensus.fasta ]; then
        echo "Medaka failed (exit code \$MEDAKA_EXIT_CODE); retrying on CPU." >&2
        rm -rf medaka_work
        CUDA_VISIBLE_DEVICES="" medaka_consensus -i ${reads} -d raw_input.fasta -o medaka_work -t ${task.cpus} -m ${params.medaka_model} -b ${params.medaka_batch_size}
        MEDAKA_EXIT_CODE=\$?
    fi
    set -e

    if [ \$MEDAKA_EXIT_CODE -eq 0 ] && [ -s medaka_work/consensus.fasta ]; then
        cp medaka_work/consensus.fasta ${meta.id}_assembly.fasta
        POLISH_STATUS="polished"
    else
        echo "Medaka polishing failed (exit code \$MEDAKA_EXIT_CODE); keeping the Flye assembly." >&2
        cp raw_input.fasta ${meta.id}_assembly.fasta
        POLISH_STATUS="failed"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$(medaka --version 2>&1 | sed 's/medaka //')
        medaka_model: ${params.medaka_model}
    END_VERSIONS
    """

    stub:
    """
    printf '>contig_1\\nACGT\\n' > ${meta.id}_assembly.fasta
    POLISH_STATUS="polished"
    printf '"%s":\\n    medaka: stub\\n' "${task.process}" > versions.yml
    """
}
