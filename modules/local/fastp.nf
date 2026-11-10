process FASTP {
    tag "${meta.id}"
    label 'process_medium'
    conda "bioconda::fastp=0.23.4"
    container "quay.io/biocontainers/fastp:0.23.4--h125f33a_5"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_R{1,2}_trimmed.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}_fastp.json"), emit: json
    tuple val(meta), path("${meta.id}_fastp.html"), emit: html
    path "versions.yml", emit: versions

    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} \\
          -o ${meta.id}_R1_trimmed.fastq.gz \\
          -O ${meta.id}_R2_trimmed.fastq.gz \\
          --json ${meta.id}_fastp.json \\
          --html ${meta.id}_fastp.html \\
          --thread ${task.cpus} \\
          --qualified_quality_phred 20 \\
          --unqualified_percent_limit 30 \\
          --length_required 50 \\
          --detect_adapter_for_pe \\
          --cut_front \\
          --cut_tail \\
          --cut_window_size 4 \\
          --cut_mean_quality 20

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed 's/fastp //')
    END_VERSIONS
    """

    stub:
    """
    echo | gzip > ${meta.id}_R1_trimmed.fastq.gz
    echo | gzip > ${meta.id}_R2_trimmed.fastq.gz
    echo '{}' > ${meta.id}_fastp.json
    touch ${meta.id}_fastp.html
    printf '"%s":\\n    fastp: stub\\n' "${task.process}" > versions.yml
    """
}
