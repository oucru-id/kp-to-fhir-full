process NANOPLOT {
    tag "${meta.id}"
    label 'process_medium'
    conda "bioconda::nanoplot=1.48.0"
    container "quay.io/biocontainers/nanoplot:1.48.0--pyhdfd78af_1"

    input:
    tuple val(meta), path(reads)
    val genome_bp

    output:
    tuple val(meta), path("${meta.id}_NanoStats.txt"), path("${meta.id}_coverage_mqc.tsv"), emit: qc
    tuple val(meta), env('COVERAGE'), emit: coverage
    tuple val(meta), path("${meta.id}_NanoPlot-report.html"), emit: html
    path "versions.yml", emit: versions

    script:
    """
    NanoPlot \\
        --fastq ${reads} \\
        --threads ${task.cpus} \\
        --outdir . \\
        --prefix ${meta.id}_ \\
        --no_static \\
        --N50

    TOTAL_BASES=\$(awk -F':' '/^Total bases/ {gsub(/[ ,]/, "", \$2); printf "%d", \$2}' ${meta.id}_NanoStats.txt)
    if [ -n "\$TOTAL_BASES" ]; then
        COVERAGE=\$(awk -v bases="\$TOTAL_BASES" -v size=${genome_bp} 'BEGIN {printf "%.1f", bases / size}')
        STATUS=\$(awk -v cov="\$COVERAGE" -v min=${params.min_nanopore_coverage} 'BEGIN {print ((cov < min) ? "LOW" : "PASS")}')
    else
        TOTAL_BASES="NA"; COVERAGE="NA"; STATUS="UNKNOWN"
    fi
    echo "Estimated coverage after trimming: \${COVERAGE}x (\$STATUS)"

    {
        echo "# id: 'nanopore_coverage'"
        echo "# section_name: 'Nanopore coverage'"
        echo "# description: 'Reads after chopper: total bases / ${genome_bp} bp. LOW = below ${params.min_nanopore_coverage}x.'"
        echo "# plot_type: 'table'"
        printf 'Sample\\tTotal bases\\tEstimated coverage (x)\\tStatus\\n'
        printf '%s\\t%s\\t%s\\t%s\\n' "${meta.id}" "\$TOTAL_BASES" "\$COVERAGE" "\$STATUS"
    } > ${meta.id}_coverage_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$(NanoPlot --version 2>&1 | sed 's/NanoPlot //')
    END_VERSIONS
    """

    stub:
    """
    printf 'Total bases:\\t165,000,000\\n' > ${meta.id}_NanoStats.txt
    printf 'Sample\\tTotal bases\\tEstimated coverage (x)\\tStatus\\n${meta.id}\\t165000000\\t30.0\\tPASS\\n' > ${meta.id}_coverage_mqc.tsv
    touch ${meta.id}_NanoPlot-report.html
    COVERAGE=30.0
    printf '"%s":\\n    nanoplot: stub\\n' "${task.process}" > versions.yml
    """
}
