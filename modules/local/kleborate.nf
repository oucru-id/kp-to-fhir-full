process KLEBORATE {
    tag "${meta.id}"
    label 'process_medium'
    conda "bioconda::kleborate=3.2.4"
    container "quay.io/biocontainers/kleborate:3.2.4--pyhdfd78af_1"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.id}_kleborate_full.txt"), emit: results
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p kleborate_output

    kleborate -a ${assembly} \\
              --modules klebsiella_pneumo_complex__mlst,klebsiella_pneumo_complex__amr,klebsiella_pneumo_complex__resistance_score,klebsiella_pneumo_complex__resistance_gene_count,klebsiella_pneumo_complex__resistance_class_count,klebsiella__ybst,klebsiella__cbst,klebsiella__abst,klebsiella__smst,klebsiella__rmst,klebsiella__rmpa2,klebsiella_pneumo_complex__virulence_score,klebsiella_pneumo_complex__kaptive,klebsiella_pneumo_complex__wzi \\
              -o kleborate_output

    if [ -f kleborate_output/results.txt ]; then
        cp kleborate_output/results.txt ${meta.id}_kleborate_full.txt
    elif [ -f kleborate_output/kleborate_results.txt ]; then
        cp kleborate_output/kleborate_results.txt ${meta.id}_kleborate_full.txt
    else
        output_file=\$(find kleborate_output -name "*.txt" | head -1)
        if [ -z "\$output_file" ]; then
            ls -la kleborate_output/ >&2
            echo "ERROR: Kleborate wrote no results table" >&2
            exit 1
        fi
        cp "\$output_file" ${meta.id}_kleborate_full.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kleborate: \$(kleborate --version 2>&1 | sed 's/Kleborate v//')
        kaptive: \$(kaptive --version 2>&1 | sed 's/[^0-9.]//g')
    END_VERSIONS
    """

    stub:
    """
    printf 'strain\\tklebsiella_pneumo_complex__mlst__ST\\n${meta.id}\\tST11\\n' > ${meta.id}_kleborate_full.txt
    printf '"%s":\\n    kleborate: stub\\n' "${task.process}" > versions.yml
    """
}
