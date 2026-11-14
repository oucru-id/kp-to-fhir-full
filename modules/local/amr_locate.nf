process AMR_LOCATE {
    tag "${meta.id}"
    label 'process_low'
    conda "bioconda::minimap2=2.26"
    container "quay.io/biocontainers/minimap2:2.26--he4a0461_2"

    input:
    tuple val(meta), path(assembly)
    path card_fasta

    output:
    tuple val(meta), path("amr_hits.paf"), emit: paf
    path "versions.yml", emit: versions

    script:
    """
    minimap2 --end-bonus=10 --eqx -c -x map-ont ${assembly} ${card_fasta} > amr_hits.paf 2> minimap2.log || : > amr_hits.paf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
    END_VERSIONS
    """

    stub:
    """
    : > amr_hits.paf
    printf '"%s":\\n    minimap2: stub\\n' "${task.process}" > versions.yml
    """
}
