process PARSE_PLASMIDS {
    tag "${meta.id}"
    label 'process_single'
    conda "conda-forge::python=3.12"
    container "quay.io/biocontainers/python:3.12.12"

    input:
    tuple val(meta), path(assembly), path(contig_report), path(mobtyper), path(recon_status), path(paf), path(typing_json)
    path card_fasta

    output:
    tuple val(meta), path("${meta.id}_plasmid.json"), emit: json
    path "versions.yml", emit: versions

    script:
    """
    parse_plasmids.py \\
        --sample_id '${meta.id}' \\
        --typing_json ${typing_json} \\
        --assembly ${assembly} \\
        --assembly_type ${meta.platform} \\
        --contig_report ${contig_report} \\
        --mobtyper ${mobtyper} \\
        --paf ${paf} \\
        --recon_status ${recon_status} \\
        --card_fasta ${card_fasta} \\
        --min_identity ${params.amr_locate_min_identity} \\
        --min_coverage ${params.amr_locate_min_coverage} \\
        --chrom_min_contig ${params.plasmid_chrom_min_contig} \\
        --output ${meta.id}_plasmid.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo '{"sample_id": "${meta.id}", "status": "no_plasmids", "plasmids": [], "gene_locations": {}}' > ${meta.id}_plasmid.json
    printf '"%s":\\n    python: stub\\n' "${task.process}" > versions.yml
    """
}
