process PLASMID_SHARING {
    label 'process_single'
    conda "conda-forge::python=3.12"
    container "quay.io/biocontainers/python:3.12.12"

    input:
    path plasmid_jsons

    output:
    path "plasmid_sharing.tsv", emit: tsv
    path "versions.yml", emit: versions

    script:
    """
    plasmid_sharing.py --inputs ${plasmid_jsons} --output plasmid_sharing.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf 'cluster_id\\n' > plasmid_sharing.tsv
    printf '"%s":\\n    python: stub\\n' "${task.process}" > versions.yml
    """
}
