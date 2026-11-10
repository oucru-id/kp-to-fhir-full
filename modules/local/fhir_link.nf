process FHIR_LINK {
    label 'process_single'
    conda "conda-forge::python=3.12"
    container "quay.io/biocontainers/python:3.12.12"

    input:
    path genomic_bundles, stageAs: 'genomic/*'
    path clinical_dir
    path external_refs, stageAs: 'external/*'

    output:
    path "bundle/*.bundle.fhir.json", emit: bundles, optional: true
    path "bundle/link_summary.tsv", emit: summary
    path "versions.yml", emit: versions

    script:
    def args = []
    if (!external_refs.toString().endsWith('NO_FILE')) args << "--external_refs ${external_refs}"
    if (params.require_clinical.toString().toBoolean()) args << "--require_clinical"
    if (params.identifier_system_base) args << "--identifier_system_base '${params.identifier_system_base}'"
    """
    link_fhir_bundle.py \\
        --genomic genomic/*.json \\
        --clinical_dir ${clinical_dir} \\
        --fhir_profile ${params.fhir_profile} \\
        --bundle_request_method ${params.bundle_request_method} \\
        --outdir bundle \\
        ${args.join(' ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir bundle
    printf 'fhir_id\\tstatus\\tbundle\\n' > bundle/link_summary.tsv
    printf '"%s":\\n    python: stub\\n' "${task.process}" > versions.yml
    """
}
