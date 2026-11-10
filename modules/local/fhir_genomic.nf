process FHIR_GENOMIC {
    tag "${meta.id}"
    label 'process_single'
    conda "conda-forge::python=3.12"
    container "quay.io/biocontainers/python:3.12.12"

    input:
    tuple val(meta), path(typing_json), path(plasmid_json, stageAs: 'plasmid/*'), path(qc_json, stageAs: 'qc/*')
    path cgmlst_alleles, stageAs: 'cgmlst/*'

    output:
    tuple val(meta), path("${meta.id}.genomic.fhir.json"), emit: json
    path "versions.yml", emit: versions

    script:
    def args = []
    if (plasmid_json.name != 'NO_FILE') args << "--plasmid_json ${plasmid_json}"
    if (qc_json.name != 'NO_FILE') args << "--qc_json ${qc_json}"
    if (cgmlst_alleles.name != 'NO_FILE') args << "--cgmlst_file ${cgmlst_alleles}"
    if (params.performer_org_id) args << "--performer_org_id '${params.performer_org_id}'"
    if (params.performer_org_name) args << "--performer_org_name '${params.performer_org_name}'"
    if (params.identifier_system_base) args << "--identifier_system_base '${params.identifier_system_base}'"
    """
    kpfhir_genomic.py \\
        --typing_json ${typing_json} \\
        --sample_id '${meta.id}' \\
        --output ${meta.id}.genomic.fhir.json \\
        --fhir_profile ${params.fhir_profile} \\
        --pipeline_version '${workflow.manifest.version}' \\
        ${args.join(' ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        fhir_profile: ${params.fhir_profile}
    END_VERSIONS
    """

    stub:
    """
    echo '{"resourceType": "Bundle", "type": "collection", "entry": []}' > ${meta.id}.genomic.fhir.json
    printf '"%s":\\n    python: stub\\n' "${task.process}" > versions.yml
    """
}
