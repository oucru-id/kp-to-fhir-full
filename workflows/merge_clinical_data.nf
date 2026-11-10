#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process MERGE_CLINICAL_FHIR {
    publishDir "${params.results_dir}/fhir_merged", mode: 'copy'

    input:
    path(fhir_bundle)
    each path(clinical_metadata)
    each path(org_metadata)
    each path(practitioner_metadata)


    output:
    path "*.merged.fhir.json", emit: merged_fhir

    script:
    def prefix = fhir_bundle.simpleName
    """
    
    merge_clinical_fhir.py \\
        --input ${fhir_bundle} \\
        --output ${prefix}.merged.fhir.json \\
        --patient_metadata ${clinical_metadata} \\
        --organization_metadata ${org_metadata} \\
        --practitioner_metadata ${practitioner_metadata}
    """
}

workflow MERGE_CLINICAL_DATA {
    take:
    fhir_ch
    clinical_ch
    org_ch
    practitioner_ch

    main:
    MERGE_CLINICAL_FHIR(fhir_ch, clinical_ch, org_ch, practitioner_ch)

    emit:
    merged_fhir = MERGE_CLINICAL_FHIR.out.merged_fhir
}