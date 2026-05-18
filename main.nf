#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ILLUMINA }         from './workflows/illumina.nf'
include { NANOPORE }         from './workflows/nanopore.nf'
include { TYPING_ANALYSIS }  from './workflows/typing.nf'
include { CGMLST }           from './workflows/cgmlst.nf'
include { GENERATE_REPORT }  from './workflows/report.nf'
include { FHIR }             from './workflows/fhir.nf'
include { VALIDATE }         from './workflows/validate_fhir.nf'
include { MERGE_CLINICAL_DATA } from './workflows/merge_clinical_data.nf'
include { UPLOAD_FHIR }      from './workflows/upload_fhir.nf'
include { VERSIONS }         from './workflows/utils.nf'

workflow {
    
log.info """
    Klebsiella pneumoniae Resistance Gene Analysis Pipeline (v${params.version})
    Developed by SPHERES-OUCRU ID
    Documentation: https://kp-pipeline-docs.readthedocs.io/
"""

    illumina_reads_ch = Channel
        .fromFilePairs("${params.reads_dir}/*_{1,2}_illumina.fastq.gz")
        .map { id, files -> tuple(id, files) }

    nanopore_reads_ch = Channel
        .fromPath("${params.reads_dir}/*_ont.fastq.gz", checkIfExists: true)
        .map { file -> tuple(file.simpleName.replaceFirst(/_ont$/, ''), file) }

    illumina_out = ILLUMINA(illumina_reads_ch)
    nanopore_out = NANOPORE(nanopore_reads_ch)

    all_assemblies = illumina_out.assembly
        .mix(nanopore_out.assembly)

    typing_results = TYPING_ANALYSIS(all_assemblies)
    
    cgmlst_results = Channel.empty()
    schema_dir = file(params.cgmlst_schema)
    
    if(schema_dir.exists() && schema_dir.isDirectory()) {
        cgmlst_out = CGMLST(all_assemblies, schema_dir)
        cgmlst_results = cgmlst_out.matrix
    } else {
        cgmlst_results = Channel.value(file("NO_CGMLST"))
    }

    all_qc = illumina_out.qc_report
        .mix(nanopore_out.qc_report)
        .map { it -> it[1] }
        .flatten()
        .filter { it.toString().endsWith('_fastqc.zip') }
        .collect()

    GENERATE_REPORT(
        all_qc, 
        typing_results.typing_report,
        typing_results.lineage_report
    )

    org_metadata_ch = Channel.fromPath(params.organization_metadata, checkIfExists: false)
        .first()

    fhir_results = FHIR(
        typing_results.typing_report,
        typing_results.lineage_report,
        cgmlst_results,
        org_metadata_ch
    )

    clinical_metadata_ch = Channel.fromPath(params.patient_metadata, checkIfExists: false)
        .first()

    practitioner_metadata_ch = Channel.fromPath(params.practitioner_metadata, checkIfExists: false)
        .first()

    merged_fhir = MERGE_CLINICAL_DATA(fhir_results.fhir_output.map { it -> it[1] }, clinical_metadata_ch, org_metadata_ch, practitioner_metadata_ch)

    validated_fhir = VALIDATE(merged_fhir.merged_fhir)

   UPLOAD_FHIR(validated_fhir.validated_fhir)

    VERSIONS()
}
