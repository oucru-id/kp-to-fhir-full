#!/usr/bin/env nextflow

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
include { KPFHIR } from './workflows/kpfhir'

workflow {
    main:
    log.info """\
        KPtoFHIR ${workflow.manifest.version}: Klebsiella pneumoniae genomics to HL7 FHIR
        ${workflow.manifest.homePage}
        """.stripIndent()
    if (params.results_dir) {
        error "--results_dir was renamed to --outdir in 2.0.0"
    }
    validateParameters()
    log.info paramsSummaryLog(workflow)
    KPFHIR()
}
