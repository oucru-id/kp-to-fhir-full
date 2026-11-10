include { CHOPPER       } from '../../modules/local/chopper'
include { NANOPLOT      } from '../../modules/local/nanoplot'
include { FLYE          } from '../../modules/local/flye'
include { MEDAKA_POLISH } from '../../modules/local/medaka_polish'
include { ASSEMBLY_QC   } from '../../modules/local/assembly_qc'

workflow ASSEMBLY_NANOPORE {
    take:
    reads
    genome_bp

    main:
    def min_coverage = params.min_nanopore_coverage as double
    CHOPPER(reads)
    NANOPLOT(CHOPPER.out.reads, genome_bp)
    FLYE(CHOPPER.out.reads)
    MEDAKA_POLISH(FLYE.out.assembly.join(CHOPPER.out.reads, failOnMismatch: true))
    ASSEMBLY_QC(NANOPLOT.out.coverage.join(MEDAKA_POLISH.out.polish_status, failOnMismatch: true))

    MEDAKA_POLISH.out.polish_status
        .filter { _meta, status -> status != 'polished' }
        .subscribe { meta, _status ->
            log.warn "${meta.id}: medaka polishing failed, so the unpolished Flye assembly is used"
        }
    NANOPLOT.out.coverage
        .filter { _meta, coverage -> coverage == 'NA' || (coverage as double) < min_coverage }
        .subscribe { meta, coverage ->
            log.warn "${meta.id}: estimated Nanopore coverage after trimming is ${coverage}x " +
                     "(minimum ${params.min_nanopore_coverage}x); assembly, MLST/cgMLST and resistance " +
                     "calls for this sample may be unreliable"
        }

    emit:
    assembly    = MEDAKA_POLISH.out.assembly
    assembly_qc = ASSEMBLY_QC.out.json
    qc          = NANOPLOT.out.qc.map { meta, stats, coverage_table -> [meta, [stats, coverage_table]] }
    versions    = CHOPPER.out.versions.mix(NANOPLOT.out.versions, FLYE.out.versions, MEDAKA_POLISH.out.versions,
                                           ASSEMBLY_QC.out.versions)
}
