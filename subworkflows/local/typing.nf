include { KLEBORATE       } from '../../modules/local/kleborate'
include { PARSE_KLEBORATE } from '../../modules/local/parse_kleborate'

workflow TYPING {
    take:
    assemblies

    main:
    KLEBORATE(assemblies)
    PARSE_KLEBORATE(KLEBORATE.out.results)

    emit:
    typing   = PARSE_KLEBORATE.out.typing
    lineage  = PARSE_KLEBORATE.out.lineage
    versions = KLEBORATE.out.versions.mix(PARSE_KLEBORATE.out.versions)
}
