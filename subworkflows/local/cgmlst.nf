include { CHEWBBACA } from '../../modules/local/chewbbaca'

workflow CGMLST {
    take:
    assemblies
    schema

    main:
    CHEWBBACA(assemblies.map { _meta, assembly -> assembly }.collect(), schema)

    emit:
    alleles  = CHEWBBACA.out.alleles
    versions = CHEWBBACA.out.versions
}
