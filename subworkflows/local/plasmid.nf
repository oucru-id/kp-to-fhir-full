include { KLEBORATE_CARD  } from '../../modules/local/kleborate_card'
include { MOB_INIT        } from '../../modules/local/mob_init'
include { MOB_RECON       } from '../../modules/local/mob_recon'
include { AMR_LOCATE      } from '../../modules/local/amr_locate'
include { PARSE_PLASMIDS  } from '../../modules/local/parse_plasmids'
include { PLASMID_SHARING } from '../../modules/local/plasmid_sharing'

workflow PLASMID {
    take:
    assemblies
    typing

    main:
    KLEBORATE_CARD()
    card = KLEBORATE_CARD.out.fasta
    ch_database = channel.empty()
    if (params.mobsuite_db) {
        ch_database = channel.value(file(params.mobsuite_db, checkIfExists: true))
    } else {
        MOB_INIT()
        ch_database = MOB_INIT.out.db
    }
    MOB_RECON(assemblies, ch_database)
    AMR_LOCATE(assemblies, card)

    ch_inputs = assemblies
        .map { meta, assembly -> [meta.id, meta, assembly] }
        .join(MOB_RECON.out.recon.map { meta, report, mobtyper, status -> [meta.id, report, mobtyper, status] },
              failOnMismatch: true)
        .join(AMR_LOCATE.out.paf.map { meta, paf -> [meta.id, paf] }, failOnMismatch: true)
        .join(typing.map { meta, json -> [meta.id, json] }, failOnMismatch: true)
        .map { _id, meta, assembly, report, mobtyper, status, paf, json ->
            [meta, assembly, report, mobtyper, status, paf, json] }
    PARSE_PLASMIDS(ch_inputs, card)
    PLASMID_SHARING(PARSE_PLASMIDS.out.json.map { _meta, json -> json }.collect())

    emit:
    plasmid  = PARSE_PLASMIDS.out.json
    sharing  = PLASMID_SHARING.out.tsv
    versions = KLEBORATE_CARD.out.versions.mix(MOB_RECON.out.versions, AMR_LOCATE.out.versions,
                                               PARSE_PLASMIDS.out.versions, PLASMID_SHARING.out.versions)
}
