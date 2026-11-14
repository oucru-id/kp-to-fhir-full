process MOB_INIT {
    label 'process_low'
    conda "bioconda::mob_suite=3.1.9"
    container "quay.io/biocontainers/mob_suite:3.1.9--pyhdfd78af_1"
    storeDir "${params.mobsuite_db_cache}"

    output:
    path "mob_suite_db", emit: db

    script:
    """
    mob_init -d mob_suite_db
    """

    stub:
    """
    mkdir mob_suite_db
    touch mob_suite_db/status.txt
    """
}
