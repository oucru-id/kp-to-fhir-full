#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process PLASMID_RECON {
    publishDir path: { "${params.results_dir}/plasmid/${sample_id}" }, mode: 'copy', pattern: '*.{txt,paf,log}'

    input:
    tuple val(sample_id), path(assembly), val(asm_type)

    output:
    tuple val(sample_id), val(asm_type), path(assembly), path("contig_report.txt"), path("mobtyper_results.txt"), path("amr_hits.paf"), path("recon_status.txt"), emit: recon
    path "versions.yml", emit: versions

    script:
    def mob_env = params.mobsuite_env ?: '$HOME/miniforge3/envs/mobsuite'
    def db_arg  = params.mobsuite_db ? "-d ${params.mobsuite_db}" : ""
    """
    CARD_FASTA=\$(parse_plasmids.py --find_card)
    MOB_VERSION=\$(PATH="${mob_env}/bin:\$PATH" mob_recon --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)

    set +e
    PATH="${mob_env}/bin:\$PATH" mob_recon \\
        --infile ${assembly} \\
        --outdir mob_out \\
        --num_threads ${task.cpus} \\
        ${db_arg} > mob_recon.log 2>&1
    MOB_EXIT_CODE=\$?
    set -e

    if [ \$MOB_EXIT_CODE -eq 0 ] && [ -s mob_out/contig_report.txt ]; then
        echo "status=ok" > recon_status.txt
        cp mob_out/contig_report.txt contig_report.txt
    else
        echo "MOB-recon failed (exit code: \$MOB_EXIT_CODE), see mob_recon.log" >&2
        echo "status=failed" > recon_status.txt
        : > contig_report.txt
    fi

    if [ -s mob_out/mobtyper_results.txt ]; then
        cp mob_out/mobtyper_results.txt mobtyper_results.txt
    else
        : > mobtyper_results.txt
    fi

    echo "mob_suite=\${MOB_VERSION:-unknown}" >> recon_status.txt
    echo "card_fasta=\$CARD_FASTA" >> recon_status.txt

    minimap2 --end-bonus=10 --eqx -c -x map-ont ${assembly} "\$CARD_FASTA" > amr_hits.paf 2> minimap2.log || : > amr_hits.paf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mob_suite: \${MOB_VERSION:-unknown}
        minimap2: \$(minimap2 --version)
    END_VERSIONS
    """
}

process PARSE_PLASMIDS {
    publishDir "${params.results_dir}/plasmid", mode: 'copy'

    input:
    tuple val(sample_id), val(asm_type), path(assembly), path(contig_report), path(mobtyper), path(paf), path(recon_status), path(typing_json)

    output:
    tuple val(sample_id), path("${sample_id}_plasmid.json"), emit: plasmid_report

    script:
    """
    parse_plasmids.py \\
        --sample_id ${sample_id} \\
        --typing_json ${typing_json} \\
        --assembly ${assembly} \\
        --assembly_type ${asm_type} \\
        --contig_report ${contig_report} \\
        --mobtyper ${mobtyper} \\
        --paf ${paf} \\
        --recon_status ${recon_status} \\
        --min_identity ${params.amr_locate_min_identity} \\
        --min_coverage ${params.amr_locate_min_coverage} \\
        --chrom_min_contig ${params.plasmid_chrom_min_contig} \\
        --output ${sample_id}_plasmid.json
    """
}

process PLASMID_SHARING {
    publishDir "${params.results_dir}/plasmid", mode: 'copy'

    input:
    path(plasmid_jsons)

    output:
    path "plasmid_sharing.tsv", emit: sharing

    script:
    """
    plasmid_sharing.py \\
        --inputs ${plasmid_jsons} \\
        --output plasmid_sharing.tsv
    """
}

workflow PLASMID {
    take:
    assemblies_ch
    typing_ch

    main:
    PLASMID_RECON(assemblies_ch)

    PARSE_PLASMIDS(PLASMID_RECON.out.recon.join(typing_ch))

    PLASMID_SHARING(PARSE_PLASMIDS.out.plasmid_report.map { _sample_id, plasmid_json -> plasmid_json }.collect())

    emit:
    plasmid_report = PARSE_PLASMIDS.out.plasmid_report
    sharing        = PLASMID_SHARING.out.sharing
    versions       = PLASMID_RECON.out.versions
}
