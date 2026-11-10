process MOB_RECON {
    tag "${meta.id}"
    label 'process_medium'
    conda "bioconda::mob_suite=3.1.9"
    container "quay.io/biocontainers/mob_suite:3.1.9--pyhdfd78af_1"

    input:
    tuple val(meta), path(assembly)
    path database

    output:
    tuple val(meta), path("contig_report.txt"), path("mobtyper_results.txt"), path("recon_status.txt"), emit: recon
    tuple val(meta), path("mob_recon.log"), emit: log
    path "versions.yml", emit: versions

    script:
    """
    MOB_VERSION=\$(mob_recon --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)

    set +e
    mob_recon \\
        --infile ${assembly} \\
        --outdir mob_out \\
        --num_threads ${task.cpus} \\
        -d ${database} > mob_recon.log 2>&1
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mob_suite: \${MOB_VERSION:-unknown}
    END_VERSIONS
    """

    stub:
    """
    printf 'status=ok\\nmob_suite=stub\\n' > recon_status.txt
    : > contig_report.txt
    : > mobtyper_results.txt
    touch mob_recon.log
    printf '"%s":\\n    mob_suite: stub\\n' "${task.process}" > versions.yml
    """
}
