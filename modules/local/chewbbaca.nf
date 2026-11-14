process CHEWBBACA {
    label 'process_high'
    conda "bioconda::chewbbaca=3.4.2"
    container "quay.io/biocontainers/chewbbaca:3.4.2--pyhdfd78af_0"

    input:
    path assemblies
    path schema_dir

    output:
    path "results_alleles.tsv", emit: alleles
    path "cgmlst_statistics.txt", emit: stats
    path "versions.yml", emit: versions

    script:
    """
    chewBBACA.py AlleleCall \\
        -i . \\
        -g ${schema_dir} \\
        -o chewbbaca_out \\
        --cpu ${task.cpus} \\
        --force-continue

    find chewbbaca_out -name "results_alleles.tsv" -exec mv {} . \\;
    find chewbbaca_out -name "results_statistics.tsv" -exec mv {} cgmlst_statistics.txt \\;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chewbbaca: \$(chewBBACA.py --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
        cgmlst_schema: ${schema_dir.name}
    END_VERSIONS
    """

    stub:
    """
    printf 'FILE\\tKP1_RS00005\\n' > results_alleles.tsv
    touch cgmlst_statistics.txt
    printf '"%s":\\n    chewbbaca: stub\\n' "${task.process}" > versions.yml
    """
}
