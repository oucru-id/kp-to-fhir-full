process KLEBORATE_CARD {
    label 'process_single'
    conda "bioconda::kleborate=3.2.4"
    container "quay.io/biocontainers/kleborate:3.2.4--pyhdfd78af_1"

    output:
    path "CARD_v*.fasta", emit: fasta
    path "versions.yml", emit: versions

    script:
    """
    CARD_FASTA=\$(parse_plasmids.py --find_card)
    cp "\$CARD_FASTA" .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kleborate_card: \$(basename "\$CARD_FASTA" .fasta)
    END_VERSIONS
    """

    stub:
    """
    printf '>1__KPC_Bla__KPC-2__1\\nACGT\\n' > CARD_v0.0.0.fasta
    printf '"%s":\\n    kleborate_card: stub\\n' "${task.process}" > versions.yml
    """
}
