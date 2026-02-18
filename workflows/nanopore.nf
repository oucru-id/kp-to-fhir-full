#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process fastqc {
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_fastqc.zip"), emit: qc_report  

    script:
    """
    export _JAVA_OPTIONS="-Xmx3g"
    fastqc ${reads} --outdir .
    """
}

process chopper {
    publishDir "${params.results_dir}/trimming", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz"), emit: trimmed_reads

    script:
    """
        zcat ${reads} | chopper -q ${params.nanopore_min_q} -l ${params.nanopore_min_l} | gzip > ${sample_id}_trimmed.fastq.gz
    """
}

process flye {
    publishDir "${params.results_dir}/assemblies", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_assembly.fasta")

    script:
    """
    GENOME_SIZE="5.5m"
    
    flye --nano-raw ${reads} \
         --out-dir flye_output \
         --genome-size \$GENOME_SIZE \
         --threads ${task.cpus} \
         --iterations 2
    
    if [ -f "flye_output/assembly.fasta" ]; then
        cp flye_output/assembly.fasta ${sample_id}_assembly.fasta
    else
        echo "Error: Flye assembly failed."
        exit 1
    fi
    
    echo "Flye assembly statistics:"
    awk '/^>/ {if (seqlen) {print name, seqlen}; name=\$0; seqlen=0; next} {seqlen += length(\$0)} END {print name, seqlen}' ${sample_id}_assembly.fasta | head -10
    
    echo "Flye assembly info:"
    if [ -f "flye_output/assembly_info.txt" ]; then
        cat flye_output/assembly_info.txt
    fi
    """
}

process medaka_polish {
    publishDir "${params.results_dir}/assemblies", mode: 'copy'
    
    input:
    tuple val(sample_id), path(raw_assembly), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_assembly.fasta")

    script:
    """
    mv ${raw_assembly} raw_input.fasta

    mkdir -p medaka_work
    MODEL="r941_e81_sup_variant_g514"
    
    echo "Running Medaka with model: \$MODEL"
    
    set +e
    medaka_consensus -i ${reads} -d raw_input.fasta -o medaka_work -t ${task.cpus} -m \$MODEL
    MEDAKA_EXIT_CODE=\$?
    set -e
    
    if [ \$MEDAKA_EXIT_CODE -eq 0 ] && [ -s medaka_work/consensus.fasta ]; then
        echo "Medaka polishing successful."
        cp medaka_work/consensus.fasta ${sample_id}_assembly.fasta
    else
        echo "Medaka polishing failed (Exit code: \$MEDAKA_EXIT_CODE)."
        cp raw_input.fasta ${sample_id}_assembly.fasta
    fi
    
    echo "Final assembly statistics:"
    awk '/^>/ {if (seqlen) {print name, seqlen}; name=\$0; seqlen=0; next} {seqlen += length(\$0)} END {print name, seqlen}' ${sample_id}_assembly.fasta | head -10
    """
}

workflow NANOPORE {
    take:
    reads_ch
    
    main:
    fastqc_out = fastqc(reads_ch)
    trimmed_ch = chopper(reads_ch)
    raw_assembly = flye(trimmed_ch)
    
    assembly_out = medaka_polish(raw_assembly.combine(trimmed_ch, by: 0))
    
    emit:
    qc_report = fastqc_out.qc_report
    assembly = assembly_out
}