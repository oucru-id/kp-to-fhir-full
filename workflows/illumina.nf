#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process fastqc {
    input:
    tuple val(sample_id), path(reads)  

    output:
    tuple val(sample_id), path("*_fastqc.zip"), emit: qc_report

    script:
    """
    fastqc ${reads.join(' ')} --outdir .
    """
}

process fastp {
    publishDir "${params.results_dir}/trimming", mode: 'copy', pattern: "*.{json,html}"
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_R{1,2}_trimmed.fastq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("${sample_id}_fastp.json"), emit: fastp_json
    tuple val(sample_id), path("${sample_id}_fastp.html"), emit: fastp_html

    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} \
          -o ${sample_id}_R1_trimmed.fastq.gz \
          -O ${sample_id}_R2_trimmed.fastq.gz \
          --json ${sample_id}_fastp.json \
          --html ${sample_id}_fastp.html \
          --thread ${task.cpus} \
          --qualified_quality_phred 20 \
          --unqualified_percent_limit 30 \
          --length_required 50 \
          --detect_adapter_for_pe \
          --cut_front \
          --cut_tail \
          --cut_window_size 4 \
          --cut_mean_quality 20
    """
}

process megahit {
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_raw_contigs.fasta"), path(reads)

    script:
    """
    megahit -1 ${reads[0]} -2 ${reads[1]} \
        -o megahit_output \
        --min-contig-len 500 \
        -t ${task.cpus} \
        --memory 0.5
    
    cp megahit_output/final.contigs.fa ${sample_id}_raw_contigs.fasta
    """
}

process pilon_polish {
    publishDir "${params.results_dir}/assemblies", mode: 'copy'
    
    input:
    tuple val(sample_id), path(raw_assembly), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_contigs.fasta")

    script:
    """
    bwa index ${raw_assembly}
    
    bwa mem -t ${task.cpus} ${raw_assembly} ${reads[0]} ${reads[1]} | \
    samtools sort -@ ${task.cpus} -o aligned.bam
    samtools index aligned.bam
    
    java -jar /opt/pilon/pilon-1.24.jar --genome ${raw_assembly} --frags aligned.bam --output polished --changes
    
    cp polished.fasta ${sample_id}_contigs.fasta
    """
}

workflow ILLUMINA {
    take:
    reads_ch
    
    main:
    fastqc_out = fastqc(reads_ch)
    trimmed = fastp(reads_ch)
    raw_assembly = megahit(trimmed.trimmed_reads)
    // raw_assembly = spades(trimmed.trimmed_reads)
    polished_assembly = pilon_polish(raw_assembly)
    
    emit:
    qc_report = fastqc_out.qc_report
    assembly = polished_assembly
    fastp_json = trimmed.fastp_json
    fastp_html = trimmed.fastp_html
}