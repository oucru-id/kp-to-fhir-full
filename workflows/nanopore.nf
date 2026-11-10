#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process chopper {
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz"), emit: trimmed_reads

    script:
    """
        zcat ${reads} | chopper -q ${params.nanopore_min_q} -l ${params.nanopore_min_l} | gzip > ${sample_id}_trimmed.fastq.gz
    """
}

process nanoplot {
    publishDir "${params.results_dir}/qc/nanoplot", mode: 'copy', pattern: "*_{NanoPlot-report.html,NanoStats.txt,coverage_mqc.tsv}"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_NanoStats.txt"), path("${sample_id}_coverage_mqc.tsv"), emit: qc_report
    tuple val(sample_id), env('COVERAGE'), emit: coverage
    path "${sample_id}_NanoPlot-report.html", emit: report

    script:
    def nanoplot_env = params.nanoplot_env ?: '$HOME/miniforge3/envs/nanoplot'
    """
    PATH="${nanoplot_env}/bin:\$PATH" NanoPlot \\
        --fastq ${reads} \\
        --threads ${task.cpus} \\
        --outdir . \\
        --prefix ${sample_id}_ \\
        --no_static \\
        --N50

    # Coverage of the reads Flye and medaka use = total bases after chopper / expected genome size
    TOTAL_BASES=\$(awk -F':' '/^Total bases/ {gsub(/[ ,]/, "", \$2); printf "%d", \$2}' ${sample_id}_NanoStats.txt)
    if [ -n "\$TOTAL_BASES" ]; then
        COVERAGE=\$(awk -v bases="\$TOTAL_BASES" -v size=${params.genome_size_bp} 'BEGIN {printf "%.1f", bases / size}')
        STATUS=\$(awk -v cov="\$COVERAGE" -v min=${params.min_nanopore_coverage} 'BEGIN {print ((cov < min) ? "LOW" : "PASS")}')
    else
        TOTAL_BASES="NA"; COVERAGE="NA"; STATUS="UNKNOWN"
    fi
    echo "Estimated coverage after trimming: \${COVERAGE}x (\$STATUS)"

    # MultiQC custom-content table; files sharing this id are merged into one table
    {
        echo "# id: 'nanopore_coverage'"
        echo "# section_name: 'Nanopore coverage'"
        echo "# description: 'Reads after chopper: total bases / ${params.genome_size_bp} bp. LOW = below ${params.min_nanopore_coverage}x.'"
        echo "# plot_type: 'table'"
        printf 'Sample\\tTotal bases\\tEstimated coverage (x)\\tStatus\\n'
        printf '%s\\t%s\\t%s\\t%s\\n' "${sample_id}" "\$TOTAL_BASES" "\$COVERAGE" "\$STATUS"
    } > ${sample_id}_coverage_mqc.tsv
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
    tuple val(sample_id), path("${sample_id}_assembly.fasta"), emit: assembly
    tuple val(sample_id), env('POLISH_STATUS'), emit: polish_status

    script:
    """
    mv ${raw_assembly} raw_input.fasta

    mkdir -p medaka_work
    MODEL="${params.medaka_model}"
    
    echo "Running Medaka with model: \$MODEL"
    
    set +e
    medaka_consensus -i ${reads} -d raw_input.fasta -o medaka_work -t ${task.cpus} -m \$MODEL -b ${params.medaka_batch_size}
    MEDAKA_EXIT_CODE=\$?
    if [ \$MEDAKA_EXIT_CODE -ne 0 ] || [ ! -s medaka_work/consensus.fasta ]; then
        # e.g. GPU out of memory; CPU inference is slower but needs no GPU memory
        echo "Medaka failed (exit code \$MEDAKA_EXIT_CODE); retrying on CPU." >&2
        rm -rf medaka_work
        CUDA_VISIBLE_DEVICES="" medaka_consensus -i ${reads} -d raw_input.fasta -o medaka_work -t ${task.cpus} -m \$MODEL -b ${params.medaka_batch_size}
        MEDAKA_EXIT_CODE=\$?
    fi
    set -e
    
    if [ \$MEDAKA_EXIT_CODE -eq 0 ] && [ -s medaka_work/consensus.fasta ]; then
        echo "Medaka polishing successful."
        cp medaka_work/consensus.fasta ${sample_id}_assembly.fasta
        POLISH_STATUS="polished"
    else
        echo "Medaka polishing failed (Exit code: \$MEDAKA_EXIT_CODE)." >&2
        cp raw_input.fasta ${sample_id}_assembly.fasta
        POLISH_STATUS="failed"
    fi
    
    echo "Final assembly statistics:"
    awk '/^>/ {if (seqlen) {print name, seqlen}; name=\$0; seqlen=0; next} {seqlen += length(\$0)} END {print name, seqlen}' ${sample_id}_assembly.fasta | head -10
    """
}

process assembly_qc {
    publishDir "${params.results_dir}/qc/assembly", mode: 'copy'

    input:
    tuple val(sample_id), val(coverage), val(polish_status)

    output:
    tuple val(sample_id), path("${sample_id}_assembly_qc.json"), emit: qc_json

    script:
    """
    assembly_qc.py \\
        --sample_id ${sample_id} \\
        --platform ont \\
        --coverage ${coverage} \\
        --min_coverage ${params.min_nanopore_coverage} \\
        --polishing ${polish_status} \\
        --polishing_tool medaka \\
        --polishing_model ${params.medaka_model} \\
        --output ${sample_id}_assembly_qc.json
    """
}

workflow NANOPORE {
    take:
    reads_ch
    
    main:
    trimmed_ch = chopper(reads_ch)
    nanoplot_out = nanoplot(trimmed_ch)
    raw_assembly = flye(trimmed_ch)

    medaka_out = medaka_polish(raw_assembly.combine(trimmed_ch, by: 0))

    medaka_out.polish_status
        .filter { _sample_id, status -> status != 'polished' }
        .subscribe { sample_id, _status ->
            log.warn "${sample_id}: medaka polishing failed, so the unpolished Flye assembly is used. " +
                     "See .command.err of medaka_polish for this sample."
        }

    qc_out = assembly_qc(nanoplot_out.coverage.join(medaka_out.polish_status))

    nanoplot_out.coverage
        .filter { _sample_id, coverage -> coverage == 'NA' || coverage.toFloat() < params.min_nanopore_coverage }
        .subscribe { sample_id, coverage ->
            log.warn "${sample_id}: estimated Nanopore coverage after trimming is ${coverage}x (minimum ${params.min_nanopore_coverage}x). " +
                     "Assembly, MLST/cgMLST and resistance calls for this sample may be unreliable."
        }

    emit:
    qc_report = nanoplot_out.qc_report.map { sample_id, stats, coverage_table -> tuple(sample_id, [stats, coverage_table]) }
    assembly = medaka_out.assembly
    assembly_qc = qc_out.qc_json
}