#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process CHEWBBACA {
    publishDir "${params.results_dir}/cgmlst", mode: 'copy'
    
    input:
    path(assemblies)
    path(schema_dir)
    
    output:
    path("results_alleles.tsv"), emit: allele_matrix
    path("cgmlst_statistics.txt"), emit: stats
    
    script:
    """
    chewBBACA.py AlleleCall \
        -i . \
        -g ${schema_dir} \
        -o chewbbaca_out \
        --cpu ${task.cpus} \
        --force-continue
        
    find chewbbaca_out -name "results_alleles.tsv" -exec mv {} . \\;
    
    find chewbbaca_out -name "results_statistics.tsv" -exec mv {} cgmlst_statistics.txt \\;
    """
}

workflow CGMLST {
    take:
    assemblies_ch
    schema_dir
    
    main:
    all_assemblies = assemblies_ch
        .map { id, assembly -> assembly }
        .collect()
        
    CHEWBBACA(all_assemblies, schema_dir)
    
    emit:
    matrix = CHEWBBACA.out.allele_matrix
    stats = CHEWBBACA.out.stats
}