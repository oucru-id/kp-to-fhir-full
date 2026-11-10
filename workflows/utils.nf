nextflow.enable.dsl = 2

process VERSIONS {
    publishDir "${params.results_dir}", mode: 'copy'

    output:
    path "software_versions.yml"

    script:
    def mob_env = params.mobsuite_env ?: '$HOME/miniforge3/envs/mobsuite'
    def nanoplot_env = params.nanoplot_env ?: '$HOME/miniforge3/envs/nanoplot'
    """
    echo "pipeline:" > software_versions.yml
    echo "  name: KPtoFHIR" >> software_versions.yml
    echo "  version: ${params.version}" >> software_versions.yml
    echo "  nextflow: $nextflow.version" >> software_versions.yml
    
    echo "databases:" >> software_versions.yml
    echo "  cgmlst_schema: \$(basename ${params.cgmlst_schema})" >> software_versions.yml
    
    echo "processing_settings:" >> software_versions.yml
    echo "  illumina:" >> software_versions.yml
    echo "    trimming: 'fastp (auto-adapter, Q20)'" >> software_versions.yml
    echo "    assembler: 'megahit (min-contig 500)'" >> software_versions.yml
    echo "    polishing: 'pilon'" >> software_versions.yml
    
    echo "  nanopore:" >> software_versions.yml
    echo "    read_qc: 'NanoPlot after chopper (warn below ${params.min_nanopore_coverage}x of ${params.genome_size_bp} bp)'" >> software_versions.yml
    echo "    trimming: 'chopper (Q10, L500)'" >> software_versions.yml
    echo "    assembler: 'flye'" >> software_versions.yml
    echo "    polishing: 'medaka (${params.medaka_model})'" >> software_versions.yml
    
    echo "  typing:" >> software_versions.yml
    echo "    tool: 'kleborate'" >> software_versions.yml
    echo "    modules: 'mlst, resistance, virulence, kaptive'" >> software_versions.yml

    echo "  plasmid:" >> software_versions.yml
    echo "    enabled: ${params.plasmid_analysis}" >> software_versions.yml
    echo "    tool: 'mob_suite (mob_recon)'" >> software_versions.yml
    echo "    amr_location: 'minimap2 vs Kleborate CARD'" >> software_versions.yml

    export BASE_DIR="${baseDir}"
    export MOBSUITE_ENV="${mob_env}"
    export NANOPLOT_ENV="${nanoplot_env}"
    get_versions.py >> software_versions.yml
    """
}