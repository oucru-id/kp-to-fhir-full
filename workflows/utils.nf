nextflow.enable.dsl = 2

process VERSIONS {
    publishDir "${params.results_dir}", mode: 'copy'

    output:
    path "software_versions.yml"

    script:
    """
    echo "pipeline:" > software_versions.yml
    echo "  name: Klebsiella pneumoniae Resistance Gene Analysis Pipeline" >> software_versions.yml
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
    echo "    trimming: 'chopper (Q10, L500)'" >> software_versions.yml
    echo "    assembler: 'flye'" >> software_versions.yml
    echo "    polishing: 'medaka'" >> software_versions.yml
    
    echo "  typing:" >> software_versions.yml
    echo "    tool: 'kleborate'" >> software_versions.yml
    echo "    modules: 'mlst, resistance, virulence, kaptive'" >> software_versions.yml
    
    export BASE_DIR="${baseDir}"
    python3 $baseDir/scripts/get_versions.py >> software_versions.yml
    """
}