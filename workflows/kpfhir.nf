
include { samplesheetToList } from 'plugin/nf-schema'
include { ASSEMBLY_ILLUMINA } from '../subworkflows/local/assembly_illumina'
include { ASSEMBLY_NANOPORE } from '../subworkflows/local/assembly_nanopore'
include { TYPING            } from '../subworkflows/local/typing'
include { PLASMID           } from '../subworkflows/local/plasmid'
include { CGMLST            } from '../subworkflows/local/cgmlst'
include { FHIR              } from '../subworkflows/local/fhir'
include { PREPARE_ASSEMBLY  } from '../modules/local/prepare_assembly'
include { MULTIQC           } from '../modules/local/multiqc'
include { SAMPLE_REPORT     } from '../modules/local/sample_report'
include { asFlag            } from '../subworkflows/local/utils_params'

def fhirId(String id) {
    return id.replaceAll(/[^A-Za-z0-9\-]/, '-').replaceAll(/-+/, '-').replaceAll(/^-+|-+$/, '')
}

def present(value) {
    return value != null && value != [] && value.toString() != ''
}

def genomeSizeBp(String size) {
    def text = size.trim().toLowerCase()
    def match = text =~ /^(\d+(?:\.\d+)?)([kmg]?)$/
    if (!match.matches()) {
        error "--genome_size '${size}': use a number with an optional k, m or g suffix (e.g. 5.5m)"
    }
    def factor = [k: 1000G, m: 1000000G, g: 1000000000G].get(match.group(2), 1G)
    return (new BigDecimal(match.group(1)) * factor).longValue()
}

def samplesFromSheet(String sheet) {
    def samples = []
    def rows = samplesheetToList(sheet, "${projectDir}/assets/schema_input.json")
    rows.eachWithIndex { row, index ->
        def meta = row[0]
        def fastq_1 = row[1]
        def fastq_2 = row[2]
        def fasta = row[3]
        def line = index + 2
        if (meta.platform == 'illumina') {
            if (!present(fastq_1) || !present(fastq_2)) {
                error "samplesheet line ${line} (${meta.id}): illumina needs fastq_1 and fastq_2"
            }
            samples << [meta: [id: meta.id, platform: 'illumina'], files: [fastq_1, fastq_2]]
        } else if (meta.platform == 'ont') {
            if (!present(fastq_1)) {
                error "samplesheet line ${line} (${meta.id}): ont needs fastq_1"
            }
            samples << [meta: [id: meta.id, platform: 'ont'], files: fastq_1]
        } else {
            if (!present(fasta)) {
                error "samplesheet line ${line} (${meta.id}): fasta needs a fasta column"
            }
            samples << [meta: [id: meta.id, platform: 'fasta'], files: fasta]
        }
    }
    return samples
}

def samplesFromDirectory() {
    def samples = []
    if (!params.reads_dir) {
        return samples
    }
    def pairs = [:]
    files("${params.reads_dir}/*_{1,2}_illumina.fastq.gz").each { f ->
        def id = f.name.replaceFirst(/_[12]_illumina\.fastq\.gz$/, '')
        if (!pairs.containsKey(id)) {
            pairs[id] = []
        }
        pairs[id] << f
    }
    pairs.keySet().sort().each { id ->
        def reads = pairs[id].sort { f -> f.name }
        if (reads.size() != 2) {
            error "Illumina sample ${id}: expected _1_illumina and _2_illumina files, found ${reads*.name}"
        }
        samples << [meta: [id: id, platform: 'illumina'], files: reads]
    }
    files("${params.reads_dir}/*_ont.fastq.gz").sort { f -> f.name }.each { f ->
        samples << [meta: [id: f.name.replaceFirst(/_ont\.fastq\.gz$/, ''), platform: 'ont'], files: f]
    }
    return samples
}

def validateSamples(samples) {
    if (!samples) {
        error "No samples found. Give --input samplesheet.csv, or --reads_dir with files named " +
              "<id>_1_illumina.fastq.gz / <id>_2_illumina.fastq.gz or <id>_ont.fastq.gz."
    }
    def seen = [:]
    def seenFhir = [:]
    samples.each { s ->
        def id = s.meta.id
        if (!(id ==~ /^[A-Za-z0-9][A-Za-z0-9._-]{0,39}$/)) {
            error "Sample id '${id}': use 1-40 letters, digits, '.', '_' or '-', starting with a letter or digit"
        }
        if (seen.containsKey(id)) {
            error "Sample id '${id}' appears more than once"
        }
        def fid = fhirId(id)
        if (seenFhir.containsKey(fid)) {
            error "Sample ids '${seenFhir[fid]}' and '${id}' become the same FHIR id '${fid}'; rename one"
        }
        if (id ==~ /^\d{16}$/) {
            log.warn "Sample id '${id}' looks like a 16-digit national ID (NIK). Sample ids appear in " +
                     "the PHI-free genomic output; use a laboratory accession instead."
        }
        seen[id] = true
        seenFhir[fid] = id
    }
    if (params.fhir_server_url && !params.patient_metadata) {
        error "--fhir_server_url uploads linked bundles, which need --patient_metadata"
    }
    if (!params.cgmlst_schema) {
        log.info "cgMLST is skipped (no --cgmlst_schema; see tools/cgmlst/fetch_schema.sh)"
    }
}

def byId(ch) {
    return ch.map { meta, value -> [meta.id, value] }
}

workflow KPFHIR {
    main:
    def samples = params.input ? samplesFromSheet(params.input) : samplesFromDirectory()
    validateSamples(samples)
    def counts = samples.countBy { s -> s.meta.platform }
    log.info "Samples: ${samples.size()} (illumina: ${counts.illumina ?: 0}, ont: ${counts.ont ?: 0}, " +
             "fasta: ${counts.fasta ?: 0})"

    def no_file = file("${projectDir}/assets/NO_FILE", checkIfExists: true)
    def genome_bp = genomeSizeBp(params.genome_size.toString())

    ch_samples  = channel.fromList(samples).map { s -> [s.meta, s.files] }
    ch_illumina = ch_samples.filter { meta, _f -> meta.platform == 'illumina' }
    ch_ont      = ch_samples.filter { meta, _f -> meta.platform == 'ont' }
    ch_fasta    = ch_samples.filter { meta, _f -> meta.platform == 'fasta' }
    ch_versions = channel.empty()

    ASSEMBLY_ILLUMINA(ch_illumina)
    ASSEMBLY_NANOPORE(ch_ont, genome_bp)
    PREPARE_ASSEMBLY(ch_fasta)
    ch_versions = ch_versions.mix(ASSEMBLY_ILLUMINA.out.versions, ASSEMBLY_NANOPORE.out.versions,
                                  PREPARE_ASSEMBLY.out.versions)
    ch_assemblies = ASSEMBLY_ILLUMINA.out.assembly.mix(ASSEMBLY_NANOPORE.out.assembly, PREPARE_ASSEMBLY.out.assembly)

    TYPING(ch_assemblies)
    ch_versions = ch_versions.mix(TYPING.out.versions)

    ch_plasmid = channel.empty()
    if (asFlag(params.plasmid_analysis)) {
        PLASMID(ch_assemblies, TYPING.out.typing)
        ch_plasmid = PLASMID.out.plasmid
        ch_versions = ch_versions.mix(PLASMID.out.versions)
    }

    ch_cgmlst = channel.value(no_file)
    if (params.cgmlst_schema) {
        CGMLST(ch_assemblies, file(params.cgmlst_schema, checkIfExists: true))
        ch_cgmlst = CGMLST.out.alleles
        ch_versions = ch_versions.mix(CGMLST.out.versions)
    }

    if ((counts.illumina ?: 0) + (counts.ont ?: 0) > 0) {
        ch_qc = ASSEMBLY_ILLUMINA.out.qc.mix(ASSEMBLY_NANOPORE.out.qc).flatMap { _meta, files -> [files].flatten() }
        MULTIQC(ch_qc.collect())
        ch_versions = ch_versions.mix(MULTIQC.out.versions)
    }

    ch_extras = byId(TYPING.out.typing)
        .join(byId(ch_plasmid), remainder: true)
        .join(byId(ASSEMBLY_NANOPORE.out.assembly_qc), remainder: true)
        .filter { row -> row[1] != null }
    ch_meta = ch_assemblies.map { meta, _f -> [meta.id, meta] }

    SAMPLE_REPORT(
        ch_meta.join(ch_extras, failOnMismatch: true)
            .join(byId(TYPING.out.lineage), failOnMismatch: true)
            .map { row -> [row[1], row[2], row[5], row[3] ?: no_file, row[4] ?: no_file] })
    ch_versions = ch_versions.mix(SAMPLE_REPORT.out.versions)

    FHIR(
        ch_meta.join(ch_extras, failOnMismatch: true)
            .map { row -> [row[1], row[2], row[3] ?: no_file, row[4] ?: no_file] },
        ch_cgmlst,
        no_file)
    ch_versions = ch_versions.mix(FHIR.out.versions)

    def versions_header = ("\"Workflow\":\n    KPtoFHIR: ${workflow.manifest.version}\n" +
                           "    Nextflow: ${nextflow.version}\n" +
                           "    fhir_profile: ${params.fhir_profile}\n" +
                           "    plasmid_analysis: ${params.plasmid_analysis}\n" +
                           "    cgmlst_schema: ${params.cgmlst_schema ? file(params.cgmlst_schema).name : 'none'}\n" +
                           "    genome_size: ${params.genome_size}\n").toString()
    ch_versions
        .map { f -> f.text.trim() }
        .unique()
        .collectFile(name: 'versions_collected.yml', sort: true, newLine: true)
        .map { collected -> versions_header + collected.text }
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info")

    emit:
    genomic = FHIR.out.genomic
    bundles = FHIR.out.bundles
}
