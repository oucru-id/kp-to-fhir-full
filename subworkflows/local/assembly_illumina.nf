include { FASTQC        } from '../../modules/local/fastqc'
include { FASTP         } from '../../modules/local/fastp'
include { MEGAHIT       } from '../../modules/local/megahit'
include { BWA_MEM       } from '../../modules/local/bwa_mem'
include { SAMTOOLS_SORT } from '../../modules/local/samtools_sort'
include { PILON         } from '../../modules/local/pilon'

workflow ASSEMBLY_ILLUMINA {
    take:
    reads

    main:
    FASTQC(reads)
    FASTP(reads)
    MEGAHIT(FASTP.out.reads)
    BWA_MEM(MEGAHIT.out.contigs.join(FASTP.out.reads, failOnMismatch: true))
    SAMTOOLS_SORT(BWA_MEM.out.sam)
    PILON(MEGAHIT.out.contigs.join(SAMTOOLS_SORT.out.bam, failOnMismatch: true))

    emit:
    assembly = PILON.out.assembly
    qc       = FASTQC.out.zip
    versions = FASTQC.out.versions.mix(FASTP.out.versions, MEGAHIT.out.versions, BWA_MEM.out.versions,
                                       SAMTOOLS_SORT.out.versions, PILON.out.versions)
}
