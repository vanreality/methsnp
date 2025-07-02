//
// Subworkflow to preprocess bam files with samtools and revelio
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMTOOLS_CONVERT } from '../../../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_INTERSECT } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_REVELIO } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FAIDX } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_VIEW } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_INTERSECT } from '../../../modules/local/samtools/intersect/main'
include { REVELIO } from '../../../modules/local/revelio/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO PREPROCESS BAM FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BAM_PREPROCESS {
    take:
    ch_samplesheet  // channel: [meta:[id:sample], bam_file_path]

    main:
    ch_fasta        = Channel.empty()
    ch_fai          = Channel.empty()
    ch_indexed_bam  = Channel.empty()
    ch_cram         = Channel.empty()
    ch_versions     = Channel.empty()

    // Index the FASTA file if FAI file does not exist
    def fasta_file_path = params.reference
    def fai_file_path   = "${params.reference}.fai"
    ch_fasta            = Channel.value([[:], file(fasta_file_path)])

    if (!file(fai_file_path).exists()) {
        SAMTOOLS_FAIDX(ch_fasta, [[:],[]])
        ch_fai  = SAMTOOLS_FAIDX.out.fai
        ch_versions     = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    } else {
        ch_fai  = Channel.value([[:], file(fai_file_path)])
    }

    // Extract a subset of the BAM file if the region BED file is provided
    params.region = params.region ?: ''
    if (!params.region) {
        ch_bam = ch_samplesheet
    } else {
        region_file = file(params.region)
        SAMTOOLS_INTERSECT(
            ch_samplesheet,
            Channel.value([[:], region_file])
        )
        ch_bam      = SAMTOOLS_INTERSECT.out.bam
        ch_versions = ch_versions.mix(SAMTOOLS_INTERSECT.out.versions)
    }

    // Index the BAM file if BAI file does not exist
    // ch_bai = ch_bam.map { meta, bam_file_path ->
    //     def bam_file_string = bam_file_path.toString()           // Convert to String
    //     def bai_file1 = bam_file_string.replace('.bam', '.bai')  // {prefix}.bai
    //     def bai_file2 = bam_file_string + '.bai'                 // {prefix}.bam.bai

    //     def index_file = file(bai_file1).exists() ? bai_file1 :
    //                      file(bai_file2).exists() ? bai_file2 : null

    //     if (!index_file) {
    //         SAMTOOLS_INDEX(meta: meta, bam_file_path: bam_file_path)
    //         ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
    //         return SAMTOOLS_INDEX.out.bai
    //     } else {
    //         return [meta, file(index_file)]
    //     }
    // }
    SAMTOOLS_INDEX_INTERSECT(ch_bam)
    ch_bai = SAMTOOLS_INDEX_INTERSECT.out.bai

    // Revilio preprocess: Mask CtoT and GtoA possible false positive variants
    REVELIO(ch_bam.join(ch_bai), ch_fasta, ch_fai)
    ch_bam = REVELIO.out.bam

    SAMTOOLS_INDEX_REVELIO(ch_bam)
    ch_bai = SAMTOOLS_INDEX_REVELIO.out.bai

    ch_bam_bai = ch_bam.join(ch_bai)

    // Convert BAM to CRAM
    SAMTOOLS_CONVERT(
        ch_bam_bai,
        ch_fasta,
        ch_fai
    )

    ch_cram_crai = SAMTOOLS_CONVERT.out.cram.join(SAMTOOLS_CONVERT.out.crai)

    emit:
    bam        = ch_bam_bai
    cram       = ch_cram_crai
    fasta      = ch_fasta
    fai        = ch_fai
    versions   = ch_versions
}
