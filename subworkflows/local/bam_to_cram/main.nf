//
// Subworkflow to convert BAM files to CRAM format
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMTOOLS_CONVERT } from '../../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO CONVERT BAM TO CRAM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BAM_TO_CRAM {
    take:
    ch_samplesheet  // channel: [meta:[id:sample], bam_file_path]

    main:
    ch_fasta        = Channel.empty()
    ch_fasta_index  = Channel.empty()
    ch_indexed_bam  = Channel.empty()
    ch_cram         = Channel.empty()
    ch_versions     = Channel.empty()

    // Get the reference FASTA file from params
    def fasta_file_path = params.reference
    def fai_file_path   = "${params.reference}.fai"
    ch_fasta            = Channel.value([[:], file(fasta_file_path)])

    // Index the FASTA file if FAI file does not exist
    if (!file(fai_file_path).exists()) {
        SAMTOOLS_FAIDX(ch_fasta, [[:],[]])
        ch_fasta_index  = SAMTOOLS_FAIDX.out.fai
        ch_versions     = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    } else {
        ch_fasta_index  = Channel.value([[:], file(fai_file_path)])
    }

    // Index the BAM file if BAI file does not exist
    ch_bam_index = ch_samplesheet.map { meta, bam_file_path ->
        def bai_file1 = bam_file_path.replace('.bam', '.bai')  // {prefix}.bai
        def bai_file2 = bam_file_path + '.bai'                 // {prefix}.bam.bai

        def index_file = file(bai_file1).exists() ? bai_file1 :
                         file(bai_file2).exists() ? bai_file2 : null

        if (!index_file) {
            SAMTOOLS_INDEX(meta: meta, bam_file_path: bam_file_path)
            ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
            return SAMTOOLS_INDEX.out.bai
        } else {
            return [meta, file(index_file)]
        }
    }

    // Convert BAM to CRAM
    SAMTOOLS_CONVERT(
        ch_samplesheet.join(ch_bam_index),
        ch_fasta,
        ch_fasta_index
    )

    emit:
    cram        = SAMTOOLS_CONVERT.out.cram
    versions    = ch_versions
}
