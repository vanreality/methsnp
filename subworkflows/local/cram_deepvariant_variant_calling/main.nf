//
// Subworkflow to call variants using Deepvariant from CRAM files
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DEEPVARIANT_RUNDEEPVARIANT } from '../../../modules/nf-core/deepvariant/rundeepvariant/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO CALL VARIANTS USING DEEPVARIANT FROM CRAM FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CRAM_DEEPVARIANT_VARIANT_CALLING {
    take:
    ch_cram      // channel: [meta, cram, crai]
    ch_fasta     // channel: [meta2, fasta]
    ch_fai       // channel: [meta3, fai]

    main:
    ch_versions = Channel.empty()

    if (!params.region) {
        ch_cram = ch_cram.map{ meta, cram, crai ->
            tuple(meta, cram, crai, [])
        }
    } else {
        ch_cram = ch_cram.map{ meta, cram, crai ->
            tuple(meta, cram, crai, file(params.region))
        }
    }

    DEEPVARIANT_RUNDEEPVARIANT(
        ch_cram,
        ch_fasta,
        ch_fai,
        [[:], []],
        [[:], []]
    )
    ch_vcf = DEEPVARIANT_RUNDEEPVARIANT.out.vcf
    ch_tbi = DEEPVARIANT_RUNDEEPVARIANT.out.vcf_tbi

    ch_versions = ch_versions.mix(DEEPVARIANT_RUNDEEPVARIANT.out.versions)

    emit:
    vcf = ch_vcf
    tbi = ch_tbi
    versions = ch_versions
}
