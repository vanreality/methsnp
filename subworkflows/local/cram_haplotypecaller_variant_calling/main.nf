//
// Subworkflow to call variants using GATK HaplotypeCaller from CRAM files
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GATK4_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_HAPLOTYPECALLER } from '../../../modules/nf-core/gatk4/haplotypecaller/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO CALL VARIANTS USING GATK HAPLOTYPECALLER FROM CRAM FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CRAM_HAPLOTYPECALLER_VARIANT_CALLING {
    take:
    ch_cram      // channel: [meta, cram, crai]
    ch_fasta     // channel: [meta2, fasta]
    ch_fai       // channel: [meta3, fai]

    main:
    ch_versions = Channel.empty()

    GATK4_CREATESEQUENCEDICTIONARY(ch_fasta)
    ch_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict
    ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)

    if (!params.region) {
        ch_cram = ch_cram.map{ meta, cram, crai ->
            tuple(meta, cram, crai, [], [])
        }
    } else {
        ch_cram = ch_cram.map{ meta, cram, crai ->
            tuple(meta, cram, crai, file(params.region), [])
        }
    }

    GATK4_HAPLOTYPECALLER(
        ch_cram,
        ch_fasta,
        ch_fai,
        ch_dict,
        [[:], []],
        [[:], []]
    )
    ch_vcf = GATK4_HAPLOTYPECALLER.out.vcf
    ch_tbi = GATK4_HAPLOTYPECALLER.out.tbi

    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions)

    emit:
    vcf = ch_vcf
    tbi = ch_tbi
    versions = ch_versions
}
