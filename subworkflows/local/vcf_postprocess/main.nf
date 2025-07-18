//
// Subworkflow to filtering variants
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BCFTOOLS_MERGE } from '../../../modules/nf-core/bcftools/merge/main'
include { BCFTOOLS_FILTER } from '../../../modules/nf-core/bcftools/filter/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO FILTERING VARIANTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VCF_POSTPROCESS {
    take:
    ch_vcf_tbi    // [meta, [gatk4_vcf, deepvariant_vcf], [gatk4_tbi, deepvariant_tbi]]
    ch_fasta      // [meta, fasta]
    ch_fai        // [meta, fai]

    main:
    ch_versions = Channel.empty()

    if (!params.region) {
        BCFTOOLS_MERGE(ch_vcf_tbi, ch_fasta, ch_fai, [[:], []])
    } else {
        BCFTOOLS_MERGE(ch_vcf_tbi, ch_fasta, ch_fai, [[:], file(params.region)])
    }

    ch_merged_vcf_tbi = BCFTOOLS_MERGE.out.vcf.join(BCFTOOLS_MERGE.out.index)
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions)

    BCFTOOLS_FILTER(ch_merged_vcf_tbi)
    ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions)

    emit:
    vcf = BCFTOOLS_FILTER.out.vcf
    tbi = BCFTOOLS_FILTER.out.tbi

    versions = ch_versions
}
