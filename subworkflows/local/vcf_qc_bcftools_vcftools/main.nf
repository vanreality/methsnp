//
// Subworkflow to QC VCF files using bcftools and vcftools
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BCFTOOLS_STATS                  } from '../../../modules/nf-core/bcftools/stats/main'
include { VCFTOOLS as VCFTOOLS_SUMMARY    } from '../../../modules/nf-core/vcftools/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO QC VCF FILES USING BCFTOOLS AND VCFTOOLS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VCF_QC_BCFTOOLS_VCFTOOLS {
    take:
    ch_vcf
    ch_tbi

    main:
    ch_versions = Channel.empty()

    // bcftools
    BCFTOOLS_STATS(ch_vcf.map{ meta, vcf -> [ meta, vcf, [] ] }, [[:],[]], [[:],[]], [[:],[]], [[:],[]], [[:],[]])

    // vcftools
    if (params.region) {
        VCFTOOLS_SUMMARY(ch_vcf, file(params.region), [])
        ch_vcftools_filter_summary = VCFTOOLS_SUMMARY.out.filter_summary
        ch_versions = ch_versions.mix(VCFTOOLS_SUMMARY.out.versions)
    } else {
        ch_vcftools_filter_summary = Channel.empty()
    }

    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions)

    emit:
    bcftools_stats          = BCFTOOLS_STATS.out.stats
    vcftools_filter_summary = ch_vcftools_filter_summary

    versions                = ch_versions
}
