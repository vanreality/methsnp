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
include { VCFTOOLS as VCFTOOLS_TSTV_COUNT } from '../../../modules/nf-core/vcftools/main'
include { VCFTOOLS as VCFTOOLS_TSTV_QUAL  } from '../../../modules/nf-core/vcftools/main'

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

    ch_vcftools_tstv_counts    = Channel.empty()
    ch_vcftools_tstv_qual      = Channel.empty()
    ch_vcftools_filter_summary = Channel.empty()

    // bcftools
    BCFTOOLS_STATS(ch_vcf.map{ meta, vcf -> [ meta, vcf, [] ] }, [[:],[]], [[:],[]], [[:],[]], [[:],[]], [[:],[]])

    // vcftools
    def region_file_path = params.region
    if (file(region_file_path).exists()) {
        VCFTOOLS_TSTV_COUNT(ch_vcf, file(region_file_path), [])
        VCFTOOLS_TSTV_QUAL(ch_vcf, file(region_file_path), [])
        VCFTOOLS_SUMMARY(ch_vcf, file(region_file_path), [])

        ch_vcftools_tstv_counts = VCFTOOLS_TSTV_COUNT.out.tstv_count
        ch_vcftools_tstv_qual = VCFTOOLS_TSTV_QUAL.out.tstv_qual
        ch_vcftools_filter_summary = VCFTOOLS_SUMMARY.out.filter_summary
    }

    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions)
    ch_versions = ch_versions.mix(VCFTOOLS_TSTV_COUNT.out.versions)

    emit:
    bcftools_stats          = BCFTOOLS_STATS.out.stats
    vcftools_tstv_counts    = ch_vcftools_tstv_counts
    vcftools_tstv_qual      = ch_vcftools_tstv_qual
    vcftools_filter_summary = ch_vcftools_filter_summary

    versions                = ch_versions
}
