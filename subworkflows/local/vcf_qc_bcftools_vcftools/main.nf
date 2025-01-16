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

    def region_file_path = params.region
    ch_region = Channel.fromPath(region_file_path)

    BCFTOOLS_STATS(ch_vcf.map{ meta, vcf -> [ meta, vcf, [] ] }, [[:],[]], [[:],[]], [[:],[]], [[:],[]], [[:],[]])
    VCFTOOLS_TSTV_COUNT(ch_vcf, ch_region, [])
    VCFTOOLS_TSTV_QUAL(ch_vcf, ch_region, [])
    VCFTOOLS_SUMMARY(ch_vcf, ch_region, [])

    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions)
    ch_versions = ch_versions.mix(VCFTOOLS_TSTV_COUNT.out.versions)

    emit:
    bcftools_stats          = BCFTOOLS_STATS.out.stats
    vcftools_tstv_counts    = VCFTOOLS_TSTV_COUNT.out.tstv_count
    vcftools_tstv_qual      = VCFTOOLS_TSTV_QUAL.out.tstv_qual
    vcftools_filter_summary = VCFTOOLS_SUMMARY.out.filter_summary

    versions                = ch_versions
}
