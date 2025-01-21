//
// Subworkflow to filtering variants
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TABIX_TABIX as TABIX_TABIX_DB } from '../../../modules/nf-core/tabix/main'
include { TABIX_TABIX as TABIX_TABIX_VCF } from '../../../modules/nf-core/tabix/main'
include { SNPSIFT_ANNOTATE } from '../../../modules/nf-core/snpsif/annotate/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO FILTERING VARIANTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VCF_ANNOTATE_SNPSIFT {
    take:
    ch_vcf_tbi    // [meta, vcf, tbi]

    main:
    ch_versions = Channel.empty()

    def snpsift_db_path = params.snpsift_db
    def snpsift_db_tbi_path = snpsift_db_path.replace('.gz', '.gz.tbi')
    if (file(snpsift_db_tbi_path).exists()) {
        SNPSIFT_ANNOTATE(ch_vcf_tbi, [[:], file(snpsift_db_path), file(snpsift_db_tbi_path)])
    } else {
        TABIX_TABIX_DB(snpsift_db_path)
        ch_versions = ch_versions.mix(TABIX_TABIX_DB.out.versions)

        SNPSIFT_ANNOTATE(ch_vcf_tbi, TABIX_TABIX_DB.out.gz_tbi)
    }

    ch_vcf = SNPSIFT_ANNOTATE.out.vcf
    ch_versions = ch_versions.mix(SNPSIFT_ANNOTATE.out.versions)

    TABIX_TABIX_VCF(ch_vcf)
    ch_vcf_tbi = TABIX_TABIX_VCF.out.gz_tbi
    ch_versions = ch_versions.mix(TABIX_TABIX_VCF.out.versions)

    emit:
    vcf = ch_vcf_tbi

    versions = ch_versions
}
