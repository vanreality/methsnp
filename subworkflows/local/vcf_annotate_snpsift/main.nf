//
// Subworkflow to filtering variants
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TABIX_TABIX } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_BGZIPTABIX } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { SNPSIFT_ANNOTATE } from '../../../modules/nf-core/snpsift/annotate/main'

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
        TABIX_TABIX([[:], file(snpsift_db_path)])
        ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

        ch_db_tbi = TABIX_TABIX.out.tbi.map {
            meta, tbi ->
            return [meta, file(snpsift_db_path), tbi]
        }

        SNPSIFT_ANNOTATE(ch_vcf_tbi, ch_db_tbi)
    }

    ch_vcf = SNPSIFT_ANNOTATE.out.vcf
    ch_versions = ch_versions.mix(SNPSIFT_ANNOTATE.out.versions)

    TABIX_BGZIPTABIX(ch_vcf)
    ch_vcf_tbi = TABIX_BGZIPTABIX.out.gz_tbi
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)

    emit:
    vcf = ch_vcf_tbi

    versions = ch_versions
}
