/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_methsnp_pipeline'

// nf-core subworkflows

// local subworkflows
include { BAM_PREPROCESS } from '../subworkflows/local/bam_preprocess/main'
include { CRAM_HAPLOTYPECALLER_VARIANT_CALLING } from '../subworkflows/local/cram_haplotypecaller_variant_calling/main'
include { CRAM_DEEPVARIANT_VARIANT_CALLING } from '../subworkflows/local/cram_deepvariant_variant_calling/main'
//include { VCF_POSTPROCESS } from '../subworkflows/local/vcf_postprocess/main'
include { VCF_QC_BCFTOOLS_VCFTOOLS } from '../subworkflows/local/vcf_qc_bcftools_vcftools/main'
include { VCF_ANNOTATE_SNPEFF } from '../subworkflows/nf-core/vcf_annotate_snpeff/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow METHSNP {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    ch_versions = Channel.empty()

    // BAM preprocessing
    BAM_PREPROCESS(ch_samplesheet)
    ch_cram     = BAM_PREPROCESS.out.cram
    ch_fasta    = BAM_PREPROCESS.out.fasta
    ch_fai      = BAM_PREPROCESS.out.fai
    ch_versions = ch_versions.mix(BAM_PREPROCESS.out.versions)

    // GATK HaplotypeCaller
    CRAM_HAPLOTYPECALLER_VARIANT_CALLING(ch_cram, ch_fasta, ch_fai)
    ch_gatk_vcf = CRAM_HAPLOTYPECALLER_VARIANT_CALLING.out.vcf
    ch_gatk_tbi = CRAM_HAPLOTYPECALLER_VARIANT_CALLING.out.tbi
    ch_versions = ch_versions.mix(CRAM_HAPLOTYPECALLER_VARIANT_CALLING.out.versions)

    // DeepVariant
    // CRAM_DEEPVARIANT_VARIANT_CALLING(ch_cram, ch_fasta, ch_fai)
    // ch_dv_vcf   = CRAM_DEEPVARIANT_VARIANT_CALLING.out.vcf
    // ch_dv_tbi   = CRAM_DEEPVARIANT_VARIANT_CALLING.out.tbi
    // ch_versions = ch_versions.mix(CRAM_DEEPVARIANT_VARIANT_CALLING.out.versions)

    // VCF postprocessing
    // VCF_POSTPROCESS(ch_)

    // QC
    // VCF_QC_BCFTOOLS_VCFTOOLS(ch_vcf, ch_tbi)

    // Variant Annotation
    def snpeff_db = params.snpeff_db
    ch_snpeff_cache = Channel.fromPath(file("${snpeff_cache}"), checkIfExists: true)
                        .collect()
                        .map( cache -> [ [ id:"${snpeff_db}" ], cache ] )

    VCF_ANNOTATE_SNPEFF(
        ch_gatk_vcf,
        params.snpeff_db,
        ch_snpeff_cache
    )
    // VCF_ANNOTATE_SNPEFF(ch_dv_vcf, params.snpeff_db, file(params.snpeff_cache))

    // TODO: Methylation Extraction

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'methsnp_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
