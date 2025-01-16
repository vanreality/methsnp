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
include { BAM_TO_CRAM } from '../subworkflows/local/bam_to_cram/main'
include { CRAM_HAPLOTYPECALLER_VARIANT_CALLING } from '../subworkflows/local/cram_haplotypecaller_variant_calling/main'
// include { VCF_FILTERING } from '../subworkflows/local/vcf_filtering/main'
include { VCF_QC_BCFTOOLS_VCFTOOLS } from '../subworkflows/local/vcf_qc_bcftools_vcftools/main'

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

    // BAM to CRAM
    BAM_TO_CRAM(ch_samplesheet)
    ch_cram  = BAM_TO_CRAM.out.cram
    ch_fasta = BAM_TO_CRAM.out.fasta
    ch_fai   = BAM_TO_CRAM.out.fai
    ch_versions = ch_versions.mix(BAM_TO_CRAM.out.versions)

    // GATK HaplotypeCaller
    CRAM_HAPLOTYPECALLER_VARIANT_CALLING(ch_cram, ch_fasta, ch_fai)
    ch_vcf = CRAM_HAPLOTYPECALLER_VARIANT_CALLING.out.vcf
    ch_tbi = CRAM_HAPLOTYPECALLER_VARIANT_CALLING.out.tbi
    ch_versions = ch_versions.mix(CRAM_HAPLOTYPECALLER_VARIANT_CALLING.out.versions)

    // Filtering
    // VCF_FILTERING(ch_vcf)

    // QC
    VCF_QC_BCFTOOLS_VCFTOOLS(ch_vcf, ch_tbi)

    // TODO: Variant Annotation

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
