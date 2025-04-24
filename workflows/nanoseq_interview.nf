/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Check input path parameters to see if they exist
checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters (missing protocol or profile will exit the run.)
if (params.input) {
    ch_input = file(params.input)
} else {
    exit 1, 'Input samplesheet not specified!'
}

// demultiplexing not needed in this version of the pipeline
ch_input_path = 'not_needed'

// Function to check if running offline
def isOffline() {
    try {
        return NXF_OFFLINE as Boolean
    }
    catch( Exception e ) {
        return false
    }
}

if (params.protocol != 'DNA' && params.protocol != 'cDNA' && params.protocol != 'directRNA') {
    exit 1, "Invalid protocol option: ${params.protocol}. Valid options: 'DNA', 'cDNA', 'directRNA'"
}


if (!params.skip_alignment) {
    if (params.aligner != 'minimap2' && params.aligner != 'graphmap2') {
        exit 1, "Invalid aligner option: ${params.aligner}. Valid options: 'minimap2', 'graphmap2'"
    }
    if (params.protocol != 'DNA' && params.protocol != 'cDNA' && params.protocol != 'directRNA') {
        exit 1, "Invalid protocol option: ${params.protocol}. Valid options: 'DNA', 'cDNA', 'directRNA'"
    }
}


if (!params.skip_quantification) {
    if (params.quantification_method != 'bambu' && params.quantification_method != 'stringtie2') {
        exit 1, "Invalid transcript quantification option: ${params.quantification_method}. Valid options: 'bambu', 'stringtie2'"
    }
    if (params.protocol != 'cDNA' && params.protocol != 'directRNA') {
        exit 1, "Invalid protocol option if performing quantification: ${params.protocol}. Valid options: 'cDNA', 'directRNA'"
    }
}

////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

ch_multiqc_config        = file("$baseDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

include { GET_TEST_DATA         } from '../modules/local/get_test_data'
include { GET_NANOLYSE_FASTA    } from '../modules/local/get_nanolyse_fasta'
include { QCAT                  } from '../modules/local/qcat'
include { BAM_RENAME            } from '../modules/local/bam_rename'
include { BAMBU                 } from '../modules/local/bambu'
include { MULTIQC               } from '../modules/local/multiqc'

/*
 * SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
 */

include { INPUT_CHECK                      } from '../subworkflows/local/input_check'
include { PREPARE_GENOME                   } from '../subworkflows/local/prepare_genome'
include { QCFASTQ_NANOPLOT_FASTQC          } from '../subworkflows/local/qcfastq_nanoplot_fastqc'
include { ALIGN_GRAPHMAP2                  } from '../subworkflows/local/align_graphmap2'
include { ALIGN_MINIMAP2                   } from '../subworkflows/local/align_minimap2'
include { BAM_SORT_INDEX_SAMTOOLS          } from '../subworkflows/local/bam_sort_index_samtools'
include { STRUCTURAL_VARIANT_CALLING       } from '../subworkflows/local/structural_variant_calling'
include { BEDTOOLS_UCSC_BIGWIG             } from '../subworkflows/local/bedtools_ucsc_bigwig'
include { BEDTOOLS_UCSC_BIGBED             } from '../subworkflows/local/bedtools_ucsc_bigbed'
include { QUANTIFY_STRINGTIE_FEATURECOUNTS } from '../subworkflows/local/quantify_stringtie_featurecounts'
include { DIFFERENTIAL_DESEQ2              } from '../subworkflows/local/differential_deseq2'
include { RNA_MODIFICATION_XPORE_M6ANET    } from '../subworkflows/local/rna_modifications_xpore_m6anet'
include { RNA_FUSIONS_JAFFAL               } from '../subworkflows/local/rna_fusions_jaffal'

////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////

/*
 * MODULE: Installed directly from nf-core/modules
 */
include { NANOLYSE                    } from '../modules/nf-core/nanolyse/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
 * SUBWORKFLOW: Consisting entirely of nf-core/modules
 */

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

// Info required for completion email and summary
def multiqc_report      = []

workflow NANOSEQ_INTERVIEW{
    /*
     * Check parameters: make sure required parameters (input sample sheet path & RNA Seq data type e.g. cDNA or direct RNA) are provided
     */

    if (!params.input) {
        exit 1, "Input sample sheet not specified! Please provide the --input parameter."
    }

    if (!params.protocol || !(params.protocol == 'cDNA' || params.protocol == 'directRNA')) {
        exit 1, "Protocol must be specified as either 'cDNA' or 'directRNA'! Please use the --protocol parameter."
    }

    /*
     * Create empty software versions channel to collect information on software versions as the pipeline runs
     */
    ch_software_versions = Channel.empty()

    /*
     * SUBWORKFLOW: Read in samplesheet, validate and stage input files
     */
    INPUT_CHECK ( ch_input, ch_input_path )
        .set { ch_sample }

    // Extract FASTQ files from the sample channel
    ch_sample
        .map { it -> if (it[6].toString().endsWith('.gz')) [ it[0], it[6], it[2], it[1], it[4], it[5] ] }
        .set { ch_fastq }
    
    /*
     * SUBWORKFLOW: Fastq QC with Nanoplot and fastqc
     */
    
    if (!params.skip_qc) {
        QCFASTQ_NANOPLOT_FASTQC ( ch_fastq, params.skip_nanoplot, params.skip_fastqc)
        ch_software_versions = ch_software_versions.mix(QCFASTQ_NANOPLOT_FASTQC.out.fastqc_version.first().ifEmpty(null))
        ch_fastqc_multiqc    = QCFASTQ_NANOPLOT_FASTQC.out.fastqc_multiqc.ifEmpty([])
    }

    /*
     * SUBWORKFLOW: Make chromosome size file and covert GTF to BED12
     */
    PREPARE_GENOME ( ch_fastq )
    ch_fasta_index = PREPARE_GENOME.out.ch_fasta_index
    ch_gtf_bed     = PREPARE_GENOME.out.ch_gtf_bed
    ch_fasta       = PREPARE_GENOME.out.ch_fasta
    ch_fai         = PREPARE_GENOME.out.ch_fai
    ch_software_versions = ch_software_versions.mix(PREPARE_GENOME.out.samtools_version.first().ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(PREPARE_GENOME.out.gtf2bed_version.first().ifEmpty(null))

    /*
     * SUBWORKFLOW: Align fastq files with minimap2 and sort bam files
     */
    ALIGN_MINIMAP2 ( ch_fasta_index, ch_fastq )
    ch_align_sam = ALIGN_MINIMAP2.out.ch_align_sam
    ch_index = ALIGN_MINIMAP2.out.ch_index
    ch_software_versions = ch_software_versions.mix(ALIGN_MINIMAP2.out.minimap2_version.first().ifEmpty(null))

    /*
     * SUBWORKFLOW: View, then  sort, and index bam files
     */
    BAM_SORT_INDEX_SAMTOOLS ( ch_align_sam, params.call_variants, ch_fasta )
    ch_view_sortbam = BAM_SORT_INDEX_SAMTOOLS.out.sortbam
    ch_software_versions = ch_software_versions.mix(BAM_SORT_INDEX_SAMTOOLS.out.samtools_versions.first().ifEmpty(null))
    ch_samtools_multiqc  = BAM_SORT_INDEX_SAMTOOLS.out.sortbam_stats_multiqc.ifEmpty([])

    /*
     * Transcript discovery and quantification with Bambu
     */

    // Prepare channels for bambu
    ch_view_sortbam
        .map { it -> [ it[0], it[3] ] }
        .set { ch_sortbam }

    // Check that reference genome and annotation are the same for all samples if perfoming quantification
    // Check if we have replicates and multiple conditions in the input samplesheet
    REPLICATES_EXIST    = false
    MULTIPLE_CONDITIONS = false
    ch_sample.map{ it[2] }.unique().toList().set { fastas }
    ch_sample.map{ it[3] }.unique().toList().set { gtfs }

    ch_r_version = Channel.empty()

    ch_sample
        .map { it -> [ it[2], it[3] ]}
        .unique()
        .set { ch_sample_annotation }
    
    BAMBU ( ch_sample_annotation, ch_sortbam.collect{ it [1] } )
    ch_gene_counts       = BAMBU.out.ch_gene_counts
    ch_transcript_counts = BAMBU.out.ch_transcript_counts
    ch_software_versions = ch_software_versions.mix(BAMBU.out.versions.first().ifEmpty(null))

    /*
     * Differential gene expression analysis with DESeq2
     */

    DIFFERENTIAL_DESEQ2(ch_gene_counts)
    ch_software_versions = ch_software_versions.mix(DIFFERENTIAL_DESEQ2.out.deseq2_version.first().ifEmpty(null))

    /*
     * MODULE: Parse software version numbers
     */
    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_software_versions.unique().collectFile()
    )
    
    workflow_summary = WorkflowNanoseq.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    MULTIQC(
        ch_multiqc_config,
        ch_multiqc_custom_config.collect().ifEmpty([]),
        [],  // No FastQC outputs
        ch_samtools_multiqc.collect().ifEmpty([]),
        [],  // No featureCounts gene multiqc
        [],  // No featureCounts transcript multiqc
        CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )
    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
