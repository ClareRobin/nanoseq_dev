/*
 * Differential Expression Analysis with DESeq2
 */

include { DESEQ2      } from '../../modules/local/deseq2'

workflow DIFFERENTIAL_DESEQ2 {
    take:
    ch_gene_counts

    main:
    /*
     * DESeq2 differential expression of genes
     */
    DESEQ2 ( ch_gene_counts )
    ch_deseq2_txt  = DESEQ2.out.deseq2_txt
    deseq2_version = DESEQ2.out.versions

    emit:
    ch_deseq2_txt
    deseq2_version
}
