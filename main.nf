#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules and subworkflows
include { BISMARK_GENOME_PREPARATION_WF } from './subworkflows/bismark_genome_preparation'
// include { BISMARK_ALIGN                 } from './modules/local/bismark_align'
// include { BISMARK_DEDUPLICATE           } from './modules/local/bismark_deduplicate'
// include { BISMARK_METHYLATION_EXTRACT   } from './modules/local/bismark_methylation_extract'
// include { BISMARK_BEDGRAPH              } from './modules/local/bismark_bedgraph'
// include { BISMARK_REPORT                } from './modules/local/bismark_report'
// include { BISMARK_SUMMARY               } from './modules/local/bismark_summary'
// include { BISMARK_CYTOSINE_REPORT       } from './modules/local/bismark_cytosine_report'


workflow {
    BISMARK_GENOME_PREPARATION_WF(
        params.fasta,
        params.aligner
    )
}
