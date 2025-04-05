#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules and subworkflows
include { BISMARK_GENOME_PREPARATION_WF } from './subworkflows/bismark_genome_preparation'
include { BISMARK_ALIGNMENT_NATIVE      } from './subworkflows/local/bismark_alignment_native'
// include { BISMARK_DEDUPLICATE           } from './modules/local/bismark_deduplicate'
// include { BISMARK_METHYLATION_EXTRACT   } from './modules/local/bismark_methylation_extract'
// include { BISMARK_BEDGRAPH              } from './modules/local/bismark_bedgraph'
// include { BISMARK_REPORT                } from './modules/local/bismark_report'
// include { BISMARK_SUMMARY               } from './modules/local/bismark_summary'
// include { BISMARK_CYTOSINE_REPORT       } from './modules/local/bismark_cytosine_report'

workflow {
    // Create a channel for input reads
    input_reads = Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .map { group_id, files -> [[id: group_id, single_end: files.size() == 1], files] }

    // Run the genome preparation workflow
    BISMARK_GENOME_PREPARATION_WF(
        params.fasta,
        params.aligner,
    )

    // Run the native implementation
    BISMARK_ALIGNMENT_NATIVE(
        input_reads,
        BISMARK_GENOME_PREPARATION_WF.out.bisulfite_genome,
        params.fasta,
        params.directional,
    )

    // Store output channels for downstream processes
    aligned_bams = BISMARK_ALIGNMENT_NATIVE.out.bam
    alignment_reports = BISMARK_ALIGNMENT_NATIVE.out.report
}
