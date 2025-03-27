#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules and subworkflows
include { BISMARK_GENOME_PREPARATION_WF } from './subworkflows/bismark_genome_preparation'
include { BISMARK_ALIGN                 } from './modules/local/bismark_align'
include { BISMARK_DEDUPLICATE           } from './modules/local/bismark_deduplicate'
include { BISMARK_METHYLATION_EXTRACT   } from './modules/local/bismark_methylation_extract'
include { BISMARK_BEDGRAPH              } from './modules/local/bismark_bedgraph'
include { BISMARK_REPORT                } from './modules/local/bismark_report'
include { BISMARK_SUMMARY               } from './modules/local/bismark_summary'
include { BISMARK_CYTOSINE_REPORT       } from './modules/local/bismark_cytosine_report'

// Parameter definitions
params.genome                  = false
params.reads                   = false
params.outdir                  = "results"
params.aligner                 = "bowtie2"
params.single_end              = false
params.comprehensive           = false
params.skip_deduplication      = false
params.skip_M_bias             = false
params.cytosine_report         = false
params.genome_preparation_only = false
params.parallel                = 1
params.path_to_aligner         = false
params.single_fasta            = false
params.large_index             = false
params.genomic_composition     = false
params.slam                    = false

workflow {
    // Check required inputs
    if (!params.genome) {
        error("Genome directory path not specified! Please provide a path with --genome")
    }

    if (!params.reads && !params.genome_preparation_only) {
        error("Input reads not specified! Please provide a path with --reads")
    }

    // Create channels for inputs
    genome_dir = Channel.fromPath(params.genome, checkIfExists: true, type: 'dir')

    if (!params.genome_preparation_only) {
        reads = Channel.fromPath(params.reads, checkIfExists: true)
    }

    // Genome preparation step
    BISMARK_GENOME_PREPARATION_WF(
        genome_dir,
        params.aligner,
    )

    // Exit here if only genome preparation is required
    if (params.genome_preparation_only) {
        return null
    }

    // Alignment
    BISMARK_ALIGN(
        reads,
        BISMARK_GENOME_PREPARATION_WF.out.bisulfite_genome,
    )

    // Deduplication (optional)
    if (!params.skip_deduplication) {
        BISMARK_DEDUPLICATE(
            BISMARK_ALIGN.out.bam
        )
        bam_for_methylation = BISMARK_DEDUPLICATE.out.bam
    }
    else {
        bam_for_methylation = BISMARK_ALIGN.out.bam
    }

    // Methylation extraction
    BISMARK_METHYLATION_EXTRACT(
        bam_for_methylation,
        BISMARK_GENOME_PREPARATION_WF.out.bisulfite_genome,
        params.comprehensive,
    )

    // Generate bedGraph files
    BISMARK_BEDGRAPH(
        BISMARK_METHYLATION_EXTRACT.out.methylation_calls
    )

    // Generate cytosine report (optional)
    if (params.cytosine_report) {
        BISMARK_CYTOSINE_REPORT(
            BISMARK_BEDGRAPH.out.bedgraph,
            BISMARK_GENOME_PREPARATION_WF.out.bisulfite_genome,
        )
    }

    // Collect all reports for summary
    reports_for_summary = Channel.empty()
    reports_for_summary = reports_for_summary.mix(
        BISMARK_ALIGN.out.report
    )

    if (!params.skip_deduplication) {
        reports_for_summary = reports_for_summary.mix(
            BISMARK_DEDUPLICATE.out.report
        )
    }

    reports_for_summary = reports_for_summary.mix(
        BISMARK_METHYLATION_EXTRACT.out.report
    )

    // Generate HTML reports
    BISMARK_REPORT(
        reports_for_summary.collect()
    )

    // Generate HTML summary
    BISMARK_SUMMARY(
        reports_for_summary.collect()
    )
}
