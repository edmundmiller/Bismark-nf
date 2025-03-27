#!/usr/bin/env nextflow

// Import the process from the module
include { BISMARK_GENOME_PREPARATION } from '../modules/local/bismark_genome_preparation'

// Define the workflow
workflow BISMARK_GENOME_PREPARATION_WF {
    take:
    genome_dir  // Directory containing genome FASTA files
    aligner     // Aligner to use (bowtie2, hisat2, minimap2)

    main:
    // Run the bismark genome preparation process
    BISMARK_GENOME_PREPARATION(
        genome_dir,
        aligner
    )

    emit:
    bisulfite_genome = BISMARK_GENOME_PREPARATION.out.bisulfite_genome
} 