#!/usr/bin/env nextflow

// Import modules
include { BISMARK_PREPARE_FOLDERS } from '../modules/local/bismark_prepare_folders'
include { BISMARK_CONVERT_GENOME } from '../modules/local/bismark_convert_genome'
include { BISMARK_INDEX_GENOME } from '../modules/local/bismark_index_genome'

// Define the workflow
workflow BISMARK_GENOME_PREPARATION_WF {
    take:
    genome_dir  // Directory containing genome FASTA files
    aligner     // Aligner to use (bowtie2, hisat2, minimap2)

    main:
    // Step 1: Create necessary folders
    BISMARK_PREPARE_FOLDERS(genome_dir)
    
    // Step 2: Create the conversion script
    CREATE_CONVERSION_SCRIPT()
    
    // Step 3: Convert genomes (C->T and G->A)
    BISMARK_CONVERT_GENOME(
        genome_dir,
        BISMARK_PREPARE_FOLDERS.out.bisulfite_dir,
        !params.single_fasta,  // multi_fasta is the inverse of single_fasta
        params.slam,
    )
    
    // Step 4: Index the converted genomes
    BISMARK_INDEX_GENOME(
        BISMARK_CONVERT_GENOME.out.bisulfite_dir,
        aligner,
        params.parallel,
        params.path_to_aligner,
        params.large_index
    )

    emit:
    bisulfite_genome = BISMARK_INDEX_GENOME.out.bisulfite_genome
} 