// Import the processes
include { CONVERT_READS } from '../../modules/local/bismark_native/convert_reads'
include { ALIGN_C2T } from '../../modules/local/bismark_native/align_c2t'
include { ALIGN_G2A } from '../../modules/local/bismark_native/align_g2a'
include { MERGE_ALIGNMENTS } from '../../modules/local/bismark_native/merge_alignments'

workflow BISMARK_ALIGNMENT_NATIVE {
    take:
    reads           // channel: [ val(meta), [ reads ] ]
    bisulfite_index // channel: /path/to/bisulfite_index
    fasta           // channel: /path/to/genome.fasta
    directional     // value: true/false for directional library
    
    main:
    ch_versions = Channel.empty()
    
    // 1. Convert reads in silico (C->T and G->A)
    CONVERT_READS ( reads )
    
    // 2. Align converted reads to respective genome indices
    ALIGN_C2T ( 
        CONVERT_READS.out.c2t_reads,
        bisulfite_index
    )
    
    ALIGN_G2A ( 
        CONVERT_READS.out.g2a_reads,
        bisulfite_index
    )
    
    // 3. Merge alignment results and perform methylation calling
    MERGE_ALIGNMENTS (
        ALIGN_C2T.out.sam,
        ALIGN_G2A.out.sam,
        fasta,
        directional
    )
    
    emit:
    bam      = MERGE_ALIGNMENTS.out.bam       // channel: [ val(meta), [ bam ] ]
    report   = MERGE_ALIGNMENTS.out.report    // channel: [ val(meta), [ txt ] ]
    versions = ch_versions                   // channel: [ versions.yml ]
} 