process MERGE_ALIGNMENTS {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}/bismark_alignments", mode: 'copy'
    
    conda "bioconda::samtools=1.15.1 conda-forge::python=3.9.5 bioconda::pysam=0.18.0"
    
    input:
    tuple val(meta), path(c2t_sam)
    path g2a_sam
    path genome_fasta
    val directional
    
    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*_report.txt"), emit: report
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def directional_flag = directional ? '--directional' : ''
    
    """
    # Run the merge alignments script
    merge_alignments.py \\
        --c2t_sam ${c2t_sam} \\
        --g2a_sam ${g2a_sam} \\
        --genome_fasta ${genome_fasta} \\
        --output_prefix ${prefix} \\
        ${directional_flag}
    
    # Name-sort for downstream deduplicate_bismark (requires paired reads adjacent)
    samtools sort -n -o ${prefix}.sorted.bam ${prefix}.bam
    mv ${prefix}.sorted.bam ${prefix}.bam
    
    # Report completion
    echo "Bismark native alignment completed successfully"
    """
} 