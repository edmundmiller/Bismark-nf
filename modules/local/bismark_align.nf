process BISMARK_ALIGN {
    tag "${meta.id}"
    
    publishDir "${params.outdir}/bismark_alignments", mode: 'copy'
    
    input:
    tuple val(meta), path(reads)
    path bismark_index
    
    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "*.txt", emit: report
    
    script:
    def args = ""
    def prefix = meta.id
    
    // Determine if single-end or paired-end
    def input_files = meta.single_end ? "-i ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    
    // Add aligner-specific options
    if (params.aligner == "bowtie2") {
        args += " --bowtie2"
    } else if (params.aligner == "hisat2") {
        args += " --hisat2"
    } else if (params.aligner == "minimap2" || params.aligner == "mm2") {
        args += " --minimap2"
    }
    
    // Add other common options
    if (params.path_to_aligner) {
        args += " --path_to_aligner ${params.path_to_aligner}"
    }
    
    if (params.parallel && params.parallel > 1) {
        args += " --parallel ${params.parallel}"
    }
    
    if (!params.directional) {
        args += " --non_directional"
    }
    
    if (params.pbat) {
        args += " --pbat"
    }
    
    if (params.unmapped) {
        args += " --unmapped"
    }
    
    if (params.ambiguous) {
        args += " --ambiguous"
    }
    
    // Add output format options
    args += " --bam"
    
    // Add prefix if specified
    if (prefix) {
        args += " --prefix ${prefix}"
    }
    
    """
    # Run the bismark alignment command
    bismark ${args} --genome ${bismark_index} ${input_files}
    
    # Check the exit status
    if [ \$? -ne 0 ]; then
        echo "Error: bismark alignment failed"
        exit 1
    fi
    
    echo "Bismark alignment completed successfully"
    """
} 