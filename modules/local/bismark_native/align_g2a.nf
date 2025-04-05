process ALIGN_G2A {
    tag "${meta.id}"
    label 'process_high'
    
    conda "bioconda::bowtie2=2.4.2 bioconda::samtools=1.15.1"
    
    input:
    tuple val(meta), path(reads)
    path ga_index_dir
    
    output:
    tuple val(meta), path("*_g2a.sam"), emit: sam
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def endedness = meta.single_end ? "-U ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    def index_base = ga_index_dir.toString() + "/GA_conversion/BS_GA"
    
    // Flag to use only forward or reverse complement depending on strand being aligned
    // For G->A conversion, reverse strand (OB)
    def strand_args = "--nofw"  // Only allow alignment to reverse strand

    """
    bowtie2 \\
        -x ${index_base} \\
        ${endedness} \\
        ${strand_args} \\
        ${args} \\
        --sam-nohead \\
        -p ${task.cpus} \\
        > ${prefix}_g2a.sam
    """
} 