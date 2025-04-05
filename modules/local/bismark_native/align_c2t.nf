process ALIGN_C2T {
    tag "${meta.id}"
    label 'process_high'
    
    conda "bioconda::bowtie2=2.4.2 bioconda::samtools=1.15.1"
    
    input:
    tuple val(meta), path(reads)
    path ct_index_dir
    
    output:
    tuple val(meta), path("*_c2t.sam"), emit: sam
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def endedness = meta.single_end ? "-U ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
    def index_base = ct_index_dir + "/CT_conversion/BS_CT"
    
    // Flag to use only forward or reverse complement depending on strand being aligned
    // For C->T conversion, forward strands (OT)
    def strand_args = "--norc"  // Only allow alignment to forward strand

    """
    bowtie2 \\
        -x ${index_base} \\
        ${endedness} \\
        ${strand_args} \\
        ${args} \\
        --sam-nohead \\
        -p ${task.cpus} \\
        > ${prefix}_c2t.sam
    """
} 