process BISMARK_GENOME_PREPARATION {
    tag "Bisulfite genome preparation"
    
    publishDir "${params.outdir}/bismark_genome", mode: 'copy'
    
    input:
    path genome_dir
    val aligner
    
    output:
    path "Bisulfite_Genome", emit: bisulfite_genome
    
    script:
    def args = ""
    
    // Add parallel option if specified
    if (params.parallel && params.parallel > 1) {
        args += " --parallel ${params.parallel}"
    }
    
    // Add aligner-specific options
    if (aligner == "bowtie2") {
        args += " --bowtie2"
    } else if (aligner == "hisat2") {
        args += " --hisat2"
    } else if (aligner == "minimap2" || aligner == "mm2") {
        args += " --minimap2"
    }
    
    // Add other options if specified
    if (params.path_to_aligner) {
        args += " --path_to_aligner ${params.path_to_aligner}"
    }
    
    if (params.single_fasta) {
        args += " --single_fasta"
    }
    
    if (params.large_index) {
        args += " --large-index"
    }
    
    if (params.genomic_composition) {
        args += " --genomic_composition"
    }
    
    if (params.slam) {
        args += " --slam"
    }
    
    """
    # Create a temporary directory for the output if it doesn't exist
    mkdir -p Bisulfite_Genome
    
    # Run the bismark_genome_preparation command
    bismark_genome_preparation ${args} ${genome_dir}
    
    # Check the exit status
    if [ \$? -ne 0 ]; then
        echo "Error: bismark_genome_preparation failed"
        exit 1
    fi
    
    echo "Bisulfite genome preparation completed successfully"
    """
} 