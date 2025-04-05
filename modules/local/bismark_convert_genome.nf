process BISMARK_CONVERT_GENOME {
    tag "${fasta.getName()}"

    input:
    path fasta
    val multi_fasta
    val slam

    output:
    path "Bisulfite_Genome", emit: bisulfite_dir

    script:
    def use_multi_fasta = multi_fasta ? "true" : "false"
    def is_slam = slam ? "true" : "false"

    """
    # Create the output directory
    mkdir -p Bisulfite_Genome
    
    # Run the conversion script with appropriate file path handling
    bismark_genome_conversion.pl \\
        --genome_folder="\$PWD/${fasta}" \\
        --output_dir="Bisulfite_Genome" \\
        --multi_fasta=${use_multi_fasta} \\
        --slam=${is_slam}
    
    echo "Bisulfite genome conversion completed successfully"
    """
}
