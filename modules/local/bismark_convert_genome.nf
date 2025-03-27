process BISMARK_CONVERT_GENOME {
    tag "Convert ${genome_dir.getName()}"

    input:
    path genome_dir
    path bisulfite_dir
    val multi_fasta
    val slam

    output:
    path "Bisulfite_Genome", emit: bisulfite_dir

    script:
    def use_multi_fasta = multi_fasta ? "true" : "false"
    def is_slam = slam ? "true" : "false"

    """
    # Copy the bisulfite directory structure
    cp -r ${bisulfite_dir}/* .
    
    # Run the conversion script
    bismark_genome_conversion.pl \\
        --genome_folder ${genome_dir} \\
        --output_dir Bisulfite_Genome \\
        --multi_fasta ${use_multi_fasta} \\
        --slam ${is_slam}
    
    echo "Bisulfite genome conversion completed successfully"
    """
}
