process BISMARK_PREPARE_FOLDERS {
    tag "${fasta.getName()}"
    
    input:
    path fasta
    
    output:
    path "Bisulfite_Genome", emit: bisulfite_dir

    script:
    """
    # Create the main Bisulfite_Genome directory
    mkdir -p Bisulfite_Genome/CT_conversion
    mkdir -p Bisulfite_Genome/GA_conversion
    
    # Check if the genome directory contains FASTA files
    if [ ! \$(find ${genome_dir} -name "*.fa" -o -name "*.fa.gz" -o -name "*.fasta" -o -name "*.fasta.gz" | wc -l) -gt 0 ]; then
        echo "Error: No FASTA files found in ${genome_dir}"
        exit 1
    fi
    
    echo "Bisulfite genome folders created successfully"
    """
} 