process BISMARK_INDEX_GENOME {
    tag "Index ${bisulfite_dir.getName()} with ${aligner}"
    
    input:
    path bisulfite_dir
    val aligner
    val parallel
    val path_to_aligner
    val large_index
    
    output:
    path "Bisulfite_Genome", emit: bisulfite_genome
    
    script:
    def cores = parallel > 1 ? parallel : 1
    def aligner_path = path_to_aligner ? "--path_to_aligner ${path_to_aligner}" : ""
    def index_args = large_index ? "--large-index" : ""
    
    // Choose indexer based on aligner
    def indexer_cmd = ""
    if (aligner == "bowtie2") {
        if (path_to_aligner) {
            indexer_cmd = "${path_to_aligner}/bowtie2-build"
        } else {
            indexer_cmd = "bowtie2-build"
        }
        index_args += " --threads ${cores}"
    } else if (aligner == "hisat2") {
        if (path_to_aligner) {
            indexer_cmd = "${path_to_aligner}/hisat2-build"
        } else {
            indexer_cmd = "hisat2-build"
        }
        index_args += " --threads ${cores}"
    } else if (aligner == "minimap2" || aligner == "mm2") {
        if (path_to_aligner) {
            indexer_cmd = "${path_to_aligner}/minimap2"
        } else {
            indexer_cmd = "minimap2"
        }
        index_args = "-k 20 -t ${cores}"
    }
    
    """
    # Copy the bisulfite directory structure
    cp -r ${bisulfite_dir}/* .
    
    # Index CT converted genome
    cd Bisulfite_Genome/CT_conversion
    FASTA_FILES=\$(ls *.fa)
    
    if [ "${aligner}" == "minimap2" ] || [ "${aligner}" == "mm2" ]; then
        ${indexer_cmd} ${index_args} -d BS_CT.mmi \$FASTA_FILES
    else
        ${indexer_cmd} ${index_args} -f \$FASTA_FILES BS_CT
    fi
    
    # Index GA converted genome
    cd ../GA_conversion
    FASTA_FILES=\$(ls *.fa)
    
    if [ "${aligner}" == "minimap2" ] || [ "${aligner}" == "mm2" ]; then
        ${indexer_cmd} ${index_args} -d BS_GA.mmi \$FASTA_FILES
    else
        ${indexer_cmd} ${index_args} -f \$FASTA_FILES BS_GA
    fi
    
    cd ../..
    
    echo "Bisulfite genome indexing completed successfully"
    """
} 