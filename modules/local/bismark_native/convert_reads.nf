process CONVERT_READS {
    tag "${meta.id}"
    label 'process_medium'
    
    conda "conda-forge::perl=5.26.2"
    
    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("*c2t.fastq"), emit: c2t_reads
    tuple val(meta), path("*g2a.fastq"), emit: g2a_reads
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    if (meta.single_end) {
        """
        # Convert C to T for forward alignment
        cat ${reads[0]} | awk '
        BEGIN { offset=0; }
        {
            if (NR % 4 == 2) {
                gsub(/C/,"T");
            }
            print \$0;
        }' > ${prefix}_c2t.fastq
        
        # Convert G to A for reverse alignment
        cat ${reads[0]} | awk '
        BEGIN { offset=0; }
        {
            if (NR % 4 == 2) {
                gsub(/G/,"A");
            }
            print \$0;
        }' > ${prefix}_g2a.fastq
        """
    } else {
        """
        # Convert read 1 (C to T)
        cat ${reads[0]} | awk '
        BEGIN { offset=0; }
        {
            if (NR % 4 == 2) {
                gsub(/C/,"T");
            }
            print \$0;
        }' > ${prefix}_1_c2t.fastq
        
        # Convert read 2 (G to A) for complementary strand
        cat ${reads[1]} | awk '
        BEGIN { offset=0; }
        {
            if (NR % 4 == 2) {
                gsub(/G/,"A");
            }
            print \$0;
        }' > ${prefix}_2_g2a.fastq
        
        # Also create reverse conversions for comprehensive alignment
        # Convert read 1 (G to A)
        cat ${reads[0]} | awk '
        BEGIN { offset=0; }
        {
            if (NR % 4 == 2) {
                gsub(/G/,"A");
            }
            print \$0;
        }' > ${prefix}_1_g2a.fastq
        
        # Convert read 2 (C to T)
        cat ${reads[1]} | awk '
        BEGIN { offset=0; }
        {
            if (NR % 4 == 2) {
                gsub(/C/,"T");
            }
            print \$0;
        }' > ${prefix}_2_c2t.fastq
        """
    }
} 