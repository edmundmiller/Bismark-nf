#!/usr/bin/env python3
import sys
import re
import pysam
import os
import gzip
import argparse

def parse_sam(sam_file, conversion_type):
    # Parse alignment file and extract information
    alignments = {}
    with open(sam_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue
            fields = line.strip().split('\t')
            read_id = fields[0]
            flag = int(fields[1])
            chrom = fields[2]
            pos = int(fields[3])
            mapq = int(fields[4])
            cigar = fields[5]
            seq = fields[9]
            
            # Skip unmapped reads
            if chrom == '*' or flag & 4:
                continue
                
            # Store alignment info
            if read_id not in alignments:
                alignments[read_id] = {}
            
            alignments[read_id][conversion_type] = {
                'flag': flag,
                'chrom': chrom,
                'pos': pos,
                'mapq': mapq,
                'cigar': cigar,
                'seq': seq,
                'line': line
            }
    return alignments

def extract_genomic_sequence(chrom, start, cigar, genome_dict):
    # Extract reference sequence
    if chrom not in genome_dict:
        return None
        
    # Parse CIGAR string to determine reference sequence
    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    
    genomic_seq = ""
    ref_pos = start - 1  # 0-based
    
    for count, op in ops:
        count = int(count)
        if op in 'M=X':  # Match or mismatch
            genomic_seq += genome_dict[chrom][ref_pos:ref_pos+count]
            ref_pos += count
        elif op == 'D':  # Deletion
            ref_pos += count
        elif op == 'I':  # Insertion
            # No reference bases consumed
            pass
        elif op == 'S':  # Soft clip
            # No reference bases consumed
            pass
        elif op == 'H':  # Hard clip
            # No reference bases consumed
            pass
        elif op == 'N':  # Skipped region
            ref_pos += count
    
    return genomic_seq

def call_methylation(read_seq, genomic_seq, conversion_type):
    # Call methylation status
    methylation_string = []
    
    if len(read_seq) != len(genomic_seq):
        return '.' * len(read_seq)
        
    for i in range(len(read_seq)):
        if conversion_type == 'c2t':
            if genomic_seq[i].upper() == 'C':
                if i+2 < len(genomic_seq) and genomic_seq[i+1].upper() == 'G':
                    # CpG context
                    if read_seq[i].upper() == 'C':
                        methylation_string.append('Z')  # Methylated CpG
                    else:
                        methylation_string.append('z')  # Unmethylated CpG
                elif i+2 < len(genomic_seq) and genomic_seq[i+1].upper() in ['A', 'C', 'T']:
                    # CHG context
                    if read_seq[i].upper() == 'C':
                        methylation_string.append('X')  # Methylated CHG
                    else:
                        methylation_string.append('x')  # Unmethylated CHG
                else:
                    # CHH context
                    if read_seq[i].upper() == 'C':
                        methylation_string.append('H')  # Methylated CHH
                    else:
                        methylation_string.append('h')  # Unmethylated CHH
            else:
                methylation_string.append('.')  # Not a C position
        elif conversion_type == 'g2a':
            if genomic_seq[i].upper() == 'G':
                if i > 0 and genomic_seq[i-1].upper() == 'C':
                    # CpG context (reverse strand)
                    if read_seq[i].upper() == 'G':
                        methylation_string.append('Z')  # Methylated CpG
                    else:
                        methylation_string.append('z')  # Unmethylated CpG
                elif i > 0 and genomic_seq[i-1].upper() in ['A', 'C', 'T']:
                    # CHG context (reverse strand)
                    if read_seq[i].upper() == 'G':
                        methylation_string.append('X')  # Methylated CHG
                    else:
                        methylation_string.append('x')  # Unmethylated CHG
                else:
                    # CHH context (reverse strand)
                    if read_seq[i].upper() == 'G':
                        methylation_string.append('H')  # Methylated CHH
                    else:
                        methylation_string.append('h')  # Unmethylated CHH
            else:
                methylation_string.append('.')  # Not a G position
    
    return ''.join(methylation_string)

def load_genome(genome_fasta):
    # Load genome into memory
    genome_dict = {}
    current_chrom = None
    
    # Check if file is gzipped
    is_gzipped = genome_fasta.endswith('.gz')
    opener = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'
    
    with opener(genome_fasta, mode) as f:
        for line in f:
            if line.startswith('>'):
                current_chrom = line.strip().split()[0][1:]
                genome_dict[current_chrom] = ''
            else:
                genome_dict[current_chrom] += line.strip()
    
    return genome_dict

def parse_args():
    parser = argparse.ArgumentParser(description='Merge Bismark alignments and call methylation')
    parser.add_argument('--c2t_sam', required=True, help='C->T converted alignments SAM file')
    parser.add_argument('--g2a_sam', required=True, help='G->A converted alignments SAM file')
    parser.add_argument('--genome_fasta', required=True, help='Reference genome FASTA file')
    parser.add_argument('--output_prefix', required=True, help='Prefix for output files')
    parser.add_argument('--directional', action='store_true', help='Whether library is directional')
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Load genome
    genome_dict = load_genome(args.genome_fasta)
    
    # Parse alignment files
    c2t_alignments = parse_sam(args.c2t_sam, 'c2t')
    g2a_alignments = parse_sam(args.g2a_sam, 'g2a')
    
    # Merge alignments
    all_read_ids = set(list(c2t_alignments.keys()) + list(g2a_alignments.keys()))
    
    # Output BAM file
    outfile = pysam.AlignmentFile(f"{args.output_prefix}.bam", "wb", reference_names=list(genome_dict.keys()), reference_lengths=[len(seq) for seq in genome_dict.values()])
    
    # Statistics for report
    total_reads = len(all_read_ids)
    mapped_reads = 0
    unique_mapped_reads = 0
    ambiguous_reads = 0
    
    for read_id in all_read_ids:
        best_alignment = None
        best_mapq = -1
        
        # Check C->T alignments
        if read_id in c2t_alignments and 'c2t' in c2t_alignments[read_id]:
            c2t_mapq = c2t_alignments[read_id]['c2t']['mapq']
            if c2t_mapq > best_mapq:
                best_alignment = ('c2t', c2t_alignments[read_id]['c2t'])
                best_mapq = c2t_mapq
        
        # Check G->A alignments if non-directional or no C->T alignment
        if not args.directional or best_mapq == -1:
            if read_id in g2a_alignments and 'g2a' in g2a_alignments[read_id]:
                g2a_mapq = g2a_alignments[read_id]['g2a']['mapq']
                if g2a_mapq > best_mapq:
                    best_alignment = ('g2a', g2a_alignments[read_id]['g2a'])
                    best_mapq = g2a_mapq
        
        if best_alignment:
            mapped_reads += 1
            if best_mapq > 0:
                unique_mapped_reads += 1
                
                # Extract alignment details
                conv_type, aln = best_alignment
                flag = aln['flag']
                chrom = aln['chrom']
                pos = aln['pos']
                mapq = aln['mapq']
                cigar = aln['cigar']
                seq = aln['seq']
                
                # Extract genomic sequence for methylation calling
                genomic_seq = extract_genomic_sequence(chrom, pos, cigar, genome_dict)
                
                if genomic_seq:
                    # Call methylation
                    methyl_string = call_methylation(seq, genomic_seq, conv_type)
                    
                    # Create BAM record
                    a = pysam.AlignedSegment()
                    a.query_name = read_id
                    a.flag = flag
                    a.reference_id = list(genome_dict.keys()).index(chrom)
                    a.reference_start = pos - 1  # 0-based
                    a.mapping_quality = mapq
                    a.cigar = [(int(count), op) for count, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar)]
                    a.query_sequence = seq
                    a.tags = [("XM", methyl_string), ("XR", conv_type)]
                    
                    outfile.write(a)
            else:
                ambiguous_reads += 1
    
    outfile.close()
    
    # Write report
    with open(f"{args.output_prefix}_report.txt", "w") as report:
        report.write(f"Bismark Native Alignment Report for {args.output_prefix}\n")
        report.write(f"Total reads processed: {total_reads}\n")
        if total_reads > 0:
            report.write(f"Mapped reads: {mapped_reads} ({mapped_reads/total_reads*100:.1f}%)\n")
            report.write(f"Uniquely mapped reads: {unique_mapped_reads} ({unique_mapped_reads/total_reads*100:.1f}%)\n")
            report.write(f"Ambiguously mapped reads: {ambiguous_reads} ({ambiguous_reads/total_reads*100:.1f}%)\n")
        else:
            report.write("No reads processed\n")

if __name__ == "__main__":
    main() 