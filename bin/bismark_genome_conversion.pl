#!/usr/bin/env perl
use strict;
use warnings;

# Parse command line arguments
my $genome_folder = "";
my $output_dir = "";
my $multi_fasta = "true";
my $slam = "false";

foreach my $arg (@ARGV) {
    if ($arg =~ /--genome_folder=(.+)/) {
        $genome_folder = $1;
    }
    elsif ($arg =~ /--output_dir=(.+)/) {
        $output_dir = $1;
    }
    elsif ($arg =~ /--multi_fasta=(.+)/) {
        $multi_fasta = $1;
    }
    elsif ($arg =~ /--slam=(.+)/) {
        $slam = $1;
    }
}

# Ensure genome folder path ends with slash
unless ($genome_folder =~ /\\/$/){
    $genome_folder .= "/";
}

# CT and GA conversion directories
my $CT_dir = "${output_dir}/CT_conversion/";
my $GA_dir = "${output_dir}/GA_conversion/";

# Check for FASTA files
chdir $genome_folder or die "Could not move to directory $genome_folder: $!";
my @filenames = <*.fa *.fasta *.fa.gz *.fasta.gz>;

# Setup output files based on multi_fasta setting
my ($CT_fh, $GA_fh);
if ($multi_fasta eq "true") {
    open($CT_fh, ">", "${CT_dir}/genome_mfa.CT_conversion.fa") 
        or die "Cannot write to ${CT_dir}/genome_mfa.CT_conversion.fa: $!";
    open($GA_fh, ">", "${GA_dir}/genome_mfa.GA_conversion.fa") 
        or die "Cannot write to ${GA_dir}/genome_mfa.GA_conversion.fa: $!";
}

# Process each input FASTA file
my %chr_names;
foreach my $filename (@filenames) {
    my $in_fh;
    if ($filename =~ /\\.gz$/) {
        open($in_fh, "gunzip -c $filename |") or die "Cannot read $filename: $!";
    } else {
        open($in_fh, "<", $filename) or die "Cannot read $filename: $!";
    }
    
    # Process each sequence
    my $line = <$in_fh>;
    chomp $line;
    if ($line !~ /^>/) {
        die "File $filename does not appear to be in FASTA format";
    }
    
    # Extract chromosome name
    my $chr = $line;
    $chr =~ s/^>//;
    $chr = (split(/\\s+/, $chr))[0];
    
    # Check for duplicates
    if (exists $chr_names{$chr}) {
        die "Duplicate chromosome name: $chr";
    }
    $chr_names{$chr} = 1;
    
    # Open new files if needed
    if ($multi_fasta eq "false") {
        close($CT_fh) if defined $CT_fh;
        close($GA_fh) if defined $GA_fh;
        open($CT_fh, ">", "${CT_dir}/${chr}.CT_conversion.fa") 
            or die "Cannot write to ${CT_dir}/${chr}.CT_conversion.fa: $!";
        open($GA_fh, ">", "${GA_dir}/${chr}.GA_conversion.fa") 
            or die "Cannot write to ${GA_dir}/${chr}.GA_conversion.fa: $!";
    }
    
    print $CT_fh ">${chr}_CT_converted\\n";
    print $GA_fh ">${chr}_GA_converted\\n";
    
    # Process the sequence
    while (my $seq = <$in_fh>) {
        chomp $seq;
        if ($seq =~ /^>/) {
            # New sequence header
            $chr = $seq;
            $chr =~ s/^>//;
            $chr = (split(/\\s+/, $chr))[0];
            
            if (exists $chr_names{$chr}) {
                die "Duplicate chromosome name: $chr";
            }
            $chr_names{$chr} = 1;
            
            if ($multi_fasta eq "false") {
                close($CT_fh);
                close($GA_fh);
                open($CT_fh, ">", "${CT_dir}/${chr}.CT_conversion.fa") 
                    or die "Cannot write to ${CT_dir}/${chr}.CT_conversion.fa: $!";
                open($GA_fh, ">", "${GA_dir}/${chr}.GA_conversion.fa") 
                    or die "Cannot write to ${GA_dir}/${chr}.GA_conversion.fa: $!";
            }
            
            print $CT_fh ">${chr}_CT_converted\\n";
            print $GA_fh ">${chr}_GA_converted\\n";
        } else {
            # Sequence data
            $seq = uc($seq);
            $seq =~ s/[^ATCGN\\n\\r]/N/g;
            
            my $CT_seq = $seq;
            my $GA_seq = $seq;
            
            if ($slam eq "true") {
                $CT_seq =~ tr/T/C/; # SLAM-seq: T->C
                $GA_seq =~ tr/A/G/; # SLAM-seq: A->G
            } else {
                $CT_seq =~ tr/C/T/; # Bisulfite: C->T
                $GA_seq =~ tr/G/A/; # Bisulfite: G->A
            }
            
            print $CT_fh $CT_seq;
            print $GA_fh $GA_seq;
        }
    }
    
    close($in_fh);
}

close($CT_fh) if defined $CT_fh;
close($GA_fh) if defined $GA_fh;

print "Genome conversion complete\\n";