#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;
use File::Spec;
use File::Path qw(make_path);

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

# Debug output
print "DEBUG: genome_folder = '$genome_folder'\n";
print "DEBUG: output_dir = '$output_dir'\n";
print "DEBUG: multi_fasta = '$multi_fasta'\n";
print "DEBUG: slam = '$slam'\n";

# Validate required parameters
die "Error: --genome_folder parameter is required\n" unless $genome_folder;
die "Error: --output_dir parameter is required\n" unless $output_dir;

# Create output directories 
my $CT_dir = File::Spec->catdir($output_dir, "CT_conversion");
my $GA_dir = File::Spec->catdir($output_dir, "GA_conversion");

# Use make_path which doesn't fail if directories already exist
make_path($output_dir, $CT_dir, $GA_dir);

# Get list of input files
my @filenames;
if (-f $genome_folder) {
    # Single file input
    @filenames = ($genome_folder);
    print "DEBUG: Processing single file: $genome_folder\n";
} elsif (-d $genome_folder) {
    # Directory input
    opendir(my $dh, $genome_folder) or die "Cannot open directory $genome_folder: $!";
    @filenames = map { File::Spec->catfile($genome_folder, $_) } 
                 grep { /\.(fa|fasta)(\.gz)?$/i } readdir($dh);
    closedir($dh);
    print "DEBUG: Found " . scalar(@filenames) . " FASTA files in directory\n";
} else {
    die "Error: $genome_folder is neither a file nor a directory\n";
}

die "Error: No input FASTA files found in $genome_folder\n" unless @filenames;

# Setup output files based on multi_fasta setting
my ($CT_fh, $GA_fh);
if ($multi_fasta eq "true") {
    open($CT_fh, ">", File::Spec->catfile($CT_dir, "genome_mfa.CT_conversion.fa")) 
        or die "Cannot write to CT conversion file: $!";
    open($GA_fh, ">", File::Spec->catfile($GA_dir, "genome_mfa.GA_conversion.fa")) 
        or die "Cannot write to GA conversion file: $!";
}

# Process each input FASTA file
my %chr_names;
foreach my $filename (@filenames) {
    print "DEBUG: Processing file: $filename\n";
    my $in_fh;
    if ($filename =~ /\.gz$/i) {
        open($in_fh, "gunzip -c $filename |") or die "Cannot read $filename: $!";
    } else {
        open($in_fh, "<", $filename) or die "Cannot read $filename: $!";
    }
    
    # Process each sequence
    my $line = <$in_fh>;
    chomp $line if defined $line;
    if (!defined $line || $line !~ /^>/) {
        die "File $filename does not appear to be in FASTA format";
    }
    
    # Extract chromosome name
    my $chr = $line;
    $chr =~ s/^>//;
    $chr = (split(/\s+/, $chr))[0];
    
    # Check for duplicates
    if (exists $chr_names{$chr}) {
        die "Duplicate chromosome name: $chr";
    }
    $chr_names{$chr} = 1;
    
    # Open new files if needed
    if ($multi_fasta eq "false") {
        close($CT_fh) if defined $CT_fh;
        close($GA_fh) if defined $GA_fh;
        open($CT_fh, ">", File::Spec->catfile($CT_dir, "${chr}.CT_conversion.fa")) 
            or die "Cannot write to CT conversion file for $chr: $!";
        open($GA_fh, ">", File::Spec->catfile($GA_dir, "${chr}.GA_conversion.fa")) 
            or die "Cannot write to GA conversion file for $chr: $!";
    }
    
    print $CT_fh ">${chr}_CT_converted\n";
    print $GA_fh ">${chr}_GA_converted\n";
    
    # Process the sequence
    while (my $seq = <$in_fh>) {
        chomp $seq;
        if ($seq =~ /^>/) {
            # New sequence header
            $chr = $seq;
            $chr =~ s/^>//;
            $chr = (split(/\s+/, $chr))[0];
            
            if (exists $chr_names{$chr}) {
                die "Duplicate chromosome name: $chr";
            }
            $chr_names{$chr} = 1;
            
            if ($multi_fasta eq "false") {
                close($CT_fh);
                close($GA_fh);
                open($CT_fh, ">", File::Spec->catfile($CT_dir, "${chr}.CT_conversion.fa")) 
                    or die "Cannot write to CT conversion file for $chr: $!";
                open($GA_fh, ">", File::Spec->catfile($GA_dir, "${chr}.GA_conversion.fa")) 
                    or die "Cannot write to GA conversion file for $chr: $!";
            }
            
            print $CT_fh ">${chr}_CT_converted\n";
            print $GA_fh ">${chr}_GA_converted\n";
        } else {
            # Sequence data
            $seq = uc($seq);
            $seq =~ s/[^ATCGN\n\r]/N/g;
            
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

print "Genome conversion complete\n";