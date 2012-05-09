##!/usr/bin/perl -w
# 
# SJ Riesenfeld
# Created: 2009
# Last updated: Feb 2012
#
# Given a file of oligomers (or constructs), find the set of
# high-frequency (HF) k-mers. Return data on their frequencies, as
# well as a set of oligomers that covers all of the HF k-mers.

use strict;
use warnings;

use File::Spec;
use Getopt::Std;

use DNA_Manip qw (get_self_rev_comps generate_all_k_mers );
use DNA_DeBruijnGraph_Utils qw( print_construct_paths print_constructs check_constructs
                                print_array remove_sublist build_extra_constructs 
                                add_flanking_seqs_to_constructs );

################### 
#### These parameters are set currently for an application-specific
#### use and can be reset by the user. They could easily be made into
#### input options. Note that they would typically be set the same for
#### all the related programs.
my $k=6; 
  # the left (right) flanking sequence is assumed not to end (begin) with 2 bases that are the same
my $flank_l = 'CTCGAG';
my $flank_r = 'AGATCT';
my $oligo_len = 15;   # oligomer length
my $num_trials_build_extras = 3; # local variable for how many times to try a heuristic, randomized process
my $freq_threshold=20; # how frequent is High-Frequency?
################### 

my $usage = qq{
Given a file of constructs,the $0 program designs constructs containing the most frequently represented kmers.\n
Usage:  $0 <options (see below)>

   -i <name of file of constructs, in format output by ./testmers.pl>
   -s <name of output solution file, that contains the new constructs>
   -r <name of file to contain the output histogram data for the
         repeated kmers in the construct design (analyzed with the
         flanking sequences on both sides); analyzed **without** the
         input constructs, only with the newly constructed oligomers>
   -a <name of file to contain the output histogram data for the
         repeated kmers in the construct design (analyzed with the
         flanking sequences on both sides); analyzed **with** the
         input constructs and the newly constructed oligomers>
   -o <name of file to contain other output of the program>
         if 0/undefined, STDOUT is used>
   -h [no value taken; this message is printed; see the program file 
         for details]
};

my %opts;
if ( !getopts('i:s:r:o:a:h', \%opts) or $opts{'h'} ) {
    die "$usage\n";
}
my @flank_seqs = ($flank_l, $flank_r);
my $path_len = $oligo_len - ($k-1);

my $insolnfile = $opts{'i'};
my $outsolnfile = $opts{'s'};
my $outfile = $opts{'o'};
my $hist_hf_file = $opts{'r'};
my $hist_all_file = $opts{'a'};
if (!$insolnfile or !$outsolnfile or !$hist_hf_file or !$hist_all_file) {
    die "Please provide arguments for options '-i', '-s', '-r', and '-a'.\n";
}

my $db_str_l = substr($flank_l, ($k-2), 2);
my $db_str_r = substr($flank_r, 0, 2);

open (INSOLN, "$insolnfile") or die "Cannot open $insolnfile: $!\n";
my $out_fh;
if ($outfile) {
    open (OUT, ">$outfile") or die "Cannot open $outfile: $!\n";
    $out_fh=\*OUT;
} else {
    $out_fh = \*STDOUT;
}

my (@constructs, @construct_paths);
my @kmers=generate_all_k_mers($k);   # all the edges in the original graph

my $line = <INSOLN>;  # skip first line
while ($line = <INSOLN>) {        
    my ($construct, $index);
    my $pattern = qr{([\d]+):[\s]+([AGCT]+)};
    if ($line =~ $pattern) {
	$index = $1;
	$construct = $2;
    }
    if ($construct) {
	if (length($construct) != $oligo_len) {
	    die "Error in construct length of construct $construct at index $index.\n";
	}
	push(@constructs, $construct);
    } 
}
close(INSOLN);
print $out_fh "Total of ".scalar(@constructs)." constructs read from $insolnfile.\n";

for my $i (0..$#constructs) {
    my $construct = $constructs[$i];
    print $out_fh "Construct $i: $construct\n";
    my @path;
    while (length($construct)>=$k) {	
	my $edge = substr($construct,0, $k);	
	$construct = substr($construct, 1);
	push(@path, $edge);
    }
    if (scalar(@path) != $path_len) {
	die "Error in path length of construct $i.\n";
    }
    $construct_paths[$i] = \@path;
}

print_construct_paths(\@construct_paths, $out_fh);

my $all_constructs_wf_ar =
    add_flanking_seqs_to_constructs($k, $oligo_len, \@construct_paths, $flank_l, $flank_r);

my ($bad_flag, $edges_by_construct_wf_hr, $repeat_edges_wf_hr) =
    check_constructs($k, $all_constructs_wf_ar, $oligo_len + 12, $out_fh, 0, 1);
if ($bad_flag) {    
    die "Input constructs with flanking sequences added do not constitute a valid solution!\n";
} else {
    print $out_fh "\nInput construct design with flanking sequences added has been validated.\n";
}

# build a few constructs from highest frequency kmers
print $out_fh "\nFinding the high-frequency kmers.\n";
my @repeated_edges = keys(%{$repeat_edges_wf_hr});
remove_sublist(\@flank_seqs, \@repeated_edges);
my @hfreq_kmers = grep {scalar(@{$repeat_edges_wf_hr->{$_}}) >= $freq_threshold} @repeated_edges;
@hfreq_kmers = sort { scalar(@{$repeat_edges_wf_hr->{$b}}) <=> scalar(@{$repeat_edges_wf_hr->{$a}}) } @hfreq_kmers;

print $out_fh "Building constructs from the following ".scalar(@hfreq_kmers)." high-frequency (i.e., frequency >= $freq_threshold) kmers:\n";
print $out_fh "kmer\tfrequency\n";
foreach  my $kmer (@hfreq_kmers) {
    print $out_fh "$kmer\t".scalar(@{$repeat_edges_wf_hr->{$kmer}})."\n"
}
    
my @hfreq_paths;
foreach my $kmer (@hfreq_kmers) {
    push(@hfreq_paths, [$kmer]);
}
my ($hf_constructs_ar, $hf_added_edges_ar) = 
    build_extra_constructs($k, \@hfreq_paths, $oligo_len, \@kmers, $db_str_l, $db_str_r, 
    $num_trials_build_extras);

print $out_fh "\nBuilt ".scalar(@{$hf_constructs_ar})." constructs, using ".scalar(@{$hf_added_edges_ar})." additional edges.\n";
print "Printing the construct paths:\n";
print_construct_paths($hf_constructs_ar, $out_fh);

print "\nChecking the design of just the new constructs.\n";
my ($hf_edges_by_construct_hr, $hf_repeat_edges_hr);
($bad_flag, $hf_edges_by_construct_hr, $hf_repeat_edges_hr)=
    check_constructs($k, $hf_constructs_ar, $oligo_len, $out_fh, 0,
		     1, $hf_added_edges_ar, 0, $db_str_l, $db_str_r, 1);
if ($bad_flag) {    
    die "New constructs are not valid!\n";
} else {
    print $out_fh "\nDesign of new constructs has been validated.\n";
}
print "Adding flanking sequence to the new constructs.\n";
my $hf_constructs_wf_ar =
    add_flanking_seqs_to_constructs($k, $oligo_len, $hf_constructs_ar, $flank_l, $flank_r); 

open(HIST_HF, ">$hist_hf_file") or die "Cannot open $hist_hf_file: $!\n";
print $out_fh "Checking new constructs with flanking sequences added.\n".
    "Printing histogram data for new contructs only in $hist_hf_file.\n";
my ($hf_edges_by_construct_wf_hr, $hf_repeat_edges_wf_hr);
($bad_flag, $hf_edges_by_construct_wf_hr, $hf_repeat_edges_wf_hr) =
    check_constructs($k, $hf_constructs_wf_ar, $oligo_len + 12, $out_fh, \*HIST_HF, 
		     1, undef, undef, undef, undef, 1);
if ($bad_flag) {    
    die "New constructs with flanking sequences added are not valid.\n";
} else {
    print $out_fh "\nInput construct design with flanking sequences added has been validated.\n";
}
close(HIST_HF);

open(HIST_ALL, ">$hist_all_file") or die "Cannot open $hist_all_file: $!\n";

print $out_fh "Checking all constructs together (with flanking sequences).\n".
    "Printing histogram data for all contructs (input and newly created) in $hist_all_file.\n";
push(@{$all_constructs_wf_ar}, @{$hf_constructs_wf_ar});
($bad_flag, $edges_by_construct_wf_hr, $repeat_edges_wf_hr) =
    check_constructs($k, $all_constructs_wf_ar, $oligo_len + 12, $out_fh, \*HIST_ALL, 1);

if ($bad_flag) {    
    die "Constructs with flanking sequences added do not constitute a valid solution!\n";
} else {
    print $out_fh "\nInput construct design with flanking sequences added has been validated.\n";
}
close(HIST_ALL);

open (OUTSOLN, ">$outsolnfile") or die "Cannot open $outsolnfile: $!\n";
print $out_fh "Printing high-frequency-kmer constructs to $outsolnfile.\n";
print_constructs($k, $hf_constructs_ar, \*OUTSOLN, $oligo_len, $db_str_l, $db_str_r);
