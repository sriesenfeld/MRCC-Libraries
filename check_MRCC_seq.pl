#!/usr/bin/perl -w
#
# SJ Riesenfeld
# Created: 2009
# Last updated: Feb 2012
#
# Check if an input sequence is an MRCC sequence (MRCC library of size
# 1) for a given value of $k; that is, every $k-mer is represented in
# the sequence either by itself or by its reverse complement (not
# both), and that representative appears exactly one time in the
# sequence.

use strict;
use warnings;
use DNA_Manip qw ( generate_all_k_mers get_self_rev_comps rev_comp );
use Seq_Utils qw (is_prefix get_matches remove);

print join(", ", @ARGV);
my ($k, $input_seq) = @ARGV;
unless ($k and $input_seq) {die "Please supply the value of k and the input sequence.\n";}

my @kmers = generate_all_k_mers($k);
my @palindromes = get_self_rev_comps(\@kmers);
my $num = (scalar(@kmers) + scalar(@self_rcs))/2;
print "\nTotal number of $k-mers: ".scalar(@kmers)." and total number palindromes (self-reverse-complementary): "
    .scalar(@palindromes)."\n";
print "Total number without reverse-complements except for palindromes: $num.\n";

my $len = length($input_seq);
print "Sequence length: ".$len."\n";

my $num_kmers_given=$len-$k+1;
if (($num_kmers_given) != $num) {
    print "Total number of $k-mers is ". (($num_kmers_given>0) ? $num_kmers_given : 0).", but it should be $num!\n";
}
my $last_start = $len - $k;
my $start;
my $last_kmer;
for ($start=0; $start <= $last_start; $start++) {
    my $cur_kmer = substr($input_seq, $start, $k);
    my $cur_kmer_rc = rev_comp($cur_kmer);
    if ($start>0) {
	my $prefix = substr($last_kmer, 1, $k-1);
	if (! is_prefix($prefix, $cur_kmer) ) {
	    die "Not a de-Bruijn-like sequence!  Error at index $start at $k-mer: $cur_kmer.\n";
	}
    }
    my @matches = get_matches($cur_kmer, 0, \@kmers);
    my @matches_rc = get_matches($cur_kmer_rc, 0, \@kmers);
    if ( (scalar(@matches)>1) or (scalar(@matches_rc)>1) ) {
	die "Error at index $start: Too many matches for $cur_kmer or $cur_kmer_rc!\n";
    }
    if ((!@matches) or (!@matches_rc)) {
	if ((!@matches) and (!@matches_rc)) {
	    die "Error at index $start: $k-mers $cur_kmer and its rev comp $cur_kmer_rc have already been represented in the string!\n";
	} elsif (@matches_rc) {
	    die "Error at index $start: $k-mer $cur_kmer has already been represented in the string!\n".
		"For some reason, it appears that $k-mer $cur_kmer_rc has not been represented yet!\n";
	} elsif (@matches) {
	    die "Error at index $start: $k-mer $cur_kmer_rc has already been represented in the string!\n".
		"For some reason, it appears that $k-mer $cur_kmer has not been represented yet!\n";
	}
    }
    remove(shift(@matches), \@kmers);
    remove(shift(@matches_rc), \@kmers);
    $last_kmer = $cur_kmer;
}
if (@kmers) {
    die "There are unrepresented $k-mers: \n".
	join(' ', @kmers);
}
