#!/usr/bin/perl 
#
# SJ Riesenfeld
# Created: 2009
# Last updated: Feb 2012
#
# Functions for some simple DNA-specific manipulations of DNA
# sequences (as strings).
package DNA_Manip;

use strict;
use warnings;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw( rev_comp is_self_rev_comp rev_comp_list
                     generate_all_k_mers get_self_rev_comps );

# Take reverse complement of DNA sequence
sub rev_comp($) {
    my $seq = shift(@_);
    $seq =~ tr/ACGTacgt/TGCAtgca/;
    $seq = reverse($seq);  # need to call reverse in separate
			   # statement from return; otherwise
			   # sometimes sequences does not get reversed
    return ($seq);
}

# Return 1 if sequence is the reverse complement of itself; 0 otherwise.
sub is_self_rev_comp($) {
    my $seq = shift(@_);
    my $revcomp_seq = &rev_comp($seq);
    if ($seq eq $revcomp_seq) {
	return 1;
    } else {
	return 0;
    }
}	

# Takes as an argument a reference to an array. Returns a reference to
# an array that has as its elements the reverse complement of each
# element in the given array. 
sub rev_comp_list($) {
    my $list_ref = shift();
    my @rc_list=();
    foreach my $seq (@$list_ref) {
	my $rc_seq = rev_comp($seq);
	# print "seq: $seq, rc_seq: $rc_seq\n";
	push(@rc_list, $rc_seq);
    }
    return \@rc_list;
}

# Return an array of the seqs that are reverse complements of
# themselves that appear in the given array (given by reference).
sub get_self_rev_comps($) {
    my $list_ref= shift(@_);
    my @rcs=();
    foreach my $seq (@$list_ref) {
	if (is_self_rev_comp($seq)) {
	    push(@rcs, $seq);
	}
    }
    return @rcs;
}	

# Generate all DNA sequences of a given length and return as an array.
sub generate_all_k_mers ($) {
    my $k = shift(@_);
    my @list_k_mers;
    my $total = 4**$k;
    my @chars = qw(A C G T); 
    my $curr_k_mer='';
    for (my $i=0; $i< $total; $i++) {
	for (my $pos = 0; $pos <$k; $pos++) { 
	    $curr_k_mer = ($chars[($i/(4**$pos))% 4]).$curr_k_mer;
	}
	$list_k_mers[$i]=$curr_k_mer;
	$curr_k_mer='';
    }
    return @list_k_mers;
}

1;
