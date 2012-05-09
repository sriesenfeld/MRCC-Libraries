#!/usr/bin/perl 
#
# SJ Riesenfeld
# Created: 2009
# Last updated: Feb 2012
#
# Basic functions for working with sequences/strings. Used in
# particular by DNA_DeBruijnGraph_Utils.pm and Palindrome_Paths.pm.

package Seq_Utils;

use strict;
use warnings;

require Exporter;
our @ISA=qw(Exporter);
our @EXPORT_OK = qw( choose_random_elt get_matches get_head get_tail
                     in_list sublist_in_list get_index remove remove_sublist
                     is_prefix is_suffix suffix_prefix_match );

# Choose a random element from an array.
sub choose_random_elt($) {
    my $array_ref = shift();  # ref to array
    my $index=int(rand(scalar(@$array_ref)));
    return ($array_ref->[$index]);
}

# Get matches in the given array to a given sequence, where the
# matches are either exact or have the sequence as a prefix, or have
# the sequence as a suffix (depending on the second arg).
#
# Returns an array of matches.
sub get_matches ($$$) {
    my ($seq,  # string of at most $k characters
	$match_type, # integer as specified: 0 = exact match; 1 =
		     # prefix; 2 = suffix;
	$seq_list_ref # a reference to an array of k-mers.
	) = @_;
    my $pattern;
    if ($match_type == 1) {
	$pattern = '^'.$seq . '.+';
    } elsif ($match_type == 2 ) {
	$pattern = '.+'.$seq .'$'; 
    }  elsif ($match_type==0) { 
	$pattern = '^'.$seq.'$';
    } else {
	return undef;
    }
    grep($_ =~ $pattern, @$seq_list_ref);
}

# Given a string of length k, return the suffix of length k-1 (or the
# prefix of length k-1, depending on the optional $dir arg). The
# returned value corresponds to the "head" of the directed edge
# labeled by that string in the de Bruijn graph ($dir can flip the
# direction of the edge).
sub get_head ($;$) {
    my ($str, $dir) = @_;
    if ($dir) {
	return (substr($str, 0, length($str)-1));
    } else {
	return (substr($str, 1, length($str)));
    }
}

# Given a string of length k, return the prefix of length k-1 (or the
# suffix of length k-1, depending on the optional $dir arg). The
# returned value corresponds to the "tail" of the directed edge
# labeled by that string in the de Bruijn graph ($dir can flip the
# direction of the edge).
sub get_tail ($;$) {
    my ($str, $dir) = @_;
    if ($dir) {
	return (substr($str, 1, length($str)));
    } else {
	return (substr($str, 0, length($str)-1));
    }       	
}

# Is the element in the list (given by reference)?
sub in_list {
    my ($elt, $list_ref) =@_;
    if (grep($_ eq $elt, @$list_ref)) {
	return 1; }
    else { return 0;}
} 
	
# Are the elements in the sublist (given by reference) contained in
# the list (given by reference)?
sub sublist_in_list {
    my ($sublist_ref, $list_ref) = @_;
    foreach my $elt (@$sublist_ref) {
	if (! in_list($elt, $list_ref) ) {
	    return 0;
	}
    }
    return 1;
}

# What is the index of the given element in the given array (given by
# reference)?
sub get_index($$) {
    my ($elt, $array_ref) = @_;
    my @array = @$array_ref;
    my ($index) = grep $array[$_] eq $elt, 0 .. $#array;
    return $index;
}

# Remove the given element from the given array (given by reference).
sub remove {
    my ($elt, $list_ref) = @_;
    @$list_ref= grep($_ ne $elt, @$list_ref);
}

# Remove the first list (passed by reference) from the second list
# (passed by reference).
sub remove_sublist {
    my ($r_list_ref, $list_ref) = @_;
    @$list_ref = grep( (!in_list($_, $r_list_ref)), @$list_ref);
}

# Is the first sequence a prefix of the second sequence?
sub is_prefix {
    my ($pot_prefix, $seq)= @_;
    unless ($pot_prefix and $seq) { return 0; }
    my $pattern = '^'.$pot_prefix . '.+';
    if ($seq =~ $pattern) { return 1;}
    else { return 0;}
}

# Is the first sequence a suffix of the second sequence?
sub is_suffix {
    my ($pot_suffix, $seq)= @_;
    unless ($pot_suffix and $seq) { return 0;}
    my $pattern = '.+'.$pot_suffix .'$';
    if ($seq =~ $pattern) { return 1;}
    else { 
	return 0;}
}

# Given two strings, find the longest suffix of the first one that is
# a prefix of the second. Returns the length of the suffix.
sub suffix_prefix_match($$) {
    my ($seq1, $seq2) = @_;
    my $len1 = length($seq1);
    my $best = 0;
    for (my $i=1; $i <= $len1; $i++) {
	my $suffix = substr($seq1, $len1-$i, $i);
	if (($suffix eq $seq2) or (is_prefix ($suffix, $seq2))) {
	    $best = $i;
	} 
    }
    return $best;
}

1;
