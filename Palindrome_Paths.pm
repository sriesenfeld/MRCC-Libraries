#!/usr/bin/perl -w
# 
# SJ Riesenfeld
# Created: 2009
# Last updated: Feb 2012
#
# Contains functions for computing matchings between pairs of DNA
# palindromes (sequences that are their own reverse complements), and
# paths in the de Bruijn graph between pairs of palindromes in the
# matching. Also contains functions for processing the paths and
# attempting to insert them into a given cycle containing other
# k-mers.
#
# Designed to work with the DNA_DeBruijnGraph_Utils package.

package Palindrome_Paths;

use strict;
use warnings;
# use List::Permutor;
use DNA_Manip qw( rev_comp is_self_rev_comp rev_comp_list );
use Seq_Utils qw (choose_random_elt get_matches in_list sublist_in_list
                  remove remove_sublist is_prefix is_suffix 
                  suffix_prefix_match get_head get_tail );                   
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw ( get_ppaths get_pdists print_insertion_points 
                      find_insertion_points 
                      get_short_path_palindromes print_ppaths
                      insert_short_ppaths compute_pconstructs );

my $verbose = 0;
my $in=1;  # direction of an edge with respect to a vertex (either $in or $out)
my $out=2; 

# The function looks for insertion points in a given cycle where 
# cyclical path between two palindromes can be spliced in.
#
# A path is 'forward' is it goes from palindrome to
# reverse(palindrome); an RC (reverse complement) path is from
# reverse(palindrome) to palindrome.  It is possible to insert the
# 2-edge cyclical paths between some pairs of palindromes (in the $in
# direction) if the last (k-1)-mer in it (the head of the path) is a
# prefix of a k-mer in the big cycle or (in the $out direction) if the
# first (k-1)-mer in the path (the tail of the path) is a suffix of a
# k-mer in the cycle.
#
# Returns an array of two references to hashes that contain
# information about possible points in the cycle where the cycle can
# be broken and the palindromic paths inserted. The keys for the first
# hash are the indices in the given cycle, and the values are the
# actual palindrome that can be inserted at that index, and the
# direction, which is $in if the path feeds into the cycle and $out if
# the cycle feeds into the (head of the) path. The second returned
# hash has as the palindromes as keys, and the value for a specific
# palindrome is a reference to a hash containing the index of a k-mer
# in the input cycle before or after which that palindrome can be
# inserted; the second value gives the direction, which is $in if the
# path feeds into the k-mer and $out if the k-mer feeds into the path.

sub find_insertion_points($$$;$) {
    my ($palindromes_ar, # ref to a list of palindromes which should
			 # be a subset of the keys of the hash
			 # referenced by the third argument;
	$big_cycle_ar, # ref to a cyclic de-Bruin-type list of k-mers
		       # (excluding palindromes or k-mers in
		       # palindrome paths);
	$matching_hr, # ref to a hash that gives a matching of
		      # palindromes by key, value;
	$out_fh) = @_; # (optional) file handle for printed output.
    my %cycle_insertion_pts;
    my %p_insertion_pts;
    my $len = scalar(@$big_cycle_ar);
    if (!$out_fh) {
	$out_fh = \*STDOUT;
    }
    my $count=0;
    foreach my $palindrome (@$palindromes_ar) {
	my $target_palindrome = $matching_hr->{$palindrome};
	if (in_list($target_palindrome, $palindromes_ar)) {
	    warn "\nNot expecting match of palindrome $palindrome to be in the list!\n";
	}
	my $path_tail = get_tail($palindrome);
	my $path_head = get_head($target_palindrome);
	my $bw_path_tail = get_tail($target_palindrome);
	my $bw_path_head = get_head($palindrome);
	for (my $i=0; $i<$len; $i++) {	    
	    my $cur_kmer = $big_cycle_ar->[$i];
	    if ( is_prefix($path_head, $cur_kmer) or is_prefix($bw_path_head, $cur_kmer) ) { 
		push ( @{ $cycle_insertion_pts{$i} }, [$palindrome, $in] );
		push( @{ $p_insertion_pts{$palindrome} }, [$i, $in]);
	    }
	    if ( is_suffix($path_tail, $cur_kmer) or is_suffix($bw_path_tail, $cur_kmer) ) {
		push (@{ $cycle_insertion_pts{$i} }, [$palindrome, $out] );
		push( @{ $p_insertion_pts{$palindrome} }, [$i, $out]);
	    }
	}
	$count++;
    }
    print $out_fh "\nInsertion points found for $count palindrome paths.\n";
    return (\%cycle_insertion_pts, \%p_insertion_pts);
}

# Finds the pairs of palindromes that are in 2-edge cycles with each
# other.  The arg is a ref to a hash of shortest distances in the de
# Bruijn graph between all palindromes.  
#
# Returns a list of refs to hashes; the first hash is a matching of
# palindromes such that the key and value have distance 1 (i.e., path
# length 2); the second hash has the same keys but has as values
# references to arrays containing the 2-kmer paths.
sub get_short_path_palindromes($) {
    my ($ppdists_hr) = @_;
    my %short_ppaths;
    my %short_pmatching;
    my @rem_palindromes = keys(%{$ppdists_hr});
    my @used_palindromes;
    for my $palindrome (@rem_palindromes) {
	if (in_list($palindrome, \@used_palindromes)) {
	    next;
	}
	for my $target (@rem_palindromes) {
	    if (in_list($target, \@used_palindromes)) {
		next;
	    }
	    if ($ppdists_hr->{$palindrome}{$target} <= 1 ) {
		$short_ppaths{$palindrome}=[$palindrome, $target];
		$short_pmatching{$palindrome}=$target;
		push(@used_palindromes, $palindrome, $target);
	    }	    
	}
    }
    return (\%short_pmatching, \%short_ppaths);
}

# Inserts 2-edge cycles containing palindromes into a larger given
# cycle.
#
# Returns a reference to an array that contains the original cycle
# plus the inserted paths.
sub insert_short_ppaths($$$$;$) {
    my ($big_cycle_ar, # ref to cycle of k-mers (stored as an array);
	$matching_hr, #  ref to a hash that gives a matching of
		      #  palindromes at distance 1;
	$p_insert_pts_sp_hr, # ref to a hash that contains the
			     # possible insertion points for the
			     # palindromes in short cycles;
	$ppaths_hr, # ref to a hash that gives the paths between those
		    # palindromes;
	$out_fh # (optional) output file handle
	) =@_;
    if (!$out_fh) {
	$out_fh = \*STDOUT;
    }
    # paths of length two can simply be inserted into big cycle without any problem    
    my %assignments;
    my @assigned_indices;

    for my $palindrome (keys (%{$p_insert_pts_sp_hr})) {
	print $out_fh "For path beginning with palindrome $palindrome:\n";
	my $insert_pts_ar = $p_insert_pts_sp_hr->{$palindrome};
	my $num_insert_pts = scalar(@{$insert_pts_ar});
	my $random_num = int(rand($num_insert_pts));
	my $insert_ar = $insert_pts_ar->[$random_num];
	my $index = $insert_ar->[0];
	my $dir = $insert_ar->[1];
	print  $out_fh "   Out of $num_insert_pts possible sites, choosing insertion site just ".
	    ($dir == $in ? "before" : "after"). " edge ".$big_cycle_ar->[$index].".\n";
	if ( in_list($index, \@assigned_indices) ) {
	    die "Should not be assigning same index twice.\n"
	}
	$assignments{$index} = [$palindrome, $dir];
	push(@assigned_indices, $index);
    }
    
    my @bigger_cycle;
    @assigned_indices = sort {$a <=> $b} (keys(%assignments));
    if (! scalar(@assigned_indices) ) {
	warn "No short palindromic path inserted!";
	print $out_fh "No short palindromic path inserted!\n";
	return $big_cycle_ar;	    
    } else {
	print $out_fh "\nAssigned indices: ". join(' ', @assigned_indices). "\n";
    }

    my $cur_assigned_index = shift(@assigned_indices);
    for (my $i=0; $i < scalar(@{$big_cycle_ar}); $i++) {
	my $cur_kmer = $big_cycle_ar->[$i];
	if (defined($cur_assigned_index) and ($i == $cur_assigned_index)) {
	    my $palindrome = $assignments{$i}->[0];
	    my $dir = $assignments{$i}->[1];
	    my @ppath = get_matching_ppath( $cur_kmer, $palindrome, $matching_hr->{$palindrome}, $dir, $ppaths_hr );
	    if ($dir == $in) {
		push(@bigger_cycle, @ppath);
		push(@bigger_cycle, $cur_kmer);
	    } else {
		push(@bigger_cycle, $cur_kmer);
		push(@bigger_cycle, @ppath);
	    }
	    $cur_assigned_index = shift(@assigned_indices);
	} else {
	    push (@bigger_cycle, $cur_kmer);
	}
    }    
    return \@bigger_cycle;
}

# Given a k-mer and a pair of palindromes, return the forward or the
# reverse path between them, depending on an input argument.
#
# Returns either the forward or the reverse palindrome path, depending
# on the fourth argument $dir, i.e., whether the palindrome path should
# come out of or go into the given $kmer.
sub get_matching_ppath($$$$$) {
    my ($kmer, # string of length $k
	$palindrome, # self-rev-complementary $k-mer
	$target_palindrome, # match for that palindrome, according to
			    # $ppaths_hr
	$dir, # integer: $in or $out
	$ppaths_hr # ref to a hash that has palindromes as keys and
		   # the value for a key p1 is the path from p1 to
		   # some other palindrome p2 to which p1 is matched.
	)=@_;
    my $path_tail = get_tail($palindrome);
    my $path_head = get_head($target_palindrome);
    my $bw_path_tail = get_tail($target_palindrome);
    my $bw_path_head = get_head($palindrome);
    my @ppath = @{$ppaths_hr->{$palindrome}};
    my @bw_path = get_rev_ppath(@ppath);
    
    if ( ( ($dir == $in) and is_prefix($path_head, $kmer) )
	 or ( ($dir == $out) and is_suffix($path_tail, $kmer) ) ) {
	return @ppath;
    } elsif ( ( ($dir == $in) and is_prefix($bw_path_head, $kmer) )
	      or ( ($dir == $out) and is_suffix($bw_path_tail, $kmer) ) ) {
	return  @bw_path;
    } else {
	die "Error: Cannot find matching ppath!\n";
    }
}

# First arg is a ref to an array containing a path from one palindrome
# to another.
# Returns an array containing the reverse-complementary path from the
# endpoint to the start point.
sub get_rev_ppath(@) {
    my @ppath=@_;
    my $palindrome = shift(@ppath);
    my $target_palindrome = pop(@ppath);
    my @bw_ppath = ( $target_palindrome );
    foreach my $kmer ( reverse(@ppath) ) {
	my $rc_kmer = rev_comp($kmer);
	push(@bw_ppath, $rc_kmer);	
    }
    push(@bw_ppath, $palindrome);
    return @bw_ppath;
}


sub print_insertion_points($$$$;$) {
    my ($cycle_ar, $cycle_insert_pts_hr, $p_insert_pts_hr, $ppaths_hr, $out_fh) = @_;
    
    if (!$out_fh) {
	$out_fh = \*STDOUT;
    }
    my $num = scalar(@$cycle_ar);
    print $out_fh "Printing cycle insertion points...\n";
    foreach (my $i=0; $i< $num; $i++) {	
	if ($cycle_insert_pts_hr->{$i}) {
	    print $out_fh "At kmer ". $cycle_ar->[$i].", index $i: ";
	    foreach my $insert_ar (@{$cycle_insert_pts_hr->{$i}}) {
		print $out_fh "  Dir " . (($insert_ar->[1] == $in) ? "\"in\"" : "\"out\"").", for path to/from palindrome ".
		    $insert_ar->[0]." of length ".scalar(@{$ppaths_hr->{$insert_ar->[0]}})."\n";
	    }	
	}	
    }
    print $out_fh "Printing palindrome insertion points...\n";
    foreach my $palindrome (keys (%$p_insert_pts_hr)) { 
	print $out_fh "Insertion points for palindrome paths from/to $palindrome: \n";
	my $i=0;
	foreach my $ar (@{ $p_insert_pts_hr->{$palindrome}}) {	    
	    my $index = $ar->[0];
	    print $out_fh "  Index $index, kmer ". $cycle_ar->[$index].', dir '. ( ($ar->[1] == $in) ? "\"in\"" : "\"out\"").";";
	    if ($i==2) {
		print $out_fh "\n";
	    }
	    $i++;
	}
	print $out_fh "\n";
    }
}			 

# Construct a set of oligomers (constructs) such that each oligomer
# corresponds to a path from one palindrome to another so that none of
# the forward paths intersect with any of the reverse-complement
# paths.
#
# Returns, if possible, an array of references to 4 arrays.  The first
# array is another array of references to arrays, such that each array
# is a path of kmers connecting 2 palindromes (one each at the start
# and end), and the total length of the path (including start and end)
# is the length given for the construct minus ($k-1); the second is an
# array containing all non-palindromic edges used in the paths; the
# third is an array containing the reverse-complements of all the
# non-palindromic edges used in the paths; the fourth is the array of
# palindromes which appear on the constructs.  If the eighth arg is
# given, then the path length is extended by two if the first
# palindrome starts with $db_str_l and by one if it starts with the
# second base of the string. If the ninth arg is given, the analogous
# (symmetric) thing is done with for $db_str_r. 
sub compute_pconstructs($$$$ $$;$$$$$) {
    my ($k, # integer
	$construct_len, # length of (number of bases in) construct;
	$num_trials, # number of times the whole randomized process of
		     # pairing and computing paths should be repeated
		     # in case of failure;
	$num_trials_longppaths, #  max number of times that finding a
				#  path of the desired length between
				#  a particular pair should attempted;
	$palindromes_ar, # ref to an array of all palindromes that
			 # should be included in the constructs
			 # created;
	$kmers_ar, # ref to an array of possible edges (kmers) that
		   # can be used in the paths;
	$used_edges_ar, # (optional) a ref to an array of edges that
			# should not be used in the paths computed;
	$db_str_l, # (optional) 2-base nucleotide string (left
		   # flanking sequence);
	$db_str_r, # (optional) 2-base nucleotide string (right
		   # flanking sequence).
	$out_fh, # (optional) file handle to pring output to
	$num_oligos_short # (optional) maximum number of oligomers
			  # that can be of the given construct length
			  # (the rest must be 1bp longer)
	) = @_;    
    if (!$out_fh) {
	$out_fh = \*STDOUT;
    }
    my ($sb_l, $fb_r);
    if ($db_str_l) {
	$sb_l = substr($db_str_l, 1, 1);
    }
    if ($db_str_r) {
	$fb_r = substr($db_str_r, 0, 1);
    }
    my $pathlen = $construct_len - ($k-1);
    my @pconstructs;
    my @used_palindromes;
    my @used_fw_edges;
    my @used_rc_edges;    
    my $countc;
    print $out_fh "\nIn compute_pconstructs, path len is $pathlen, and "
	. scalar(@$palindromes_ar)." palindromes to be used.\n";
    for (my $trial=0; $trial< $num_trials; $trial++ ) {	
	my @rem_palindromes = @{$palindromes_ar};
	if ($num_oligos_short == scalar(@pconstructs)) {
	    print $out_fh "Constructed all oligomers of length $construct_len.\n";
	    $construct_len++;
	    $pathlen++;
	    print $out_fh "Target construct length now increased to $construct_len.\n";
	}
	$countc=0;
	@pconstructs=();
	@used_palindromes=();
	@used_fw_edges=();
	@used_rc_edges=();
	print $out_fh "\nTrial number $trial.\n";
	while (scalar(@rem_palindromes)) {
	    my $palindrome=choose_random_elt(\@rem_palindromes);
	    remove($palindrome, \@rem_palindromes);
	    print $out_fh "\nCurrent palindrome: $palindrome\n";
	    push(@used_palindromes, $palindrome);
	    while (scalar(@rem_palindromes)) {
		my $target = choose_random_elt(\@rem_palindromes);
		print $out_fh "Current target palindrome: $target\n";
		my $swap_order=0;
		if (is_prefix($db_str_l, $target) or is_suffix($db_str_r, $palindrome)) {
		    print $out_fh "Switching palindrome pair ordering.\n";
		    $swap_order=1;
		}
		my @temp_used_edges = (@used_fw_edges, @used_rc_edges);
		if ($used_edges_ar) {
		    push(@temp_used_edges, @{$used_edges_ar});
		}
		my ($ppath_ar, $used_fw_edges_ar, $used_rc_edges_ar);
		for (my $path_trial =0; $path_trial < $num_trials_longppaths; $path_trial++) { 
		    print $out_fh "\nTrying to get a long palindrome path: trial $path_trial.\n";
		    if (!$swap_order) {
			($ppath_ar, $used_fw_edges_ar, $used_rc_edges_ar) =
			    get_ppath($k, $palindrome, $target, \@temp_used_edges, $pathlen, $kmers_ar, $db_str_l, $db_str_r, $out_fh);
		    } else {
			($ppath_ar, $used_fw_edges_ar, $used_rc_edges_ar) =
			    get_ppath($k, $target, $palindrome, \@temp_used_edges, $pathlen, $kmers_ar, $db_str_l, $db_str_r, $out_fh);
		    }
		    if (scalar(@{$ppath_ar})) {
			last;
		    }
		}
		if (!(scalar(@{$ppath_ar}))) {
		    print $out_fh "In compute_pconstructs:  Could not find a long path from ".
			(!$swap_order ? "$palindrome to target $target" : "$target to target $palindrome")."!\n";
		    next;
		} else {
		    print $out_fh "Path of length ". scalar(@$ppath_ar)." found: ".join(' ', @$ppath_ar)."\n";
		    $pconstructs[$countc]= $ppath_ar;
		    push(@used_palindromes, $target);
		    remove($target, \@rem_palindromes);
		    push(@used_fw_edges, @{$used_fw_edges_ar});
		    push(@used_rc_edges, @{$used_rc_edges_ar});
		    last;
		}
	    }  # end loop over ending palindromes
	    if ($pconstructs[$countc]) {
		$countc++;
	    } elsif ($trial == ($num_trials-1)) {
		die "Could not build a construct starting with palindrome $palindrome.\n";
	    } else {
		last;
	    }
	} # end loop over starting palindromes
    } # end loop over trials
    print $out_fh "\n$countc constructs completed, using ".scalar(@used_palindromes). " palindromes.\n";
    return (\@pconstructs, \@used_fw_edges, \@used_rc_edges, \@used_palindromes);
}

# Computes a path of a given length (if one exists) from one given
# palindrome to another, such that the path includes no reverse
# complements of edges in the path and no other palindromes.
#
# Returns a list of two references, the first to an array that
# contains the path from the second arg to the third, the second to an
# array of the non-palindromic edges used in the path.  
# If the seventh arg is given, then the path length is extended by two
# if the first palindrome starts with $db_str_l and by one if it
# starts with the second base of the string.  The eigth argument, if
# given, is treated symmetrically.  
#
# If a desired length is not given, a short, direct path via shifting
# is computed; otherwise the desired length should be >= ($k+1) to
# make it easier to try to compute such a path.
sub get_ppath($$$;$$$ $$$) {
    my ($k, # integer
	$p1, # palindrome of length $k;
	$p2, # palindrome of length $k;
	$used_edges_ar, # (optional) ref to an array of kmers (edges)
			# that are not to be used in computing the
			# path;
	$total_path_len, # (optional) integer: desired total length of
			 # path, including palindromes at start and
			 # end.

	$edges_ar, # (optional) ref to an array of possible edges --
		   # this arg must be specified if the fifth arg is
		   # specified;
	$db_str_l, # (optional) 2-base nucleotide string (left
		   # flanking sequence); assumes string is not
		   # repeated base.
	$db_str_r, # (optional) 2-base nucleotide string (right
		   # flanking sequence); assumes string is not
		   # repeated base.
	$out_fh  # file handle where to write output.
	) =@_;

    if (!$out_fh) {
	$out_fh = \*STDOUT;
    }
    # print "In get_ppath: looking for path from $p1 to $p2.\n";
    if ($used_edges_ar) {
    } 
    my ($sb_l, $fb_r);
    if ($db_str_l) {
	$sb_l = substr($db_str_l, 1, 1);
    }
    if ($db_str_r) {
	$fb_r = substr($db_str_r, 0, 1);
    }
    
    my @ppath;    
    my @used_fw_edges;
    my @used_rc_edges;
    my $extra_len;
    if (!$total_path_len) {
	$extra_len=0;
    } else {
	if (!$edges_ar) {
	    die "Specify available edges for computing longer ppaths!\n";
	}
	if ($total_path_len < ($k+1)) {
	    die "Desired length for ppath should be at least ". ($k+1). "!\n";
	}
	$extra_len = $total_path_len-($k+1);
    }
    if ($db_str_l) {
	if (is_prefix($db_str_l, $p1)) {
	    $extra_len += 2;
	} elsif (is_prefix($sb_l, $p1)) {
	    $extra_len += 1;
	}
    }
    if ($db_str_r) {
	if (is_suffix($db_str_r, $p2)) {
	    $extra_len += 2;
	} elsif (is_suffix($fb_r, $p2)) {
	    $extra_len += 1;
	}
    }
    push(@ppath, $p1);
    my $head = get_head($p1);
    my $cur_edge;
    for (my $i=0; $i < $extra_len; $i++) {
	my @matches=get_matches($head, 1, $edges_ar);
	if (!(scalar(@matches))) {
	    warn "In get_ppath:  no edges out of $head!\n";
	    @ppath=();
	    return \@ppath;
	}
	$cur_edge=choose_random_elt(\@matches);
	my @temp_used_edges = (@used_fw_edges, @used_rc_edges);
	if ($used_edges_ar) {
	    push(@temp_used_edges, @{$used_edges_ar});
	}
	while (is_self_rev_comp($cur_edge) or in_list($cur_edge, \@temp_used_edges)) {
	    remove($cur_edge, \@matches);
	    if (scalar(@matches)) {
		$cur_edge=choose_random_elt(\@matches);
	    } else {
		warn "In get_ppath:  Cannot find unused, non-palindromic edge out of $head.\n";
		@ppath=();
		return \@ppath;
	    }
	}
	push(@ppath, $cur_edge);
	push(@used_fw_edges, $cur_edge);
	my $rc_cur_edge = rev_comp($cur_edge);
	push(@used_rc_edges, $rc_cur_edge);
	$head = get_head($cur_edge);
    }    
    if ($cur_edge) {
	for (my $i=1; $i <$k; $i++) {
	    my $next_str = substr($cur_edge, $i, $k-$i).substr($p2, 0, $i);
	    my @temp_used_edges = (@used_fw_edges, @used_rc_edges);
	    if ($used_edges_ar) {
		push(@temp_used_edges, @{$used_edges_ar});
	    }
	    if (is_self_rev_comp($next_str)) {
		print $out_fh "Trying to use palindrome $next_str as an edge in long palindromic path from $p1 to $p2, via $cur_edge!\n";
		@ppath=();
		return \@ppath;
	    } elsif (in_list($next_str, \@temp_used_edges)) {
		print $out_fh "Tried to use used edge $next_str in the long palindromic path from $p1 to $p2, via $cur_edge!\n"; 
		@ppath=();
		return \@ppath;
	    }	    
	    push(@ppath, $next_str);
	    push(@used_fw_edges, $next_str);
	    my $rc_next_str = rev_comp($next_str);
	    push(@used_rc_edges, $rc_next_str);
	}
    } else {
	my $len_suffix = suffix_prefix_match($p1, $p2);
	my $rem_len = $k-$len_suffix;
	for (my $i=1; $i<$rem_len; $i++) {
	    my $next_str = substr($p1, $i, $k-$i).substr($p2, $len_suffix, $i);
	    my @temp_used_edges = (@used_fw_edges, @used_rc_edges);
	    if ($used_edges_ar) {
		push(@temp_used_edges, @{$used_edges_ar});
	    }
	    if (is_self_rev_comp($next_str)) {
		warn "Trying to use palindrome $next_str as an edge in the path from $p1 to $p2!\n";
		@ppath=();
		return \@ppath;
	    } elsif (in_list($next_str, \@temp_used_edges)) {
		warn "Trying to use used edge $next_str in palindromic path from $p1 to $p2.\n";
		@ppath=();
		return \@ppath;
	    }
	    push(@ppath, $next_str);
	    push(@used_fw_edges, $next_str);
	    my $rc_next_str = rev_comp($next_str);
	    push(@used_rc_edges, $rc_next_str);
	}
    }
    push(@ppath, $p2);
    print $out_fh "Found ppath of length ".scalar(@ppath)." from $p1 to $p2: ".join(' ', @ppath)."\n";
    return(\@ppath, \@used_fw_edges, \@used_rc_edges);    
}

# Compute shortest distances between palindromes.
#
# Builds hash of hashes such the keys are palindromes and the value
# for a given pair of palindromes is the shortest distance from the
# first to the second in the de Bruijn graph of order k-1. Also
# computes the closest neighboring palindrome for each palindrome.
sub get_pdists($$) {
    my ($k, # integer 
	$palin_ref # ref to array of palindromes of length $k
	)=@_;
    my %ppdists;
    my @palindromes = @$palin_ref;
    my $num = scalar(@palindromes);
    my $bigval=$k*$num;
    my (%bestp, %bestd);
    for(my $i=0; $i< $num; $i++) {
	for (my $j=0; $j<$num; $j++) {
	    my $palindrome1 = $palindromes[$i];
	    my $palindrome2 = $palindromes[$j];
	    my ($distfw, $distbw);
	    if ($palindrome1 eq $palindrome2) {
		$distfw=$bigval;
		$distbw=$bigval;
	    } else {
		$distfw = $k-suffix_prefix_match($palindrome1, $palindrome2);
		$distbw = $k-suffix_prefix_match($palindrome2, $palindrome1);
	    }
	    $ppdists{$palindrome1}{$palindrome2}=$distfw;
	    $ppdists{$palindrome2}{$palindrome1}=$distbw;
	    if ($distfw != $distbw) {
		die "Forward and reverse distances are not equal!\n";
	    }
	    if ( !($bestd{$palindrome1}) or ($bestd{$palindrome1}>$distfw)) {
		$bestd{$palindrome1}=$distfw;
		$bestp{$palindrome1}=$palindrome2;
	    }
	}
    }            
    return (\%ppdists, \%bestp, \%bestd);
}

1;
