#!/usr/bin/perl 
#
# SJ Riesenfeld
# Created: 2009
# Last updated: Feb 2012
#
# A library of functions for computing MRCC libraries and other
# functions of the de Bruijn graph (with the DNA alphabet).
#
# This package is designed to work with the Palindrome_Paths package.
# The intended use of these functions is illustrated in
# design_oligomers.pl.
#
##### 
# This package is somewhat large and complicated. However, the
# functions for computing a simple MRCC library are easy to use and
# not too difficult to understand. 
#
# The complexity is introduced when flanking sequences that may be
# required in experimental applications are considered in the design
# (modifying the output so it is no longer an MRCC library). The
# functions here implement a heuristic, randomized approach for
# reducing the repetitions that occur in the k-mers bridging the
# junctions between oligomers (or constructs) and flanking sequences.
# 
# To simply compute an MRCC library, the only high-level functions
# needed are:
#
# find_eulerian_cycle, find_eulerian_pair, check_connected, dfs,
# reorder_cycle, check_degrees, start_new_cycle,
# break_cycle_into_constructs.
#####
package DNA_DeBruijnGraph_Utils;

use strict; 
use warnings;
use DNA_Manip qw( rev_comp rev_comp_list generate_all_k_mers );
use Seq_Utils qw( choose_random_elt get_matches get_head get_tail 
                  in_list remove is_prefix is_suffix suffix_prefix_match );
use POSIX qw( ceil floor modf );

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw ( find_eulerian_cycle find_eulerian_pair
                      check_connected dfs reorder_cycle
                      check_degrees start_new_cycle
                      break_cycle_into_constructs
                      build_extra_constructs build_runs_table
                      check_constructs print_construct_paths
                      print_constructs print_array
                      find_cycles build_cycle_cover
                      get_cycles_with_node add_flanking_seqs_to_constructs
                      compute_lengths
);

my $verbose = 0;

# Check subgraph is strongly connected.
#
# Returns 1 if set of edges (second arg) forms a strongly connected
# graph (on incidental ($k-1)-mers); returns 0 otherwise.
sub check_connected($$;$$) {
    my ($k,  # integer
	$kmers_ref,  # ref to array of $k-mers
	$km1mers_ref, # (optional) ref to an array of ($k-1)-mers
	$out_fh # (optional) file handle for printed output
	) = @_;
    if (!$out_fh) {
	$out_fh = \*STDOUT;
    }
    my $total_nodes;
    if ($km1mers_ref) {
	$total_nodes= scalar(@$kmers_ref);
    } else {
	$total_nodes = 4**($k-1);
    } 
    my $start = get_head(@$kmers_ref[0]);
    # print "Starting node: $start.\n";
    my $seen_fw_ref = dfs($kmers_ref, $start, 0, $total_nodes, $out_fh);
    my $seen_bw_ref = dfs($kmers_ref, $start, 1, $total_nodes, $out_fh);
    foreach my $node (@$seen_fw_ref) {
	if (!in_list($node, $seen_bw_ref)) {
	    warn "\nNot connected:  missing node $node!\n";
	    return 0;
	} else {
	    remove($node, $seen_bw_ref)
	}
    }
    if (!(scalar(@$seen_bw_ref))) {
	print $out_fh "\nStrongly connected!\n";
	return 1;
    } else {
	warn "\nNot connected:  remaining nodes: ".join(' ', @$seen_bw_ref)."\n";
	return 0;
    }
}

# DFS = Depth First Search
#
# Returns array of (k-1)-mers reached via DFS started at tail of first
# edge in array.
sub dfs($$;$$$) {
    my ($kmers_ref, # reference to an array of k-mers;
	$start_node, # start node, i.e., (k-1)-mer;
	$dir, # (optional) direction of dfs:
              # 1 to follow edges in reverse direction
              # 0 to follow edges in forward (normal) direction;
	$total_nodes, # total number of nodes (e.g., (k-1)-mers) that
		      # may be incident on edges in the array.
	$out_fh # file handle for printed output
	) = @_;
    if (!$dir) {
	$dir=0;
    }
    if (!$out_fh) {
	$out_fh = \*STDOUT;
    }
    my @seen = ($start_node);
    my @edges=get_matches($start_node, $dir+1, $kmers_ref);
    while (@edges) {
	my $edge = shift(@edges);
	my $head = get_head($edge, $dir);
	if (!(in_list($head, \@seen))) {
	    push(@seen, $head);
	    my @new_edges = get_matches($head, $dir+1, $kmers_ref);
	    unshift(@edges, @new_edges);
	    if ($total_nodes and (scalar(@seen) == $total_nodes)) {
		print $out_fh "\nDFS: All nodes visited!\n";
		last;
	    }
	} else {
	    # print "Node $head already seen.\n";
	}
    }
    return \@seen;
}   

# Check if degrees of all vertices are balanced.
#
# Returns 1 if in- and out- degree is equal for every ($k-1)-mer 
#      in the induced de Bruijn graph on the given $k-mers;
# returns 0 otherwise.
sub check_degrees ($ $;$){
    my ($k, #  integer $k;
	$seqs_ref, # a reference to an array of $k-mers.
	$out_fh # (optional) file handle for printed output.
	) = @_; 
    if (!$out_fh) {
	$out_fh = \*STDOUT;
    }
    my @km1mers = generate_all_k_mers($k-1);
    my $equal_degs = 1;
    foreach my $km1mer (@km1mers) {
	my $out_deg = scalar(get_matches($km1mer, 1, $seqs_ref));
	my $in_deg = scalar(get_matches($km1mer, 2, $seqs_ref));
	# print ''. ($k-1)."-mer $km1mer in-deg: $in_deg, out-deg: $out_deg\n";
	if ($out_deg != $in_deg) {
	    warn "\nNode $km1mer has unequal in- and out- degree: in-deg $in_deg, out-deg $out_deg.\n";
	    $equal_degs=0;
	    last;
	}	
    }
    if ($equal_degs) {
	    print $out_fh "\nEvery node has its in-degree equal to its out-degree.\n";
    }
    return $equal_degs;
}

# Find (if possible) a partitioning of the edges in the input graph
# such that one cycle is the reverse complement of the other.
#
# Returns a reference to an array of references to two cycles, 
# one containing the edges that are the reverse complements of the other,
# such that together they cover all edges in the initial array exactly once. 
# (Such a solution should exist when this is called.)
sub find_eulerian_pair($$;$$$$) {
    my ($k, # integer
	$edges_ref, # a ref to array of edges, i.e., kmers;
	$period, # (optional) integer giving periodicity with which
		 # the third optional arg should appear as the initial
		 # part of the edge label, if possible, in the forward
		 # cycle (if the desired construct length is B, not
		 # including flanking sequence, and the left flanking
		 # sequence ends with the sequence given as the third
		 # argument, then the period is probably B+2); if
		 # zero/null, this arg and the fourth arg are ignored;
		 # if the fifth arg is also provided, then the period
		 # is lengthened by 2 if both can appear within the
		 # same segment;
	$db_str_l, # (optional) two-base sequence of DNA (the left flanking sequence);
	$db_str_r, #  (optional) two-base sequence of DNA (the right flanking sequence);
	$out_fh # (optional) file handle for printed output.
	) = @_;
    if ($period and !$db_str_l) {
	die "In find_eulerian_pair: If a period is specified, a di-nucleotide string must also be specified.\n";
    }
    if (!$out_fh) {
	$out_fh = \*STDOUT;
    }
    my ($cycles_ref, $rc_cycles_ref) =build_cycle_cover($k, $edges_ref, $period, $db_str_l, $db_str_r, $out_fh);
    print $out_fh "\n".scalar(@$cycles_ref). " forward cycles in cover:\n";
    # print_cycles($cycles_ref);
    print $out_fh ''.scalar(@$rc_cycles_ref). " RC cycles in cover:\n";
    if (!check_connected($k, $cycles_ref, undef, $out_fh) or !check_connected($k, $rc_cycles_ref, undef, $out_fh)) {
	die "Have not found connected subgraphs!\n";
    } else {
	print $out_fh "Both forward and reverse-complemtary subgraphs of cycle covers are strongly connected.\n";
    }
    my $big_cycle_ref= find_eulerian_cycle($cycles_ref, $edges_ref, $period, $db_str_l, $out_fh);
    if (!$big_cycle_ref) {
	die "Could not find Eulerian cycle!\n";
    }
    print $out_fh "\nForward cycle has ". scalar(@$big_cycle_ref). " edges. \n";
    # print_array($big_cycle_ref);
    my $rc_big_cycle_ref = rev_comp_list($big_cycle_ref);
    @$rc_big_cycle_ref=reverse(@$rc_big_cycle_ref);
    print $out_fh "RC cycle has ". scalar(@$rc_big_cycle_ref). " edges. \n";
    # print_array($rc_big_cycle_ref);
    return ($big_cycle_ref, $rc_big_cycle_ref);
}


# Given a partitioning of a subset of edges into cycles, compute a
# single cycle (if possible) that contains all the edges.
#
# Returns a reference to a single cycle which covers all the edges once.
sub find_eulerian_cycle($$;$$$) {
    my ($cycles_ref, # ref to array of cycles of edges (subgraph
		     # containing all these edges must be Eulerian)
		     # such that each edge is covered once;
	$edges_ref, # ref to an array of the edges covered in the
		    # above cycles;
	# see the comments for "find_eulerian_pair" for info on the optional final three arguments.
	$period, $db_str_l, $out_fh) = @_;
    my $num_tries = 10;
    if (!$out_fh) {
	$out_fh = \*STDOUT;
    }
    print $out_fh "\nBuilding Eulerian cycle from cycle cover.\n";
    my ($start_node, $cycles_stack_ref);
    for (my $i=0; $i < $num_tries; $i++) {
	my $random_edge;
	if ($period) {
	    my @poss_edges = get_matches($db_str_l, 1, $edges_ref);
	    if (scalar(@poss_edges)) {
		$random_edge = choose_random_elt(\@poss_edges);
	    }
	}
	if (!$random_edge) {
	    $random_edge = choose_random_elt($edges_ref);
	} 
	$start_node = get_tail($random_edge);
	$cycles_stack_ref = get_cycles_with_node($cycles_ref, $start_node);
	print $out_fh "Found ". scalar(@$cycles_stack_ref)." cycles containing start node $start_node.\n";
	if (!scalar(@$cycles_stack_ref)) {
	    print $out_fh "Cannot find cycle with start node: $start_node!\n";
	    if ($i == ($num_tries - 1)) {
		return ();
	    } else {
		next;
	    }	
	} 
	last;
    }
    my $cur_cycle_ref = shift(@$cycles_stack_ref);
    my ($cur_edge, $cur_node);
    my @big_cycle;
    print $out_fh "Start node: $start_node.\n";
    while(1) {
	# take one step in current cycle
	$cur_edge = shift(@$cur_cycle_ref);
	if ($cur_edge) {
	    # print "Current edge: $cur_edge.\n";
	}
	if (!$cur_edge) {
	    print $out_fh "Finished cycle.\n";
	    # finished current cycle
	    if ($cur_cycle_ref=shift(@$cycles_stack_ref)) {
		print $out_fh "Getting next cycle off stack.\n";
		$cur_edge = shift(@$cur_cycle_ref);
	    } else {
		# no more cycles in stack
		last;
	    }
	}
	# save this step
	push (@big_cycle, $cur_edge);
	$cur_node = get_head($cur_edge);
	# look for unseen cycle containing current node
	my $new_cycles_ref = get_cycles_with_node($cycles_ref, $cur_node);
	if (scalar(@$new_cycles_ref)) {
	    print $out_fh  "Found ". scalar(@$new_cycles_ref)." new cycles containing node $cur_node.\n";	    
	    print $out_fh "Saving rest of current cycle to stack.\n";
	    unshift (@$cycles_stack_ref, $cur_cycle_ref);	
    	    print $out_fh "Saving new cycles to stack.\n";
	    unshift (@$cycles_stack_ref, @$new_cycles_ref);	
	    $cur_cycle_ref=shift(@$cycles_stack_ref);
	    print $out_fh "Taking new cycle off stack.\n";
	}
    }
    return \@big_cycle;
}

# Given a reference to an array where each entry is a reference to an
# array of edges composing a cycle, print out the cycles.
sub print_cycles($;$) {
    my $cycles_ref=shift();
    my $out_fh = shift;
    if (!$out_fh) {
	$out_fh = \*STDOUT;
    }
    foreach my $cycle_ref (@$cycles_ref) {
	print $out_fh "cycle: \n";
	print_array($cycle_ref, $out_fh);
    }
}		 

# Rotate the current "start" of the cycle to the given node (if possible).
#
# first arg is ref to cycle (array of kmers);
# second arg is node, i.e., a ($k-1)-mer;
# Returns a ref to a cycle that is a rotation of the given cycle 
#    so that the returned cycle starts with that node (if it occurs in the cycle).
sub reorder_cycle($$) {
    my ($cycle_ref, $node) = @_;
    if (!$cycle_ref) {
	return $cycle_ref;
    }
    my @cycle = @$cycle_ref;
    my @reordered_cycle;
    my $last_index= scalar(@cycle)-1;
    for (my $i=0; $i<= $last_index; $i++) {
	if (get_tail($cycle[$i]) eq $node) {
	    @reordered_cycle = @cycle[$i..$last_index];
	    push(@reordered_cycle, @cycle[0..$i-1]);
	}
    }
    return \@reordered_cycle;
}

# Given a reference to an array of cycles, find cycles that start with
# a given node.
#
# Returns ref to array of reordered cycles that contain the node so
# they start with the node.  deletes references to cycles containing
# the node from the original array of references -- destructive!
sub get_cycles_with_node($$) {
    my ($cycles_ref, # ref to array of references to cycles (arrays of
		     # kmers);
	$node # a node, i.e., a (k-1)-mer.
	) = @_;
    my @matched_cycles=();
    foreach my $cycle_ref (@$cycles_ref) {
	if (get_matches($node, 1, $cycle_ref)) {
	    # print "Found cycle containing node $node.\n";
	    remove ($cycle_ref, $cycles_ref);
	    my $reo_cycle_ref = reorder_cycle($cycle_ref, $node);
	    my $index = int(rand(scalar(@matched_cycles)));	
	    splice (@matched_cycles, $index, 0, $reo_cycle_ref);
	}
    }     
    return \@matched_cycles;
}

# Partition the edges of the graph into two sets of cycles: such that
# each cycle in the first set has a reverse-complementary cycle in the
# second set.
#
# Returns array of two references, one to an array of references to
# forward cycles, the other to an array of reference to the
# reverse-complementary cycles, such that all edges are covered once.
sub build_cycle_cover($$;$$$$) {
    my ($k, # integer
	$edges_ref, # reference to array of edges, i.e., kmers
		    # (subgraph containing these edges should be
		    # Eulerian);
    # see "find_eulerian_pair" for info on optional final four args.
	$period, $db_str_l, $db_str_r, $out_fh) = @_;
    if (!$out_fh) {
	$out_fh = \*STDOUT;
    }
    print $out_fh "\nBuilding cycle cover...\n";
    my @edges = @$edges_ref;
    my @cycles;
    my @rc_cycles;
    my @total_hits=(0,0,0,0);  # in order left-2, left-1, right-2, right-1
    if (!$out_fh) {
	$out_fh = \*STDOUT;
    }
    while (@edges) {
	print $out_fh "There are ".scalar(@edges)." edges left to cover.\n";
	my ($cycle_ref, $rc_cycle_ref, $hit_counts_ar)=find_cycles($k, \@edges, $period, $db_str_l, $db_str_r, 
								   $out_fh);
	if ($cycle_ref and $rc_cycle_ref) {
	    print $out_fh  "Found forward cycle.\n"; #.join(' ', @$cycle_ref)."\n";
	    print $out_fh "Found rc cycle.\n"; #.join(' ', @$rc_cycle_ref)."\n";
	    if ($period) {
		foreach my $i (0..3) {
		    $total_hits[$i] += $hit_counts_ar->[$i];
		}
	    }
	    push(@cycles, $cycle_ref);
	    push(@rc_cycles, $rc_cycle_ref);
	} else {
	    die "Cannot find cycle!\n";
	}
    }
    print $out_fh "\nFinished building cycle cover.\n";
    if ($period) {
	print $out_fh "Total number of times left flanking seq $db_str_l was inserted at correct periodicity: ".$total_hits[0]."\n";
	print $out_fh "Total number of times just the second base of left flanking seq was inserted at correct periodicity: ".$total_hits[1]."\n";
	print $out_fh "Total number of times right flanking seq $db_str_r was inserted at correct periodicity: ".$total_hits[2]."\n";
	print $out_fh "Total number of times just the second base of right flanking seq was inserted at correct periodicity: ".$total_hits[3]."\n";

    }
    return (\@cycles, \@rc_cycles);
}

# Find two cycles in the graph such that one is the reverse-complement
# of the other.
#
# Returns an array whose first two elements are references to a pair of (somewhat randomly chosen) cycles 
#    containing a subset of the given edges, such that one cycle contains the reverse complement of the edges in the other
#    (if such a pair of cycles exists);
### may change the following:!!!
# If the third and fourth args are non-null/non-zero, then the returned array contains two additional elements:
#     the third element is the total number of times the sequence specified by the fourth arg was inserted at 
#     the correct periodicity;
#     the fourth element is the total number of times just the second base of the sequence was inserted at the correct periodicity.
# Destructive to the given array of edges!
sub find_cycles($$;$$$$) {
    my ($k, # integer
	$edges_ref, # ref to array of edges, i.e., kmers, (subgraph of
		    # edges should be Eulerian);
        # see "find_eulerian_pair" for info on the optional four args.
	$period, $db_str_l, $db_str_r, $out_fh) = @_;
    if (!$out_fh) {
	$out_fh = \*STDOUT;
    }
    if ($period and !$db_str_l) {
	die "In find_cycles: If a period is specified, a di-nucleotide left flanking sequence must also be specified.\n";
    }
    
    my ($sb_l, $fb_r, $fb_l, $sb_r);
    if ($db_str_l) {
	$fb_l = substr($db_str_l, 0, 1);
	$sb_l = substr($db_str_l, 1, 1);
    }
    if ($db_str_r) {	
	$fb_r = substr($db_str_r, 0, 1);
	$sb_r = substr($db_str_r, 1, 1);
    }    

    if (!scalar(@$edges_ref)) {
	return;
    }
    my ($path_period, $check_time, $check_later, $check_latest);
    if ($period) {
	$path_period = $period - $k +1;
	$check_time = $path_period -$k+1;
	$check_later = $check_time+1;
	$check_latest = $check_time + 2;
	if ($check_time < 0) {
	    die "Period $period given is too short!\n";
	}
    }
    print $out_fh "\nFinding pair of fw and rc cycles...\n";
    my @poss_start_edges;
    my ($hit1_count, $hit2_count)=(0,0);
    my ($hit1_r_count, $hit2_r_count)=(0,0);
    my ($hit1_r, $hit2_r)=(0,0);
    my ($hit1_later, $hit2_later)=(0,0);
    my ($hit1_latest, $hit2_latest)=(0,0);
    my ($hit1_prev, $hit1_cur, $hit2_prev, $hit2_cur) = (0,0,0,0);
    my ($partial_hit, $partial_hit_later, $partial_hit_latest)=(0,0,0);

    my $start_edge;    
    my $phase = 0;  # where are we in the current stretch of sequence, relative to the period

    if ($period) {
	@poss_start_edges = get_matches($db_str_l, 1, $edges_ref);
	if (scalar(@poss_start_edges)) {
	    $hit2_prev = 1;
	    $start_edge = choose_random_elt(\@poss_start_edges);
	    print $out_fh "Using 2-base left flanking sequece $db_str_l to choose starting edge $start_edge of cycle.\n";
	} else {
	    print $out_fh "No edges starting with 2-base flanking sequence $db_str_l!\n";
	    @poss_start_edges = get_matches($sb_l, 1, $edges_ref);
	    if (scalar(@poss_start_edges)) {
		$phase = 1;  # push the phase one base ahead since we want the period to be one base shorter
		$hit1_prev=1;
		$start_edge = choose_random_elt(\@poss_start_edges);
		print $out_fh "Using second base $sb_l of flanking sequence to choose starting edge $start_edge of cycle.\n";		
	    } else {
		print $out_fh "No edges starting with second base $sb_l of flanking sequence.\n";
	    }
	}
    }
    
    if (!$start_edge) {
	$start_edge=choose_random_elt($edges_ref);
    }
    my $rc_start_edge = rev_comp($start_edge);
    
    remove ($start_edge, $edges_ref);
    remove ($rc_start_edge, $edges_ref);
    my @cycle_edges = ($start_edge);
    my @rc_cycle_edges = ($rc_start_edge);
    my $end_node = get_tail($start_edge);    
    my $head = get_head($start_edge);

    my $over=0;  
    if ($period) {
	print $out_fh "Initially:  path_period: $path_period, phase: $phase.\n";
    }
    while ($head ne $end_node) {
	my @matches=get_matches($head, 1, $edges_ref);
	# print "Matches: ".join(' ', @matches)."\n";
	if (!(scalar(@matches))) {
	    die "In find_cycles:  No edges out of $head!\n";
	}
	my $cur_edge;
	if ( $period ) {
	    $phase = ($phase + 1) % $path_period;
	    print $out_fh "Phase: $phase. ";
	    if ( $phase==0 ) {  # finish a stretch of path length $path_period if it begins with $db_str_l;
		                # if stretch instead begins with $sb_l, its length is $path_period-1;
		                # if it begins with a different sequence, its length is $path_period-2.
		print $out_fh "Finishing a stretch -- checking to see if we should try to add right flanking seq $db_str_r...\n";
		print $out_fh "   Hit2_prev: $hit2_prev, hit1_prev: $hit1_prev.\n";
		print $out_fh "   Hit2_cur: $hit2_cur, hit1_cur: $hit1_cur, over: $over.\n";
		print $out_fh "   Hit1_later: $hit1_later, hit2_later: $hit2_later.\n";
		print $out_fh "   Hit1_latest: $hit1_latest, hit2_latest: $hit2_latest.\n";
		if (!$hit2_cur and (!$hit1_cur or $hit1_later or $hit2_later or $hit1_latest or $hit2_latest)) {
		    print $out_fh "Trying to add right flanking sequence at end of stretch\n";
		    # don't want phase to increase from 0 while we add flanking sequence on the right
		    my @fb_r_matches = get_matches($fb_r, 2, \@matches);
		    if (scalar(@fb_r_matches)) {			
			$cur_edge = shift(@fb_r_matches);
			print $out_fh "Found match $cur_edge for first base $fb_r of right flanking seq.\n";
			$head = process_edge ($cur_edge, $edges_ref, \@cycle_edges, \@rc_cycle_edges);
			$hit1_r=1;
			$hit1_r_count++;
			$over=0;
			$hit1_cur=0;   # the start is now later in the phase
			print $out_fh  "   Setting: over: $over.\n";  # .", hit1_cur: $hit1_cur.\n";
			if ($head eq $end_node) {   # cycle will end 
			    $hit1_r_count++;
			    print $out_fh  "Hit1_r_count increased: $hit1_r_count.\n";
			    last;
			} elsif (!$hit1_later and !$hit2_later) {  
			    print $out_fh "Checking for a match to the second base of the right flanking seq.\n";
			    my $hit2_str_r = $head.$sb_r;
			    my @sb_r_matches = get_matches($hit2_str_r,0, $edges_ref);
			    if (scalar(@sb_r_matches)) {
				$cur_edge = shift(@sb_r_matches);
				$head = process_edge ($cur_edge, $edges_ref, \@cycle_edges, \@rc_cycle_edges);
				$hit2_r=1;
				$hit1_r=0;
				$hit2_r_count++;       
				$hit1_r_count--;
				print $out_fh "Found match $cur_edge for second base $sb_r of right flanking seq.\n";
				print $out_fh "Hit2_r_count increased: $hit2_r_count.\n";
				if ($head eq $end_node) {
				    last;
				}
			    } else {
				print $out_fh "No match found for second base $sb_r of right flanking seq.\n";
				print $out_fh "Hit1_r_count increased: $hit1_r_count.\n";
			    }
			} else {
			    print $out_fh "Skipping check for match to second base of the right flanking seq:\n".
				"   hit1_later: $hit1_later, hit2_later: $hit2_later\n";
			}
			@matches=get_matches($head, 1, $edges_ref);
			if (!(scalar(@matches))) {
			    die "In find_cycles:  No edges out of $head!\n";
			}			
			$cur_edge=undef;
			print $out_fh "Updated values:   Hit2_r: $hit2_r, hit1_r: $hit1_r\n";
		    }
		}
		# since the phase got lengthened, these later hits turn into real current hits
		if ($hit2_r) {
		    if ($hit2_latest) {
			$hit2_cur=1;
		    } elsif ($hit1_latest) {
			$hit1_cur=1;
		    }
		} elsif ($hit1_r) {
		    if ($hit2_later) {
			$hit2_cur=1;
		    } elsif ($hit1_later) {
			$hit1_cur=1;
		    }
		}
		    
		# set the current phase correctly according to how this stretch began (with a 2-, 1-, or 0-base hit)		
		if ($hit1_cur) {
		    print $out_fh "\nSkipping ahead one in the phase.\n";
		    $phase = ($phase+1) % $path_period;
		} elsif (!$hit2_cur) {   
		    print $out_fh  "Skipping ahead two in the phase.\n";
		    $phase = ($phase + 2) % $path_period;
		    if (!$hit1_prev and !$hit2_prev and !$hit1_r and !$hit2_r) {
			# keep track: if we have just had one stretch of length $path_period-2 that began randomly,
                        # then we don't need another; will end current stretch when we find a good next starting edge.
			$over=1; 			
			print $out_fh "Setting over to $over.\n"
		    }
		}
		if ( $hit2_prev ) {
		    $hit2_count++;
		    print $out_fh  "Hit2_count increased: $hit2_count.\n";
		} elsif ($hit1_prev) {
		    $hit1_count++;
		    print $out_fh  "Hit1_count increased: $hit1_count.\n";
		}
		$hit1_prev=$hit1_cur;
		$hit2_prev=$hit2_cur;
		($hit1_later, $hit2_later)=(0,0);
		($hit1_latest, $hit2_latest)=(0,0);
		($hit1_r, $hit2_r) = (0,0);
		($hit1_cur, $hit2_cur)=(0,0);
		print $out_fh "Values reset to:\n";
		print $out_fh "   Hit1_r: $hit1_r, hit2_r: $hit2_r, hit2_prev: $hit2_prev, hit1_prev: $hit1_prev.\n";
		print $out_fh "   Hit2_cur: $hit2_cur, hit1_cur: $hit1_cur, over: $over.\n";
		print $out_fh "   Hit1_later: $hit1_later, hit2_later: $hit2_later.\n";
		print $out_fh "   Hit1_latest: $hit1_latest, hit2_latest: $hit2_latest.\n";
	    }
	    
	    if ( $over or ($phase==$check_time) or $partial_hit 
		 or $partial_hit_later or $partial_hit_latest
		 or (!$hit2_cur and (($phase==$check_later) 
				     or (!$hit1_later and ($phase==$check_latest)))) ) {
		
		if (!$partial_hit and !$partial_hit_later and !$partial_hit_latest) {
		    print $out_fh "Looking for match to first base of left flanking sequence.\n";
		    my @fb_matches = get_matches($fb_l, 2, \@matches);
		    if (scalar(@fb_matches)) {
			my $fb_match = shift(@fb_matches);
			my $poss_head = get_head($fb_match);
			print $out_fh "Found possible match for first base of left flanking sequence -- checking if it leads to second base match.\n";
			if ($poss_head eq $end_node) {   # cycle will end 
			    $cur_edge = $fb_match;
			} else {   # check to see if the first base match will lead to a match for $db_str_l
			    my $hit2_str = $poss_head.$sb_l;
			    my @next_matches = get_matches($hit2_str,0, $edges_ref);
			    if (scalar(@next_matches)) {
				$cur_edge = $fb_match;
				print $out_fh "Found match $cur_edge for first base $fb_l of left flanking seq $db_str_l that should lead ".
				    "to match for second base $sb_l.\n";
				if ($phase==$check_latest) {
				    $partial_hit_latest=1;
				} elsif ($phase==$check_later) {
				    $partial_hit_later=1;
				} else {				
				    $partial_hit=1;
				}
			    }
			}
		    }
		}
		if (!$cur_edge  # either $partial_hit, $partial_hit_later, or $partial_hit_latest is true or we didn't find a match for $db_str_l
		    and (($phase != $check_time) or $over)) {     # don't bother to find second base match in phase == $check_time
		    my @sb_matches = get_matches($sb_l, 2, \@matches);
		    print $out_fh "Looking for match to second base of left flanking sequence.\n";
		    if (scalar(@sb_matches)) {
			$cur_edge = shift(@sb_matches);
			print $out_fh "Found match $cur_edge for second base $sb_l.\n";
			if ($partial_hit) {
			    print $out_fh "Finished partial hit.\n";
			    $partial_hit=0;
			    $hit2_cur=1;
			    $over=0;
			} elsif ($partial_hit_later) {
			    print $out_fh "Finished later partial hit.\n";
			    $partial_hit_later=0;
			    $hit2_later=1;
			} elsif ($partial_hit_latest) {
			    print $out_fh "Finished latest partial hit.\n";
			    $partial_hit_latest=0;
			    $hit2_latest=1;
			} else {
			    if ($phase==$check_latest) {
				$hit1_latest=1;
			    } elsif ($phase==$check_later) {
				$hit1_later=1;
			    } else {
				$over=0;							
				$hit1_cur=1;
			    }
			}
		    }
		}
	    }
	}
	if (!$cur_edge) {
	    $cur_edge=choose_random_elt(\@matches);
	    print $out_fh "Choosing random edge: $cur_edge.\n";
	    $partial_hit=0;
	}	    
	
	$head = process_edge ($cur_edge, $edges_ref, \@cycle_edges, \@rc_cycle_edges);
	
    } # end while loop over unfinished cycle
    print $out_fh "Finished cycle.\n";
    my @hit_counts;
    if ($period) {  # final update of hits
	$phase = ($phase+1)%$path_period;
	if ( !$phase ) {
	    if($hit1_prev ) {
		$hit1_count++;
	    } elsif ($hit2_prev) {
		$hit2_count++;
	    }
	}
	print $out_fh "Number of times left flanking seq $db_str_l was inserted at the correct periodicity: $hit2_count.\n";
	print $out_fh "Number of times just the second base $sb_l was inserted at the correct periodicity: $hit1_count.\n";
	print $out_fh "Number of times right flanking seq $db_str_l was inserted at the correct periodicity: $hit2_r_count.\n";
	print $out_fh "Number of times just the first base $fb_r was inserted at the correct periodicity: $hit1_r_count.\n";
	@hit_counts = ($hit2_count, $hit1_count, $hit2_r_count, $hit1_r_count);
    }
    return (\@cycle_edges, \@rc_cycle_edges, \@hit_counts);
}
	
# Process a given edge by removing it and its reverse complement from the
# given edge set and adding them to the corresponding given cycles.
sub process_edge($$$$) {
    my ($cur_edge, $edges_ar, $cycle_edges_ar, $rc_cycle_edges_ar)=@_;
    my $rc_cur_edge = rev_comp($cur_edge);
    if (!in_list($rc_cur_edge, $edges_ar)) {
	die "In process_edge:  Rev comp of edge $cur_edge is missing!\n";
    }
    remove($cur_edge, $edges_ar);
    remove ($rc_cur_edge, $edges_ar);
    push(@{$cycle_edges_ar}, $cur_edge);
    unshift(@{$rc_cycle_edges_ar}, $rc_cur_edge);
    
    my $head = get_head($cur_edge);
    return $head;		 
}

# This function assesses exhaustively how well a cycle periodically
# incorporates the appropriate flanking sequence bases. It builds a
# table containing the best "run" value starting at each index in the
# given cycle, where a run is a sequence of appropriately spaced
# incorporations. (The appearance of a particular base may or may not
# count as an incorporation of the flanking sequence in the optimal
# run; this decision is a branching point in the exhaustive search.)
#
# Builds a table of where the values are run values and the keys are indices in the cycle of kmers;
#     the run value for each index is value for the run starting at that index,
#     where one segment of the run is a sequence of $path_len+2 k-mers
#     such that the first kmer in the sequence starts with the di-nucleotide string given
#     or, if the flag is not set, a segment may also be a sequence of $path_len+1 k-mers
#     such that the first kmer in the sequence starts with the second base of the string given.
#     If the fourth arg is provided, a segment may similarly take into account the right flanking sequence.
#     The "value" for each segment in the run is number of base pairs longer than $path_len of that segment
#     (see the sub update_run_val() ).
# Returns an array of 3 values:  
#     the first is a ref to the table of runs values built; the second is a ref to the table of
#     runs paths built, the third is the max value over the values of all runs stored.
sub build_runs_table($$$$;$$$) {
    my ($cycle_ar, # ref to an array containing a cycle of kmers;
	$path_len,  # second arg is the base length of a path, i.e. a
		    # sequence of kmers, without the kmers starting
		    # (or ending) in $str_l (or $str_r);
	$str_l,  # third arg is a string, assumed to be a 2-base DNA
		 # sequence (left flanking sequence);
	$str_r, #  a 2-base DNA sequence (for the right flanking
		#  sequence)
	$flag_2bp, # (optional) flag; if zero/undefined, then all
		   # segments which save at least 1 base pair count in
		   # a run; otherwise, only segments which save at
		   # least 2bp count.
	$out_fh # (optional) file handle for printed feedback
	)= @_;
    if (scalar(@{$cycle_ar}) < ($path_len+2)) {
	die "In build_runs_table:  the length of the cycle should be at least 2+ the length of the path.\n";
    }
    if (!$out_fh) {
	$out_fh = \*STDOUT;
    }
    my ($second_base_l, $first_base_r);
    if (length($str_l) !=2) {
	die "The left string given to build_runs_table should have length 2!\n";
	if ($str_r and length($str_r !=2)) {
	    die "The right string given to build_runs_table should have length 2!\n";
	}
    }
    my (%runs_vals, %runs_paths);
    my $max_val=0;
    my $max_val_index;
    
    print $out_fh "\nBuilding runs table for cycle length: ". scalar(@{$cycle_ar}). 
	", minimum path length: $path_len,\n".
	"   left flanking seq: $str_l".($str_r ? ", and right flanking seq: $str_r":'').".\n";
    
    foreach (my $i=0; $i < scalar(@{$cycle_ar}); $i++) {
	set_best_run($i, $cycle_ar, \%runs_vals, \%runs_paths, $path_len, $str_l, $str_r, $flag_2bp, $out_fh);
	if ($runs_vals{$i} > $max_val) {
	    $max_val = $runs_vals{$i};
	    print $out_fh "Found new best run value $max_val at index $i.\n";
	    $max_val_index=$i;
	}
    }
    print $out_fh "\nPrinting runs table:\n";
    my @runs_keys = sort { $a <=> $b } keys(%runs_vals);
    if (scalar(@runs_keys) != scalar(@{$cycle_ar})) {
	die "Badly formed runs table!\n";
    }
    foreach (my $i=0; $i<scalar(@{$cycle_ar}); $i++)  {
	if ($i != $runs_keys[$i]) {
	    die "Badly formed runs table!\n";
	}
    }
    my @single_runs = grep { $runs_vals{$_} == 1 } @runs_keys;
    my @mult_runs = grep {$runs_vals{$_} > 1} @runs_keys;
    @mult_runs = sort { $runs_vals{$b} <=> $runs_vals{$a} } @mult_runs;
    print $out_fh "   The following indices have run value 1:\n";
    print_array (\@single_runs, $out_fh);
    print $out_fh "   The following indices have run values greater than 1:\n";
    foreach my $index (@mult_runs) {
	print $out_fh "Start index: $index,  run value: ". $runs_vals{$index}."\n";
    }
    print $out_fh "   All other run values are 0\n";
    print $out_fh "A max run value index: $max_val_index, with run value: ".$runs_vals{$max_val_index}."\n".
	"Path indices: ".join(' ', @{$runs_paths{$max_val_index}})."\n";
    return (\%runs_vals, \%runs_paths, $max_val);
}
    
# Does a depth first search to set the best run value starting at a
# given index.  Parents and children refer to branching points in the
# search. They are represented by indices that are structured so the
# indices actually fall out of the range of $cycle_ar, but the indices
# actually written in the runs tables are corrected to match those of
# $cycle_ar.  The index given as an argument is expected to fall in
# the range of $cycle_ar indices.
sub set_best_run($$$$$$;$$$) {
    my ($index, $cycle_ar, $runs_vals_hr, $runs_paths_hr, $path_len, $db_str_l, 
	$db_str_r, $flag_2bp, $out_fh)=@_;
    if (!$out_fh) {
	$out_fh = \*STDOUT;
    }
    print $out_fh "\nSetting the best run value for index $index, corresponding to kmer ".$cycle_ar->[$index].".\n";

    if ($runs_vals_hr and defined($runs_vals_hr->{$index})) {
	print $out_fh "Value already set.\n";
	return;
    }
    my $cycle_len = scalar(@{$cycle_ar});
    my $max_poss_index = ($index + $cycle_len) - 1;
    
    # find children of current index
    my $indices_to_check_ar = get_indices_to_check($index, $cycle_ar, $path_len, $db_str_l, $db_str_r, $flag_2bp);
    my @indices_to_check = @{$indices_to_check_ar};

    if (!scalar(@indices_to_check)) {
	$runs_vals_hr->{$index}=0;
	$runs_paths_hr->{$index}=[];
	"No children, so setting value to ".$runs_vals_hr->{$index}.".\n";
	return;
    }
    # initial setup of children edges of this index
    my (%children, %parents);
    $children{$index}=$indices_to_check_ar;
    print $out_fh "The children of index $index: ".join(' ', @{$children{$index}})."\n";

    # push the children on a stack and then check each child for the best run value starting from that child;
    # when the values for all children have been computed, then compute the value for the parent.
    foreach my $child (@indices_to_check) {
	push(@{$parents{$child}},$index);    
	print $out_fh "The parents of index $child: ".join(' ', @{$parents{$child}})."\n";
    }
    my @cur_path =($index);        
    while(scalar(@indices_to_check)) {	
	print $out_fh "Current path: ".join(' ', @cur_path)."\n";	
	print $out_fh "Current list of indices to check: ".join(' ', @indices_to_check)."\n";
	my $cur_index = shift(@indices_to_check);	
	print $out_fh "Current index: $cur_index, corresponding kmer: ".$cycle_ar->[$cur_index % $cycle_len]."\n";
	if ($cur_index > $max_poss_index) {
	    die "In set_best_run, looped over entire cycle!  Current index; $cur_index, start index: $index\n"
	}	
	my $next_indices_ar;
	
	if (!( $runs_vals_hr and defined($runs_vals_hr->{$cur_index % $cycle_len}) )) {
	    # the best run value has not already been computed for this index
	    $next_indices_ar =get_indices_to_check($cur_index, $cycle_ar, $path_len, $db_str_l, $db_str_r, $flag_2bp);
	    
	    $children{$cur_index}=$next_indices_ar;
	    if (!scalar(@{$next_indices_ar})) {
		# cur_index is a leaf, i.e., no run starts from it
		my $real_cur_index = $cur_index % $cycle_len;
		$runs_vals_hr->{$real_cur_index}=0;
		$runs_paths_hr->{$real_cur_index}=[];
		print $out_fh "$cur_index has no children; setting real cur index $real_cur_index runs value to 0\n";
	    } 
	} 
	if (!$next_indices_ar or !scalar(@{$next_indices_ar})) {
	    # the current index is done being processed, either because it has been previously processed or
	    # because it is a leaf;
	    print $out_fh "Current index has already been processed.\n";
	    while(scalar(@cur_path)) {
		# check if all children of most recent index in path are done being processed; 
		if (!scalar(@indices_to_check) or !in_list($indices_to_check[0], $children{$cur_path[0]})) {
		    print $out_fh "Updating the value of the top index ".$cur_path[0]." in the current path stack.\n";
		    # if so, update that parents value and pop parent off the cur_path stack		    
		    update_run_val($cur_path[0], \%children, $runs_vals_hr, $runs_paths_hr, $cycle_ar, $path_len,
			$out_fh);
		    shift(@cur_path);
		} else {
		    last;
		}
	    }
	} else {	   
	    print $out_fh "Current index has children: ".join(' ', @{$children{$cur_index}})."\n";
	    # the current index has children which must be explored
	    foreach my $child (@{$next_indices_ar}) {
		push(@{$parents{$child}}, $cur_index);
		print $out_fh "The parents of child $child are: ".join(' ', @{$parents{$child}})."\n";
	    }
	    # put children of current index on indices_to_check stack and put cur_index on cur_path stack
	    unshift(@cur_path, $cur_index);
	    unshift(@indices_to_check, @{$next_indices_ar});
	}
    }   
}

# Updates the best run value for a parent based on the values set for the children.
sub update_run_val($$$$$$;$) {
    my ($parent, $children_hr, $runs_vals_hr, $runs_paths_hr, $cycle_ar, $path_len, $out_fh) = @_;
    if (!$out_fh) {
	$out_fh = \*STDOUT;
    }
    print $out_fh "Updating the run values for index $parent.\n";
    my $best_val=0;
    my @best_path;    
    foreach my $child (@{$children_hr->{$parent}}) {
	my $real_child = $child % scalar(@{$cycle_ar});
	my $dist = ($real_child - $parent) % @{$cycle_ar};
	my $incr = $dist - $path_len;
	my $val = $incr + $runs_vals_hr->{$real_child};
	if ($val == 0) {
	    die "In update_run_val:  should not have value of run equal to 0!\n";
	}
	print $out_fh "Child: $child, real child: $real_child, has distance $dist, incremental value $incr, total value $val, path: "
	    .join(' ', @{$runs_paths_hr->{$real_child}})."\n";
	if ($val > $best_val) {
	    $best_val=$val;
	    @best_path = @{$runs_paths_hr->{$real_child}};
	    unshift(@best_path, $real_child);	    
	}
    }
    
    my $real_parent=$parent%scalar(@{$cycle_ar});
    $runs_vals_hr->{$real_parent} = $best_val;
    $runs_paths_hr->{$real_parent} = \@best_path;
    print $out_fh "Setting the values for parent $parent, real parent $real_parent, to run val: ".
	$runs_vals_hr->{$real_parent}." and path: ".
	join(' ', @{$runs_paths_hr->{$real_parent}})."\n";
    
}

# Figure out the children for the given index (i.e., branch point).
# This function assumes that
#    $str_l does not contain repeated bases
#    (and likewise for $str_r).
# Finds potential starting indices for next segment of a run.
sub get_indices_to_check($$$$;$$) {
    my ($index, $cycle_ar, $path_len, $db_str_l, 
	$db_str_r, $flag_2bp)=@_;
    my ($sb_l, $fb_r);
    if ($db_str_l) {
	$sb_l = substr($db_str_l, 1, 1);
    }
    if ($db_str_r) {
	$fb_r = substr($db_str_r, 0, 1);
    }    
    my $cycle_len =scalar(@{$cycle_ar});
    my $index0 = $index % $cycle_len;
    my $index_pp = ($index + $path_len) % $cycle_len;
    my $index_pp1 = ($index_pp +1) % $cycle_len;
    my $index_pp2 = ($index_pp1+1) % $cycle_len;
    my $index_pp3 = ($index_pp2+1) % $cycle_len;
    my $index_pp4 = ($index_pp3+1) % $cycle_len;

    my @indices_to_check=();
    if (is_prefix($db_str_l, $cycle_ar->[$index0])) {
	# start with left flanking seq, 2bp
	push(@indices_to_check, $index_pp2);
	if($db_str_r) {
	    if (is_suffix($db_str_r, $cycle_ar->[$index_pp3])) {
		# 2bp, right
		push(@indices_to_check, $index_pp4);
	    } elsif (is_suffix($fb_r, $cycle_ar->[$index_pp2])) {
		#1bp, right
		push(@indices_to_check, $index_pp3);
	    }
	}
    } elsif (is_prefix($sb_l, $cycle_ar->[$index0])) {  
	# 1bp, left flanking seq
	if (!$flag_2bp) {
	    push(@indices_to_check, $index_pp1);
	}
	if ($db_str_r) {
	    if (is_suffix($db_str_r, $cycle_ar->[$index_pp2])) {
		# 2bp, right flanking seq
		push(@indices_to_check, $index_pp3);
	    } elsif (is_suffix($fb_r, $cycle_ar->[$index_pp1])) {
		# this saves 2bp so use this even if !$flag_2bp
		# 1bp, right flanking seq
		push(@indices_to_check, $index_pp2);
	    }	
	} 
    } else {  # no left flanking seq
	if ($db_str_r) {
	    if (is_suffix($db_str_r, $cycle_ar->[$index_pp1])) {
		# 2bp, right flanking seq
		push(@indices_to_check, $index_pp2);
	    } elsif (!$flag_2bp and is_suffix($fb_r, $cycle_ar->[$index_pp])) {
		# 1bp, right flanking seq
		push(@indices_to_check, $index_pp1);
	    }
	}
    }    
    # my @corrected_indices;
    # while (my $index = pop(@indices_to_check)) {
    # my $cor_index = check_flip($index, scalar(@{$cycle_ar}), $flipped_sr);
    # unshift(@corrected_indices, $cor_index);
    # }
    # return(\@corrected_indices);
    return(\@indices_to_check);
}


# returns 1 if not over-looped
sub check_index_bnds($$$) {
    my ($index, $start, $flipped) =@_;
    if ((!$flipped and ($index >= $start)) or ($flipped and ($index < $start))) {
	return 1;
    } else {
	return 0;
    }		     
}

# check whether the index has extended beyond the cycle length
sub check_flip($$;$) {
    my ($index, $len, $flipped_sr) =@_;
    if ($index >= $len) {
	if ($flipped_sr) {
	    ${$flipped_sr}=1;
	}
    } elsif ($index < 0) {
	if ($flipped_sr) {
	    ${$flipped_sr}=0;
	}
    }
    $index = $index % $len;
    return $index;
}

sub split_final_segment($$$$;$) {
    my ($path_ar, $default_path_len, $db_str_l, $db_str_r, $out_fh) = @_;
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
    my @rem_path;
    my @ret_path;
    # my $special_flag=0;
    my @path= @{$path_ar};
    my $path_len = scalar(@path);
    print $out_fh "In split_final_segment: input path of length $path_len: ".join(' ', @path)."\n";

    my $split=1;
    if ($path_len < $default_path_len) {
	$split=0;
	@rem_path=@path;	
    } elsif ($path_len == $default_path_len) {
	$split=0;
	@ret_path=@path;
    } elsif ($path_len > $default_path_len) {	
	if (is_prefix($db_str_l, $path[0])) {
	    if (($path_len == ($default_path_len+1)) and !is_suffix($fb_r, $path[$#path])) {
		@rem_path=@path[$#path..$#path];
		@ret_path=@path[0..($default_path_len-1)];
	    } elsif (($path_len == ($default_path_len+3)) and !is_suffix($fb_r, $path[$#path])) {
		@rem_path=@path[$#path..$#path];
		@ret_path=@path[0..($default_path_len+1)];	    
	    } 
	} elsif (is_suffix($db_str_r, $path[$#path])) {
	    if (($path_len == ($default_path_len+1)) and !is_prefix($sb_l, $path[0])) {
		@rem_path=@path[0..0];
		@ret_path=@path[1..$#path];
	    } elsif (($path_len == ($default_path_len+3)) and !is_prefix($sb_l, $path[0])) {
		@rem_path=@path[0..0];
		@ret_path=@path[1..$#path];	    
	    } 
	} elsif ($path_len == ($default_path_len+2)) {
	    if (is_prefix($sb_l, $path[0]) and !is_suffix($fb_r, $path[$#path])) {
		@rem_path = @path[$#path..$#path];
		@ret_path=@path[0..($default_path_len)];
	    } elsif (is_suffix($fb_r, $path[$#path]) and !is_prefix($sb_l, $path[0])) {
		@rem_path=@path[0..0];
		@ret_path=@path[1..$#path];
	    } elsif ( !is_prefix($sb_l, $path[0]) and !is_suffix($fb_r, $path[$#path])) {
		@rem_path=@path[$default_path_len..$#path];
		@ret_path=@path[0..($default_path_len-1)];
	    }
	} elsif ( (($path_len == ($default_path_len+1)) and !is_prefix($sb_l, $path[0])
		   and !is_suffix($fb_r, $path[$#path]))
		  or ($path_len > ($default_path_len+1)) ) {
	    @rem_path=@path[$default_path_len..$#path];
	    @ret_path=@path[0..($default_path_len-1)];
	} 
	if (!scalar(@ret_path)) {
	    $split=0;
	    @ret_path=@path;
	}
    }
    if ($split) {
	print $out_fh "Split paths as follows:\n";
    } else {
	print $out_fh "Did not split paths:\n";
    }
    if (scalar(@ret_path)) {
	print $out_fh "ret_path: ".join(' ', @ret_path)."\n";
    } else {
	print $out_fh "ret_path not set.\n";
    }
    if (scalar(@rem_path)) {
	print $out_fh "rem_path: ".join(' ', @rem_path)."\n";
    } else {
	print $out_fh "rem_path not set.\n";
    }
    
    return (\@ret_path, \@rem_path);    
}

# Compute a connected path in the cycle from the $beg_index to the $end_index.
# $start_index marks the beginning of the cycle.
# (Needed in case the path wraps around.)
#   assumes indices are both > 0 and <= $#cycle_ar
#   assumes that it is ok if $beg_index == $start_index but not if $end_index == $start_index
sub get_conn_path($$$$;$) {
    my ($cycle_ar, $beg_index, $end_index, $start_index,
	$out_fh) = @_;
    if (!$out_fh) {
	$out_fh = \*STDOUT;
    }
    my @cur_path;
    print $out_fh "In get_conn_path, getting path in cycle from index $beg_index (".$cycle_ar->[$beg_index].
	") to index $end_index (".$cycle_ar->[$end_index].");  start_index: $start_index.\n";
    my $cycle_len = scalar(@{$cycle_ar});
    my $short_flag = 0;
    if ($beg_index == $end_index) {
	if ($end_index == $start_index) {
	    die "In get_conn_path:  cannot have all three given indices equal (they are all $end_index)!\n";
	} else {
	    @cur_path=$cycle_ar->[$beg_index];
	}
    } elsif ($beg_index < $end_index) {
	if (($beg_index < $start_index) and ($end_index >= $start_index)) {
	    $end_index = $start_index-1;
	    $short_flag=1;
	    print $out_fh "In get_conn_path:  Resetting end index to $end_index, and setting short flag.\n";
	}
	@cur_path = @{$cycle_ar}[$beg_index..$end_index];
    } else {
	@cur_path = @{$cycle_ar}[$beg_index..($cycle_len-1)];
	if ($end_index >= $start_index) {
	    $short_flag=1;
	    $end_index = $start_index-1;
	    print $out_fh "In get_conn_path:  Resetting end index to $end_index, and setting short flag.\n";
	}	    
	if ($end_index >= 0) {
	    push(@cur_path, @{$cycle_ar}[0..$end_index]);	    	    
	}
    }
    return (\@cur_path, $short_flag);
}

# Check the trade-off in starting a new run at some point after the
# given $beg_index or just making a construct starting from that index.
#
# assumes $beg_index <= $max_poss_index
sub check_for_run($$$$$;$) {
    my ($runs_vals_hr, $beg_index, $path_len, $cycle_len, $max_poss_index, 
	$out_fh)=@_;    
    if (!$out_fh) {
	$out_fh = \*STDOUT;
    }
    # if a run actually begins at the beginning of the path, this is good; otherwise, there is cost trade-off
    print $out_fh "In check_for_run, looking in the range $beg_index to ".($beg_index+$path_len-1)."\n";
    print $out_fh "   real indices: ".($beg_index % $cycle_len)." to ".(($beg_index+$path_len-1)%$cycle_len)."\n";
    my $max_net_val=$runs_vals_hr->{$beg_index % $cycle_len};    
    my $max_val_index = $beg_index;
    print $out_fh "Max val initialized to $max_net_val\n";
    # in particular, for now, if there is only room left for two normal constructs, don't bother checking for a run
    # may want to change!!
    if ($max_poss_index <= ( $beg_index+ 2*$path_len )) {
	return ($max_val_index, $max_net_val);
    }
    foreach (my $index=$beg_index; $index < ($beg_index+$path_len); $index++) {
	# this may not be the right cost value to use!!
	my $cost = $path_len - ($index-$beg_index);  
	my $net_val = ($runs_vals_hr->{$index % $cycle_len}) - $cost;
	print $out_fh "Current index: $index, val: ".($runs_vals_hr->{$index%$cycle_len}).", cost: $cost, net value: $net_val\n";
	if ($net_val > $max_net_val) {
	    $max_net_val = $net_val;
	    $max_val_index = $index;
	    print $out_fh "Max val updated to $max_net_val, for run starting at real index ".($index % $cycle_len)."\n";
	}
    }
    return ($max_val_index, $max_net_val);
}
    
    
# Break the given cycle into a set of paths (constructs). If flanking
# sequences are being taken into consideration, then try to do it in a
# way that optimizes the frequency that the correct bases appear on
# each end of the constructs; this may require cutting out small
# pieces of the cycle, which do not end up in constructs at this point
# and are returned.
#
# Returns a ref to  array of references to constructs, where each
#     construct is an array of kmers that follow a de-bruijn-like
#     pattern.
sub break_cycle_into_constructs($$$;$$$$$) {
    # assumes single construct length
    my ($k, # integer
	$construct_len, # integer, target construct length;
	$cycle_ar, # ref to an array containing kmers that form a big
		   # cycle;
	$db_str_l, # (optional) string, assumed to be a two-base
		   # sequence of DNA (the left flanking seq) -- if
		   # non-zero/non-null, an attempt is made to break
		   # the cycle into constructs so that as many as
		   # possible begin with a kmer that begins with the
		   # two-base sequence given as the fifth arg (or, at
		   # least, the second base of that sequence);
		   # constructs are made longer by 2 base pairs (or,
		   # at least, 1 base pair), if they start with the
		   # specified sequence (the second base of the
		   # specified sequence); this attempt results in
		   # leftover excised partial paths which may be sewn
		   # together into "extra" constructs (that may
		   # contain repeats of previously used edges);
	$db_str_r, # (optional) string, assumed to be a two-base
		   # sequence of DNA (the right flanking sequence),
		   # which is treated in a symmetric way to the fourth
		   # arg;
	$num_trials,  # number of trials to try to create fewer extra
		      # paths;
	$out_fh, # file handle for printing output
	$num_oligos_short # max number of oligos that can be of the
			  # given construct length (the rest must be
			  # 1bp longer)
	)=@_;
    if (!$out_fh) {
	$out_fh = \*STDOUT;
    }
    print $out_fh "\nBreaking cycle into constructs...\n";
    my $path_len = $construct_len - ($k-1);
    my @best_constructs;
    my @best_extra_bits;
    my $best_total_ex_edges;
    my $check_run=0;  # 0 if regular construct, 1 if at beginning of a run, 2 if should check for a run
    my $all_reg=0;
    if (!$db_str_l and !$db_str_r) {
	$all_reg=1;
    }
    
    if (!$num_trials) {
	$num_trials=1;
    }

    my ($start_index, $runs_vals_hr, $runs_paths_hr, $max_run_val, @best_keys, @good_keys);
    if ($db_str_l) {
	print $out_fh "In break_cycle_into_constructs:  building runs table.\n";
	($runs_vals_hr, $runs_paths_hr, $max_run_val) = 
	    build_runs_table($cycle_ar, $path_len, $db_str_l, $db_str_r, undef, $out_fh);
	@best_keys = grep { $runs_vals_hr->{$_} == $max_run_val } keys(%{$runs_vals_hr});	
	if (!scalar(@best_keys) or ($max_run_val != $runs_vals_hr->{$best_keys[0]})) {
	    die "\nIn break_cycle_into_constructs:  Error getting best keys for runs table.\n";
	}
	@good_keys =grep { $runs_vals_hr->{$_} >= ( (($max_run_val-5) > 0) ? ($max_run_val-5) : 1) } keys(%{$runs_vals_hr});
	if (!scalar(@good_keys)) {
	    die "\nIn break_cycle_into_constructs:  No positive runs values in table?\n";
	}
	print $out_fh "\nKeys with max value $max_run_val in runs table: ".join(' ', @best_keys)."\n";
	print $out_fh "Keys with value at least max(".($max_run_val-5).", 1) in runs table: ".join(' ', @good_keys)."\n";	

	}

    for (my $trial=0; $trial < $num_trials; $trial++) {
	my @constructs;
	if ($db_str_l) {
	    if (!$trial) {
		print $out_fh "\nChoosing start index from best keys:  trial $trial.\n";
		$start_index = choose_random_elt(\@best_keys);	
		remove($start_index, \@good_keys);
	    } else {
		if (!scalar(@good_keys)) {
		    last;
		}
		print $out_fh "\nChoosing start index from good keys:  trial $trial.\n";
		$start_index = choose_random_elt(\@good_keys);	
		remove($start_index, \@good_keys);
	    }
	    $check_run=1;
	} else {
	    print $out_fh "\nChoosing start edge randomly.\n";
	    $start_index = int(rand(scalar(@{$cycle_ar})));
	}
	my $cycle_len = scalar(@{$cycle_ar});
	my $max_poss_index = ($start_index + $cycle_len) - 1;
	
	print $out_fh "Breaking cycle into constructs, starting at index: $start_index,  start edge: ".
	    $cycle_ar->[$start_index].",\n".
	    "cycle len: $cycle_len, max_poss_index: $max_poss_index, which corresponds to real index ".
	    ($max_poss_index % $cycle_len).".\n";
	
	my $cur_index=$start_index;
	my @extra_bits;
	
	my $num_specialc=0;
	my $end_flag=0;
	
	while( $cur_index <= $max_poss_index ) {
	    if (defined($num_oligos_short) and ($num_oligos_short == scalar(@constructs))) {
		print $out_fh "Constructed all oligomers of length $construct_len.\n";
		$construct_len++;
		$path_len++;
		print $out_fh "Target construct length now increased to $construct_len.\n";
	    }
	    my ($sub_path_ar, $rem_path_ar);
	    print $out_fh "\nCurrent index: $cur_index, corresponds to real index: ".($cur_index % $cycle_len)."\n";
	    if ($check_run == 2) {
		# check for a run with a good benefit/cost trade-off here
		my ($max_val_index, $max_net_val) =
		    check_for_run( $runs_vals_hr, $cur_index, $path_len, $cycle_len, $max_poss_index, $out_fh);   
		if (!$max_net_val) {
		    print $out_fh "No good restart found from index $cur_index, real index ".($cur_index%$cycle_len).".\n";
		    $check_run=0;
		} else {
		    $check_run=1;		    
		    # put piece of cycle onto partial paths list
		    if ($max_val_index != $cur_index) {		    
			($sub_path_ar, $end_flag) = get_conn_path($cycle_ar, $cur_index % $cycle_len, 
								  ($max_val_index-1) % $cycle_len, $start_index, 
								  $out_fh);		    
			push(@extra_bits, $sub_path_ar);
			print $out_fh "Restart found from index $cur_index, real index ".($cur_index%$cycle_len).".\n".
			    "stored partial path up to ".($max_val_index-1)% $cycle_len.", at index ".
			    $#extra_bits.".\n";
			if ($end_flag) {
			    die "In break_cycle_into_constructs: Unexpectedly ended loop through cycle!\n";
			}
			$cur_index = $max_val_index;
		    }
		}
	    }	
	    my $beg_index = $cur_index % $cycle_len;
	    my $end_index;
	    if ($check_run == 0) {
		print $out_fh "Trying to make a normal construct.\n";
		# make a normal construct, if there's enough of the cycle left
		$end_index = ($cur_index+$path_len-1) % $cycle_len;	    
		($sub_path_ar, $end_flag) = get_conn_path ($cycle_ar, $beg_index, $end_index, $start_index, $out_fh);
		if ($end_flag) {
		    # short path
		    if (scalar(@{$sub_path_ar})) {
			$end_index = $cur_index + scalar(@{$sub_path_ar})-1;
			push(@extra_bits, $sub_path_ar);
			print $out_fh "Stored a partial path instead, up to ".$end_index.", real index ".
			    ($end_index % $cycle_len).", at index ".$#extra_bits.".\n";
			if(!($db_str_l or $db_str_r)) {
			    die "\nCycle does not break evenly into constructs of length $construct_len.  Last construct is short!\n";
			}
		    } 
		    last;
		} else {
		    print $out_fh "Stored normal construct through index $end_index.\n";
		    # normal construct
		    push(@constructs, $sub_path_ar);
		    $cur_index += $path_len;
		    if (!$all_reg) {
			$check_run=2;
		    }
		    next;
		}
	    }
	    if ($check_run == 1) {
		# at the beginning of a run
		# $runs_vals_hr->{$cur_index} should be set and positive
		print $out_fh "\nStarting a run at index $cur_index, real index $beg_index.\n";
		my @cur_path_indices = @{$runs_paths_hr->{$beg_index}};			
		print $out_fh "Run path in table: ".join(' ', @cur_path_indices)."\n";
		$end_index = $cur_path_indices[$#cur_path_indices];
		while (scalar(@cur_path_indices)) {
		    my $next_index = shift(@cur_path_indices);
		    ($sub_path_ar, $end_flag) = get_conn_path ($cycle_ar, $beg_index, $next_index-1, $start_index, 
							       $out_fh);
		    if (!$end_flag) {
			push(@constructs, $sub_path_ar);
			$num_specialc++;
			$cur_index += (($next_index-$beg_index) % $cycle_len);
			$beg_index = $next_index;
			print $out_fh "Stored a special construct, through real index ". ($next_index-1).".\n";
			# print $out_fh "Updated cur_index to $cur_index, beg_index to $beg_index.\n";
			if ($cur_index > $max_poss_index) {
			    $end_flag=1;
			    last;
			}
		    } else {
			# last special construct ended early, may need to split it up
			if (scalar(@{$sub_path_ar})) {
			    print $out_fh "Short flag was set.  Checking returned path: ".join(' ', @{$sub_path_ar})."\n";
			    $end_index = $cur_index + scalar(@{$sub_path_ar})-1;
			    ($sub_path_ar, $rem_path_ar)= 
				split_final_segment ($sub_path_ar, $path_len, $db_str_l, $db_str_r, $out_fh);
			    if ($sub_path_ar) {
				print $out_fh "Path returned by split_final_segment: ".join(' ', @{$sub_path_ar})."\n";
			    } 
			    if ($rem_path_ar) {
				print $out_fh "Partial path returned by split_final_segment: ".join(' ', @{$rem_path_ar})."\n";
			    }
			    if (!$sub_path_ar and !$rem_path_ar) {
				print $out_fh "Nothing returned by split_final_segment!\n";
			    }
			    my $sub_end_index;
			    if ($sub_path_ar and scalar(@{$sub_path_ar})) {	
				if ($rem_path_ar and scalar(@{$rem_path_ar})) { 
				    $sub_end_index = $cur_index + scalar(@{$sub_path_ar})-1;
				} else {
				    $sub_end_index = $end_index;
				}
				print $out_fh "Set (non-real) end index (for path returned by split_final_segment) to $sub_end_index.\n"
			    }
			    if ( $sub_path_ar and scalar(@{$sub_path_ar}) ) {
				push(@constructs, $sub_path_ar);
				if (scalar(@{$sub_path_ar}) == $path_len) {
				    print $out_fh "Stored a normal (final) path instead, up to ".$sub_end_index.", real index ".
					($sub_end_index % $cycle_len)."\n";
				} elsif (scalar(@{$sub_path_ar}) > $path_len) {
				    $num_specialc++;
				    print $out_fh "Stored a (final) special path, through index "
					.$sub_end_index.", real index ".($sub_end_index % $cycle_len)."\n";
				} else {
				    die "Path returned by split_final_segment is unexpected length: ".
					scalar(@{$sub_path_ar})."!\n";
				}
			    }
			    if ( $rem_path_ar and scalar(@{$rem_path_ar})) {
				push(@extra_bits, $rem_path_ar);
				print $out_fh "Stored a partial (final) path instead,".
				    ( defined($sub_end_index) ? " from ".($sub_end_index+1) : '').
				    " up to ".$end_index.", real index ".($end_index % $cycle_len).
				    ", at index ". $#extra_bits.".\n";	
			    }
			}
			last;
		    }
		}
		$check_run=2;
		if ($end_flag) {
		    last;
		}
	    }
	}
	
	print $out_fh "\nFor trial $trial:  Constructs built: ".scalar(@constructs)." total, ".
	    "$num_specialc special.\n";
	my $total_ex_edges=0;
	if (scalar(@extra_bits)) {
	    for (my $i=0; $i < scalar(@extra_bits); $i++) {
		my $ex_path_ar = $extra_bits[$i];
		$total_ex_edges += scalar(@{$ex_path_ar});
	    }
	    print $out_fh "There are a total of $total_ex_edges edges contained in ".scalar(@extra_bits)." extra paths.\n";
	}
	if (!@best_constructs or (scalar(@extra_bits) < scalar (@best_extra_bits))
	    or ( (scalar(@extra_bits) == scalar(@best_extra_bits))
		 and ($total_ex_edges < $best_total_ex_edges) ) ) {
	    print $out_fh "Found improvement for breaking up the cycle.\n";
	    @best_constructs = @constructs;
	    @best_extra_bits = @extra_bits;
	    $best_total_ex_edges = $total_ex_edges;
	} else {
	    print $out_fh "Did not find improvement this trial for breaking up the cycle.\n";
	}
    }
    if (scalar(@best_extra_bits)) {
	print $out_fh "\nThere are ".scalar(@best_extra_bits). " extra partial paths that have not been made into constructs.\n";
	print $out_fh "Printing remaining partial paths: \n";
	for (my $i=0; $i < scalar(@best_extra_bits); $i++) {
	    my $ex_path_ar = $best_extra_bits[$i];
	    print $out_fh "Partial path index $i:\t".join(' ', @{$ex_path_ar})."\n";
	}
	print $out_fh "There are a total of $best_total_ex_edges edges contained in these extra paths.\n";
	return (\@best_constructs, \@best_extra_bits, $best_total_ex_edges);
    } else {
	return (\@best_constructs);
    }
}

# Sews the extra pieces of paths into constructs as compactly as
# possible, may add kmers if necessary to put different partial paths
# onto the same construct.
#
# Returns an array of references to two arrays; the first contains the constructs, and the second contains
#        the list of edges added during this construction process.
sub build_extra_constructs($$$$; $$$$$) {
    my ($k, # integer
	$extra_paths_ar, # array ref to an array of references to
			 # partial paths, i.e., de-bruijn-like
			 # sequences of kmers stored in arrays;
	$construct_len, # an integer, the desired length of the
			# construct;
	$kmers_ar, # ref to an array of edges that can be used to pad the paths;
	$db_str_l, # (optional) string, a two-base sequence of DNA (left flanking sequence);
	$db_str_r, # (optional) string, a two-base sequence of DNA
		   # (right flanking sequence); -- see the discussion
		   # of the use of the last two args in the comments
		   # for sub break_cycle_into_constructs (they are
		   # ignored if zero/undefined);
	$num_trials, # (optional) number of trials (if greater than 1)
		     # to try to build the smallest number of paths;
	$out_fh, # (optional) file handle for printed output 
	$num_oligos_short  # maximum number of oligomers that can be
			   # of the given construct length (the rest
			   # must be 1bp longer).
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
    my $default_path_len = $construct_len - ($k-1);
    my @sorted_extra_paths = sort { scalar(@{$a}) <=> scalar(@{$b}) } @{$extra_paths_ar};
    my $biggest_path_ar = $sorted_extra_paths[$#sorted_extra_paths];
    if ( (!defined($num_oligos_short) and (scalar(@{$biggest_path_ar}) >= $default_path_len))
	 or (defined($num_oligos_short) and (scalar(@{$biggest_path_ar}) >= ($default_path_len+1)))) {
	die "Unexpectedly long (>= $default_path_len) paths in extra paths array: at least one of length ".
	    scalar(@{$biggest_path_ar})."!\n";
    }
	
    my @best_added_edges;
    my @best_constructs;
    if (!$num_trials) {
	$num_trials=1;
    }    
    for (my $trial=0; $trial < $num_trials; $trial++) {
	print $out_fh "\nBuilding constructs out of partial paths: trial $trial.\n\n";
	my @extra_paths = @sorted_extra_paths;
	my @constructs;
	my @added_edges;
	my @cur_path;	
	while( scalar(@extra_paths) ) {	    
	    if (defined($num_oligos_short) and ($num_oligos_short == scalar(@constructs))) {
		print $out_fh "Constructed all oligomers of length $construct_len.\n";
		$construct_len++;
		$default_path_len++;
		print $out_fh "Target construct length now increased to $construct_len.\n";
	    }
	    print $out_fh ''. scalar(@extra_paths) ." extra paths remain.\n";
	    my $cur_path_ar;
	    if (!$trial) {
		# on first trial, just do them in increasing order
		print $out_fh "\nTaking smallest extra path remaining (first trial).\n";
		$cur_path_ar = shift(@extra_paths);
	    } else {
		# on future trials, try randomly 
		print $out_fh "\nChoosing from remaining extra paths randomly.\n";
		$cur_path_ar = choose_random_elt(\@extra_paths);
		remove($cur_path_ar, \@extra_paths);
	    }
	    @cur_path = @{$cur_path_ar};
	    print $out_fh "Current path: ".join(' ', @cur_path)."\n";
	    my ($construct_ar, $add_path_ar) = 
		get_best_extra_path_construct($k, \@cur_path, \@extra_paths, $default_path_len, $kmers_ar, 
					      $db_str_l, $db_str_r, $out_fh);
	    if (!$trial ) {
		# re-sort since the extra paths may have changed
		@extra_paths = sort { scalar(@{$a}) <=> scalar(@{$b}) } @extra_paths;
	    }
	    if (scalar(@{$construct_ar})) {
		push(@constructs, $construct_ar);
		print $out_fh "Created construct at index ".$#constructs.":\n".
		    join(' ', @{$construct_ar})."\n";		
	    }
	    if (scalar(@{$add_path_ar})) {
		print $out_fh "Used ".scalar(@{$add_path_ar})." additional edges: ".join(' ', @{$add_path_ar})."\n";		
		push(@added_edges, @{$add_path_ar});
	    }
	}
	print $out_fh "\n".scalar(@constructs)." constructs built from extra paths;\n".
	    scalar(@added_edges)." edges used that have already been represented.\n";
	if (!scalar(@best_constructs) or (scalar(@constructs) < scalar(@best_constructs))
	    or ( (scalar(@constructs) == scalar(@best_constructs)) 
		 and (scalar(@added_edges) < scalar(@best_added_edges)) ) ){
	    print $out_fh "Found improvement in extra constructs: ".scalar(@constructs)." extra constructs\n";
	    if (scalar(@best_constructs)) {
		print $out_fh "   (down from ".scalar(@best_constructs)." constructs).\n";
		if (scalar(@added_edges) < scalar(@best_added_edges)) {
		    print $out_fh ''.scalar(@added_edges). " added edges (down from ".scalar(@best_added_edges).").\n";
		}
	    }
	    @best_constructs = @constructs;
	    @best_added_edges = @added_edges;
	} else {
	    print $out_fh "No improvement found in number of extra constructs: ".scalar(@constructs).", up from best ".
		scalar(@best_constructs).".\n";
	}
    }
    return (\@best_constructs, \@best_added_edges);    
}

# Build a new constrct by connecting the given path in $path_ar in the
# most efficient greedy way possible with path(s) in $extra_paths_ar
# (if possible). Assumes that $extra_paths_ar does not point to an
# empty array.
sub get_best_extra_path_construct ($$$$$; $$$) {
    my ($k, $path_ar, $extra_paths_ar, $default_path_len, $kmers_ar, 
	$db_str_l, $db_str_r, $out_fh) = @_;
    my ($sb_l, $fb_r);
    if (!$out_fh) {
	$out_fh = \*STDOUT;
    }
    if ($db_str_l) {
	$sb_l = substr($db_str_l, 1, 1);
    }
    if ($db_str_r) {
	$fb_r = substr($db_str_r, 0, 1);
    }
    my @cur_path=@{$path_ar};
    print $out_fh "In get_best_extra_path_construct: cur_path is ".join(' ', @cur_path)."\n";
    if (scalar(@cur_path) > $default_path_len) {
	die "In get_best_extra_path_construct:  cur_path has unexpected length ".scalar(@cur_path)."!\n";
    }
    my (@construct, @add_path);  # return vars
    my ($construct_ar, $add_path_ar);
    
    if (!$extra_paths_ar or !scalar(@{$extra_paths_ar})) {
	print $out_fh "No paths in extra_path array.\n";
	($construct_ar, $add_path_ar) = randomly_finish_path(\@cur_path, $kmers_ar, $default_path_len, undef,
							     $out_fh);
	return ($construct_ar, $add_path_ar);
    }

    my @target_path_len = ($default_path_len, $default_path_len);
    if (is_prefix($db_str_l,$cur_path[0])) {
	# @cur_path is first
	$target_path_len[0] += 2;
    } elsif (is_prefix($sb_l, $cur_path[0])) {
	$target_path_len[0] += 1;
    }
    if (is_suffix($db_str_r, $cur_path[$#cur_path])) {
	# @cur_path is last
	$target_path_len[1] +=2;
    } elsif (is_suffix($fb_r, $cur_path[$#cur_path])) {
	$target_path_len[1] +=1;
    }
    print $out_fh "Looking for path to put before or after current path on construct.\n";
    print $out_fh "Cur path: ".join(' ', @cur_path)."\n";
    print $out_fh "Target path len (if cur_path starts or ends construct): ".join(' or ', @target_path_len)."\n";
    my $pot_next_path_index=0;	
    
    my @max_len_suffix;
    my @match_path_index;
    my @len_suffix;
    for (my $i=0; $i< scalar(@{$extra_paths_ar}); $i++) {	    
	my @check_path =@{$extra_paths_ar->[$i]};
	$len_suffix[0]=suffix_prefix_match($cur_path[$#cur_path], $check_path[0]);
	$len_suffix[1]=suffix_prefix_match($check_path[$#check_path], $cur_path[0]);	
	print $out_fh "Checked index $i, len_suffix: ".join(', ', @len_suffix).", path: ".join(' ', @check_path)."\n";
	for (my $order=0; $order<2; $order++) {
	    if (!@max_len_suffix or !$max_len_suffix[$order] or ($len_suffix[$order] > $max_len_suffix[$order])) {
		$max_len_suffix[$order] = $len_suffix[$order];
		$match_path_index[$order] = $i;
		print $out_fh "Updating max_len_suffix for order $order to ".$max_len_suffix[$order].".\n";
	    }
	}
    }
    my @full_len;   
    my @match_paths;
    print $out_fh "Indices in match_path_index: ".join(', ', @match_path_index)."\n";
    $match_paths[0]=$extra_paths_ar->[$match_path_index[0]];
    $match_paths[1]=$extra_paths_ar->[$match_path_index[1]];    
    print $out_fh "The corresponding paths:\n".
	join(' ', @{$match_paths[0]})."\n".
	join(' ', @{$match_paths[1]})."\n";
    my $last_edge_index_0 = scalar(@{$match_paths[0]})-1;
    
    # too complicated (in case actual path len is too long) to allow extra len at both ends at the same time
    if ( is_suffix($db_str_r, $match_paths[0]->[$last_edge_index_0]) ) {
	$target_path_len[0] = $default_path_len+2;
    } elsif ( is_suffix($fb_r, $match_paths[0]->[$last_edge_index_0])  
	      and ($target_path_len[0] != ($default_path_len+2)) ) {
	$target_path_len[0] = $default_path_len+1;
    }
    if ( is_prefix($db_str_l, $match_paths[1]->[0]) ) {
	$target_path_len[1] = $default_path_len + 2;
    } elsif ( is_prefix($sb_l, $match_paths[1]->[0]) 
	      and ($target_path_len[1] != ($default_path_len+2)) ) {
	$target_path_len[1] = $default_path_len+1;
    }

    print $out_fh "Target path lens adjusted (if cur_path starts or ends construct): ".join(' or ', @target_path_len)."\n";
			  
    $full_len[0] = scalar(@cur_path) + scalar(@{$extra_paths_ar->[$match_path_index[0]]}) + ($k-$max_len_suffix[0]-1);
    $full_len[1] = scalar(@cur_path) + scalar(@{$extra_paths_ar->[$match_path_index[1]]}) + ($k-$max_len_suffix[1]-1);

    print $out_fh "Full path lengths computed: ".join(', ', @full_len)."\n";
    my $best_order;
    if ($full_len[0] <= $target_path_len[0]) {
	$best_order=0;	    
    } elsif ($full_len[1] <= $target_path_len[1]) {
	$best_order=1;
    } else {
	$best_order = ( ($max_len_suffix[0] > $max_len_suffix[1]) ? 0 : 1);
    }

    my $best_index=$match_path_index[$best_order];
    my $add_len = $k-$max_len_suffix[$best_order]-1;
    print $out_fh "Best order $best_order chosen is with cur_path ".(!$best_order ? "starting" : "ending").".\n".
	"best index: $best_index, with add_len: $add_len.\n";
    if ( ( !$best_order
	   and ( (scalar(@cur_path) + $add_len) >= $target_path_len[$best_order] ) )
	 or
	 ( $best_order
	   and ( (scalar(@{$extra_paths_ar->[$best_index]}) + $add_len) >= $target_path_len[$best_order]) ) ) {
	# should only put the cur_path on a construct alone and not bother with a match
	print $out_fh "Can only fit current path on a construct alone -- randomly finishing construct.\n";
	($construct_ar, $add_path_ar) = randomly_finish_path(\@cur_path, $kmers_ar, $default_path_len, undef, 
							     $out_fh);
	return ($construct_ar, $add_path_ar);
    } 

    my @best_match_path=@{$extra_paths_ar->[$best_index]};
    remove($extra_paths_ar->[$best_index], $extra_paths_ar);       
    my @longer_path;
    if (!$best_order) {
	@add_path = get_inner_path($k,$cur_path[$#cur_path], $best_match_path[0]);
	@longer_path=(@cur_path, @add_path, @best_match_path);
    } else {
	@add_path = get_inner_path($k,$best_match_path[$#best_match_path], $cur_path[0]);
	@longer_path=(@best_match_path, @add_path, @cur_path);
    }   
    print $out_fh "Created longer path: ".join(' ', @longer_path)."\n";
    if (scalar(@longer_path) > $target_path_len[$best_order]) {
	my @rem_path;
	if ( ($target_path_len[$best_order] > $default_path_len) 
	     and ( is_suffix($db_str_r, $longer_path[$#longer_path]) 
		   or ( !is_prefix($db_str_l, $longer_path[0]) and is_suffix($fb_r, $longer_path[$#longer_path]) ) ) ) {
	    my $init_index = scalar(@longer_path)- $target_path_len[$best_order];
	    @construct = @longer_path[$init_index..$#longer_path];
	    @rem_path=@longer_path[0..($init_index-1)];
	} else {
	    my $last_index = $target_path_len[$best_order]-1;
	    @construct = @longer_path[0..$last_index];
	    @rem_path = @longer_path[($last_index+1)..$#longer_path];
	}
	print $out_fh "Breaking up longer path and putting remaining path back on extra_paths list.\n".
	    "remaining path: ".join(' ', @rem_path)."\n";
	push(@{$extra_paths_ar}, \@rem_path);
    } elsif (scalar(@longer_path) < $default_path_len) {
	print $out_fh "Longer path is shorter than default len, putting back on extra_paths list.\n";
	push(@{$extra_paths_ar}, \@longer_path);	    
    } elsif ((scalar(@longer_path)==$default_path_len) 
	     or (scalar(@longer_path) == $target_path_len[$best_order])) {
	@construct = @longer_path;
    } else {
	# may need to pad if we are in between acceptable lengths
	$add_path_ar=0;
	if (is_prefix($db_str_l, $longer_path[0])) {
	    ($construct_ar, $add_path_ar) = randomly_finish_path(\@longer_path, $kmers_ar, $target_path_len[$best_order], undef, $out_fh);
	} elsif (is_suffix($db_str_r, $longer_path[$#longer_path])) {
	    # pad at the front
	    ($construct_ar, $add_path_ar) = randomly_finish_path(\@longer_path, $kmers_ar, $target_path_len[$best_order], 1, $out_fh);
	} else {
	    $construct_ar = \@construct;
	}
	if ($add_path_ar and scalar(@{$add_path_ar})) {
	    push(@add_path, @{$add_path_ar});
	}
	return ($construct_ar, \@add_path);
    }
    return (\@construct, \@add_path);
}

# Randomly finish a given path so it has the right length.
#
# Adds edges onto the path referred to by first arg;
# returns an array or two refs, one to the array containing the new whole path
# and the other to an array of the edges added.
# optional fourth arg is a flag; if set, path is padded at the front instead of the end.
sub randomly_finish_path($$$;$$) {
    my ($path_ar, $kmers_ar, $target_path_len, 
	$dir, $out_fh) = @_;
    if (!$out_fh) {
	$out_fh = \*STDOUT;
    }
    my @path = @{$path_ar};
    if (!scalar(@path)) {
	die "In randomly_finish_path: Expecting a path of at least length 1.\n"
    }
    print $out_fh "Randomly finishing path: ".join(' ', @path)."\n";
    my @add_edges;
    if ($dir) {
	print $out_fh "Padding at the front of the path.\n";
    }
    for (my $i=scalar(@path); $i< $target_path_len; $i++) {
	my $edge;
	if (!$dir) {
	    my $prev_edge = $path[$i-1];
	    my $head = get_head($prev_edge);
	    my @poss_next_edges = get_matches($head, 1, $kmers_ar);
	    $edge=choose_random_elt(\@poss_next_edges);
	    $path[$i]=$edge;
	} else {
	    my $prev_edge = $path[0];
	    my $tail = get_tail($prev_edge);
	    my @poss_next_edges = get_matches($tail, 2, $kmers_ar);
	    $edge = choose_random_elt(\@poss_next_edges);
	    unshift(@path, $edge);
	}
	print $out_fh "Adding $edge to ".(!$dir ? "end" : "front")." of path.\n";
	push(@add_edges, $edge);
    }
    return(\@path, \@add_edges);			 
}
    

# Convert the paths representing the constructs into paths where
# the flanking sequences are represented on each end.
sub add_flanking_seqs_to_constructs($$$$$;$$$) {
    my ($k, $construct_len, $constructs_ar, $flank_l, $flank_r, 
	$out_fh, $num_oligos_short, $no_path_len_check) = @_;
    if (!$out_fh) {
	$out_fh = \*STDOUT;
    }
    my $path_len = $construct_len - ($k-1);
    my (@path_l, @path_r);
    my $db_str_l = substr($flank_l, $k-2, 2);
    my $lastb_l = substr($flank_l, $k-1, 1);
    my $db_str_r = substr($flank_r, 0, 2);
    my $firstb_r = substr($flank_r, 0, 1);
    my @constructs_wf;
    my $short_path_len = $path_len;
    foreach my $construct_ar (@{$constructs_ar}) {	
	if (defined($num_oligos_short) and (scalar(@constructs_wf) == $num_oligos_short)) {
	    print $out_fh "Added flanking sequences to all oligomers of length $construct_len.\n";
	    $construct_len++;
	    $path_len++;
	    print $out_fh "Target construct length now increased to $construct_len.\n";
	}

	my @construct_wf = @{$construct_ar};
	my $first_edge = $construct_wf[0];
	my $last_edge = $construct_wf[$#construct_wf];
	my $cur_path_len = scalar(@construct_wf);
	my $temp_path_len = $cur_path_len;
	my @add_path_l;
	my @add_path_r;
	if ( is_prefix($db_str_l, $first_edge) and ($temp_path_len >= ($path_len + 2)) ) {	    
	    @add_path_l = get_inner_path($k, $flank_l, $first_edge, 2);
	    $temp_path_len -= 2;
	} 
	if ( is_suffix($db_str_r, $last_edge) and ($temp_path_len >= ($path_len + 2)) ) {	    
	    @add_path_r = get_inner_path($k, $last_edge, $flank_r, 2);
	    $temp_path_len -= 2;
	} 
	if ( !@add_path_l and is_prefix($lastb_l, $first_edge) and ($temp_path_len >= ($path_len + 1)) ) {	    
	    @add_path_l = get_inner_path($k, $flank_l, $first_edge, 1);
	    $temp_path_len--;
	}
	if ( !@add_path_r and is_suffix($firstb_r, $last_edge) and ($temp_path_len >= ($path_len + 1)) ) {	    
	    @add_path_r = get_inner_path($k, $last_edge, $flank_r, 1);
	    $temp_path_len--;
	}
	if (!@add_path_l) {
	    @add_path_l = get_inner_path($k, $flank_l, $first_edge, 0);
	}
	if (!@add_path_r) {
	    @add_path_r = get_inner_path($k, $last_edge, $flank_r, 0);
	}
	if ((!$no_path_len_check) and ($temp_path_len != $path_len)) {
	    die "In add_flanking_sequences_to_constructs:  Unexpected path length: $cur_path_len! Anticipated path length: $path_len\n";
	}
	unshift(@add_path_l, $flank_l);
	unshift(@construct_wf, @add_path_l);

	push(@add_path_r, $flank_r);
	push(@construct_wf, @add_path_r);
	
	if (!$no_path_len_check and (scalar(@construct_wf) != ($path_len + 2*$k))){
	    die "In add_flanking_sequences_to_constructs:  Unexpected construct length with flanking seqs: \n".
		scalar(@construct_wf)." instead of ".($path_len + 2*$k)."!\n";
	}
	push(@constructs_wf, \@construct_wf);
    }
    return \@constructs_wf;
}

# Check whether the given constructs have the desired form of either
# an MRCC library or an MRCC library modified to avoid some of the
# repetitions of k-mers caused by flanking sequences.
#
# Returns an array containing:  an integer that is 0 if there are no errors or unexpected repeats, 1 otherwise;
#      a ref to a hash table such that the keys are edges and the value for an edge/key is the
#      index of the first construct in which that edge appears;
#      a ref to a hash table such that the keys are repeated edges and the value for an edge/key is the number
#      of times that edge appears in the given contruct design.
sub check_constructs ($$; $$$ $$$ $$$$) {
    my ($k, # integer
	$constructs_ar, # ref to an array of "constructs", i.e. each
			# element in the array is a ref to an array of
			# kmers;
        ## all args below are optional
	$construct_len, # (optional) integer, the length each
			# construct should be; lengths are not checked
			# if this arg is 0.
	$out_fh, # (optional) file handle for writing output; if zero/undefined, set to STDOUT;
	$hist_out_fh, # (optional) file handle for writing the
		      # histogram data for the repeats, if a separate
		      # file for this data is desired; if
		      # zero/undefined, it is set to the fourth arg
		      # value;
	$repeats_flag, # (optional) flag; if non-zero/non-null,
		       # repeated edges are expected in some
		       # constructs and do not force a bad return
		       # value;
	$expected_repeats_ar, # (optional) ref to an array of kmers
			      # that are expected to be seen more than
			      # once in the constructs -- will check;
			      # if 0/undefined and repeats are
			      # expected, then they can be any of the
			      # kmers.
	$start_index_extras, # starting index in the constructs array
			     # of the constructs in which repeats are
			     # expected; if 0 or undefined and repeats
			     # are expected, then they can occur in
			     # any construct.
	$db_str_l, # di-nucleotide sequence (left flanking sequence)
		   # -- constructs that begin with a kmer that begins
		   # with this sequence must have path length that is
		   # longer by 2 than normal; constructs that begin
		   # with a kmer that begins with the second base of
		   # this sequence must have path length that is
		   # longer by 1;
	$db_str_r, # di-nucleotide sequence (right flanking sequence)
		   # -- treated symmetrically relative to the left 
	$not_all_flag, # flag that indicates not all $k-mers are
		       # expected to appear in the constructs;
	$num_oligos_short # maximum number of oligomers that can be of
			  # the given construct length (the rest must
			  # be 1bp longer).
	) = @_;
    if (!$out_fh) {
	$out_fh = \*STDOUT;
    }
    if (!$hist_out_fh) {
	$hist_out_fh = $out_fh;
    }
    print $out_fh "\nChecking validity of construct design...\n";
    my $bad_flag = 0;
    my $path_len;
    my $short_path_len;
    if ($construct_len) {
	$path_len = $construct_len - ($k-1);
	$short_path_len = $path_len;
    }
    my @rem_kmers = generate_all_k_mers($k); # keeps track of which edges have not been seen yet;
                                             # an edge counts as having been seen if its reverse complement is seen

    my ($sb_l, $fb_r);
    if ($db_str_l) {
	$sb_l = substr($db_str_l, 1, 1);
    }
    if ($db_str_r) {
	$fb_r = substr($db_str_r, 0, 1);
    }
    my %repeat_edges;  # keeps track of all constructs in which a repeated edge appears
    my %edges_by_construct;  # keeps track of first construct in which an edge appears

    for (my $cindex=0; $cindex < scalar(@{$constructs_ar}); $cindex++) {
	## assumes the shorter ones are at the beginning
	if (defined($num_oligos_short) and (defined ($path_len)) and ($cindex == $num_oligos_short)) {
	    print $out_fh "Checked all $cindex oligomers of length $construct_len.\n";
	    $construct_len++;
	    $path_len++;
	    print $out_fh "Target construct length now increased to $construct_len.\n";
	}       
	my @path = @{$constructs_ar->[$cindex]};	
	if (defined($path_len)) {
	    my $cur_path_len = scalar(@path);
	    if ($cur_path_len != $path_len) {
		if (!$db_str_l and !$db_str_r) {
		    print $out_fh "In check_constructs:  ".
			"Construct $cindex is the wrong length:  $cur_path_len (instead of $path_len)!\n";	    
		    $bad_flag = 1;
		} else {
		    # prioritize the double base prefix or suffixes since the single bases can be used or ignored
		    # as needed
		    my $temp_path_len = $cur_path_len;
		    if ( ($temp_path_len > ($short_path_len+1)) and (is_prefix($db_str_l, $path[0])) ) {
			$temp_path_len -=2;
		    }
		    if ( ($temp_path_len > ($short_path_len+1)) and (is_suffix($db_str_r, $path[$#path])) ) {
			$temp_path_len-=2;
		    }
		    if ( ($temp_path_len > $short_path_len) and (is_prefix($sb_l, $path[0])) ) {
			$temp_path_len--;			
		    }
		    if ( ($temp_path_len > $short_path_len) and (is_suffix($fb_r, $path[$#path])) ) {
			$temp_path_len--;
		    }
		    if (($temp_path_len != $path_len) and ($temp_path_len != $short_path_len)) {
			$bad_flag=1;
			print $out_fh "Construct $cindex has unexpected path length: $cur_path_len!\n";
		    }	    
		}
	    }
	}
	for (my $eindex=0; $eindex < scalar(@path); $eindex++) {
	    my $edge = $path[$eindex];
	    my $rc_edge = rev_comp($edge);
	    
	    if (defined($edges_by_construct{$edge}) or defined($edges_by_construct{$rc_edge})) {
		my $first_cindex = (defined($edges_by_construct{$edge}) ? 
				    $edges_by_construct{$edge} : $edges_by_construct{$rc_edge});
		# the edge or its rev comp has been seen before
		my $print_err=0;
		if (!$repeats_flag) {
		    print $out_fh "Edges are not expected to be represented more than once!\n";
		    $bad_flag = 1;
		    $print_err=1;
		} elsif (defined($start_index_extras) and ($cindex < $start_index_extras)) {
		    print $out_fh "Repeated edge $edge unexpected in construct with index $cindex ".
			"(less than $start_index_extras)!\n";
		    $bad_flag = 1;
		    $print_err=1;
		} elsif (defined ($expected_repeats_ar) and !(in_list($edge, $expected_repeats_ar))
			 and !(in_list($rc_edge, $expected_repeats_ar))) {
		    print $out_fh "Edge $edge is not in list of expected repeats!\n";
		    $print_err=1;
		    $bad_flag=1;
		}
		if ($print_err) {
		    print $out_fh "Edge $edge, in construct $cindex, index $eindex in path, has already been represented,\n".
			"   first in construct with index $first_cindex.\n";		    
		}
		# save info about in which constructs the repeated edge or its rev comp occurs
		if ( defined($repeat_edges{$edge}) ) {
		    push(@{$repeat_edges{$edge}}, $cindex);
		} elsif ( defined($repeat_edges{$rc_edge}) ) {
		    push(@{$repeat_edges{$rc_edge}}, $cindex);
		} else {
		    if (defined($edges_by_construct{$edge})) {
			$repeat_edges{$edge} = [$edges_by_construct{$edge}];
			push(@{$repeat_edges{$edge}}, $cindex);
		    } else {
			$repeat_edges{$rc_edge} = [$edges_by_construct{$rc_edge}];
			push(@{$repeat_edges{$rc_edge}}, $cindex);
		    }
		}		    
	    } else {
		# first time seeing the edge
		$edges_by_construct{$edge} = $cindex;
		remove($edge, \@rem_kmers);
		if ($edge ne $rc_edge) {
		    remove($rc_edge, \@rem_kmers);
		}
	    }
	    if ($eindex < $#path) {
		if (! is_prefix ( get_head( $edge ), $path[$eindex+1]) ) {
		    print $out_fh "Not de-Bruijn-like sequence in construct with index $cindex, at edge with index $eindex!\n";
		    $bad_flag = 1;
		}
	    }
	} # end loop over edges in path	
    } # end loop over all constructs
    my @edges_appearing = keys(%edges_by_construct);
    print $out_fh ''.scalar(@edges_appearing)." edges represented at least once in constructs.\n";
    if (scalar(@rem_kmers) and !$not_all_flag) {
	print $out_fh "Not all edges (or their reverse complements) are used in the constructs!  ". scalar(@rem_kmers)." remain:\n";
	print_array(\@rem_kmers, $out_fh);
	print $out_fh "\n";
	$bad_flag = 1;
    }
    # if (scalar(@all_fw_edges) != (scalar(@rc_edges) + scalar(@palindromes)) ) {
	# print $out_fh "Edges in constructs don't add up -- not sure why!\n";
    # return ();
    # }
    my @repeated_edges = sort { scalar(@{$repeat_edges{$b}}) <=> scalar(@{$repeat_edges{$a}}) } keys(%repeat_edges);
    my @multip_repeats = grep { scalar(@{$repeat_edges{$_}}) > 2 } @repeated_edges;
    my @single_repeats = grep {scalar(@{$repeat_edges{$_}}) == 2 } @repeated_edges;
    print $out_fh ''.scalar(@repeated_edges)." edges represented repeatedly:  ".
	scalar(@single_repeats)." represented exactly twice; ".scalar(@multip_repeats)." are repeated more times.\n";
    if ( (scalar(@multip_repeats) + scalar(@single_repeats)) != scalar(@repeated_edges) ) {
	print $out_fh "Something in the repeat info doesn't add up.\n";
	$bad_flag = 1;
    }
    print $out_fh "\nPrinting histogram data for repeated $k-mers ($k-mer and number of occurrences).\n\n";
    print $hist_out_fh "kmer\tfrequency\n";
    for my $repeat (@repeated_edges) {	
	print $hist_out_fh "$repeat\t".scalar(@{$repeat_edges{$repeat}})."\n";
    }
    print $out_fh "\nConstruct info for repeated $k-mers ($k-mer and indices of constructs containing it):\n\n";
    for my $repeat (@repeated_edges) {
	print $out_fh "$repeat\t".join(', ', @{$repeat_edges{$repeat}})."\n";
    }
    
    if (!$bad_flag) {
	print $out_fh "\nConstructs look good!\n";
    }    
    return ($bad_flag, \%edges_by_construct, \%repeat_edges);
}

# print the paths in the de Bruijn graph that represent the constructs
sub print_construct_paths($;$) {    
    my ($constructs_ar, $out_fh) = @_;
    if (!$out_fh) {
	$out_fh = \*STDOUT;
    }
    print $out_fh "\nTotal of ". scalar(@{$constructs_ar}). " constructs.\n";
    print $out_fh "Printing the paths of edges represented in the constructs:\n\n";
    for (my $i=0; $i < scalar(@{$constructs_ar}); $i++) {
    	print $out_fh "Construct $i: \t". join(' ', @{$constructs_ar->[$i]}) ."\n";
    }
    print $out_fh "\n";
}

# Print the actual constructs, incorporating (if appropriate) the end
# of the left flanking sequence and the beginning of the right
# sequence at the beginnings and ends, respectively, of constructs.
sub print_constructs($$;$$$$$$) {
    my ($k, $constructs_ar, 
	$out_fh, $construct_len, $db_str_l, $db_str_r, $num_oligos_short, $no_path_len_check) = @_;
    if (!$out_fh) {
	$out_fh = \*STDOUT;
    }
    my $default_path_len;
    if ($construct_len) {
	$default_path_len= $construct_len - ($k-1);
    }
    my ($sb_l, $fb_r);
    if ($db_str_l) {
	$sb_l = substr($db_str_l, 1, 1);
    }
    if ($db_str_r) {
	$fb_r = substr($db_str_r, 0, 1);
    }

    print $out_fh "Total of ". scalar(@{$constructs_ar}). " constructs:\n";
    for (my $i=0; $i < scalar(@{$constructs_ar}); $i++) {
	# assume all short oligos are at the beginning
	if (defined($num_oligos_short) and ($i== $num_oligos_short)) {
	    $construct_len++;
	    $default_path_len++;
	}
	my @path = @{$constructs_ar->[$i]};
    	print $out_fh "Construct $i: \t";
	my $pathlen = scalar(@path);
	my $first_edge = $path[0];
	my $last_edge = $path[$#path];
	if ($default_path_len) {
	    if ( $db_str_l and is_prefix($db_str_l, $first_edge) and ($pathlen >= ($default_path_len + 2)) ) {
		shift(@path);
		shift(@path);
		$pathlen -= 2;
	    } 
	    if ( $db_str_r and is_suffix($db_str_r, $last_edge) and ($pathlen >= ($default_path_len + 2)) ) {
		pop(@path);
		pop(@path);
		$pathlen -=2;
	    }
	    if ($db_str_l and is_prefix($sb_l, $first_edge) and ($pathlen >= ($default_path_len+1)) ) {
		shift(@path);
		$pathlen--;
	    } 	    
	    if ($db_str_r and is_suffix($fb_r, $last_edge) and ($pathlen >= ($default_path_len+1)) ) {
		pop(@path);
		$pathlen--;
	    } 	    
	}  
	my $construct = shift(@path);
	while (my $edge = shift(@path)) {
	    $construct = $construct.chop($edge);
	}
	if (!$no_path_len_check && ($construct_len and (length($construct) != $construct_len))){	    
	    die "In print_constructs:  construct $i has incorrect construct length:  ".
		length($construct)." instead of $construct_len!\n";
	}
	print $out_fh "$construct\n";
    }
    print $out_fh "\n";
}

sub print_array($;$) {
    my ($cycle_ar, $out_fh) = @_;
    if (!$out_fh) {
	$out_fh = \*STDOUT;
    }
    my $num_per_line = 10;
    for (my $i=0; $i < scalar(@{$cycle_ar}); $i++) {
	print $out_fh ' '. $cycle_ar->[$i];
	if (! (($i+1) % $num_per_line) ) {
	    print $out_fh "\n";
	}
    }
    print $out_fh "\n";
}

# Computes the number of oligomers that are either all of the same
# size or all of one of two sizes (1 bp difference) so that constructs
# of those lengths can be designed to be an MRCC library.
#		
# The segments must overlap by $k-1 bases, and the cycle contains no
# reverse complements except for palindromes. Returns an array of 2
# or 3 values; the first value $len_floor is the floor of the optimal,
# possibly fractional length of the oligo; the second value is the
# number of pieces that should have this smaller length; if $len_floor
# is not integral, then a third value is returned which is the number
# of oligos which should have length $len_floor+1.
sub compute_lengths ($$;$) {
    my ($k, # integer
	$num_pieces, # integer number of oligos to create from cycle
	$out_fh # file hand for printed output
	)=@_;
    if (!$out_fh) {
	$out_fh = \*STDOUT;
    }
    my $num_kmers = (4**$k + 2**$k)/2;
    my $len = ($num_kmers/$num_pieces) + $k-1;
    my $len_floor = floor($len);
    my $len_ceil = ceil($len);
    my @returnvals;
    if ($len_floor == $len_ceil) {	
	# length is integral
	return ($len_floor, $num_pieces);
    }
    print $out_fh "Shorter length: $len_floor, Longer length: ".($len_floor+1)."\n";
    # length is fractional
    my ($frac_big, $num_big, $num_small);
    $frac_big = $len-$len_floor;
    # print $out_fh "Fraction of oligos of longer length: $frac_big.\n";
    $num_big = floor($frac_big*$num_pieces);
    $num_small = $num_pieces-$num_big;
    print $out_fh "Initial values: Num longer oligos: $num_big, Num shorter oligos: $num_small.\n";
    my $total =  ($num_big*$len_ceil) + ($num_small*$len_floor);
    my $expected_total = $num_kmers + ($k-1)*$num_pieces;
    print $out_fh "Initial check on computation:  Total num bases: $total, Anticipated total: $expected_total\n";
    if (($total < $expected_total)) {
	print $out_fh " -- Rounding problem, updating values...\n";
	$num_big++;
	$num_small--;
	$total =  ($num_big*$len_ceil) + ($num_small*$len_floor);
	print $out_fh "Updated: Num longer oligos: $num_big, Num shorter oligos: $num_small.\n";
	print $out_fh "Final check on computation: Total num bases: $total, Anticipated total: $expected_total\n";
	if ($total != $expected_total) {
	    die "Cannot find a good solution to the rounding problem!"
	}
    }
    return ($len_floor, $num_small, $num_big);
}

# Computes a short path between two sequences by shifting as few times
# as needed
sub get_inner_path($$$;$) {
    my ($k,$str1, $str2, $len_suffix) = @_;
    if (! defined($len_suffix)) {
	$len_suffix = suffix_prefix_match($str1, $str2);
    } elsif (($len_suffix > $k) or ($len_suffix < 0)){
	die "Length of suffix to be used in get_inner_path must be at least 0 and at most $k!\n";
    }
    my $rem_len = $k-$len_suffix;
    my @inner_path;
    for (my $i=1; $i<$rem_len; $i++) {
	my $next_str = substr($str1, $i, $k-$i).substr($str2, $len_suffix, $i);
	push(@inner_path, $next_str);
    }
    return @inner_path;
}

1;
