#!/usr/bin/perl -w
# 
# SJ Riesenfeld
# Created: 2009
# Last updated: Feb 2012
#
# Designs an MRCC library of order $k (set below) and size $num_oligos
# (below), i.e., a set of oligomers (chosen to be as even in length as
# possible) such that each sequence of length k (k-mer) is represented
# exactly once in the set, either as itself or as its reverse
# complement.
#
# It is an option (see option -f below) to design oligomers with
# awareness of the flanking sequences that will eventually be adjacent
# to each end of an oligomer in experimental applications. Repetitions
# of k-mers in this case are unavoidable, due to the k-mers that
# bridge junctions with the flanking sequences. So the program designs
# a set of oligomers that reduces these repetitions by attempting to
# build slightly longer constructs that begin (and/or end) with the
# final (and/or initial) bases in the left (and/or right) flanking
# sequences. These bases are then eliminated from the design, without
# losing coverage of any k-mers. However, the output design is
# intended to work with the specific flanking sequences, and it may
# not cover every k-mer when the flanking sequences are not
# present. The output in this case is not a true MRCC library.
####
# 
# To build an example of a true MRCC library for k=6, of size 208
# oligomers, each 15bp in length, run at the command-line:
#
#  perl design_oligomers.pl 
#
# To change the values of parameters such as k, the number of
# oligomers, and the flanking sequences, see below.
#
####

use strict;
use warnings;

use File::Spec;
use Getopt::Std;

use DNA_Manip qw (get_self_rev_comps generate_all_k_mers );
use Seq_Utils qw (remove_sublist);
use DNA_DeBruijnGraph_Utils qw(find_eulerian_pair check_connected check_degrees 
                      print_construct_paths print_constructs check_constructs
                      break_cycle_into_constructs print_array 
                      build_extra_constructs add_flanking_seqs_to_constructs 
                      compute_lengths);                       
use Palindrome_Paths qw(find_insertion_points compute_pconstructs 
                        get_pdists get_short_path_palindromes 
                        print_insertion_points insert_short_ppaths ); 
                        

################### 
#### These parameters are set currently for an application-specific
#### use and can be reset by the user. They could easily be made into
#### input options. Note that they would typically be set the same for
#### all the related programs.
my $k=6;  # length of sequences to be covered by output oligomers
my $num_oligos = 208;   # number of oligomers for a true MRCC library,
			# i.e., if flanking sequences are *not*
			# incorporated in design

# the left (right) flanking sequence is assumed not to end (begin)
# with 2 bases that are the same

# my $flank_l = 'CTCGAG'; # if option -f is not selected: these
my $flank_l='';		  # flanking sequences are ignored to compute an
			  # MRCC library; then, if the seqs are defined 
# my $flank_r = 'AGATCT'; # and not empty, they are concatenated to the
my $flank_r='';		  # oligomers afterward for the purpose of counting
			  # repeated $k-mers; the design validation
			  # will not work once they are concatenated
			  # if the flanking seqs do not have length $k
			  # (also see option -f below)

################### 
my $solnfile = 'oligomer_library_design.txt';
my $histfile = 'design_histogram_info.txt';
my $outfile = 'design_process.log';

my $usage = qq{
The $0 program designs a library of oligomers, all of roughly equal
   length, that includes exactly one copy of a single representative
   for each k-mer and reverse complement.

Usage:  $0 <options (see below)>

   -f [no value taken; if set, 2 bases of the left flanking sequence
       are incorporated into the design; otherwise, the flanking
       sequence is not used in the design]
   -d <name for output directory, to be made and concatenated onto the
       beginning of the names of any output files>
   -s <name of solution file, the file that contains the actual
       construct design; if 0/undefined, the argument for the '-o'
       option is used>
   -a <name of file to contain the histogram data for the repeated
       kmers in the construct design (analyzed with the flanking
       sequences on both sides); if 0/undefined, the argument for the
       '-o' option is used>
   -o <name of file to contain other output of the program, including
       the steps of the algorithm, the cycoe built from the graph, the
       design validation, etc; if 0/undefined, STDOUT is used (not
       recommended)>
   -k <integer value of k>
   -n <integer target number of oligomers (enforced if 'f' is not set;
       approximate otherwise)>
   -l <string: left flanking sequence>
   -r <string: right flanking sequence>
   -h [no value taken; this message is printed]

  To build an example of an MRCC library for k=6, of 208 oligomers,
  each 15bp in length, run:

     perl design_oligomers.pl

  For more information, see the program file.

};

my %opts;
if ( !getopts('fd:s:a:o:k:n:l:r:h', \%opts) or $opts{'h'} ) {
    die "$usage\n";
}
if ($opts{'k'}) { $k = $opts{'k'};}
if ($opts{'n'}) { $num_oligos = $opts{'n'};}
if ($opts{'a'}) { $histfile = $opts{'a'};}
if ($opts{'s'}) { $solnfile = $opts{'s'};}
if ($opts{'o'}) { $outfile = $opts{'o'};}
######## These variables affect now many times heuristic, randomized
######## processes are attempted (related to flanking-sequence
######## incorporation).
my $num_trials_longppaths = 6;
my $num_trials_pconstructs = 3;
my $num_trials_cycle = 5;
my $num_trials_break_cycle = 3;
# can increase this
my $num_trials_build_extras = 3;
#####
my ($flank_seq_flag, $outdir);
$flank_seq_flag = $opts{'f'};
if ($opts{'l'}) { $flank_l = $opts{'l'};}
if ($opts{'r'}) { $flank_r = $opts{'r'};}
my $curdir = File::Spec->curdir();
my $k_even;
if ( ($k %2) == 0) {
    $k_even=1;
} 
my $n_palindromes = (($k_even) ? (2**($k-1)) : 0);
if (($k_even) && ($num_oligos <= $n_palindromes)) {
    die ("Option -n requires a value larger than ".$n_palindromes.", the number of palindromes!")
}
if ($outdir=$opts{'d'}) {
    unless ( (-d $outdir) and (-w $outdir)) {
	mkdir($outdir) or die "Cannot find, write, or make directory $outdir.\n";
    }
}
if (!$outdir) { $outdir = $curdir; }

# Last two DNA bases of left flanking sequence and first two DNA bases
# of right flanking sequence
my ($db_str_l, $db_str_r);
if ($flank_seq_flag) {
    if (!($flank_l)) {
	die "Left flanking sequence is not defined or is empty";
    } else {
	$db_str_l = substr($flank_l, ($k-2), 2);
    }
    if (!($flank_r)) {
	die "Right flanking sequence is not defined or is empty";
    } else {
	$db_str_r = substr($flank_r, 0, 2);
    }
}
my ($out_fh, $soln_fh, $hist_out_fh);
$outfile = File::Spec->catfile($outdir, $outfile);
open ($out_fh, ">$outfile") or die "Cannot open $outfile for writing: $!\n";
print "\nInformation about this design process is being written in file $outfile.\n";

# if the target number of oligos does not result in an integral oligo length, 
# two lengths are used; the longer length is 1+shorter length.
print $out_fh "Computing the correct oligomer length(s).\n";
# desired length, num of oligomers (from $num_oligos)
my ($oligo_len_short, $num_oligos_short, $num_oligos_long) = 
    compute_lengths($k, $num_oligos, $out_fh); 
my $num_oligos_cur_limit = $num_oligos_short;
my $oligo_len = $oligo_len_short;

print $out_fh "\nThis program is designing constructs that contain each $k-mer or its reverse complement.\n";
if ($flank_seq_flag) {
    print $out_fh "Incorporating the flanking sequence into the design.  Some $k-mers may be represented repeatedly.\n";
}
my @kmers=generate_all_k_mers($k);   # all the edges in the original graph
print $out_fh "\nBuilding a de Bruijn graph that contains an edge for each $k-mer.  There are ". @kmers." $k-mers.\n";
my @s_rev_comps = get_self_rev_comps(\@kmers);
print $out_fh "   There are ". @s_rev_comps." self-reverse-complementary $k-mers (a.k.a. palindromes).\n";
my $total= (scalar(@kmers)+scalar(@s_rev_comps))/2;  # total num. edges wanted (#edges in MRCC library)
print $out_fh "   There are $total $k-mers, excluding reverse complements (except palindromes, which are included).\n"; 
my @km1mers = generate_all_k_mers($k-1);

# first build paths between palindromes and remove these edges from the graph
# get all shortest distances between palindromes
my ($ppdists_hr, $bestp_hr, $bestd_hr) = get_pdists($k, \@s_rev_comps);
# get matching of palindromes that have path length 2 (including endpoints);
# these can easily be inserted in the big cycle later 
my ($short_pmatching_hr, $short_ppaths_hr) = get_short_path_palindromes($ppdists_hr);
my @short_p_palindromes = %{$short_pmatching_hr};
my @short_ppkeys = keys(%{$short_pmatching_hr});
remove_sublist(\@short_p_palindromes, \@s_rev_comps);
remove_sublist(\@short_p_palindromes, \@kmers);
if ($k_even) {
    print $out_fh "\nThere are ".scalar(@short_p_palindromes)." palindromes\n".
	"   (".join(' ',@short_p_palindromes).")\n".
	"   that can be paired so that pairs have a connecting path of length 2\n".
	"   (including start and end edges).\n".
	"   These paths will be inserted into the big forward cycle computed later.\n";
}
# build constructs containing remaining palindromes:
# find longer paths that take up an entire construct between pairs of remaining palindromes.
if ($k_even) {
    print $out_fh "\nBuilding the constructs that contain the remaining palindromes...\n";
}
my ($pconstructs_ar, $pconstr_fw_edges_ar, $pconstr_rc_edges_ar, $pconstr_palindromes_ar);
if ($flank_seq_flag) {
    ($pconstructs_ar, $pconstr_fw_edges_ar, $pconstr_rc_edges_ar, $pconstr_palindromes_ar)= 
	compute_pconstructs ($k, $oligo_len, $num_trials_pconstructs, $num_trials_longppaths, 
			     \@s_rev_comps, \@kmers, 0, $db_str_l, $db_str_r, $out_fh,
			     $num_oligos_cur_limit);
} else {
    ($pconstructs_ar, $pconstr_fw_edges_ar, $pconstr_rc_edges_ar, $pconstr_palindromes_ar)= 
	compute_pconstructs ($k, $oligo_len, $num_trials_pconstructs, $num_trials_longppaths, 
			     \@s_rev_comps, \@kmers, undef, undef, undef, $out_fh,
			     $num_oligos_cur_limit);
}
if (($oligo_len == $oligo_len_short) and  (scalar(@{$pconstructs_ar}) >= $num_oligos_cur_limit)) {
    $num_oligos_cur_limit = $num_oligos_long;
    $oligo_len++;
} else {
    $num_oligos_cur_limit -= scalar(@{$pconstructs_ar});
}
if ($k_even) {
    print $out_fh ''. scalar(@$pconstructs_ar). " constructs built, containing ".
	scalar(@{$pconstr_palindromes_ar})." palindromes.\n";
    print $out_fh "These constructs also contain ".scalar(@$pconstr_fw_edges_ar)." non-palindromic edges.\n";
    if ($out_fh != \*STDOUT) {
	print "\n". scalar(@$pconstructs_ar). " oligomers designed that correspond to paths between pairs among the ".
	    scalar(@{$pconstr_palindromes_ar})." isolated palindrome edges.\n";
    }
    print_construct_paths($pconstructs_ar, $out_fh);
}
if ( scalar(@$pconstr_fw_edges_ar) != scalar(@$pconstr_rc_edges_ar) ) {
    die "Expecting same number of (non-palindromic) forward and reverse-complement edges to be used in building constructs!\n";
}
my $expected_num;
remove_sublist($pconstr_fw_edges_ar, \@kmers);
remove_sublist($pconstr_rc_edges_ar, \@kmers);
remove_sublist($pconstr_palindromes_ar, \@kmers);
remove_sublist($pconstr_palindromes_ar, \@s_rev_comps);
if (scalar(@s_rev_comps)) {
    die "Should be no more palindromes remaining, but there are ".scalar(@s_rev_comps)."!\n";
}


my $num_remaining = scalar(@kmers);
if ($k_even) {
    print $out_fh "After removing the edges in these constructs and their reverse complements from the graph,\n".
	"   there are $num_remaining edges remaining in graph.\n";
# check that remaining graph is Eulerian
    print $out_fh "\nChecking that remaining graph is Eulerian...\n";
}
my $conn=check_connected($k, \@kmers, undef, $out_fh);
my $deg_good=check_degrees($k, \@kmers, $out_fh);
if (!$conn or !$deg_good) {
    die "Remaining graph is not Eulerian.\n";
} else {
    print $out_fh "Checked: Graph is Eulerian.\n";
}
print "Computing a partitioning of ".
    ($k_even ? "the rest of ":''). "the graph into two reverse-complementary cycles.\n";
my ($big_cycle_ar, $rc_big_cycle_ar);    
for (my $trial=0; $trial < $num_trials_cycle; $trial++) {
    print $out_fh "\nFinding an \"Eulerian pair\" of cycles:  trial $trial.\n";
    # find forward and reverse-complementary cycles that together cover every edge in graph exactly once.
    my	$period = $oligo_len+2;  # periodicity for incorporating the flanking sequence bases
    if ($flank_seq_flag) {
	($big_cycle_ar, $rc_big_cycle_ar) = find_eulerian_pair($k,\@kmers, $period, $db_str_l, $db_str_r, $out_fh);
    } else {
	($big_cycle_ar, $rc_big_cycle_ar) = find_eulerian_pair($k,\@kmers, undef, undef, undef, $out_fh);
    }
    print $out_fh "\nForward Eulerian cycle has ". scalar(@$big_cycle_ar)." $k-mers: \n\n";
    $expected_num = $num_remaining/2;
    if (scalar(@$big_cycle_ar)!=$expected_num) {
	print "Error: should have $expected_num rather than ".scalar(@$big_cycle_ar). " $k-mers in cycle!\n";
	print "Trying again...\n";
	if ($trial == ($num_trials_cycle-1)) {
	    die "Cannot find an Eulerian pair of cycles!\n";
	} else {
	    next;
	}
    }
    print_array ($big_cycle_ar, $out_fh);    
    print $out_fh "\nForward cycle has ". scalar(@$big_cycle_ar)." $k-mers: \n\n";
    # print_array ($big_cycle_ar);
    print $out_fh "\nReverse-complement Eulerian-like cycle has ". scalar(@$rc_big_cycle_ar) ." $k-mers: \n";
    print_array ($rc_big_cycle_ar, $out_fh);
    if (scalar(@$rc_big_cycle_ar) != $expected_num) {
	die "Error: should have $expected_num rather than ".scalar(@$rc_big_cycle_ar). " $k-mers in rc cycle!\n";
    }
    last;
}

# find insertion points for the pairs of palindromes that have path length 2;
# pick randomly among possibilities, and insert them.
if ($k_even) {
    print $out_fh "\nInserting the ". scalar(@short_ppkeys). " palindrome paths of length 2 into the cycle...\n";
    my ($cycle_insert_pts_sp_hr, $p_insert_pts_sp_hr) = find_insertion_points(\@short_ppkeys, 
									      $big_cycle_ar, 
									      $short_pmatching_hr, $out_fh);
    print_insertion_points($big_cycle_ar, $cycle_insert_pts_sp_hr, $p_insert_pts_sp_hr, 
			   $short_ppaths_hr, $out_fh);
    $big_cycle_ar = insert_short_ppaths($big_cycle_ar, $short_pmatching_hr, $p_insert_pts_sp_hr, 
					$short_ppaths_hr, $out_fh);
    print $out_fh "\nForward cycle with short palindrome paths inserted has ". 
	scalar(@$big_cycle_ar)." $k-mers.\n";
}
my ($reg_constructs_ar, $extra_paths_ar, $total_rem_edges);
print $out_fh "\nBreaking cycle up into constructs...\n";
if ($flank_seq_flag) {
    ($reg_constructs_ar, $extra_paths_ar, $total_rem_edges) =
	break_cycle_into_constructs($k, $oligo_len, $big_cycle_ar, $db_str_l, $db_str_r,
				    $num_trials_break_cycle, $out_fh, $num_oligos_cur_limit);
    print $out_fh "Found solution that leaves ".scalar(@{$extra_paths_ar}). " partial paths.\n";
} else {
    ($reg_constructs_ar) = 
	break_cycle_into_constructs($k, $oligo_len, $big_cycle_ar, undef, undef, undef, 
				    $out_fh, $num_oligos_cur_limit);
}

if (($oligo_len == $oligo_len_short) and (scalar(@{$reg_constructs_ar}) >= $num_oligos_cur_limit)) {
    $num_oligos_cur_limit = $num_oligos_long;
    $oligo_len++;
} else {
    $num_oligos_cur_limit -= scalar(@{$reg_constructs_ar});
}
print $out_fh ''.scalar(@{$reg_constructs_ar})." constructs created via normal construction ".
    "from the edges in the forward cycle.\n";
print $out_fh "\nPrinting the ". scalar(@{$reg_constructs_ar}). " regular construct paths created from forward cycle.\n";
print_construct_paths($reg_constructs_ar, $out_fh);

my ($extra_constructs_ar, $added_edges_ar);

my @copy_kmers=generate_all_k_mers($k);   # all the edges in the original graph
if ($flank_seq_flag and $total_rem_edges) {
    print $out_fh "There are $total_rem_edges edges remaining in the partial paths.\n";
    print $out_fh "\nBuilding constructs from the partial paths...\n";
    ($extra_constructs_ar, $added_edges_ar) = 
	build_extra_constructs($k, $extra_paths_ar, $oligo_len, \@copy_kmers, $db_str_l, $db_str_r, 
			       $num_trials_build_extras, $out_fh, $num_oligos_cur_limit);
    print $out_fh ''.scalar(@{$extra_constructs_ar})." constructs created from the partial paths.\n".
	scalar(@{$added_edges_ar})." repeat (i.e., previously represented) edges were used in these constructs.\n";
    print $out_fh "\nPrinting the ".scalar(@{$extra_constructs_ar})." construct paths from forward cycle's extra partial paths.\n";
}


my @all_constructs = (@{$pconstructs_ar}, @{$reg_constructs_ar});
if ($extra_constructs_ar) {
    push(@all_constructs, @{$extra_constructs_ar});
}
if (!$flank_seq_flag and (scalar(@all_constructs) != $num_oligos) ) {
    die "Unexpected number ".scalar(@all_constructs)." of constructs created (instead of $num_oligos)!\n";
}

my $start_index_reg_constructs = scalar(@{$pconstructs_ar});
my $start_index_extra_constructs = scalar(@{$pconstructs_ar}) + scalar(@{$reg_constructs_ar});
print $out_fh "\nPrinting paths for all constructs now:  \n".
    (scalar(@{$pconstructs_ar}) ? 
     "   0-".($start_index_reg_constructs-1)." index palindrome constructs;\n" : '').
    "   $start_index_reg_constructs-".($start_index_extra_constructs-1)." index regular constructs;\n";
if ($extra_constructs_ar) {
    print $out_fh
	"   $start_index_extra_constructs-".$#all_constructs.
	" index constructs built from partial paths (may contain repeats).\n";
}

print_construct_paths(\@all_constructs, $out_fh);

my ($bad_flag, $edges_by_construct_hr, $repeat_edges_hr);
if ($flank_seq_flag) {
    ($bad_flag, $edges_by_construct_hr, $repeat_edges_hr) =
	check_constructs($k, \@all_constructs, 
			 0, $out_fh, 0,			 
			 1, $added_edges_ar, $start_index_extra_constructs, 
			 $db_str_l, $db_str_r, undef, $num_oligos_short);
} else {
    ($bad_flag, $edges_by_construct_hr, $repeat_edges_hr) =
	check_constructs($k, \@all_constructs, 
			 $oligo_len_short, $out_fh, undef,
			 undef, undef, undef, 
			 undef, undef, undef, $num_oligos_short);			 
}
if ($bad_flag) {    
    die "Constructs do not constitute a valid solution!\n";
} else {
    print $out_fh "\nConstruct design has been validated!\n";
    warn "\nConstruct design has been validated!\n";
}

if ($histfile=$opts{'a'}) {
    $histfile = File::Spec->catfile($outdir, $histfile);
    open ($hist_out_fh, ">$histfile") or die "Cannot open $histfile for write!";
    print $out_fh "\nVerifying construct design -- \n".
	"   histogram data for repeated $k-mers (analyzed with flanking sequences, if specified) will be printed in $histfile.\n";
    print "\nVerifying construct design -- \n".
	"   histogram data for repeated $k-mers (analyzed with flanking sequences, if specified) will be printed in $histfile.\n";
}
if (!$hist_out_fh) {
    $hist_out_fh = $out_fh;
    print $out_fh "\nVerifying construct design -- \n".
	"   histogram data for repeated $k-mers (analyzed with flanking sequences) will be printed.\n";
}
my $all_constructs_wf_ar;

if ($flank_seq_flag) {
    print $out_fh "\nAdding flanking sequences to constructs in order to count repeats.\n";
    print $out_fh "Some of the paths in the design may be appropriately shortened;\n"
}
if (defined ($flank_l) && ($flank_l ne '') && defined ($flank_r) && ($flank_r ne '')) {
    $all_constructs_wf_ar =
	add_flanking_seqs_to_constructs($k, $oligo_len_short, \@all_constructs, $flank_l, $flank_r, $out_fh, 
					$num_oligos_short, 1);
    print $out_fh "\nPrinting the ".scalar(@{$all_constructs_wf_ar})." construct paths with flanking sequences.\n";
    print_construct_paths($all_constructs_wf_ar, $out_fh);
    
    print $out_fh "Validating the construction with the flanking sequences added and analyzing repeated $k-mers.\n";
    my ($edges_by_construct_wf_hr, $repeat_edges_wf_hr);
    
    if ($flank_seq_flag) {
	($bad_flag, $edges_by_construct_wf_hr, $repeat_edges_wf_hr) =
	    check_constructs($k, $all_constructs_wf_ar, 
			     0, $out_fh, $hist_out_fh, 
			     1, undef, undef, 
			     undef, undef, undef, $num_oligos_short);		     
    } else {
	($bad_flag, $edges_by_construct_wf_hr, $repeat_edges_wf_hr) =
	    check_constructs($k, $all_constructs_wf_ar, 
			     ($oligo_len_short+length($flank_l)+length($flank_r)), 
			     $out_fh, $hist_out_fh, 
			     1, undef, undef, 
			     undef, undef, undef, $num_oligos_short);		     
    }
    if ($bad_flag) {    
	warn "Warning: Constructs with flanking sequences added may not constitute a valid solution!\n";
	print $out_fh "Warning: Constructs with flanking sequences added may not constitute a valid solution!\n";
    } else {
	print $out_fh "\nConstruct design with flanking sequences added has been validated.\n";
    }
}
$solnfile = File::Spec->catfile($outdir, $solnfile);
open ($soln_fh, ">$solnfile") or die "Cannot open $solnfile for write!";
print $out_fh "\nPrinting construct design (without flanking sequence) in file $solnfile.\n";

if ($flank_seq_flag) {
    print_constructs($k, \@all_constructs, $soln_fh, $oligo_len_short, $db_str_l, $db_str_r, $num_oligos_short, 1);
} else {
    print_constructs($k, \@all_constructs, $soln_fh, $oligo_len_short, undef, undef, $num_oligos_short);
}
print "\nComplete oligomer library of ".scalar(@all_constructs)." constructs ".
    "printed (without flanking sequence) in file ".$solnfile.".\n\n";

