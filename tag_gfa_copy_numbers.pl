#!/usr/bin/perl -w

# Analyses a GFA file to estimate the coverage depth.
# From here we turn it to a node coverage estimation. (Floating point)

use strict;

my $verbose = 0;
if ($ARGV[0] eq "-v") {
    shift(@ARGV);
    $verbose++;
}

# Parse GFA
my @node_order;
my %seq;      # node sequence
my %node;     # node GFA line (S)
my @edge;     # edge GFA line (L)
my $edge_num=0;
my %edge_in;  # index into @edge above
my %edge_out; # index into @edge above
my %self_loop;
local $"="\t";
while (<>) {
    chomp();
    if (/^S\s+(\S+)/) {
	push(@node_order, $1);
        $node{$1} = $_;
        my @N = split("\t", $_);
        $seq{$N[1]} = $N[2];
    } elsif (/^L\s+(\S+)\s+(.)\s+(\S+)\s+(.)/) {
	# ASSUMPTION: nodes before edges
	if (!exists($node{$1}) || !exists($node{$3})) {
	    print STDERR "Removing edge $_\n";
	    next;
	}
        $edge[$edge_num] = $_;
	$self_loop{$1} = 1 if $1 eq $3;
#	if ($2 eq "+") {
#	    push(@{$edge_out{$1}}, $edge_num);
#	} else {
#	    push(@{$edge_in{$1}}, $edge_num);
#	}
#	if ($4 eq "+") {
#	    push(@{$edge_in{$3}},  $edge_num);
#	} else {
#	    push(@{$edge_out{$3}},  $edge_num);
#	}
        $edge_num++;
    }
}


my $min_depth = 0.1;

# An array of raw depths
my @depths = ();
my $avg_depth = 0;
my $total = 0;
my $tlen = 0;
foreach my $n (sort keys %node) {
    my ($d) = $node{$n} =~ m/SC:f:(-?\d+(\.\d+)?)/;
    if ($d > $min_depth) {
	# For initial average, exclude self loops
	next if exists($self_loop{$n});

	push(@depths, $d);
	$total += $d*(length($seq{$n}));
 	$tlen += (length($seq{$n}));
    }
}
$avg_depth = $total / $tlen;
if ($verbose) {
    print "Avg depth $avg_depth\n";
}

# When running on an trimmed sub-graph our stats can be skewed.
# Rerun on the entire file to get a better estimate of the starting average
# depth.
#
# Detrimental to minigraph based GFAs.
# Marginally helpful to alfapang graphs, but maybe not worth it.
# Relies on the print "Avg depth" above, but we're not using this now
# anyway.
#
#if ($ARGV =~ /\.edited\./) {
#    my ($base) = $ARGV =~ m/(.*)\.edited\..*/;
#    $_=`tag_gfa_copy_numbers.pl $base.gfa | grep "Avg depth"`;
#    my @F = split(/\s+/, $_);
#    print "Avg $F[-1]\n";
#    $avg_depth = $F[-1];
#}

# Alternative; try fitting depth to D, 1*D, 2*D, etc.
my $best_try=$avg_depth;
my $best_delta=1e10;
for (my $try=$avg_depth/1.3; $try<$avg_depth*1.3; $try+=0.5) {
    my $delta = 0;

    foreach my $n (sort keys %node) {
	my ($d) = $node{$n} =~ m/SC:f:(-?\d+(\.\d+)?)/;
	if ($d > $avg_depth/4) {
	    my $mult = int(0.5+$d/$try);
	    my $diff = abs($d-$mult*$try)/$try;
	    #my $diff = $d % int($try+0.5)/$try; # better with minigraph?
	    $delta += $diff*(length($seq{$n}));
	    #printf("    %s\t%.2f\t%d\t%d\t%.2f\t%8.2f\n",
	    #	   $n,$d,$mult,length($seq{$n}),$diff,
	    #	   $diff*(length($seq{$n})));
	}
    }

    # Normalise by try itself?  Smaller values inherently have smaller
    # remainders.
    #$delta /= $try;

    if ($best_delta > $delta) {
	$best_delta = $delta;
	$best_try = $try;
    }
    #print "try $try, delta=$delta\n";
}

if ($verbose) {
    print "Possibly 1-depth = $best_try\n";
}

$avg_depth = $best_try;

# Another possible test is if our copy number estimations are now mostly even
# and rarely odd, then we have half the depth.  Similarly if they go in
# multiples of 3.  There's likely too much noise for this to work however,
# and if it's not noisy then the chances are we found the correct depth. 


# Report
# TODO: maybe output a new GFA with CN:i: or CN:f: values?

if ($verbose) {
    foreach my $n (@node_order) {
	my ($d) = $node{$n} =~ m/SC:f:(-?\d+(\.\d+)?)/;
	printf("%s\t%.2f\t%.2f\t%d\n", $n, $d, $d/$avg_depth,
	       int($d/$avg_depth+0.5));
    }
} else {
    my @copy;
    foreach my $n (@node_order) {
	my ($d) = $node{$n} =~ m/SC:f:(-?\d+(\.\d+)?)/;
	push(@copy, int($d/$avg_depth+0.5));
    }
    local $"=",";
    print "@copy\n";
}

__END__
compare:

~/lustre/tmp/_ga_pf2k-t2.25_00100/seq_1050-0003-#1#1.edited.gfa  (27.9)
~/lustre/tmp/_ga_pf2k-t2.25_00100/seq_1089-0054-#1#1.edited.gfa  (113.5)

[NB: works better on the unedited.gfa to get initial depth]

[But pathfinder picks the wrong values too]


See also ../_km_pf2k-t1.25_00102/seq_1055-0054-#1#1.edited.gfa for
demonstration of graph simplification which hasn't worked.
