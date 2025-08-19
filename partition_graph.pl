#!/usr/bin/perl -w

# Partitions a graph based on observed node weights.
# Any node within DIST edges from a node with weight > WEIGHT will be
# included in that group.
# We find all such groups and generate sub-GFAs for each.

use strict;

if (scalar(@ARGV) < 3) {
    print STDERR "Usage: partition_graph.pl in.gfa dist min_weight [min_report_weight]\n";
    exit 1;
}
my $gfa = shift(@ARGV);
my $dist = shift(@ARGV);
my $weight = shift(@ARGV) + 1;
my $weight_report = scalar(@ARGV) ? shift(@ARGV) : $weight/2;

# Parse GFA
my %node;     # node GFA line (S)
my @edge;     # edge GFA line (L)
my $edge_num=0;
my %edge_in;  # index into @edge above
my %edge_out; # index into @edge above
open(my $fh, "<", $gfa) || die;
local $"="\t";
while (<$fh>) {
    chomp();
    if (/^S\s+(\S+)/) {
        $node{$1} = $_;
    } elsif (/^L\s+(\S+)\s+(.)\s+(\S+)\s+(.)/) {
        $edge[$edge_num] = $_;
        push(@{$edge_out{$1}}, $edge_num);
        push(@{$edge_in{$3}},  $edge_num);
        $edge_num++;
    }
}
close($fh);

my %visited;  # Which nodes we've visited
my $nvisited = 0;

sub partition {
    my ($n, $ldist, $depth, @node_list) = (@_);
    return if (exists($visited{$n}));
    $visited{$n}=1;
    $nvisited++;

    $node{$n} =~ m/SC:f:(-?\d+(\.\d+)?)/;
    my $node_count = $1;
    #print "  "x$depth, "Node $n, Weight $node_count, dist $ldist\n";

    if ($node_count < $weight) {
	$ldist--;
    } else {
	$ldist = $dist;
    }

    push(@node_list, $n);

    foreach my $e (@{$edge_out{$n}}) {
	my @F = split("\t", $edge[$e]);
	$edge[$e] =~ m/EC:i:(\d+)/;
	my $edge_count = $1;
	$node{$F[3]} =~ m/SC:f:(\d+(\.\d+)?)/;
	my $node_ocount = $1;

	# Either/or for node and edge weights is sufficient to include.
	my $edist = $ldist;
	$edist++ if ($node_count < $weight && ($edge_count >= $weight || $node_ocount >= $weight));
	#print "  "x$depth, "  Edge out $e: $F[1]/$node_count $F[3]/$node_ocount $ldist/$edist\n";
	push(@node_list, partition($F[3], $edist, $depth+1)) if $edist>0;
    }

    foreach my $e (@{$edge_in{$n}}) {
	my @F = split("\t", $edge[$e]);
	$edge[$e] =~ m/EC:i:(\d+)/;
	my $edge_count = $1;
	$node{$F[1]} =~ m/SC:f:(\d+(\.\d+)?)/;
	my $node_icount = $1;

	# Either/or for node and edge weights is sufficient to include.
	my $edist = $ldist;
	$edist++ if ($node_count < $weight && ($edge_count >= $weight || $node_icount >= $weight));
	#print "  "x$depth, "  Edge in $e: $F[1]/$node_icount $F[3]/$node_count $ldist/$edist\n";
	push(@node_list, partition($F[1], $edist, $depth+1)) if $edist>0;
    }

    return @node_list;
}


my $gnum =0;
sub sub_graph {
    open(my $fh, ">$gfa.sub_graph.$gnum") || die;
    $gnum++;

    print "======\n";
    my %included = ();
    foreach my $n (@_) {
	$included{$n}++;
    }

    # Nodes
    foreach my $n (@_) {
	print $fh "$node{$n}\n";
    }

    # Edges
    my %edge_done = ();
    foreach my $n (@_) {
	foreach my $e (@{$edge_in{$n}}) {
	    next if exists($edge_done{$e});
	    my @F = split("\t", $edge[$e]);
	    if (exists($included{$F[1]}) && exists($included{$F[3]})) {
		print $fh "$edge[$e]\n";
		$edge_done{$e}=1;
	    }
	}

	foreach my $e (@{$edge_out{$n}}) {
	    next if exists($edge_done{$e});
	    my @F = split("\t", $edge[$e]);
	    if (exists($included{$F[1]}) && exists($included{$F[3]})) {
		print $fh "$edge[$e]\n";
		$edge_done{$e}=1;
	    }
	}
    }
    close($fh);
}

sub nat_sort {
    my ($s) = @_;
    my ($d) = $s =~ /(\d+)/;
    return $d;
}

while ($nvisited < scalar(keys %node)) {
    foreach my $n (sort keys(%node)) {
	if (!$visited{$n}) {
	    my @node_list = partition($n, $dist, 0);
	    @node_list = sort { nat_sort($a) <=> nat_sort($b) } @node_list;

	    # Also check the node list has at least one node with minimum depth
	    my $non_zero_depth = 0;
	    foreach my $n2 (@node_list) {
		$node{$n2} =~ m/SC:f:(\d+(\.\d+)?)/;
		$non_zero_depth++ if ($1 >= $weight_report);
	    }

	    if ($non_zero_depth) {
		print "Nodes: $non_zero_depth @node_list\n";
		sub_graph @node_list;
	    }
	}
    }
}

# See /nfs/users/nfs_j/jkb/lustre/tmp/pf2-base_illumina_km_00104/seq_1089-0082-#1#1.gfa
# pathfinder: s1-16 (minus eg s2 etc), 31-33, s47, s54
#             s22-30, s37-39
#
# this (2 5): s1-18, s31-36, s41-51, s54
#             s21-30, s37-40, s44, s52-53
