#!/usr/bin/perl -w

# Removes simple tiny SNP bubbles in GFA.  We're happy to have SNPs when
# aligning against a longer node so this isn't a major issue.
#
# Note: this does not support updating CIGAR strings, or anything other than
# S and L lines.

use strict;

# GFA indexed by node name
# gfa{}{seq}    node seq content
# gfa{}{edge}   [node dir node dir cigar]
my %gfa;
my %incoming_edge_count;

my $min_size = 3;


# Parse GFA; minimally
open(my $gfa,   "<", shift(@ARGV)) || die;
while (<$gfa>) {
    chomp($_);
    my @F = split(/\s+/, $_);
    #next unless scalar(@F) && $F[0] eq "S"; # skip other fields for now
    if ($F[0] eq "S") {
	$gfa{$F[1]}{seq} = uc($F[2]);
    } elsif ($F[0] eq "L") {
	# Push edges to from node
	push(@{$gfa{$F[1]}{edge}}, [@F]);
	$incoming_edge_count{$F[3]}++;
    }
}

# Find bubbles
local $"="\t";
foreach my $node (keys %gfa) {
#    print "$node\n";
    my $bubble = 1;
    my $nedge = 0;
    foreach my $edge (@{$gfa{$node}{edge}}) {
#	print "\t@$edge\n";
	$nedge++;
	$bubble = 0 if (length($gfa{$edge->[3]}{seq}) >= $min_size ||
	    $incoming_edge_count{$edge->[3]} != 1)
    }

    if ($bubble && $nedge > 1) {
	# Burst by picking sequence from one edge and merging into prev node
	print STDERR "Burst bubble $gfa{$node}{edge}[0]->[3] $gfa{$node}{edge}[1]->[3]\n";
	$gfa{$node}{seq} .= $gfa{$gfa{$node}{edge}[0]->[3]}{seq};

	# Take all out edges from bubbles and copy over to this node
	my @new_edge = ();
	my %new_edge = ();
	foreach my $edge (@{$gfa{$node}{edge}}) {
	    my $out_node = $edge->[3];
	    my @F = @{$gfa{$out_node}{edge}};
	    foreach my $out_edge (@{$gfa{$out_node}{edge}}) {
		$out_edge->[1] = $node;
		#print "OUT EDGE $out_node @$out_edge\n";
		$new_edge{"@$out_edge"}=1;
		push(@new_edge, [@$out_edge]);
		$incoming_edge_count{$out_edge->[3]}--;
	    }
	    # Why doesn't this always work?  Ref counting?
	    delete $gfa{$out_node};
	}

	# Uniquify new outgoing edges
	$gfa{$node}{edge} = ();
	foreach (keys %new_edge) {
	    my @F = split("\t", $_);
	    push(@{$gfa{$node}{edge}}, [@F]);
	    $incoming_edge_count{$F[3]}++;
	}
    }
}

# Squash A->B->C to ABC
foreach my $node (sort keys %gfa) {
    next unless $gfa{$node}{seq};
    next unless exists($gfa{$node}{edge});

    if (scalar(@{$gfa{$node}{edge}}) == 1 &&
	$incoming_edge_count{$gfa{$node}{edge}[0]->[3]} == 1) {
	print STDERR "Merge linear $node $gfa{$node}{edge}[0]->[3]\n";

	my $onode = $gfa{$node}{edge}[0]->[3];
	$gfa{$node}{seq} .= $gfa{$onode}{seq};
	$gfa{$node}{edge} = $gfa{$onode}{edge};
	foreach my $edge (@{$gfa{$node}{edge}}) {
	    $edge->[1] = $node;
	}
	
	delete $gfa{$onode};
    }
}

# Print up new graph
foreach my $node (sort keys %gfa) {
    next unless $gfa{$node}{seq};
    print "S\t$node\t$gfa{$node}{seq}\n";
#    if (exists $incoming_edge_count{$node}) {
#	print STDERR "$node incoming $incoming_edge_count{$node}\n";
#    }
    foreach my $edge (@{$gfa{$node}{edge}}) {
	print "@$edge\n";
    }
}
