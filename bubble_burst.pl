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
@ARGV = ("/dev/stdin") if scalar(@ARGV) == 0;

open(my $gfa,   "<", shift(@ARGV)) || die;
while (<$gfa>) {
    chomp($_);
    my @F = split(/\s+/, $_);

    if ($F[0] eq "S") {
	$gfa{$F[1]}{seq} = uc($F[2]);
	$incoming_edge_count{$F[1]}=0 if !defined($incoming_edge_count{$F[1]});
    } elsif ($F[0] eq "L") {
	# Push edges to from node
	push(@{$gfa{$F[1]}{edge}}, [@F]);
	$incoming_edge_count{$F[3]}++;
    }
}

# Find bubbles
new_pass:
print STDERR "PASS\n";

my $edited = 0;

local $"="\t";
foreach my $node (keys %gfa) {
    next unless $gfa{$node}{edge};

    my $nedge = scalar @{$gfa{$node}{edge}};
    my $indel_bubble = 0;
    my $snp_bubble = 0;

    # TODO: we could in theory have a 3-way SNP.
    # Eg A->{B,C,D}->E where B, C and D are all 1bp long.
    if ($nedge == 2) {
	my $e1 = @{$gfa{$node}{edge}->[0]}[3];
	my $e2 = @{$gfa{$node}{edge}->[1]}[3];

	# B->C->D or  B---->D  (B edge 0)
	# B---->D     B->C->D  (B edge 1)
	if ($incoming_edge_count{$e1} == 1 &&
	    $incoming_edge_count{$e2} == 1 &&
	    defined $gfa{$e1}{edge} &&
	    defined $gfa{$e2}{edge} &&
	    scalar @{$gfa{$e1}{edge}} == 1 &&
	    scalar @{$gfa{$e2}{edge}} == 1 &&
	    length($gfa{$e1}{seq}) < $min_size &&
	    length($gfa{$e2}{seq}) < $min_size) {
	    # A->B->D   (A edge 0)
	    # A->C->D   (A edge 1)
	    my $dest = @{$gfa{$e1}{edge}->[0]}[3];
	    print STDERR "SNP:\t$node -> {$e1,$e2} -> $dest\n";
	    $snp_bubble = 1;
	    $edited = 1;
	}

	if ($incoming_edge_count{$e1} == 1 &&
	    $incoming_edge_count{$e2} == 2 &&
	    defined $gfa{$e1}{edge} &&
	    scalar @{$gfa{$e1}{edge}} == 1 &&
	    @{$gfa{$e1}{edge}->[0]}[3] eq $e2 &&
	    length($gfa{$e1}{seq}) < $min_size) {
	    # A->B->C   (A edge 0)
	    # A---->C   (A edge 1)
	    print STDERR "INDEL1:\t$node -> {$e1 -> $e2, $e2}\n";
	    $indel_bubble = 1;
	    $edited = 1;
	}

	if ($incoming_edge_count{$e2} == 1 &&
	    $incoming_edge_count{$e1} == 2 &&
	    defined $gfa{$e2}{edge} &&
	    scalar @{$gfa{$e2}{edge}} == 1 &&
	    @{$gfa{$e2}{edge}->[0]}[3] eq $e1 &&
	    length($gfa{$e2}{seq}) < $min_size) {
	    # A---->C   (A edge 0)
	    # A->B->C   (A edge 1)
	    print STDERR "INDEL2:\t$node -> {$e2 -> $e1, $e1}\n";
	    $indel_bubble = 2;
	    $edited = 1;
	}
    }

    if ($snp_bubble) {
	# Burst by picking sequence from one edge and merging into prev node
	#print STDERR "Burst bubble $gfa{$node}{edge}[0]->[3] $gfa{$node}{edge}[1]->[3]\n";
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

    if ($indel_bubble == 1) {
	# Merge first edge into node and remove from graph
	my $e1 = @{$gfa{$node}{edge}->[0]}[3];
	my $e2 = @{$gfa{$node}{edge}->[1]}[3];

	$gfa{$node}{seq} .= $gfa{$e1}{seq};
	delete $gfa{$e1};
	shift(@{$gfa{$node}{edge}});
	$incoming_edge_count{$e2}--;
    }

    if ($indel_bubble == 2) {
	# Merge second edge into node and remove from graph
	my $e1 = @{$gfa{$node}{edge}->[0]}[3];
	my $e2 = @{$gfa{$node}{edge}->[1]}[3];

	$gfa{$node}{seq} .= $gfa{$e2}{seq};
	delete $gfa{$e2};
	pop(@{$gfa{$node}{edge}});
	$incoming_edge_count{$e1}--;
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

goto new_pass if $edited;

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
