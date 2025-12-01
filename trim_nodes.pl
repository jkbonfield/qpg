#!/usr/bin/perl -w

# Removes nodes that aren't covered (weight <= $W) or are withing $D edges
# from such a node.

# TODO: also remove disconnected tips with weight <= $W.
# Equivalent to saying keep everything between nodes of weight > $D.

my $D = 3;
my $W = 0; # minimum depth to keep

use strict;
$"="\t";

use Getopt::Long;

GetOptions("w|min-node-weight=f" => \$W,
	   "d|node-distance=i"   => \$D)
or die("Usage: trim_edges.pl [-l|-level trim-level] [-n|-min-node-depth depth] [-e|-min-edge-depth depth]\n");

# Parse GFA
my %seq;      # node sequence
my %node;     # node GFA line (S)
my @edge;     # edge GFA line (L)
my $edge_num=0;
my %edge_in;  # index into @edge above
my %edge_out; # index into @edge above
local $"="\t";
while (<>) {
    chomp();
    if (/^S\s+(\S+)/) {
        $node{$1} = $_;
        my @N = split("\t", $_);
        $seq{$N[1]} = $N[2];
    } elsif (/^L\s+(\S+)\s+(.)\s+(\S+)\s+(.)/) {
	#print STDERR "EDGE: $edge_num $_ // $2 $4\n";
        $edge[$edge_num] = $_;
	if ($2 eq "+") {
	    push(@{$edge_out{$1}}, $edge_num);
	} else {
	    push(@{$edge_in{$1}}, $edge_num);
	}
	if ($4 eq "+") {
	    push(@{$edge_in{$3}},  $edge_num);
	} else {
	    push(@{$edge_out{$3}},  $edge_num);
	}
	#print STDERR "edge_in{$1}=@{$edge_in{$1}}\n"   if exists($edge_in{$1});
	#print STDERR "edge_in{$3}=@{$edge_in{$3}}\n"   if exists($edge_in{$3});
	#print STDERR "edge_out{$1}=@{$edge_out{$1}}\n" if exists($edge_out{$1});
	#print STDERR "edge_out{$3}=@{$edge_out{$3}}\n" if exists($edge_out{$3});

        $edge_num++;
    }
}

# Fill out %node_list for all nodes within $d distance from this node.
# (Irrespective of depth)
sub connected_nodes {
    my ($n,$d,$node_list) = @_;
    $node_list->{$n} = 1;
    foreach my $e (@{$edge_out{$n}},@{$edge_in{$n}}) {
	my ($n1,$n2,$EC) = $edge[$e] =~ m/^L\s+(\S+)\s+[+-]\s+(\S+)\s+[+-].*EC:i:(\d+)/;
	my $ne = $n1 eq $n ? $n2 : $n1;
	if ($d > 1 && !exists($node_list->{$ne})) {
	    connected_nodes($ne, $d-1, $node_list);
	}
    }
}

# For each node with zero-weight, find all others D steps away.
# Fill out %node_list with 1 for weighty and 0 otherwise.
sub disconnected_nodes {
    my ($n,$d,$node_list,$val) = @_;
    $node_list->{$n} = 0;
    foreach my $e (@{$edge_out{$n}},@{$edge_in{$n}}) {
	my ($n1,$n2,$EC) = $edge[$e] =~ m/^L\s+(\S+)\s+[+-]\s+(\S+)\s+[+-].*EC:i:(\d+)/;
	my $ne = $n1 eq $n ? $n2 : $n1;
	my ($SC) = $node{$ne} =~ m/SC:f:(\d+(\.\d+)?)/;

	#print "$ne EC=$EC SC=$SC\n";
	if ($SC>$W || $EC>$W) {
	    $node_list->{$ne} = $val;
	} else {
	    if ($d > 1 && !exists($node_list->{$ne})) {
		disconnected_nodes($ne, $d-1, $node_list, $val);
	    }
	    $node_list->{$ne} = 0;
	}
    }
}

sub disconnected_nodes_dir {
    my ($n,$d,$node_list) = @_;
    $node_list->{$n} = 0;

    foreach my $e (@{$edge_out{$n}}) {
	my ($n1,$n2,$EC) = $edge[$e] =~ m/^L\s+(\S+)\s+[+-]\s+(\S+)\s+[+-].*EC:i:(\d+)/;
	my $ne = $n1 eq $n ? $n2 : $n1;
	my ($SC) = $node{$ne} =~ m/SC:f:(\d+(\.\d+)?)/;

	#print "$ne EC=$EC SC=$SC\n";
	if ($SC>$W || $EC>$W) {
	    $node_list->{$ne} = 1;
	} else {
	    if ($d > 1 && !exists($node_list->{$ne})) {
		disconnected_nodes($ne, $d-1, $node_list, 1);
	    }
	    $node_list->{$ne} = 0;
	}
    }
    foreach my $e (@{$edge_in{$n}}) {
	my ($n1,$n2,$EC) = $edge[$e] =~ m/^L\s+(\S+)\s+[+-]\s+(\S+)\s+[+-].*EC:i:(\d+)/;
	my $ne = $n1 eq $n ? $n2 : $n1;
	my ($SC) = $node{$ne} =~ m/SC:f:(\d+(\.\d+)?)/;

	#print "$ne EC=$EC SC=$SC\n";
	if ($SC>$W || $EC>$W) {
	    $node_list->{$ne} = -1;
	} else {
	    if ($d > 1 && !exists($node_list->{$ne})) {
		disconnected_nodes($ne, $d-1, $node_list, -1);
	    }
	    $node_list->{$ne} = 0;
	}
    }
}

# For each node with weight, find all other D steps away and mark them to keep.
my %node_list = ();
my %linked_nodes = ();
foreach my $n (sort keys %node) {
    my ($d) = $node{$n} =~ m/SC:f:(\d+(\.\d+)?)/;
    if ($d <= $W) {
	# Zero depth: so see if linking two non-zero depth nodes
	my %dis_nodes = ();
	disconnected_nodes_dir($n, $D,\%dis_nodes);
	my $nw = 0; # number of weighty nodes
	my $nwf = 0; # number of weighty nodes fwd
	my $nwb = 0; # number of weighty nodes bck
	#print "$n:";
	foreach (keys(%dis_nodes)) {
	    $nw += abs($dis_nodes{$_});
	    #print " $dis_nodes{$_}:$_";
	    if ($dis_nodes{$_} > 0) {
		$nwf++;
	    } elsif ($dis_nodes{$_} < 0) {
		$nwb++;
	    } 
	}
	#$linked_nodes{$n} = $nw;
	$linked_nodes{$n} = ($nwf>0) + ($nwb>0);
	#print " = $linked_nodes{$n} +/$nwf -/$nwb\n";
    } else {
	# Non-zero depth, so track what else we link to
	connected_nodes($n, $D,\%node_list);
    }
}

# Print up new graph
foreach my $n (sort keys %node) {
    if (defined($node{$n}) && exists($node_list{$n}) &&
	(!exists($linked_nodes{$n}) || $linked_nodes{$n} > 1)) {
	print $node{$n},"\n";
    }
    #print $node{$n},"\n" if defined($node{$n}) && exists($node_list{$n});
}

foreach my $e (@edge) {
    next unless defined($e);
    my ($n1,$n2,$EC) = $e =~ m/^L\s+(\S+)\s+[+-]\s+(\S+)\s+[+-].*EC:i:(\d+)/;
    next unless exists($node_list{$n1}) && exists($node_list{$n2});
    print "$e\n";
}
