#!/usr/bin/perl -w

use strict;

# Parse GFA
my %seq;      # node sequence
my %node;     # node GFA line (S)
my @edge;     # edge GFA line (L)
my $edge_num=0;
my %edge_in;  # index into @edge above
my %edge_out; # index into @edge above
my %inv_edge_in;
my %inv_edge_out;
local $"="\t";
while (<>) {
    chomp();
    if (/^S\s+(\S+)/) {
        $node{$1} = $_;
        my @N = split("\t", $_);
        $seq{$N[1]} = $N[2];
    } elsif (/^L\s+(\S+)\s+(.)\s+(\S+)\s+(.)/) {
        $edge[$edge_num] = $_;
        push(@{$edge_out{$1}}, $edge_num);
        push(@{$edge_in{$3}},  $edge_num);

        push(@{$inv_edge_out{$3}}, $edge_num);
        push(@{$inv_edge_in{$1}}, $edge_num);

        $edge_num++;
    }
}

# Count valid routes in and out per node, where valid implies sufficient depth
# to this node and sufficient edge weight.
my %route_in;
my %route_out;
foreach my $n (sort keys %node) {
    $route_in{$n}=0;
    $route_out{$n}=0;
#    my $min_node_depth = 10;
#    my $min_edge_depth = 10;
    my $min_node_depth = 3;
    my $min_edge_depth = 1;
    foreach my $e (@{$edge_in{$n}}) {
	my ($n1,$n2,$EC) = $edge[$e] =~ m/^L\s+(\S+)\s+[+-]\s+(\S+)\s+[+-].*EC:i:(\d+)/;
	my ($d1) = $node{$n1} =~ m/SC:f:(\d+(\.\d+)?)/;
	my ($d2) = $node{$n2} =~ m/SC:f:(\d+(\.\d+)?)/;
	$route_in{$n}++ if ($d1 >= $min_node_depth &&
			    $d2 >= $min_node_depth &&
			    $EC >= $min_edge_depth);
    }

    foreach my $e (@{$edge_out{$n}}) {
	my ($n1,$n2,$EC) = $edge[$e] =~ m/^L\s+(\S+)\s+[+-]\s+(\S+)\s+[+-].*EC:i:(\d+)/;
	my ($d1) = $node{$n1} =~ m/SC:f:(\d+(\.\d+)?)/;
	my ($d2) = $node{$n2} =~ m/SC:f:(\d+(\.\d+)?)/;
	$route_out{$n}++ if ($d1 >= $min_node_depth &&
			     $d2 >= $min_node_depth &&
			     $EC >= $min_edge_depth);
    }
}

# Go through the nodes again looking for nodes with edges that go between nodes
# with a valid route in/out, but the edge has zero count.
# These are candidates for removal without breaking the graph into fragments.
foreach my $n (sort keys %node) {
    #print "$n routes: IN $route_in{$n}     OUT $route_out{$n}\n";
    next unless ($route_in{$n} && $route_out{$n});

    foreach my $e (@{$edge_in{$n}}) {
	next unless defined($edge[$e]);
	my ($n1,$n2,$EC) = $edge[$e] =~ m/^L\s+(\S+)\s+[+-]\s+(\S+)\s+[+-].*EC:i:(\d+)/;
	if ($EC == 0 &&
	    $route_in{$n1} && $route_out{$n1} &&
	    $route_in{$n2} && $route_out{$n2}) {
	    #print "\tCull $n1 to $n2\n";
	    $edge[$e] = undef;
	}
    }

    foreach my $e (@{$edge_out{$n}}) {
	next unless defined($edge[$e]);
	my ($n1,$n2,$EC) = $edge[$e] =~ m/^L\s+(\S+)\s+[+-]\s+(\S+)\s+[+-].*EC:i:(\d+)/;
	if ($EC == 0 &&
	    $route_in{$n1} && $route_out{$n1} &&
	    $route_in{$n2} && $route_out{$n2}) {
	    #print "\tCull $n1 to $n2\n";
	    $edge[$e] = undef;
	}
    }
}


# Print up new graph
foreach my $n (sort keys %node) {
    print $node{$n},"\n";
}

foreach my $n (@edge) {
    print "$n\n" if (defined($n));
}
