#!/usr/bin/perl -w

use strict;
$"="\t";

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
	# ASSUMPTION: nodes before edges
	if (!exists($node{$1}) || !exists($node{$3})) {
	    print STDERR "Removing edge $_\n";
	    next;
	}
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
    foreach my $e (@{$edge_in{$n}}) {
	my ($n1,$n2,$EC) = $edge[$e] =~ m/^L\s+(\S+)\s+[+-]\s+(\S+)\s+[+-].*EC:i:(\d+)/;
	my ($d1) = $node{$n1} =~ m/SC:f:(\d+(\.\d+)?)/;
	my ($d2) = $node{$n2} =~ m/SC:f:(\d+(\.\d+)?)/;
	$route_in{$n}++;
	#print STDERR "$n: in $edge[$e] = $d1 $d2 $EC => $route_in{$n}\n";
    }

    foreach my $e (@{$edge_out{$n}}) {
	my ($n1,$n2,$EC) = $edge[$e] =~ m/^L\s+(\S+)\s+[+-]\s+(\S+)\s+[+-].*EC:i:(\d+)/;
	my ($d1) = $node{$n1} =~ m/SC:f:(\d+(\.\d+)?)/;
	my ($d2) = $node{$n2} =~ m/SC:f:(\d+(\.\d+)?)/;
	$route_out{$n}++;
	#print STDERR "$n: out $edge[$e] = $d1 $d2 $EC => $route_out{$n}\n";
    }
}

# Collapse any runs of node A->B->C together if there are no branch points
# For now simplify only + + transitions
restart:
foreach my $n1 (sort keys %node) {
    next unless defined($node{$n1});
    my $n1_out=0;
    my $n2_in=0;
    my $n2;
    my $e_out;
    foreach my $e (@{$edge_out{$n1}}) {
	if (defined($edge[$e])) {
	    $e_out = $e;
	    $edge[$e] =~ m/^L\s+\S+\s+[+-]\s+(\S+)\s+[+-]/;
	    $n2 = $1;
	    $n1_out++;
	}
    }

    my $e_in;
    if ($n1_out == 1) {
	foreach my $e (@{$edge_in{$n2}}) {
	    if (defined($edge[$e])) {
		$e_in = $e;
		$n2_in++;
	    }
	}	    
    }

    # Also check n2-out cannot link back to n1?

    if ($n1_out == 1 && $n2_in == 1 && $e_in == $e_out) {
	if ($edge[$e_in] =~ /^L\s+\S+\s+\+\s+\S+\s+\+/) {
	    print STDERR "Merge $n1 -> $n2\n";
	    #print STDERR "edge_in{$n1}=@{$edge_in{$n1}}\n";
	    #foreach my $e (@{$edge_in{$n1}}) {
	    #	next unless (defined($edge[$e]));
	    #	print STDERR "$e\t$edge[$e]\n";
	    #}
	    #print STDERR "edge_out{$n1}=@{$edge_out{$n1}}\n";
	    #foreach my $e (@{$edge_out{$n1}}) {
	    #	next unless (defined($edge[$e]));
	    #	print STDERR "$e\t$edge[$e]\n";
	    #}
	    #print STDERR "edge_in{$n2}=@{$edge_in{$n2}}\n";
	    #foreach my $e (@{$edge_in{$n2}}) {
	    #	next unless (defined($edge[$e]));
	    #	print STDERR "$e\t$edge[$e]\n";
	    #}
	    #print STDERR "edge_out{$n2}=@{$edge_out{$n2}}\n";
	    #foreach my $e (@{$edge_out{$n2}}) {
	    #	next unless (defined($edge[$e]));
	    #	print STDERR "$e\t$edge[$e]\n";
	    #}
	    #print STDERR "node n1 = $node{$n1}\n";
	    #print STDERR "node n2 = $node{$n2}\n";

	    my (@N1) = split(/\s+/, $node{$n1});
	    my (@N2) = split(/\s+/, $node{$n2});
	    $N1[2] .= $N2[2];     # SEQ
	    my $SC_pos = 0;
	    my $LN = 0;
	    my $KC = 0;
	    for (my $i = 3; $i < scalar(@N1); $i++) {
		if ($N1[$i] =~ /^LN:/) {
		    $LN = length($N1[2]);
		    $N1[$i] = "LN:i:$LN";
		} elsif ($N1[$i] =~ /^KC:/) {
		    for (my $j = 3; $j < scalar(@N2); $j++) {
			next unless $N2[$j] =~ /^KC:/;
			($KC) = $N1[$i] =~ /KC:i:(\d+)/;
			$KC += [$N2[$j] =~ /KC:i:(\d+)/]->[0];
			#print STDERR "KC: $N1[$i]+$N2[$j] = $KC\n";
			$N1[$i] = "KC:i:$KC";
		    }
		} elsif ($N1[$i] =~ /SC:/) {
		    $SC_pos = $i;
		}
	    }

	    if ($SC_pos) {
		my $sc = $KC/$LN;
		$N1[$SC_pos] = "SC:f:$sc";
	    }
	    $node{$n1}="@N1";
	    #print STDERR "@N1\n";

	    #print STDERR "Undef $edge[$e_in]\n";
	    $edge[$e_in] = undef;
	    #print STDERR "edge_out{$n1}=@{$edge_out{$n1}}\n";
	    $edge_out{$n1} = ();
	    #print STDERR "edge_out{$n2}=@{$edge_out{$n2}}\n";
	    foreach my $e (@{$edge_out{$n2}}) {
		next unless (defined($edge[$e]));
		my @E = split(/\t/,$edge[$e]);
		if ($E[1] eq $n2) {
		    @E[1] = $n1;
		} else {
		    $E[3] = $n1;
		}
		$edge[$e] = "@E";
		#print STDERR "New edge $e => $edge[$e]\n";
		push(@{$edge_out{$n1}}, $e);
		#print STDERR "$n1 IN:  @{$edge_in{$n1}}\n";
		#print STDERR "$n1 OUT: @{$edge_out{$n1}}\n";
		#print STDERR "$n2 IN:  @{$edge_in{$n2}}\n";
		#print STDERR "$n2 OUT: @{$edge_out{$n2}}\n";
	    }
	    #print STDERR "edge_out{$n1}=@{$edge_out{$n1}}\n";
	    $node{$n2} = undef;
	    #print STDERR "UNDEF node $n2\n";

	    goto restart;
	    last;
	}
    }
}

# Print up new graph
foreach my $n (sort keys %node) {
    print $node{$n},"\n" if defined($node{$n});
}

foreach my $n (@edge) {
    print "$n\n" if (defined($n));
}
