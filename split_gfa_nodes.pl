#!/usr/bin/perl -w
use strict;

my $gfa = shift(@ARGV);
my $max_size = shift(@ARGV);

# Parse GFA
my %node;
my @edge;
my $edge_num=0;
my %edge_in;  # index into @edge above
my %edge_out; # index into @edge above
open(GFA, "<", $gfa) || die;
while (<GFA>) {
    if (/^S\s+(\S+)/) {
        $node{$1} = $_;
    } elsif (/^L\s+(\S+)\s+(.)\s+(\S+)\s+(.)\s+(.*)/) {
	$edge[$edge_num] = $_;
        push(@{$edge_out{$1}}, $edge_num);
        push(@{$edge_in{$3}},  $edge_num);
	$edge_num++;
    }
}
close(GFA);

# Split long nodes
local $"="\t";
foreach my $n (sort keys %node) {
    my @F = split("\t", $node{$n});
    next if length($F[2]) <= $max_size;
    my $len = length($F[2]);
    #print "$n $len\n";
    my $last = undef;

    my @out_edges = exists($edge_out{$n}) ? @{$edge_out{$n}} : ();
    my $last_edge = int(($len-1)/$max_size)*$max_size;
    for (my $i=0;$i<$len;$i+=$max_size) {
	# S nodes
	my $sub_len = ($i+$max_size<$len?$i+$max_size:$len)-$i;
	my $sub = substr($F[2],$i,$sub_len);
	my $sub_node = "$n.sub".sprintf("%08d",$i);
	my @sub = @F;
	$sub[1] = $sub_node;
	$sub[2] = $sub;
	# FIXME: also fix LN: length field
	$node{$sub_node} = "@sub";

	# L edges
	if (!$last) {
	    #print "$sub_node\n";
	    my @new_list = ();
	    foreach my $e (@{$edge_in{$n}}) {
		my @E = split("\t", $edge[$e]);
		if ($E[4] eq "-") {
		    $E[3] = "$n.sub".sprintf("%08d",$last_edge);
		} else {
		    $E[3] = $sub_node;
		}
		$edge[$e] = "@E";
		#print "SPLIT \t@E";
	    }
	} else {
	    #print "$sub_node\n";
	    my @E = ("L", $last, "+", $sub_node, "+", "0M\n");
	    push(@edge, "@E");
	    #print "SPLIT2 \t@E";
	}
	$last = $sub_node;
    }

    foreach (@out_edges) {
	my @E = split("\t", $edge[$_]);
	if ($E[2] eq "-") {
	    $E[1]="$n.sub".sprintf("%08d",0);
	} else {
	    $E[1]=$last;
	}
	$edge[$_] = "@E";
	#print "LAST @E";
    }

    delete $node{$n};
}

# Print graph
foreach (sort keys %node) {
    print $node{$_};
}
foreach (@edge) {
    print $_;
}
#for (my $i=0; $i<=$#edge; $i++) {
#    print "$i $edge[$i]";
#}
