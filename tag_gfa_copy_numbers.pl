#!/usr/bin/perl -w

# Analyses a GFA file to estimate the coverage depth.
# From here we turn it to a node coverage estimation. (Floating point)

use strict;
use Getopt::Long;

my $verbose = 0;
my $offset = 0;
my $float_copy = 0;
my $min_depth = 0.1;
my $depth_div = 4;
my $max_copy = 10;

GetOptions("v|verbose"    => \$verbose,
	   "offset=f"     => \$offset,
	   "f|float-copy" => \$float_copy,
           "m|min-depth:f"=> \$min_depth,
           "d|depth-div:f"=> \$depth_div,
           "C|max-copy:i" => \$max_copy);

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
	my ($SC) = $node{$N[1]} =~ m/SC:f:(-?\d+(\.\d+)?)/;
	if (!defined($SC)) {
	    my ($KC) = $node{$N[1]} =~ m/KC:i:(\d+)/;
	    if (!defined($KC)) {
		print STDERR "No SC or KC present\n";
		exit(1);
	    }
	    $SC = $KC / length($seq{$N[1]});
	    $node{$N[1]} .= "\tSC:f:$SC";
	}

    } elsif (/^L\s+(\S+)\s+(.)\s+(\S+)\s+(.)/) {
	# ASSUMPTION: nodes before edges
	if (!exists($node{$1}) || !exists($node{$3})) {
	    print STDERR "Removing edge $_\n";
	    next;
	}
        $edge[$edge_num] = $_;
	$self_loop{$1} = 1 if $1 eq $3;
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

# Find simple cycles.  This isn't particularly fast on large graphs.
# TODO: look at Johnson's algorithm, or rewrite-in / callout to Python which
# has an implementation of this in networkx:
#
# my $code = "import networkx as nx\nG=nx.DiGraph()\n";
# $code .= "G.add_nodes_from([\"" . join("\",\"",@node_order) . "\"])\n";
# my @edge_str;
# foreach (@edge) {
#     m/^L\s+(\S+)\s+(.)\s+(\S+)\s+(.)/;
#     push(@edge_str, "(\"".$1."\",\"".$3."\")");
# }
# $code .= "G.add_edges_from([" . join(",", @edge_str) . "])\n";
# $code .= "print(list(nx.simple_cycles(G)))\n";
# my $simple_cycles = `python -c '$code'`;
# $simple_cycles =~ tr/[],'/ /;
# $simple_cycles =~ s/^\s+//;
# foreach (split(/\s+/,$simple_cycles)) {
#     $loop{$_}=1;
# }

sub find_all_cycles {
    my @all_cycles;

    for my $start (@node_order) {
        my @stack = ([$start, [$start]]);
        while (@stack) {
            my ($node, $path) = @{pop @stack};
	    my @next;

	    # Do we also need to do edge_in here, or take into account orientation?
	    foreach (@{$edge_out{$node}}) {
		$edge[$_] =~ m/^L\s+(\S+)\s+(.)\s+(\S+)\s+(.)/;
		if ($1 eq $node) {
		    push(@next, $3);
		} else {
		    push(@next, $1);
		}
	    }

            for my $neighbour (sort @next) {
                if ($neighbour eq $path->[0] && @$path > 1) {
                    push @all_cycles, [@$path];
                } elsif (!grep { $_ eq $neighbour } @$path[1..$#$path]) {
                    push @stack, [$neighbour, [@$path, $neighbour]];
                }
            }
        }
    }
    return @all_cycles;
}

my @all_cycles = find_all_cycles();
my %loop;
foreach my $cycle (@all_cycles) {
    foreach (@$cycle) {
	$loop{$_}=1;
    }
}
if ($verbose) {
    print "In loops:";
    foreach (sort keys %loop) {
	print " $_";
    }
    print "\n";
}


# An array of raw depths.  Filter very low figures (unused nodes), compute
# average, and then cycle again with a revised estimate of what low-depth
# means.
my $avg_depth = 0;
my $total;
my $tlen;

my $last_avg_depth;
my $ncycles = 0;
do {
    $last_avg_depth = $avg_depth;
    $avg_depth = 0;
    $total = 0;
    $tlen = 0;

    foreach my $n (sort keys %node) {
	my ($d) = $node{$n} =~ m/SC:f:(-?\d+(\.\d+)?)/;
	#if ($verbose) {
	#    print "Node $n\tdepth ",int(100*$d)/100,"\tlen ",length($seq{$n}),"\tin-loop ",exists($loop{$n}) ? 1 : 0, "\n";
	#}
	if ($d > $min_depth) {
	    # For initial average, exclude self loops
	    next if exists($self_loop{$n});
	    next if $loop{$n};

	    #if ($verbose) {
	    #	print "Using $n len ",length($seq{$n}), " depth $d\n";
	    #}

	    $total += $d*(length($seq{$n}));
	    $tlen += (length($seq{$n}));
	}
    }

    $avg_depth = $total / $tlen;
    if ($verbose) {
	print "Avg depth $avg_depth\n";
    }
    $min_depth = $avg_depth/$depth_div;
    $ncycles++;
} while ($avg_depth != $last_avg_depth && $ncycles < 10);



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
#    $_=`tag_gfa_copy_numbers.pl -v $base.gfa | grep "1-depth"`;
#    my @F = split(/\s+/, $_);
#    #print "Avg $F[-1]\n";
#    $avg_depth = $F[-1];
#}

# Alternative; try fitting depth to D, 1*D, 2*D, etc.
my $best_try=$avg_depth;
my $best_delta=1e10;
for (my $try=$avg_depth/1.3; $try<$avg_depth*1.5; $try+=0.1) {
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
    print "try $try, delta=$delta\n" if $verbose;
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
	my $copy;
	if ($float_copy) {
	    $copy = $d/$avg_depth+$offset;
	} else {
	    $copy = int($d/$avg_depth+0.5+$offset);
	}
	$copy = $max_copy if ($copy > $max_copy);
	push(@copy, $copy);
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
