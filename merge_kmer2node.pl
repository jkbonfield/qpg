#!/usr/bin/perl -w

# Loads multiple kmer2node outputs and merges them.
# This is a simple summation, but we may wish to use weighted sums
# or a union in some manner.

use strict;

my %node; # keyed on node name
my %edge; # keyed on node1,dir1,node2,dir2

# Parse GFA edges
my $gfa = shift(@ARGV);
open(my $fh, "<", $gfa) || die;
while (<$fh>) {
    if (/^L/) {
	my ($n1,$d1,$n2,$d2) = m/^L\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
	$edge{"$n1.$d1.$n2.$d2"} = 0;
	#print "Set edge $n1.$d1.$n2.$d2\n";
    }
}
close($fh);

foreach my $file (@ARGV) {
    open(my $fh, "<", $file) || die;
    while (<$fh>) {
	chomp();
	if (/Node/) {
	    my @F = split(/\s+/, $_);
	    if (exists($node{$F[1]})) {
		$node{$F[1]}[-1] = $F[-1] if ($F[-1] > $node{$F[1]}[-1]);
		#$node{$F[1]}[-1] += $F[-1];
	    } else {
		$node{$F[1]} = \@F;
	    }
	} else {
	    if (/Edge/) {
		my ($n1,$d1,$n2,$d2,$c)
		    = m/^Edge\s+(\S+)(\S)\s+(\S+)(\S)\s+(\d+)/;
		my $d1m = $d2 eq "+" ? "-" : ($d2 eq "-" ? "+" : "b");
		my $d2m = $d1 eq "+" ? "-" : ($d1 eq "-" ? "+" : "b");
		my $e1 = "$n1.$d1.$n2.$d2";
		my $e2 = "$n2.$d1m.$n1.$d2m";
		#print "Edge\t$e1\t$c\t",exists($edge{$e1}) ? 1 : 0,"\n";
		#print "Edge\t$e2\t$c\t",exists($edge{$e2}) ? 1 : 0,"\n";
		if (exists($edge{$e1})) {
		    # Basic fwd strand case
		    $edge{$e1} += $c;
		} elsif (exists($edge{$e2})) {
		    # Basic reverse-complement case
		    $edge{$e2} += $c;
		} elsif ($d2 eq "b" && $d1 ne "b") {
		    # Second element in both orientations but we can resolve.
		    # We distribute evenly if both are found.
		    my $p = exists($edge{"$n1.$d1.$n2.+"});
		    my $m = exists($edge{"$n1.$d1.$n2.-"});
		    $edge{"$n1.$d1.$n2.+"} += $c / ($p+$m) if $p;
		    $edge{"$n1.$d1.$n2.-"} += $c / ($p+$m) if $m;

		    $p = exists($edge{"$n2.$d1m.$n1.+"});
		    $m = exists($edge{"$n2.$d1m.$n1.-"});
		    $edge{"$n2.$d1m.$n2.+"} += $c / ($p+$m) if $p;
		    $edge{"$n2.$d1m.$n1.-"} += $c / ($p+$m) if $m;
		} elsif ($d1 eq "b" && $d2 ne "b") {
		    # First element in both orientations but we can resolve.
		    # We distribute evenly if both are found.
		    my $p = exists($edge{"$n1.+.$n2.$d2"});
		    my $m = exists($edge{"$n1.-.$n2.$d2"});
		    $edge{"$n1.+.$n2.$d2"} += $c / ($p+$m) if $p;
		    $edge{"$n1.-.$n2.$d2"} += $c / ($p+$m) if $m;

		    $p = exists($edge{"$n2.+.$n1.$d2m"});
		    $m = exists($edge{"$n2.-.$n1.$d2m"});
		    $edge{"$n2.+.$n2.$d2m"} += $c / ($p+$m) if $p;
		    $edge{"$n2.-.$n1.$d2m"} += $c / ($p+$m) if $m;
		} else {
		    # Both first and second elements existing in two dirs.
		    # Not implemented yet.
		}
	    }
	}
    }
    close($fh);
}

$"="\t";
foreach (sort keys %node) {
    print "@{$node{$_}}\n";
}

foreach (sort keys %edge) {
    next unless $edge{$_} > 0;
    m/(\S+)\.(\S+)\.(\S+)\.(\S+)/;
    print "Edge\t$1$2\t$3$4\t$edge{$_}\n";
}
