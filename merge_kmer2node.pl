#!/usr/bin/perl -w

# Loads multiple kmer2node outputs and merges them.
# This is a simple summation, but we may wish to use weighted sums
# or a union in some manner.

use strict;

my %node; # keyed on node name

foreach my $file (@ARGV) {
    open(my $fh, "<", $file) || die;
    while (<$fh>) {
	chomp();
	my @F = split(/\s+/, $_);
	if (exists($node{$F[1]})) {
	    $node{$F[1]}[-1] = $F[-1] if ($F[-1] > $node{$F[1]}[-1]);
	    #$node{$F[1]}[-1] += $F[-1];
	} else {
	    $node{$F[1]} = \@F;
	}
    }
    close($fh);
}

$"="\t";
foreach (sort keys %node) {
    print "@{$node{$_}}\n";
}
