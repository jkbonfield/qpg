#!/usr/bin/perl -w

# Process a GraphAligner GAF file to build the most common path
# Usage: ga2path.pl input.gaf

use strict;

# Load GFA paths
my @paths; # total list of paths
my %nodes; # node list
while (<>) {
    chomp($_);
    my @F=split("\t", $_);

    push(@paths, $F[5]);
    foreach ($F[5] =~ m/[<>]([^<>]+)/g) {
	$nodes{$_}++;
    }
}

print "=== NODES ===\n";
foreach (sort keys %nodes) {
    print "$_ $nodes{$_}\n";
}
print "\n";

# Foreach GFA node, get the predecessor.
my %prefix = ();
foreach (@paths) {
    my $last1 = "";
    my $last2 = "";
    foreach (m/[<>][^<>]+/g) {
	my ($dir,$node) = m/(.)(.*)/;
	# last1 is the middle entry.
	if ($last2 ne "") {
	    my $prev = "";
	    if ($last1 =~ />/) {
		$prev = $last2;
	    } else {
		$prev = $_;
		$prev =~ tr/<>/></;
	    }
	    my $p = $last1;
	    $p =~ tr/<>//d;
	    $prefix{$p}{$prev}++ if ($prev ne "");
	    #if ($prev ne "") {print "PREV: $last1 $prev\n"}
	}
	$last2 = $last1;
	$last1 = $_;
    }    
}

foreach (sort keys(%prefix)) {
    my @a = sort {$prefix{$_}{$a} <=> $prefix{$_}{$b}} (keys(%{$prefix{$_}}));
    print ":", scalar(@a),"\n";
    foreach my $sub (@a) {
	print "$_ $sub = $prefix{$_}{$sub}\n";
    }
    #print "$_ $a[0]/$prefix{$_}{$a[0]} $a[-1]/$prefix{$_}{$a[-1]}\n";
}
