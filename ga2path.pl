#!/usr/bin/perl -w

# Process a GraphAligner GAF file to build the most common path
# Usage: ga2path.pl input.gaf

use strict;

# Load GFA paths
my @paths; # total list of paths
my %nodes; # node list; redundant
my %in;    # list of incoming nodes per node
my %out;   # list of outgoing nodes per node
while (<>) {
    chomp($_);
    my @F=split("\t", $_);

    push(@paths, $F[5]);
    foreach ($F[5] =~ m/[<>]([^<>]+)/g) {
	$nodes{$_}++;
    }
}

loop:
my $loop = 0;

print "=== NODES ===\n";
foreach (sort keys %nodes) {
    print "$_ $nodes{$_}\n";
}
print "\n";

# FIXME: Basically we can do sequence assembly on the node structures.


# Eg     ,-->u15->.
#       /          `v
# u3-->u13-->u16-->u17-->u7
#       ^           /
#        `<--u14<--'
#
#           vvv             vvv
# So u3 u13 U15 u17 u14 u13 U16 u17 u7
# or u3 u13 U16 u17 u14 u13 U15 u17 u7
#           ^^^             ^^^
#
# Look for paths u3 u13 u16, u3 u13 u15, u14 u13 u13 and u14 u13 u16.
# We disambiguate the circuit then.
#

# 1. Identify triplets (IN, NODE, OUT)
# 2. If <NODE then swap IN/OUT around so order aware
# 3. Skip if *,NODE,* is just 1 route. (Ie A->B->C chain)
# 4. Foreach IN->NODE
#        Replace NODE:  IN->NODE->{OUT,...} with IN->IN:NODE->{OUT}
#        where {OUT} is a reduced OUT list based only those from IN.
# ie.
#
# I     P        I--I:A--P
#  \_A_/    =>               (or I-I:A-Q,  J-J:A-P)
#  /   \
# J     Q        J--J:A--Q

# Foreach GFA node, get the predecessor.
my %prefix = ();
foreach (@paths) {
    my $last1 = "";
    my $last2 = "";
    # FIXME: match things in @nodes only
    foreach (m/[<>][^<>]+/g) {
	my ($dir,$node) = m/(.)(.*)/;
	# last1 is the middle entry.
#	if ($last1 ne "") {
#	    my $t = $_;
#	    my $s = $last1;
#	    $t=~tr/<>//d;
#	    $s=~tr/<>//d;
#	    if ($_ =~ /</) {
#		my $x = $t;
#		$t = $s;
#		$s = $x;
#	    }
#	    $in{$t}{$s}++;
#	    $out{$s}{$t}++;
#	}
	if ($last2 ne "") {
	    my $prev = "";
	    my $next = "";
	    if ($last1 =~ />/) {
		$prev = $last2;
		$next = $_;
	    } else {
		$prev = $_;
		$prev =~ tr/<>/></;
		$next = $last2;
		$next =~ tr/<>/></;
	    }
	    my $p = $last1;
	    $p =~ tr/<>//d;
	    $prefix{$p}{$prev}++ if ($prev ne "");
	    if ($prev ne "") {print "PREV: $prev $last1 $next\n"}
	}
	$last2 = $last1;
	$last1 = $_;
    }    
}

#print "=== EDGES ===\n";
#foreach (sort keys(%in)) {
#    foreach my $e (sort keys(%{$in{$_}})) {
#	print "IN $e $_\n";
#    }
#    foreach my $e (sort keys(%{$out{$_}})) {
#	print "OUT $_ $e\n";
#    }
#}
#print "\n";

foreach (sort keys(%prefix)) {
    my @a = sort {$prefix{$_}{$a} <=> $prefix{$_}{$b}} (keys(%{$prefix{$_}}));
    if (scalar(@a) == 1) {
	# Collapse neighbouring nodes together and remove the old one
	$nodes{"$_$a[0]"}++;
	delete $nodes{$_};
	$loop = 1;
    }
    print ":", scalar(@a),"\n";
    foreach my $sub (@a) {
	print "$_ $sub = $prefix{$_}{$sub}\n";
    }
    #print "$_ $a[0]/$prefix{$_}{$a[0]} $a[-1]/$prefix{$_}{$a[-1]}\n";
}

goto loop if ($loop);
