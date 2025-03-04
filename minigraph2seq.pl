#!/usr/bin/perl -w

use strict;

# Take a GFA and a full genome, minigraph to get the path, and turn that
# pathfinder compatible path
my $gfa = shift(@ARGV);
my $seq = shift(@ARGV);
open(FH, "minigraph -x lr -l 100 -d 100 --vc -m 100,100 -n 10,10 $gfa $seq 2>/dev/null |");
while (<FH>) {
    my @F = split("\t", $_);
    my $n=1;
    print "PATH node name [from minigraph]\n";
    foreach ($F[5]=~m/[<>][^<>]+/g) {
	my ($dir,$node) = ($_=~m/(.)(.*)/);
	printf("[ %3d] %3s%s\n", $n++, $node, $dir eq ">" ? "+" : "-");
    }
}
close(FH);

