#!/usr/bin/perl -w

use strict;

# Take a GFA and a shredded fasta.  Use minigraph to annotate the GFA.
# Then process that to get node coverage data.
# a better annotation for Bandage and Pathfinder.

my $gfa = shift(@ARGV);
my $seq = shift(@ARGV);
my $expected_depth = shift(@ARGV);

# Long read
#open(FH, "minigraph -x lr -j 0.01 -w 15 -l 100 -d 100 $gfa $seq -m 50,40 -n 5,5 --cov 2>/dev/null |");

# Short read
open(FH, "minigraph -x sr --cov  $gfa $seq 2>/dev/null |");

# Output is a GFA
while (<FH>) {
    next unless /^S/;
    chomp($_);

    m/dc:f:([0-9.]*)/;
    print int($1/$expected_depth + .8), " ";
    #printf("%3d ", int($1+.5));
}
print "\n";
close(FH);
