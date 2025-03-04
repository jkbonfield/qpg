#!/usr/bin/perl -w

use strict;

# Take a GFA and a shredded fasta.  Perform alignments and produce an annotated
# GFA.  Most of this is down to minigraph, with some modifications to support
# a better annotation for Bandage and Pathfinder.

# On seq nodes we create KC (kmer count) and dp (depth) stats.
# On edges we add EC (edge count).
my $gfa = shift(@ARGV);
my $seq = shift(@ARGV);
# Long read
#open(FH, "minigraph -x lr -j 0.01 -w 15 -l 100 -d 100 $gfa $seq -m 50,40 -n 5,5 --cov 2>/dev/null |");

# Short read
open(FH, "minigraph -x sr --cov  $gfa $seq 2>/dev/null |");

# Output is a GFA
while (<FH>) {
    chomp($_);

    m/dc:f:([0-9.]*)/;
    my $dc = $1;
    if (/^S/) {
	my @F = split("\t",$_);
	print "$_\tKC:i:",int($dc*length($F[2])+.5),"\tDP:i:",int($dc+.9),"\n";
    } else {
	print "$_\tEC:i:",int($dc+.9),"\n";
    }
}
close(FH);

# cov1/cov2/nbreak on tangle_sim_0000?/seq_004?#1#1
# -x lr -l 100 -d 100 $gfa $seq -m 100,100 -n 10,10   87.2744 83.2393 6.73333
# -m 50,40                                          + 87.5936 83.7799 6.84444
# -m 150,140                                          87.0502 83.9876 6.85556
# -n 3,3                                              87.1736 83.6458 6.81111
# -n 5,5                                            + 87.4017 83.4251 6.85556
# -n 7,7                                              87.1801 83.4192 6.85556
# -n 20,20                                            86.1912 82.4994 7.24444
# -l 1000 -d 1000                                     87.2744 83.2393 6.73333
# -l 10   -d 10                                       87.2744 83.2393 6.73333
# -q 40                                               87.2744 83.2393 6.73333
# -r 100,1000                                       + 87.8049 83.5514 6.84444
# -r 1000,40000                                       87.2744 83.2393 6.73333
# -j 0.2                                              84.0150 71.9070 9.74444
# -j 0.05                                             88.9823 85.0627 6.47778
# -j 0.02                                           + 88.9400 85.7458 6.47778
# -U 25,125                                           87.2744 83.2393 6.73333
# -U 100,500                                          87.2744 83.2393 6.73333
# -w 10                                               87.4530 82.9884 6.95556
# -w 12                                               87.9290 84.1211 6.62222
# -w 14                                             + 87.4967 84.4070 6.54444
# -w 16                                               87.5278 83.6670 6.67778
# -k 16                                             + 88.1528 83.4111 6.90000
# -k 18                                               87.0453 82.6343 6.78889

# -w 14 -k 16                                         87.7236 83.5770 7.00000
# -w 14 -m 50,40 -n 5,5                               88.0123 84.4770 6.53333
# -w 14 -m 50,40 -n 5,5 -j 0.02                       89.4712 86.7221 6.52222
# -w 14 -m 50,40 -n 5,5 -j 0.02 -k16                  89.5540 86.1299 6.63333
# -w 13 -m 50,40 -n 5,5 -j 0.02                       89.4757 86.7008 6.51111
# -w 15 -m 50,40 -n 5,5 -j 0.02                       89.5312 86.7742 6.43333
# -w 15 -m 30,20 -n 3,3 -j 0.02                       89.2870 86.6721 6.44444
# -w 15 -m 50,40 -n 5,5 -j 0.01                     * 89.7382 87.0002 6.50000
# " -x sr  	    	   			      86.6428 80.6556 8.35556


