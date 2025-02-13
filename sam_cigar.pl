#!/usr/bin/perl -w

# Process a SAM file POS + CIGAR to give alignment start..end for query / ref
# F[0] = name
# F[1] = flag
# F[3] = pos
# F[5] = cigar

use strict;

print "Contig_name          flag       qlen      qstart     qend       rstart     rend      delta  indelsz\n";

# Assumption: primaries come before secondaries
my $slen = 0;
while (<>) {
    chomp($_);
    next if (/^@/);
    my @F = split("\t", $_);
    my $rstart = $F[3];
    my $rend = $F[3];
    my $qstart = 0;
    $slen = length($F[9]) if (($F[1] & 0x900) == 0); # primary
    my $qend = $slen;

    # Find query extents
    $qstart += $1 if ($F[5] =~ m/^(\d+)[SH]/);
    $qend   -= $1 if ($F[5] =~ m/(\d+)[SH]$/);
    $F[5] =~ s/(\d+)S//g;

    # Find ref extents.  This needs full CIGAR processing
    my $isz=0;
    foreach ($F[5] =~ m/\d+./g) {
	m/(\d+)(.)/;
	if ($2 eq "M" || $2 eq "D") {
	    $rend += $1;
	}
	if ($2 eq "I" || $2 eq "D") {
	    $isz += $1;
	}
    }

    my $delta = ($rend-$rstart) - ($qend-$qstart);
    printf("%-20s %4d   %8d    %8d %8d     %8d %8d   %8d %8d\n",
	   $F[0], $F[1], $slen, $qstart, $qend, $rstart, $rend, $delta, $isz);
    #print "$F[0]\t$F[1]\t$slen\t$qstart\t$qend\t$rstart\t$rend\n";
}
