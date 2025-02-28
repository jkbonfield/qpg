#!/usr/bin/perl -w

# Turn a GFA + GAF alignment path string ">node>node<node..." into a
# sequence by copying the elements out of the GFA and stitching them
# together.

# Usage: path2seq.pl in.gfa path

if (scalar @ARGV != 2) {
    print STDERR "Usage: path2seq.pl in.gfa path\n";
    exit 1;
}

# Parse GFA; minimally
open(my $gfa, "<", shift(@ARGV)) || die;
while (<$gfa>) {
    chomp($_);
    my @F = split(/\s+/, $_);
    next unless scalar(@F) && $F[0] eq "S"; # skip other fields for now
    $gfa{$F[1]}{seq} = uc($F[2]);
}


# Parse the path
while (<>) {
    chomp();
    my $seq="";

    foreach ($_=~m/[<>][^<>]*/g) {
	my ($dir,$node) = ($_=~m/(.)(.*)/);
	my $gseq = $gfa{$node}{seq};
	if ($dir eq "<") {
	    $gseq =~ tr/ACGT/TGCA/;
	    $gseq = reverse($gseq);
	}
	$seq .= $gseq;
    }
    
    print "$seq\n";
}
