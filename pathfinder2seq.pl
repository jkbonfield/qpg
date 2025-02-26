#!/usr/bin/perl -w

# Turn a GFA + pathfinder PATH outputs to sequences by copying the
# elements out of the GFA and stitching them together.

# Usage: pathfinder foo.gfa | pathfinder2seq.pl foo.gfa > foo.consensus.fa

# Parse GFA; minimally
open(my $gfa, "<", shift(@ARGV)) || die;
while (<$gfa>) {
    chomp($_);
    my @F = split(/\s+/, $_);
    next unless scalar(@F) && $F[0] eq "S"; # skip other fields for now
    $gfa{$F[1]}{seq} = uc($F[2]);
}


# Parse the path
my @path = ();
my $in_path=0;
my $contig=0;
my $seq = "";
while (<>) {
    if (!/^\[/) {
	if (/^PATH/) {
	    print ">contig_$contig\n$seq\n" if ($seq);
	    $in_path = 1;
	    $contig++;
	    $seq = "";
	} else {
	    $in_path = 0;
	}
	next;
    }

    next unless $in_path;
    chomp();
    my ($node,$dir) = ($_=~m/(\S+)([-+])$/);
    push(@path, [$dir, $node]);
}

print "@path\n";

$seq = "";
foreach my $p (@path) {
    my ($dir,$node) = @{$p};
    my $gseq = $gfa{$node}{seq};
    if ($dir eq "-") {
	$gseq =~ tr/ACGT/TGCA/;
	$gseq = reverse($gseq);
    }
    $seq .= $gseq;
}

print ">contig_$contig\n$seq\n" if ($seq);
