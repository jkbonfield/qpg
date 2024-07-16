#!/usr/bin/perl -w
use strict;

# Shreds a FASTA file into a series of small reads, randomly subsampled
# many times over in either orientation.
# Errors may be introduced too.

# The purpose of this is to simulate short-read sequencing (albeit
# single ended) of a sample so we can compare kmer mapping of
# sequenced fragments vs the one complete original fasta entry and see
# whether the kmer index has sufficient power to reconstruct the path.

# eg ./shred.pl -s 0 -d 100 < s1.fa | ./kmer2node train8.fa.nodeseq - 1

use Getopt::Long;

my $length = 100;
my $err = 0;  # specify as a fraction
my $depth = 50;
my $seed = -1;
my $kmer = 14; # corrects depth for missing matches at fragment ends
my $circular = 0;
GetOptions ("length=i" => \$length,
	    "err=f"    => \$err,
	    "depth=i"  => \$depth,
            "seed=i"   => \$seed,
            "kmer=i"   => \$kmer,
            "circular" => \$circular)
or die("Usage: shred [-c] [-s seed] [-l length] [-e error_fraction] [-d depth] in.fa > out.fa\n");

srand($seed) if $seed>=0;


# Parse FASTA
my $name = "";
my $seq = "";
while (<>) {
    chomp($_);
    if (/^>/) {
        if ($seq ne "") {
	    shred($name, uc($seq));
            $seq = "";
        }
        s/^>//;
        s/\s.*//;
        $name = $_;
    } else {
        $seq .= $_;
    }
}
shred($name, uc($seq));

# Shred a single fasta entry
sub shred {
    my ($name, $seq) = @_;

    my $slen = length($seq);
    $seq .= $seq; # permits circularisation
    die if $slen < $length;

    my $rseq = reverse $seq;
    $rseq =~ tr/ACGT/TGCA/;

    # permit overlapping start/end with short reads
    my $mod_slen = $circular ? $slen : $slen+$length;
    my $nreads = int($mod_slen / ($length-$kmer+1) * $depth);

    print STDERR "shred $name, length $slen => $nreads reads\n";

    for (my $i = 0; $i < $nreads; $i++) {
	my $pos = $circular
	    ? int(rand($slen))
	    : int(rand($slen+$length)) - $length;

	# Truncate pos and/or seq for end overlaps
	my $len2 = $length;
	if ($pos < 0) {
	    $len2 += $pos;
	    $pos = 0;
	}
	if (!$circular && $pos+$length >= $slen) {
	    $len2 -= $pos+$length-$slen;
	}
	next if $len2 <= 0;

	# Produce a subsequence and r/c at random
	my $sub;
	if (rand()>=0.5) {
	    $sub = substr($seq, $pos, $len2);
	    $pos++;
	} else {
	    $sub = substr($rseq, $pos, $len2);
	    if (!$circular) {
		$pos = $slen - ($pos+$len2-1);
	    }
	}

	my $nerr = $err * $len2;
	my $ierr = int($nerr);
	$ierr++ if (rand() < $nerr-$ierr);

	for (my $e=0; $e<$ierr; $e++) {
	    my $p = int(rand($len2));
	    my $base = substr($sub, $p, 1);
	    my @bases = qw/A C G T/;
	    while ($base eq substr($sub, $p, 1)) {
		$base = $bases[int(rand(4))];
	    }
	    substr($sub, $p, 1) = $base;
	}
	print ">$name#$pos#$i\n$sub\n";
    }
}

