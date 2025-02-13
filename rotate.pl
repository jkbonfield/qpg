#!/usr/bin/perl
# Reads a FASTA file from stdin and rotates each sequence by a random amount.

use strict;

if (scalar(@ARGV)) {
    srand(shift(@ARGV));
}

sub rotate {
    my ($seq) = @_;
    my $len = length($seq);
    my $pos = int(rand($len));
    return substr($seq,$pos) . substr($seq,0,$pos);
}

my $seq="";
my $name="";
while (<>) {
    chomp();
    if (/^>/) {
	print "$name\n", rotate($seq), "\n" if ($seq);
	$name=$_;
	$seq="";
    } else {
	$seq .= $_;
    }
}
print "$name\n", rotate($seq), "\n" if ($seq);
