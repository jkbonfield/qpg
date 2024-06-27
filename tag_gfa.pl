#!/usr/bin/perl -w

# Parses the output of kmer2node and a GFA file to incorporate node weights.
# Bandage uses LN for length and KC for kmer coverage, with the ratio being
# the total base coverage SC (KC/LN).

use strict;

# Stdin = kmer2node output
# 1st argument = GFA file
# Stdout = amended GFA
my $gfa = shift(@ARGV);

# Parse kmer2node output
my %coverage;
while (<>) {
    next unless /^Node/;
    chomp($_);
    my @F = split(/\s+/, $_);
    $coverage{$F[1]} = $F[-1];
}

# Parse GFA and mark it up
open(my $fh, "<", $gfa) || die;
local $"="\t";
while (<$fh>) {
    if (/^S\s/) {
	chomp($_);
	my @F = split(/\s+/, $_);
	push(@F, "LN:i:".length($F[2]));
	push(@F, "KC:i:".int($coverage{$F[1]} * length($F[2])));
	push(@F, "SC:f:".$coverage{$F[1]});
	print "@F\n";
    } else {
	print;
    }
}
close($fh);
