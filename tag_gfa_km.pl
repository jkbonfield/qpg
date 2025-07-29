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
my %EC;
while (<>) {
    chomp($_);
    my @F = split(/\s+/, $_);
    if (/^Node/) {
       $coverage{$F[1]} = $F[-1];
    } elsif (/^Edge/) {
	m/^Edge\s+(\S+)\s+(\S+)\s+(\d+)/;
	$EC{"$1.$2"} += $3;
	#print STDERR "EC{$1.$2} += $3\n";
    }
}

# Parse GFA and mark it up
open(my $fh, "<", $gfa) || die;
local $"="\t";
while (<$fh>) {
    if (/^S\s/) {
	chomp($_);
	my @F_ = split(/\s+/, $_);
	my @F = @F_[0..2];
	push(@F, "LN:i:".length($F[2]));
	push(@F, "KC:i:".int($coverage{$F[1]} * length($F[2])));
	push(@F, "SC:f:".$coverage{$F[1]});
	print "@F\n";
    } elsif (/^L\s/) {
	chomp($_);
	my @F = split(/\s+/, $_);
	my $e = "$F[1]$F[2].$F[3]$F[4]";
	my $ec = exists($EC{$e}) ? int($EC{$e}+.99) : 0;
	print "$_\tEC:i:$ec\n";
    } else {
	print;
    }
}
close($fh);
