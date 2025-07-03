#!/usr/bin/perl -w

# Parses the output of GraphAligner and a GFA file to incorporate node weights.
# Bandage uses LN for length and KC for kmer coverage, with the ratio being
# the total base coverage SC (KC/LN).

use strict;

# Stdin = GraphAligner
# 1st argument = GFA file
# Stdout = amended GFA
my $gfa = shift(@ARGV);

# Parse GFA
my %seq;      # node sequence
my %node;     # node GFA line (S)
my @edge;     # edge GFA line (L)
my $edge_num=0;
my %edge_in;  # index into @edge above
my %edge_out; # index into @edge above
open(my $fh, "<", $gfa) || die;
local $"="\t";
while (<$fh>) {
    chomp();
    if (/^S\s+(\S+)/) {
        $node{$1} = $_;
        my @N = split("\t", $_);
        $seq{$N[1]} = $N[2];
    } elsif (/^L\s+(\S+)\s+(.)\s+(\S+)\s+(.)/) {
        $edge[$edge_num] = $_;
        push(@{$edge_out{$1}}, $edge_num);
        push(@{$edge_in{$3}},  $edge_num);
        $edge_num++;
    }
}
close($fh);

# Parse GraphAligner output.  Uses:
# F[0] for name to identify secondary alignments
# F[5] for the path ([<>]node)*
# F[10] for alignment length
# F[11] for alignment mapping quality
# F[2] and F[3] for alignment start/end
my %coverage;
my %EC;
my @lines = ();
my %ids;
while (<>) {
    chomp($_);
    push(@lines,$_);
    my @F = split(/\s+/, $_);
    $ids{$F[0]}++;
}

for (my $i=0; $i <= $#lines; $i++) {
    $_ = $lines[$i];
    my @F = split(/\s+/, $_);
    if ($ids{$F[0]} > 1) {
	# Pick best
	my $best = 0;
	my $best_line = "";
	#my $best_mq = 0;
	#my $dup = 0;
	my $j;
	for ($j = 0; $j < $ids{$F[0]}; $j++) {
	    my @F = split(/\s+/,$lines[$i+$j]);
	    if ($F[10] >= $best && $F[11] > 0) {
		#if ($F[10] > $best || $best_mq < $F[11]) {
		#    $best_mq = $F[11];
		#    $best_line = $lines[$i+$j];
		#    $best = $F[10];
		#}
		#$dup = ($F[10] > $best) ? 0 : 1;
		$best = $F[10];
		$best_line = $lines[$i+$j];
	    }
	}
	$i = $i+$j-1;
	#if ($best > 0 && !$dup) { # doesn't help
	if ($best > 0) {
	    @F = split(/\s+/, $best_line);
	} else {
	    next;
	}
    }
    my @nodes = ();
    my @dirs = (); # FIXME: unused currently
    foreach ($F[5]=~m/[<>][^<>]*/g) {
	m/(.)(.*)/;
	push (@nodes, $2);
	push (@dirs, $1);
    }
    my $len_used = 0;
    # Internal nodes, excluding the first and last
    for (my $i = 1; $i < $#nodes; $i++) {
	my $node_len = length($seq{$nodes[$i]});
	my $d2 = $nodes[$i-1] eq ">" ? "+" : "-";
	my $d1 = $nodes[$i  ] eq ">" ? "+" : "-";
	$EC{"$nodes[$i-1]$d2.$nodes[$i]$d1"}++;
	$coverage{$nodes[$i]} += $node_len;
	$len_used += $node_len;
    }
    if ($#nodes > 0) {
	my $d2 = $dirs[-2] eq ">" ? "+" : "-";
	my $d1 = $dirs[-1] eq ">" ? "+" : "-";
	$EC{"$nodes[-2]$d2.$nodes[-1]$d1"}++;
    }

    # Remaining length to distribute between first and last node
    my $align_len = $F[3]-$F[2]+1;
    my $len_used1 = ($align_len - $len_used)/2;
    $len_used1 = 1 if ($len_used1 < 1);

    my $node_len = length($seq{$nodes[0]});
    $len_used1 = $node_len if ($len_used1 > $node_len);

    $coverage{$nodes[0]} += $len_used1;
    if ($#nodes > 0) {
	my $len_used2 = $align_len - $len_used - $len_used1;
	my $node_len = length($seq{$nodes[-1]});
	$len_used2 = $node_len if ($len_used2 > $node_len);
	$coverage{$nodes[-1]} += $len_used2;
    }

#    # This sounds like a better strategy, with proportional distribution,
#    # but for some reason it doesn't work as well as the above (pathfinder).
#    my $align_len = $F[3]-$F[2]+1;
#    my $remaining = $align_len - $len_used;
#    my $node_len1 = length($seq{$nodes[0]});
#    my $node_len2 = $#nodes > 0 ? length($seq{$nodes[-1]}) : 0;
#
#    my $len_used1 = $remaining * $node_len1 / ($node_len1 + $node_len2);
#    my $len_used2 = $remaining * $node_len2 / ($node_len1 + $node_len2);
#    $len_used1 = 1 if ($len_used1 < 1);
#    $len_used2 = 1 if ($len_used2 < 1);
#    $len_used1 = $node_len1 if ($len_used1 > $node_len1);
#    $len_used2 = $node_len2 if ($len_used2 > $node_len2);
#
#    $coverage{$nodes[0]} += $len_used1;
#    if ($#nodes > 0) {
#    	$coverage{$nodes[-1]} += $len_used2;
#    }
}


# Parse GFA and mark it up
open($fh, "<", $gfa) || die;
local $"="\t";
while (<$fh>) {
    if (/^S\s/) {
	chomp($_);
	my @F_ = split(/\s+/, $_);
	my @F = @F_[0..2];
	push(@F, "LN:i:".length($F[2]));
	my $cov = exists($coverage{$F[1]}) ? int($coverage{$F[1]}+.99) : 0;
	#$cov = $cov * 30 / length($F[2]);
	push(@F, "KC:i:".$cov);
	push(@F, "SC:f:".$cov / length($F[2]));
	print "@F\n";
    } elsif (/^L\s/) {
	chomp($_);
	my @F = split(/\s+/, $_);
	my $e = "$F[1]$F[2].$F[3]$F[4]";

	my $ec = exists($EC{$e}) ? int($EC{$e}+.99) : 0;
	$F[2]=~tr/+-/-+/;
	$F[4]=~tr/+-/-+/;
	my $f = "$F[3]$F[4].$F[1]$F[2]";
	$ec += exists($EC{$f}) ? int($EC{$f}+.99) : 0;

	print "$_\tEC:i:$ec\n";
    } else {
	print;
    }
}
close($fh);
