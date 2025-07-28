#!/usr/bin/perl -w
use strict;

# Estimates the quality of a candidate solution by analysing the bam alignment stats
# Usage: estimate_qual.pl run/seq_base ...
# eg: estimate_qual.pl seq_1050-0003-#1#1

# Output is the estimated quality and the actual quality (eval_cons)

foreach my $base (@ARGV) {
    foreach my $fn (sort glob("$base.bam.*")) {
	my ($nmapped,$mapped,$supp);
	open(FH, "samtools flagstat $fn |") || die;
	while (<FH>) {
	    chomp($_);
	    my @F = split(/\s+/, $_);
	    $mapped  = $F[5] if /primary mapped/;
	    $nmapped = $F[0] if /primary mapped/;
	    $supp    = $F[0] if /supplementary/;
	}
	my $supp_perc = int($supp/($supp+$nmapped) * 100 * 1000 + 0.5)/1000;
	$mapped=~tr/%(//d;
	$fn =~ s/bam/eval_seq/;
	open(FH, "<$fn") || die "$fn";
	my ($cov,$used,$contigs,$breaks);
	while (<FH>) {
	    chomp($_);
	    next if /contig/;
	    s/%//g;
	    my @F = split("\t", $_);
	    $cov     = $F[3];
	    $used    = $F[4];
	    $contigs = $F[5];
	    $breaks  = $F[6];
	}
	my $f1 = int(2*$cov*$used/($cov+$used)*1000+.5)/1000;
	close(FH);
	print "$fn\t$mapped\t$supp_perc\t$cov\t$used\t$f1\t$contigs\t$breaks\n";
    }
}
