#!/usr/bin/perl -w

# Processes pangene.js BB output and splits output gfa
use strict;

my $min_size=50;

if ($#ARGV != 1) {
    print STDERR "Usage: k8 pangene.js call in.gfa | bubble_split.pl in.gfa out_prefix\n";
    exit 1;
}

my $gfa = shift(@ARGV);
my $prefix = shift(@ARGV);

# Parse GFA
print "Parsing GFA\n";
my $nseq = 0;
my $nlink = 0;

my %node;
my %edge;
open(GFA, "<", $gfa) || die;
while (<GFA>) {
    if (/^S\s+(\S+)/) {
	$node{$1} = $_;
	$nseq++;
    } elsif (/^L\s+(\S+)\s+(.)\s+(\S+)\s+(.)\s+(.*)/) {
	push(@{$edge{$1}},[$3, $_]);
	push(@{$edge{$3}},[$1, $_]);
	$nlink++;
    }
}
close(GFA);
print "Nseq $nseq, Nlink $nlink\n";

# Parse pangene output looking for bubbles of $min_size or more nodes
while (<>) {
    next unless /^BB/;
    chomp();
    my @F = split(/\s+/, $_);
    next unless $F[7] >= $min_size;
    $F[4]=~s/>//;
    $F[5]=~s/>//;
    my @nodes = split(",", $F[8]);

    print "Writing $prefix.$F[4].$F[5] with ",scalar(@nodes)+2, " nodes\n";
    open(BUBBLE, ">", "$prefix.$F[4].$F[5]");

    my %filter;

    $filter{$F[4]}=1;
    print BUBBLE $node{$F[4]};

    foreach (@nodes) {
	$filter{$_}=1;
	print BUBBLE $node{$_};
    }

    $filter{$F[5]}=1;
    print BUBBLE $node{$F[5]};


    foreach (@nodes) {
	foreach (@{$edge{$_}}) {
	    my ($node,$line)=@{$_};
	    next unless exists($filter{$node});
	    print BUBBLE $line;
	}
    }
    close(BUBBLE);
}
