#!/usr/bin/perl -w

# Load gaf2nodeseq but minus a graph alignment.
# Instead we use the GFA alone, but we augment this with incoming edges to work
# out the possible kmers we may wish to match against.

use strict;

# Usage: gfa2nodeseq.pl graph.gfa [kmer]

open(my $gfa,   "<", shift(@ARGV)) || die;
my $kmer = shift(@ARGV);
$kmer = 14 unless defined($kmer);


# Parse GFA
my %node;
my @edge;
my $edge_num=0;
my %edge_in;  # index into @edge above
my %edge_out; # index into @edge above
while (<$gfa>) {
    chomp();
    if (/^S\s+(\S+)/) {
        $node{$1} = $_;
    } elsif (/^L\s+(\S+)\s+(.)\s+(\S+)\s+(.)\s+(.*)/) {
	$edge[$edge_num] = $_;
        push(@{$edge_out{$1}}, $edge_num);
        push(@{$edge_in{$3}},  $edge_num);
	$edge_num++;
    }
}
close($gfa);


#-----------------------------------------------------------------------------
# Report graph sequences per node, appending incoming edges (maximum 1 node away)
foreach my $n (sort keys %node) {
    my @N = split("\t", $node{$n});
    my $nseq = $N[2];
    print "\@$n\n$nseq\n";

    # Incoming for + dir
    foreach my $e (@{$edge_in{$n}}) {
	my @E=split("\t", $edge[$e]);
	next unless $E[4] eq "+";
	my @N2 = split("\t", $node{$E[1]});
	my $pseq = $N2[2];
	if ($E[2] eq "-") {
	    $pseq =~ tr/ACGT/TGCA/;
	    $pseq = reverse($pseq);
	}
	print substr($pseq,-($kmer-1)) . $nseq, "\n";
    }

    # Outgoing for - dir
    foreach my $e (@{$edge_out{$n}}) {
	my @E=split("\t", $edge[$e]);
	next unless $E[2] eq "-";
	my @N2 = split("\t", $node{$E[3]});
	my $pseq = $N2[2];
	if ($E[4] eq "+") {
	    $pseq =~ tr/ACGT/TGCA/;
	    $pseq = reverse($pseq);
	}
	print substr($pseq,-($kmer-1)) . $nseq, "\n";
    }
}

