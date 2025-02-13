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
my %seq;      # node sequence
my %node;     # node GFA line (S)
my @edge;     # edge GFA line (L)
my $edge_num=0;
my %edge_in;  # index into @edge above
my %edge_out; # index into @edge above
while (<$gfa>) {
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
close($gfa);


#-----------------------------------------------------------------------------
# Report graph sequences per node, appending incoming edges (maximum 1 node away)

# Get the previous node sequences and return them as an array.
# This is recursive until $plen bases have been copied or we hit the end of
# the previous links.
sub prev_seq {
    my ($n, $plen) = @_;
    
    #print ">>> PREV SEQ for $n $plen\n";

    my @pseq = ();

    # Incoming for + dir
    if (exists($edge_in{$n})) {
	foreach my $e (@{$edge_in{$n}}) {
	    my @E=split("\t", $edge[$e]);
	    next unless $E[4] eq "+";
	    my @N2 = split("\t", $node{$E[1]});
	    my $pseq = $N2[2];
	    if ($E[2] eq "-") {
		$pseq =~ tr/ACGT/TGCA/;
		$pseq = reverse($pseq);
	    }

	    if (length($pseq) < $plen) {
		my $p2len = $plen-length($pseq);
		my @p = prev_seq($E[1], $p2len);
		#print "Too short $plen => @p\n";
		if (scalar(@p)) {
		    foreach (@p) {
			push(@pseq, substr($_, -$p2len) . $pseq);
		    }
		} else {
		    #print "No previous edges\n";
		    push(@pseq, substr($pseq,-$plen));
		}
	    } else {
		#print "Long enough: ",length($pseq),">=$plen\n";
		push(@pseq, substr($pseq,-$plen));
	    }
	}
    }

    # Outgoing for - dir
    if (exists($edge_out{$n})) {
	foreach my $e (@{$edge_out{$n}}) {
	    my @E=split("\t", $edge[$e]);
	    next unless $E[2] eq "-";
	    my @N2 = split("\t", $node{$E[3]});
	    my $pseq = $N2[2];
	    if ($E[4] eq "+") {
		$pseq =~ tr/ACGT/TGCA/;
		$pseq = reverse($pseq);
	    }
	    push(@pseq, substr($pseq,-$plen));
	}
    }

    #print "<<< return @pseq\n";
    return @pseq;
}

foreach my $n (sort keys %node) {
    my $nseq = $seq{$n};
    print "\@$n\n$nseq\n";

    foreach my $prev (prev_seq($n, $kmer-1)) {
	print $prev . $nseq, "\n";
    }
}

