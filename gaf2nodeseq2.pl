#!/usr/bin/perl -w

# Applies a GAF (graph alignment file) aligned against a GFA to identify
# the location that each sequence base aligns against.

# The GFA contains sequence nodes (name + SEQ)
# The GAF contains paths (>cons1>insA>cons2 etc), plus CIGAR strings
# The FASTA contains the sequence.

# We need to marry up all 3 so we can attribute actual FASTA sequence content
# to the node they were aligned against.  This means we can hash the
# aligned data rather than the node data, producing a better mapping strategy.

use strict;

# Usage: gaf2gfa.pl align.gaf in.fasta graph.gfa [kmer]

open(my $gaf,   "<", shift(@ARGV)) || die;
open(my $fasta, "<", shift(@ARGV)) || die;
open(my $gfa,   "<", shift(@ARGV)) || die;
my $kmer = shift(@ARGV);
$kmer = 14 unless defined($kmer);

# Data structures, indexed by name
# seq{}{seq}    sequence
# seq{}{path}   >node path
# seq{}{start}  first base coord in alignment
# seq{}{end}    last base coord in alignment
# seq{}{pstart} first base in path
# seq{}{pend}   last base in path
# seq{}{cigar}  CIGAR string
my %seq;

# GFA indexed by node name
# gfa{}{seq}    node seq content
# gfa{}{edge}   [node dir node dir cigar]
# gfa{}{query}  [list of query seqs mapped to node seq]
# gfa{}{=}      [list of match counts per node seq base]
# gfa{}{X}      [list of mismatch counts per node seq base]
my %gfa;

#-----------------------------------------------------------------------------
# Load data

# Parse FASTA
my $name = "";
my $seq = "";
while (<$fasta>) {
    chomp($_);
    s/\r//;
    if (/^>/) {
	if ($seq ne "") {
	    $seq{$name}{seq} = uc($seq);
	    $seq = "";
	}
	s/^>//;
	s/\s.*//;
	$name = $_;
    } else {
	$seq .= $_;
    }
}
$seq{$name}{seq} = uc($seq);

# Parse gaf
my $rep=0;
while (<$gaf>) {
    chomp($_);
    my @F = split("\t", $_);
    $F[0]=~s/\s.*//;
    if (exists($seq{$F[0]}{start})) {
	next if ($F[11]<60 or $F[10]<100);
	my $old_name = $F[0];
	$F[0]="seq_$rep";
	$rep++;
	print STDERR "Dup $old_name $F[0] mapq $F[11]\n";
	$seq{$F[0]}{seq} = $seq{$old_name}{seq};
    }
    #next if $F[11]<60;   nodeseq-2
    #next if $F[10]<100;  nodeseq-2
    $seq{$F[0]}{start} = $F[2];
    $seq{$F[0]}{end} = $F[3];
    $seq{$F[0]}{path} = $F[5];
    $seq{$F[0]}{pstart} = $F[7];
    $seq{$F[0]}{pend} = $F[8];
    m/.*cg:Z:([^\t]*)/;
    $seq{$F[0]}{cigar} = $1;
}

# Parse GFA; minimally
while (<$gfa>) {
    chomp($_);
    my @F = split(/\s+/, $_);
    next unless scalar(@F) && $F[0] eq "S"; # skip other fields for now
    $gfa{$F[1]}{seq} = uc($F[2]);
}


#-----------------------------------------------------------------------------
# Marry up seq, path and cigar
#print STDERR keys(%seq),"\n";
foreach my $name (sort keys %seq) {
    #print "Processing $name\n";
    my $seq = $seq{$name}{seq};
    my $cig = $seq{$name}{cigar};
    my $qstart = $seq{$name}{start};
    my $pstart = $seq{$name}{pstart};

    $qstart = $seq{$name}{start};
    $pstart = $seq{$name}{pstart};

    my $op;
    my $oplen = 0;

    my $ppos = $seq{$name}{pstart};

    my $last_node = "";

    while ($seq{$name}{path} =~ m/>([^>]*)/g) {
	my $seq_start = $qstart;
	my $node = $1;
#	my $seqn = $gfa{$1}{seq};
#	my $lenn = length($seqn)-$pstart;

	my $seqn = substr($gfa{$1}{seq}, $pstart);
	my $lenn = length($seqn);
	if ($ppos + $lenn > $seq{$name}{pend}) {
	    # Truncated as sequence runs out before graph ends
	    $lenn = $seq{$name}{pend}-$ppos;
	}
	$pstart = 0;

	my $subn = substr($seqn,$pstart,10);
	my $subq = substr($seq,$qstart,$lenn<10?$lenn:10);
	#print "\t$1\t",length($seq),"\t$subn\t$subq\t$lenn,$pstart,$qstart\t",length($seqn),"\n";

	my $query = "";
	my $prefix;
	if ($qstart >= $kmer-1) {
	    $prefix = substr($seq,$qstart-($kmer-1),$kmer-1);
	} else {
	    #$prefix = substr($seq, 0, $qstart);
	    $prefix = "N" x ($kmer-1 - $qstart) . substr($seq, 0, $qstart);
	}

	while ($lenn > 0) {
	    if ($oplen == 0) {
		$cig =~ m/(\d+)(.)(.*)/;
		$op = $2;
		$oplen = $1;
		$cig = $3;
	    }
	    my $x = substr($seq, $qstart,5);
	    my $y = substr($seqn,$pstart,5);
	    #print "> $op $oplen rem $lenn\t$x $y\t$ppos\n";

	    if ($op eq "=" || $op eq "X") {
		my $iend = $lenn<$oplen?$lenn:$oplen;
		for (my $i = 0; $i < $iend; $i++) {
		    #print "INCR gfa{$node}{$op}[$pstart+$i]\n";
		    @{$gfa{$node}{$op}}[$pstart+$i]++;
		}
		if ($lenn <= $oplen) {
		    $query .= substr($seq,$qstart,$lenn);
		    $qstart += $lenn;
		    $pstart += $lenn;
		    $oplen -= $lenn;
		    $lenn = 0;
		} else {
		    $query .= substr($seq,$qstart,$oplen);
		    $qstart += $oplen;
		    $pstart += $oplen;
		    $lenn -= $oplen;
		    $oplen = 0;
		}
	    } elsif ($op eq "D") {
		my $iend = $lenn<$oplen?$lenn:$oplen;
		for (my $i = 0; $i < $iend; $i++) {
		    @{$gfa{$node}{$op}}[$pstart+$i]++;
		}
		if ($lenn <= $oplen) {
		    $pstart += $lenn;
		    $oplen -= $lenn;
		    $lenn = 0;
		} else {
		    $pstart += $oplen;
		    $lenn -= $oplen;
		    $oplen = 0;
		}
	    } elsif ($op eq "I") {
		@{$gfa{$node}{$op}}[$pstart]++;
		$query .= substr($seq,$qstart,$oplen);
		$qstart += $oplen;
		$oplen = 0;
	    } else {
		die "unknown op $op";
	    }
	}

	#print "# $node $seqn\n@ $node $query\n";

	if ($last_node ne $node) {
	    push(@{$gfa{$node}{query}}, $prefix . $query);
	} else {
	    @{$gfa{$node}{query}}[-1] .= $query;
	}


	my $seq_end = $qstart;
	$ppos += $pstart;
	#print "QRY $node ",substr($seq,$seq_start,$seq_end-$seq_start),"\n";
	#print "REF $node $seqn\n";
	#print "POS $pstart,$qstart\n";

	$pstart = 0; # start of next segment
	#$qstart += length($seqn);

	$last_node = $node;
    }
}

# Report graph and query sequences per node name
foreach my $node (sort keys %gfa) {
    print "\@$node\n";
    print "$gfa{$node}{seq}\n";
    foreach (@{$gfa{$node}{query}}) {
	print "$_\n";
    }
}

# DEBUGGING below
__END__

# Report frequency of mismatch by node by position
foreach my $node (sort keys %gfa) {
    print "#node $node gfa{$node}{=}\n";
    next unless exists $gfa{$node}{"="};
    # assume "M" is last; fixme
    my $end = scalar(@{$gfa{$node}{"="}}) + 0;
    for (my $i=0; $i<$end; $i++) {
	print "# $node\t$i ", substr($gfa{$node}{seq}, $i, 1),
	    "\t",
	    exists($gfa{$node}{"="}[$i])?$gfa{$node}{"="}[$i]:0,
	    "\t",
	    exists($gfa{$node}{X}[$i])?$gfa{$node}{X}[$i]:0,
	    "\t",
	    exists($gfa{$node}{D}[$i])?$gfa{$node}{D}[$i]:0,
	    "\t",
	    exists($gfa{$node}{I}[$i])?$gfa{$node}{I}[$i]:0,
	    "\n";
    }
}
