#!/usr/bin/perl -w

# Evaluate a SAM alignment of assembly contigs to count:
# - percent genome coverage,
# - number of fragments
# - number of indels >= 10bp,
# - number of regions with 30%+ diff >= 100bp (excludes long deletions above)

use strict;

# Load reference if set
my $ref_file = pop(@ARGV) if (scalar(@ARGV) == 2);
my %rseq;
open(my $fd, "<", $ref_file) || die $!;
my $tmp_seq = "";
my $tmp_name = "";
while (<$fd>) {
    if (/^>/) {
	if ($tmp_name ne "") {
	    $rseq{$tmp_name} = $tmp_seq;
	}
	m/^>(\S+)/;
	$tmp_name = $1;
	$tmp_seq = "";
	next;
    }
    chomp($_);
    $tmp_seq .= $_;
}
if ($tmp_name ne "") {
    $rseq{$tmp_name} = $tmp_seq;
}
close($fd);

my %SQ_len;   # Sequence lines, mapping SN name to length
my %cov;      # Coverage string per ref base (0/1)
my %delta;    # Delta string per ref base (0/1)
my %match;    # y=match, n=mismatch, .=unknown
my %ndiffs;   # No. regions of >= $diff_max long deltas
my %nindels;  # No. regions of >= $indel_max long indels
my %nfrags;   # No. of alignments (minus 1 is no. of breaks)

my $indel_max=10;
my $diff_max=100;
my $diff_perc=30;

# Read header and SAM records
while (<>) {
    chomp($_);
    my @F = split("\t", $_);

    #--- header
    if (/^\@SQ/) {
	my ($rname) = ($_ =~ m/SN:(\S+)/);
	my ($ln) = ($_ =~ m/LN:(\d+)/);
	$SQ_len{$rname} = $ln;
	$cov{$rname} = "0" x $SQ_len{$rname};
	$delta{$rname} = "0" x $SQ_len{$rname};
	$match{$rname} = "." x $SQ_len{$rname};
	$ndiffs{$rname} = 0;
	$nindels{$rname} = 0;
	$nfrags{$rname} = 0;
    }
    next if (/^\@/);

    #--- alignments

    # Skip secondary alignments
    next if ($F[1] & 0x100);

    my $rname = $F[2];
    my $rpos  = $F[3]-1;
    my $qpos  = 0;
    my $cigar = $F[5];
    my $seq   = $F[9];

    $nfrags{$rname}++;

    # Process cigar
    foreach ($cigar =~ m/\d+./g) {
	m/(\d+)(.)/;
	if ($2 eq "M") {
	    substr($cov{$rname}, $rpos, $1) = "1" x $1;
	    if ($ref_file) {
		for (my $i=0;$i<$1;$i++) {
		    if (substr($rseq{$rname},$rpos+$i,1) ne
			substr($seq,$qpos+$i,1)) {
			substr($delta{$rname}, $rpos+$i, 1) = "1";
			substr($match{$rname}, $rpos+$i, 1) = "n";
		    } else {
			substr($match{$rname}, $rpos+$i, 1) = "y";
		    }
		}
	    }
	    $rpos += $1;
	    $qpos += $1;
	} elsif ($2 eq "D") {
	    $rpos += $1;
	    $nindels{$rname}++ if ($1 >= $indel_max);
	    if ($1 < $indel_max) {
		substr($delta{$rname}, $rpos, $1) = "1" x $1;
		substr($match{$rname}, $rpos, $1) = "n" x $1;
	    }
	} elsif ($2 eq "I") {
	    $qpos += $1;
	    $nindels{$rname}++ if ($1 >= $indel_max);
	    substr($delta{$rname}, $rpos, 1) = "1";
	    substr($match{$rname}, $rpos, 1) = "n";
	} elsif ($2 eq "S") {
	    $qpos += $1;
	} elsif ($2 ne "H") {
	    print STDERR "Unexpected CIGAR op $2\n";
	}
    }
}

my $total_len = 0;
my $total_cov = 0;
my $total_break = 0;
my $total_indels = 0;
my $total_diffs = 0;
print "Per contig: (rname len covered nbreaks nindel ndiff)\n";
foreach (sort keys %SQ_len) {
    # Count coverage as number of non-zeros.
    my $n = 0;
    $cov{$_} =~ tr/1//dc;

    # Turn delta into a scoring function and track start / stop locations
    # based on %age diff
    #print $delta{$_},"\n";
    my $delta_score = 0;
    my $start_pos = 0;
    my $end_pos = 0;
    my $max_score = 0;
    my $len = length($delta{$_});
    my $diff_cost = 1/($diff_perc/100);
    for (my $pos = 0; $pos < $len; $pos++) {
	my $d = substr($delta{$_}, $pos, 1);
	#print "$pos $d ",int($delta_score),"/",int($max_score),"\n";
	if ($d eq 0) {
	    if (--$delta_score < 0) {
		if ($max_score > 0 && $end_pos - $start_pos >= $diff_max) {
		    print "DELTA $start_pos $end_pos ",$end_pos-$start_pos,"\n";
		    $pos = $end_pos; # restart to find up-down-up combos.
		}
		$delta_score = 0;
		$max_score = 0;
	    }
	} else {
	    if ($delta_score == 0) {
		$start_pos = $pos;
		$end_pos = $pos+1;
	    }
	    $delta_score += $diff_cost;
	    if ($delta_score > $max_score) {
		$max_score = $delta_score;
		$end_pos = $pos+1;
	    }
	}
    }
    if ($max_score > 0 && $end_pos - $start_pos >= $diff_max) {
	print "DELTA $start_pos $end_pos ",$end_pos-$start_pos,"\n";
    }

    printf("%-20s\t%d\t%d\t%d\t%d\t%d\n", $_, $SQ_len{$_}, length($cov{$_}),
	   $nfrags{$_}-1, $nindels{$_}, $ndiffs{$_});
    $total_len    += $SQ_len{$_};
    $total_cov    += length($cov{$_});
    $total_break  += $nfrags{$_}-1;
    $total_indels += $nindels{$_};
    $total_diffs  += $ndiffs{$_};

}

print "\nTotals: (rname len covered nbreaks nindel ndiff)\n";
printf("%-20s\t%d\t%d\t%d\t%d\t%d\n", "Total", $total_len, $total_cov,
       $total_break, $total_indels, $total_diffs);

