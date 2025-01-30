#!/usr/bin/perl -w

use strict;

# Evaluate a SAM alignment of assembly contigs to count:
# - percent genome coverage,
# - number of fragments
# - number of indels >= 10bp,
# - number of regions with 30%+ diff >= 100bp (excludes long deletions above)
# - percent identity of aligned segments

my $indel_max=10;
my $diff_max=100;
my $diff_perc=30;

if (scalar @ARGV != 2) {
    print STDERR "Usage: candidate_stats.pl candidate.fa true.fa\n";
    exit 1;
}
my ($ref_file, $qry_file) = @ARGV;

# Load reference
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
my %match;    # y=match, n=mismatch, d=big deletion, .=unknown
my %ins;      # 1=ins, 0/no
my %ndiffs;   # No. regions of >= $diff_max long deltas
my %nindels;  # No. regions of >= $indel_max long indels
my %nfrags;   # No. of alignments (minus 1 is no. of breaks)

# Run minimap2 on the candidate, aligning against truth.
# NB candidate.fa may be in multiple contigs.

my $minimap2 = "minimap2"; # just find it in the PATH for now
my $minimap2_opts = "-a -O 2 -E 4 -r 100 -x lr:hq --no-long-join --secondary=no";
open(my $mm, "minimap2 $minimap2_opts @ARGV 2>/dev/null|") || die "$!";

# Read header and SAM records
while (<$mm>) {
    chomp($_);
    my @F = split("\t", $_);

    #--- header
    if (/^\@SQ/) {
	my ($rname) = ($_ =~ m/SN:(\S+)/);
	my ($ln) = ($_ =~ m/LN:(\d+)/);
	$SQ_len{$rname} = $ln;
	$cov{$rname} = "0" x $SQ_len{$rname};
	$match{$rname} = "." x $SQ_len{$rname};
	$ins{$rname} = "0" x $SQ_len{$rname};
	$ndiffs{$rname} = 0;
	$nindels{$rname} = 0;
	$nfrags{$rname} = 0;
    }
    next if (/^\@/);

    #--- alignments

    # Skip secondary alignments
    next if ($F[1] & 0x100);

    my $rname = $F[2];
    my $rpos = $F[3]-1;
    my $qpos  = 0;
    my $cigar = $F[5];
    my $seq   = $F[9];

    $nfrags{$rname}++;

    # Process cigar
    foreach ($cigar =~ m/\d+./g) {
	m/(\d+)(.)/;
	#print "$rpos ",$rpos+$1," $1$2\n";
	if ($2 eq "M") {
	    substr($cov{$rname}, $rpos, $1) = "1" x $1;

	    #print "FILL $rpos to ",$rpos+$1,"\n";
	    for (my $i=0;$i<$1;$i++) {
		#print substr($rseq{$rname},$rpos+$i,1), " ",
		#	substr($seq,$qpos+$i,1), " ";
		if (substr($rseq{$rname},$rpos+$i,1) ne
		    substr($seq,$qpos+$i,1)) {
		    substr($match{$rname}, $rpos+$i, 1) = "n";
		    #print "match n ",$rpos+$i,"\n";
		} else {
		    substr($match{$rname}, $rpos+$i, 1) = "y";
		    #print "match y ",$rpos+$i,"\n";
		}
	    }
	    $rpos += $1;
	    $qpos += $1;
	} elsif ($2 eq "D") {
	    $nindels{$rname}++ if ($1 >= $indel_max);
	    if ($1 < $indel_max) {
		substr($match{$rname}, $rpos, $1) = "n" x $1;
		#print "match N $rpos to ",$rpos+1,"\n";
	    } else {
		substr($match{$rname}, $rpos, $1) = "d" x $1;
	    }
	    $rpos += $1;
	} elsif ($2 eq "I") {
	    $nindels{$rname}++ if ($1 >= $indel_max);
	    substr($ins{$rname}, $rpos, 1) = "1";
	    $qpos += $1;
	} elsif ($2 eq "S") {
	    $qpos += $1;
	} elsif ($2 ne "H") {
	    print STDERR "Unexpected CIGAR op $2\n";
	}
    }
    #print STDERR "Read $F[0] covers $F[3] to $rpos\n";
}

my $total_len = 0;
my $total_cov = 0;
my $total_break = 0;
my $total_indels = 0;
my $total_diffs = 0;
my $total_id_yes = 0;
my $total_id_no = 0;

print "Per contig: (rname      len     covered nbreaks nindel  ndiff   %identity)\n";
foreach (sort keys %SQ_len) {
    # Count coverage as number of non-zeros.
    my $n = 0;
    $cov{$_} =~ tr/1//dc;

    # Turn $match into a scoring function and track start / stop locations
    # based on %age diff
    my $delta_score = 0;
    my $start_pos = 0;
    my $end_pos = 0;
    my $max_pos = 0;
    my $max_score = 0;
    my $len = length($match{$_});
    my $diff_cost = 1/($diff_perc/100);
    my $id_yes = 0;
    my $id_no = 0;
    for (my $pos = 0; $pos < $len; $pos++) {
	$max_pos = $pos if ($max_pos < $pos);
	my $d = substr($match{$_}, $pos, 1);
	my $i = substr($ins{$_}, $pos, 1);
	if ($d ne "." && $pos >= $max_pos) {
	    if ($d eq "y" && $i eq "0") {
		$id_yes++;
	    } else {
		$id_no++;
	    }
	}
	#print "$pos $d $i ",int($delta_score),"/",int($max_score),"\n";

	if (($d eq "y" || $d eq "d") && $i eq "0") {
	    # Match or a big deletion
	    if (--$delta_score < 0) {
		if ($max_score > 0 && $end_pos - $start_pos >= $diff_max) {
		    #print "DELTA $start_pos $end_pos ",$end_pos-$start_pos,"\n";
		    $ndiffs{$_}++;
		    $pos = $end_pos; # restart to find up-down-up combos.
		}
		$delta_score = 0;
		$max_score = 0;
	    }
	} elsif ($d ne "." || $i eq "1") {
	    if ($delta_score == 0) {
		$start_pos = $pos;
		$end_pos = $pos+1;
	    }
	    $delta_score += $diff_cost;
	    if ($delta_score > $max_score) {
		$max_score = $delta_score;
		$end_pos = $pos+1;
	    }
	} # else not covered
    }
    if ($max_score > 0 && $end_pos - $start_pos >= $diff_max) {
	$ndiffs{$_}++;
	#print "DELTA $start_pos $end_pos ",$end_pos-$start_pos,"\n";
    }


    my $identity = ($id_yes+.01)/($id_yes+$id_no+.01);
    printf("%-20s\t%d\t%d\t%d\t%d\t%d\t%.2f%%\n",
	   $_, $SQ_len{$_}, length($cov{$_}), $nfrags{$_}-1, $nindels{$_},
	   $ndiffs{$_}, 100*$identity);
    $total_len    += $SQ_len{$_};
    $total_cov    += length($cov{$_});
    $total_break  += $nfrags{$_}-1;
    $total_indels += $nindels{$_};
    $total_diffs  += $ndiffs{$_};
    $total_id_yes += $id_yes;
    $total_id_no  += $id_no;
}

my $identity = ($total_id_yes+.01)/($total_id_yes+$total_id_no+.01);
print "\nTotals:     (rname      len     covered nbreaks nindel  ndiff   %identity)\n";
printf("%-20s\t%d\t%d\t%d\t%d\t%d\t%.2f%%\n", "Total", $total_len, $total_cov,
       $total_break, $total_indels, $total_diffs, 100*$identity);

