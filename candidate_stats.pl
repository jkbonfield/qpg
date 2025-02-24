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
    print STDERR "Usage: candidate_stats.pl true.fa candidate.fa\n";
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
my %ncontigs; # No. of contigs
my %nbreaks;  # No. of alignments minus ncontigs

# Run aligner the candidate, aligning against truth.
# NB candidate.fa may be in multiple contigs.

# Should we remove --no-long-join?
#my $aligner = "minimap2";
#my $aligner_index = "/bin/true";
#my $aligner_opts = "-a -O 2 -E 4 -r 100 -x lr:hq --no-long-join --secondary=no";

my $aligner = "bwa mem";
my $aligner_index = "bwa index";
my $aligner_opts = "-B4 -O4 -E2";
#my $aligner_opts = "";
system("$aligner_index $ARGV[0] 2>/dev/null") && die "$!";
open(my $mm, "$aligner $aligner_opts @ARGV 2>/dev/null|tee /tmp/_.sam|samtools sort -O sam|") || die "$!";

# Read header and SAM records
my $last_contig = "";
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
	$ncontigs{$rname} = 0;
	$nbreaks{$rname} = 0;
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
    #print "Read $F[0] at $rpos\n";
    my $hc = 0;

    if ($last_contig eq "" || $last_contig ne $F[0]) {
	$ncontigs{$rname}++;
    } else {
	$nbreaks{$rname}++;
    }
    $last_contig = $F[0];

    # Process cigar
    foreach ($cigar =~ m/\d+./g) {
	m/(\d+)(.)/;
	#print "$rpos ",$qpos+$hc," $1$2\n";
	if ($2 eq "M") {
	    substr($cov{$rname}, $rpos, $1) = "1" x $1;

	    #print "FILL $rpos to ",$rpos+$1,"\n";
	    for (my $i=0;$i<$1;$i++) {
		#print substr($rseq{$rname},$rpos+$i,1), " ",
		#	substr($seq,$qpos+$i,1), " ";
		if (substr($rseq{$rname},$rpos+$i,1) ne
		    substr($seq,$qpos+$i,1)) {
		    substr($match{$rname}, $rpos+$i, 1) = "n";
		    #print "MIS   n ",$rpos+$i,"\n";
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
		#print "del N $rpos to ",$rpos+$1,"\n";
	    } else {
		substr($match{$rname}, $rpos, $1) = "d" x $1;
		#print "del D $rpos to ",$rpos+$1,"\n";
	    }
	    $rpos += $1;
	} elsif ($2 eq "I") {
	    $nindels{$rname}++ if ($1 >= $indel_max);
	    my $ichar = chr("0" + ($1>255-48?255:ord($1+48)));
	    substr($ins{$rname}, $rpos, 1) = $ichar;
	    #$ins += $i;
	    #print "ins $rpos to ",$rpos+$1,"\n";
	    $qpos += $1;
	} elsif ($2 eq "S") {
	    $qpos += $1;
	} elsif ($2 eq "H") {
	    $hc += $1; # used for debugging only
	} else {
	    print STDERR "Unexpected CIGAR op $2\n";
	}
    }

    #print STDERR "Read $F[0] covers $F[3] to $rpos\n";
}

my $total_len = 0;
my $total_cov = 0;
my $total_contig = 0;
my $total_break = 0;
my $total_indels = 0;
my $total_diffs = 0;
my $total_id_yes = 0;
my $total_id_no = 0;

print "Per contig: (rname      len     covered ncontig nbreaks nindel  ndiff   %identity)\n";
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
		if ($d ne "y") {
		    $id_no++;
		}
		if ($i ne "0") {
		    $id_no += ord($i)-ord("0");
		}
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
    printf("%-20s\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f%%\n",
	   $_, $SQ_len{$_}, length($cov{$_}), $ncontigs{$_}, $nbreaks{$_},
	   $nindels{$_}, $ndiffs{$_}, 100*$identity);
    $total_len    += $SQ_len{$_};
    $total_cov    += length($cov{$_});
    $total_contig += $ncontigs{$_};
    $total_break  += $nbreaks{$_};
    $total_indels += $nindels{$_};
    $total_diffs  += $ndiffs{$_};
    $total_id_yes += $id_yes;
    $total_id_no  += $id_no;
}

# Not needed anymore

# my $identity = ($total_id_yes+.01)/($total_id_yes+$total_id_no+.01);
# print "\nTotals:     (rname      len     covered ncontig nbreaks nindel  ndiff   %identity)\n";
# printf("%-20s\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f%%\n", "Total",
#        $total_len, $total_cov, $total_contig, $total_break,
#        $total_indels, $total_diffs, 100*$identity);

