#!/usr/bin/perl -w

# Compares a bunch of runs with different parameters to report which work best
# Usage: compare_runs.pl seed_start seed_end prefix1 prefix2 prefix3...
#
# Output will be seed number and then delta against the mean score for each
# prefix, where score is the F1 score of covered+used percentage (harmonic
# mean).

print "Usage: compare_runs.pl seed_start seed_end prefix1 prefix2 prefix3...\n"
    if scalar(@ARGV) < 3;

my $start = shift(@ARGV);
my $end   = shift(@ARGV);

sub score_run {
    my ($run) = @_;

    my $score = 0;
    my $count = 0;
    foreach my $seq (glob("$run/*eval_seq*")) {
	open(my $fd, "<$seq") || die $seq;
	while (<$fd>) {
	    next if /contig/;
	    my @F = split("\t", $_);
	    # Or harmonic mean F1 score?
	    $F[3]=~s/%//; $F[3]/=100;
	    $F[4]=~s/%//; $F[4]/=100;
	    my $avg = ($F[3]+$F[4])/2;
	    my $f1 = 2*($F[3]*$F[4])/($F[3]+$F[4]);
	    $count++;
	    $score += $f1;
	}
	close($fd);
    }

    return $score/$count;
}

print "#seed";
my @prefixes = ();
foreach (@ARGV) {
    my $path = sprintf("%s%05d", $_, $start);
    next if ( ! -e $path );
    push(@prefixes, $_);
    $p=$_;
    $p=~s#.*/##;
    print "\t$p";
}
print "\n";

for (my $i = $start; $i <= $end; $i++) {
    my $seed = sprintf("%05d", $i);
    my @scores = ();
    my $tscore = 0;
    my $nscore = 0;
    foreach my $p (@prefixes) {
	my $s = score_run($p . $seed);
	push(@scores, $s);
	$tscore += $s;
	$nscore++;
    }
    print "$seed";
    foreach (@scores) {
	my $delta = $_ - ($tscore/$nscore);
	printf("\t%7.3f", $delta*100);
    }
    print "\n";
}

