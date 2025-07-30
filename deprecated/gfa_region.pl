#!/usr/bin/perl -w
use strict;

# use GFA SN:Z and SO:i tags to identify the portion of a contig referred
# to by the nodes in it.

my $first = 0;
my $last = 0;
my $first_sn = "";
my $last_sn = "";
while (<>) {
    next unless (/^S/);
    m/SN:Z:(\S+)/;
    my $sn=$1;
    m/LN:i:(\d+)/;
    my $ln=$1;
    m/SO:i:(\d+)/;
    my $so=$1;
    if (!$first) {
	$first = $so;
	$first_sn = $sn;
    } else {
	$last = $so + $ln;
	$last_sn = $sn;
    }
}
if ($first_sn ne $last_sn) {
    print "ERROR: first and last SN:Z: differ: $first_sn, $last_sn\n";
    exit 1;
}
print "$first_sn:$first-$last\n";
