#!/usr/bin/perl -w

# Turn a GAF alignment output from GraphAligner + GFA file with node sequence
# into a path of appended GFA sequences.

# Ie this is the sequence to which the query was aligned, and the sequence
# for which the CIGAR string applies.

# Usage: gaf2seq in.gfa [in.gaf]

# Parse GFA; minimally
open(my $gfa, "<", shift(@ARGV)) || die;
while (<$gfa>) {
    chomp($_);
    my @F = split(/\s+/, $_);
    next unless scalar(@F) && $F[0] eq "S"; # skip other fields for now
    $gfa{$F[1]}{seq} = uc($F[2]);
}


# Parse GAF
while (<>) {
    chomp();
    my @F=split("\t", $_);
    $F[5]=join(">",reverse split("<",$F[5])) if ($F[5]=~/</);

    my $seq="";
    foreach(split(">", $F[5])) {
	next unless $_;
	$seq .= $gfa{$_}{seq};
    }
    print ">$F[0]\n$seq\n";
}
