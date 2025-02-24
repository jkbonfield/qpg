#!/usr/bin/perl -w

use strict;

sub eord {
    my ($seed,$str,$enc) = @_;
    my $x = int(rand(1<<30));
    srand($seed);
    my $s="";

    for (my $i=0; $i < 4; $i++) {
	$s .= chr(int(rand(256))) if $enc;
    }

    foreach (split("", $str)) {
	$s .= chr(ord($_) ^ int(rand(256)))
    }

    for (my $i=0; $i < 4; $i++) {
	$s .= chr(int(rand(256))) if $enc;
    }

    srand($x);

    return $enc ? $s : substr($s, 4, -4);
}

# Load a GAF file, process the paths, and turn into a DNA sequence string.
# These are always A + 1-bit encoded ([CG]) versions of the node name.
# We reverse complement the entire thing if the node starts with <

# Turn GAF paths into a fasta sequence
my $count = 0;
my %seen = ();
while (<>) {
    my @F = split(/\t/, $_);
    next if $F[5] =~ /^</; # TODO: reverse complement

    # Remove duplicates
    next if $seen{$F[5]};
    $seen{$F[5]}=1;

    print STDERR $F[5],"\n";

    my $seq = "";
    foreach ($F[5] =~ m/[<>][^<>]*/g) {
	m/(.)(.*)/;
	#$seq .= ($1 eq ">") ? "AA" : "TT";
	#$seq .= "AA";
	#$seq .= "N";
	foreach (split("", eord(1, $2, 1))) {
	#foreach (split("", $2)) {
	    my $o = ord($_);
	    my @bit = qw/C G/;
	    #print "$o\n";
	    for (my $i=7; $i>=0; $i--) {
		$seq .= $bit[$o & (1<<$i) ? 1 : 0];
		#print $o & (1<<$i) ? 1 : 0;
	    }
	    #print "\n";
	}
    }

    print ">$count\n$seq\n";
    $count++;
}

# Assemble (phrap on FASTA)

# Map FASTA back to contigs

# Count depth of contig to determine real routes from incorrect routes



