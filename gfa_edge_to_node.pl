#!/usr/bin/perl -w

# FAILS as an experiment.  It needs full graph traversal to do this.

# Reads a GFA file with edge and node weights and applies edge transformations
# to split node weights into plus and minus direction.
#
# We read edge EC:i node SC:f, and add node sp:f and sm:f tags

use strict;

# Parse GFA
my %seq;      # node sequence
my %node;     # node GFA line (S)
my @edge;     # edge GFA line (L)
my $edge_num=0;
my %edge_in;  # index into @edge above
my %edge_out; # index into @edge above
local $"="\t";
while (<>) {
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

# Count routes into a node via edges and emit the modified GFA
foreach my $n (sort keys %node) {
    print "$n \n";
    my ($SC) = $node{$n} =~ m/SC:[if]:(\d+(\.\d+)?)/;
    my %dc=qw/+ 0 - 0/;
    foreach my $e (@{$edge_in{$n}}) {
    print "edge in num: $e\n";
	my ($d1,$d2,$EC) = $edge[$e] =~ m/.*([+-]).*([+-]).*EC:i:(\d+)/;
	$dc{$d2}+=$EC;
    }
    foreach my $e (@{$edge_out{$n}}) {
    print "edge out num: $e\n";
	my ($d1,$d2,$EC) = $edge[$e] =~ m/.*([+-]).*([+-]).*EC:i:(\d+)/;
	$dc{$d1}+=$EC;
    }
    my $dn=$dc{"+"}+$dc{"-"};
    print "dn: $dn \n";
    print "dc+: ", $dc{"+"}, "\n";
    print "dc-: ", $dc{"-"}, "\n";
    print "\n";
    $dc{"+"}=($dc{"+"}+0)/($dn+1) * $SC;
    $dc{"-"}=($dc{"-"}+0)/($dn+1) * $SC;

    # printf("%-10s %5.1f %5.1f %5.1f\n",
	#    $n, $SC, $dc{"+"}, $dc{"-"});
    # print "$node{$n}\tsp:f:",$dc{"+"},"\tsm:f:",$dc{"-"},"\n";
}
# foreach my $e (@edge) {
#     print "$e","\n"
# }


__END__;

/tmp/tangle/illumina_km_00103/seq_1091-0027-#1#1.gfa

[  21] s42+
[  22] s18+
[  23] s43-
[  24] s21-
[  25] s20-
[  26] s19-
[  27] s44+
[  28] s23+

L	s16	+	s42	+	SR:i:1	L1:i:233	L2:i:56	EC:i:41
L	s17	+	s18	+	SR:i:0	L1:i:1	L2:i:24	EC:i:0
L	s18	+	s19	+	SR:i:0	L1:i:24	L2:i:20	EC:i:0
L	s18	+	s43	-	SR:i:4	L1:i:24	L2:i:12	EC:i:0
L	s19	+	s20	+	SR:i:0	L1:i:20	L2:i:76	EC:i:16
L	s19	+	s49	+	SR:i:3	L1:i:20	L2:i:76	EC:i:0
L	s19	-	s44	+	SR:i:4	L1:i:20	L2:i:45	EC:i:27
L	s20	+	s21	+	SR:i:0	L1:i:76	L2:i:120	EC:i:26
L	s21	+	s22	+	SR:i:0	L1:i:120	L2:i:168	EC:i:0
L	s21	+	s43	+	SR:i:1	L1:i:120	L2:i:12	EC:i:0
L	s22	+	s23	+	SR:i:0	L1:i:168	L2:i:905	EC:i:2
L	s23	+	s24	+	SR:i:0	L1:i:905	L2:i:492	EC:i:85
L	s23	+	s45	+	SR:i:1	L1:i:905	L2:i:565	EC:i:0
L	s42	+	s18	+	SR:i:1	L1:i:56	L2:i:24	EC:i:47
L	s43	+	s44	+	SR:i:1	L1:i:12	L2:i:45	EC:i:0
L	s44	+	s23	+	SR:i:1	L1:i:45	L2:i:905	EC:i:69
L	s49	+	s21	+	SR:i:3	L1:i:76	L2:i:120	EC:i:0

     +    +    -    -    -    -    +    +
So +42+ +18+ +43- +21+ +20+ +19+ -44+ +23+

SR:i  rank
L1:i  arc_len,  length of first seq
L2:i  arc_lw,   length of second seq
