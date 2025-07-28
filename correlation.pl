#!/usr/bin/perl -w
use strict;

# Runs a correlation analysis on estimate_qual.pl output given a set of rules.

my @x = ();
my @y = ();
my $sum_x = 0;
my $sum_y = 0;
my $sum_x2 = 0;
my $sum_y2 = 0;
my $sum_xy = 0;
while (<>) {
    chomp($_);
    my @F = split(/\t/, $_);
#    my $x = $F[1] - 3*$F[2];    #1  0.745
#    my $y = $F[5] - 0.5*$F[7];

    # F1 score vs primary-2*supp
    my $x = $F[1] - 2*$F[2];    #2  0.751   %primary_mapped - 2*%supp
    my $y = $F[5];              #           F1 score

#    my $x = $F[1] - 5*$F[2];    #3  0.712
#    my $y = $F[5] - 1.5*$F[7] - 0.5*($F[6]-1);

    # cov% mainly, but not helpful as an overall quality metric
#    my $x = $F[1] - 2*$F[2];    #4  0.803
#    my $y = $F[3] + 0.3*$F[5];

#    my $x = $F[1] - 2.5*$F[2];    #5  0.743
#    my $y = 2*($F[3]<$F[4]?$F[3]:$F[4]) + 2*$F[5] - 2.5*$F[7] - 0*($F[6]-1);

    #4 0.882667072607652 0.791658437480038 0.823728869474239 0.847899706619491
    #2 0.879090564062543 0.763626912230842 0.796498220624491 0.825303243395947
    #1 0.870630645569386 0.778913227959159 0.795614432686065 0.820286443681295
    #3 0.837021988717656 0.765224675820016 0.767925027389423 0.786678927613394

    $sum_x  += $x;
    $sum_x2 += $x*$x;
    $sum_y  += $y;
    $sum_y2 += $y*$y;
    $sum_xy += $x*$y;
    push(@x, $x);
    push(@y, $y);
}

# Means
my $n = scalar(@x);
my $mean_x = $sum_x / $n;
my $mean_y = $sum_x / $n;

my $sx = sqrt($n*$sum_x2 - $sum_x*$sum_x);
my $sy = sqrt($n*$sum_y2 - $sum_y*$sum_y);
my $r = ($n*$sum_xy - $sum_x * $sum_y) / ($sx * $sy);
print "R = $r\n";
