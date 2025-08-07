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

my $tot_x = 0;
my $tot_y = 0;
my $cnt = 0;
my $mse = 0;
my $mean_diff = 0;
while (<>) {
    chomp($_);
    my @F = split(/\t/, $_);

    # F1 score of used and coverage %age
    my $y = $F[6];             # F1 score; R=0.870, 0.898
    my $x = $F[11];             # final estimated score from estimate_qual.pl

# Slight refinements
#    $x = $F[1] - 2*$F[2] - 2*($F[9]-1) - $F[10]/100 - $F[3]/30;
#    $x = ($x**0.9) + 41.7;   # .902 16.7
#
#    $x = $F[1] - 2*$F[2]**0.7*1.2 - 2*($F[9]-1) - $F[10]/100 - $F[3]/30;
#    $x = ($x**0.9) + 41.6;   # .911 14.9
#
#    $x = $F[1] - 2*$F[2]**0.5*1.37 - 2*($F[9]-1) - $F[10]/100 - $F[3]/30;
#    $x = ($x**0.9) + 41.8;   # .913 14.8
#
#    $x = $F[1]**1.1*.634 - 2*$F[2]**0.5*1.37 - 2*($F[9]-1) - $F[10]/100 - $F[3]/30;
#    $x = ($x**0.9) + 41.8;   # .916 14.4

    $x = $F[1]**1.4*.1614 - 2*$F[2]**0.5*1.37 - 2*($F[9]-1) - $F[10]/100 - $F[3]/30;
    $x = 0.1 if $x < 0.1;
    $x = ($x**0.9) + 41;   # .919 13.7

#    $x = $F[1]**1.4*.165 - 2*$F[2]**0.5*1.30 - 2*($F[9]-1)**0.7*1.0 - $F[10]/110 - $F[3]/30;
#    $x = 0.1 if $x < 0.1;
#    $x = ($x**0.9) + 39.4;   # .919 12.7

#    # On remapped output in estimate_qual.pl (see comment)
#    $x = $F[1]**1.4*.162 - 2*$F[2]**0.5*1.25 - 2*($F[9]-1)**0.7*0.8 - $F[10]/115 - $F[3]/28;
#    $x = 0.1 if $x < 0.1;
#    $x = ($x**0.9) + 39;   # .925 11.8


    #print STDERR "M: $F[1] $F[2] $F[3] $F[9] $F[10]\n";
    # AVG             95.5  1.88  254.6 2.41  390.4 

    #$x = $F[1] - 2*$F[2] - 2*($F[9]-1) - $F[10]/100 - ($F[3]**0.8)/6;
    #$x=0.1 if ($x<0.1);
    #$x = ($x**0.8) + 60;   # .902 16.7

    #$x = $F[1] - 2*$F[2] - 2*($F[9]-1) - $F[10]/100 - ($F[3]**1.2)/100;
    #$x=0.1 if ($x<0.1);
    #$x = ($x**0.9) + 41.5;   # .895 19.5

#    # Coverage only.  Great correlation
#    $y=$F[4]; $x = $F[1] - $F[2]*2 - $F[9]*1.75 - $F[3]/35; $x *= 1.13; # cov% R=0.909
#    $x = 0.1 if $x < 0.1;
#    $x = 100 if $x > 100;
#    my $x1 = $x;

#    # Used% only.
#    # OK correlation, but too big a spread of -ve to high +ve
#    $y=$F[5]; $x = $F[1] - $F[2]*0.5 - $F[9]*3 - $F[3]/35 - $F[10]/18; # used% R=0.753
#    $x = 0.1 if $x < 0.1;
#    $x = 100 if $x > 100;
#    my $x2 = $x;

#    # Try combining the Cov only and Used only together into a unified
#    # recalibrated score.  Gain is minimal and it's really complex.
#    $y = $F[6];                  # F1 real
#    $x = (3.5*$x1+$x2)/4.5;     # 0.901 16.3
#    $x = $x**0.9 + 38;

#    $tot_x += $x; $tot_y += $y; $cnt++;
#    print STDERR "$cnt $tot_x $tot_y ",$tot_x/$tot_y,"\n" if $cnt==255;

    $x = 0.1 if $x < 0.1;
    $x = 100 if $x > 100;

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
    $mse += ($x - $y)**2;
    $mean_diff += $x - $y;

    print STDERR "PLOT $y $x\n";
}

# Means
my $n = scalar(@x);
my $mean_x = $sum_x / $n;
my $mean_y = $sum_x / $n;
$mse /= $n;
$mean_diff /= $n;

my $sx = sqrt($n*$sum_x2 - $sum_x*$sum_x);
my $sy = sqrt($n*$sum_y2 - $sum_y*$sum_y);
my $r = ($n*$sum_xy - $sum_x * $sum_y) / ($sx * $sy);
print "R = $r   mse = $mse   delta = $mean_diff\n";
