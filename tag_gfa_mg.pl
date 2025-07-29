#!/usr/bin/perl -w

use strict;
while (<>) {
    chomp($_);
    my @F = split("\t", $_);
    m/dc:f:([0-9.]*)/;
    my $dc=$1;
    if (/^S/) {
        print "$_\tKC:i:",int($dc*length($F[2])+.5),"\tSC:f:",$dc,"\n";
    } else {
        print "$_\tEC:i:",int($dc+0.5),"\n";
    }
}

