#!/bin/bash
mode=$3
. ${CONFIG:-$QDIR/sim_path_config_hifi.sh} $mode

# Adds weights using minigraph to align and amend the GFA

gfa=$1
query=$2

# Shred the FASTA file into small bits, simulating sequencing
echo === Shredding genome $query
shred_fa=$query.shred.fa
shred.pl -s 1 -l $shred_len -e $shred_err -d $shred_depth $query > $shred_fa

# Length specific minigraph options.  Long-read vs short-read
if [ $shred_len -gt 500 ]
then
    mg_opts="-x lr -w 15 -l 100 -d 100 -m 50,40 -n 5,5"
else
    mg_opts="-x sr"
fi
mg_opts="$mg_opts -j 0.01 --cov ${MINIGRAPH_OPTS}"


# Run minigraph to get a GFA
eval minigraph $mg_opts $gfa $shred_fa 2>mg.err > $query.mg

# Quantise dc:f field to kmer-count by multplying up the length.
# For edges we use it as EC direct, but in both cases it must be integer.
perl -we '
use strict;
while (<>) {
    chomp($_);
    my @F = split("\t", $_);
    m/dc:f:([0-9.]*)/;
    my $dc=$1;
    if (/^S/) {
        print "$_\tKC:i:",int($dc*length($F[2])+.5),"\n";
    } else {
        print "$_\tEC:i:",int($dc+0.5),"\n";
    }
}' $query.mg > $query.gfa
