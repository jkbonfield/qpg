#!/bin/bash

seed=$1
k1=$2
k2=$3
k3=$4
num_training=$5
# mode=$6


. ${CONFIG:-$QDIR/config_hifi_km.sh} # $mode

# Produce sequences
echo "=== Creating population $seed"
eval genome_create "$genome_opts" -s "$seed" -P1 -G0 > pop0.fa
eval genome_create "$genome_opts" -s "$seed" > pop.fa

# Use the first 30 to build a pangenome
perl -lane 'if (/>/) {s/>//;close(FH),open(FH,">$_");print FH ">$_";next} {print FH $_} END {close(FH)}' pop.fa

# Random selection of 20 from the last 30 for training
ls -1 seq_* | tail -50 | shuf --random-source=/usr/bin/emacs > fofn.all
head -40 fofn.all > fofn.train

# Use the last 10 as the test set
#ls -1 seq_* | tail -10 > fofn.test
tail -"$num_training" fofn.all > fofn.test

[ "x$use_syncasm" = "x1" ] && exit 0

echo "=== Building GFA"
#minigraph -l 1000 -d 1000 -cxggs `cat fofn.train` > pop.gfa 2>pop.minigraph.err
#minigraph  -l 1000 -d 1000 -n 5,30 -cxggs `cat fofn.train` > pop.gfa 2>pop.minigraph.err
minigraph  -l 1000 -d 10000 -n 5,20 -cxggs $(cat fofn.train) > pop.gfa 2>pop.minigraph.err
echo "Node count:" "$(grep -E -c '^S' pop.gfa)"

cat `cat fofn.train` > train.fa

# Build the nodeseq index
if expr "$annotate" : ".*km.*"
then
    echo "=== Creating nodeseq files, kmers $k1, $k2, $k3"
    /nfs/sam_scratch/jkb/conda22.old/bin/GraphAligner -g pop.gfa -f train.fa -x vg -a pop.gaf >pop.GraphAligner.out
    gaf2nodeseq2.pl pop.gaf train.fa pop.gfa "$k1" > pop.gfa.ns"$k1"
    gaf2nodeseq2.pl pop.gaf train.fa pop.gfa "$k2" > pop.gfa.ns"$k2"
    gaf2nodeseq2.pl pop.gaf train.fa pop.gfa "$k3" > pop.gfa.ns"$k3"

#    # Faster alternative using population GFA file only
#    # Oddly this is sometimes better than using GraphAligner despite fewer kmers
#    gfa2nodeseq.pl pop.gfa $k1 > pop.gfa.ns$k1
#    gfa2nodeseq.pl pop.gfa $k2 > pop.gfa.ns$k2
#    gfa2nodeseq.pl pop.gfa $k3 > pop.gfa.ns$k3
fi

exit 0
