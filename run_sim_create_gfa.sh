#!/bin/bash

seed=$1
num_training=$2

. ${CONFIG:-$QDIR/config_hifi_km.sh}

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

exit 0
