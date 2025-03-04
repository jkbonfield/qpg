#!/bin/bash

. ${CONFIG:-$QDIR/sim_path_config_hifi.sh}

seed=$1
k1=$2
k2=$3
k3=$4


# Produce sequences
echo === Creating population
eval genome_create $genome_opts -s $seed > pop.fa

# Use the first 30 to build a pangenome
echo === Building GFA
perl -lane 'if (/>/) {s/>//;close(FH),open(FH,">$_");print FH ">$_";next} {print FH $_} END {close(FH)}' pop.fa

ls -1 seq_* | head -30 > fofn.train
minigraph -l 1000 -d 1000 -cxggs `cat fofn.train` > pop.gfa 2>pop.minigraph.err
echo "Node count:" `egrep -c '^S' pop.gfa`

head -40 pop.fa > train.fa

# Build the nodeseq index
if [ "x$use_mg" != "x1" ]
then
    echo === Creating nodeseq files, kmers $k1, $k2, $k3
    /nfs/sam_scratch/jkb/conda22.old/bin/GraphAligner -g pop.gfa -f train.fa -x vg -a pop.gaf >pop.GraphAligner.out
    gaf2nodeseq2.pl pop.gaf train.fa pop.gfa $k1 > pop.gfa.ns$k1
    gaf2nodeseq2.pl pop.gaf train.fa pop.gfa $k2 > pop.gfa.ns$k2
    gaf2nodeseq2.pl pop.gaf train.fa pop.gfa $k3 > pop.gfa.ns$k3
fi

# Use the last 10 as the test set
ls -1 seq_* | tail -10 > fofn.test
