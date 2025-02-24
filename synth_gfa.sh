#!/bin/bash

qdir=${QDIR:-.}
seed=$1

test=$2

# Produce 30 sequences
#$qdir/genome_create -P 30 -s 0 -l 10000 -S 0.001 -C 0.001 -N 0.02 -n 0.01 \
#     -A 0.00005 -L 0.00002 -T 0.0002 -s $seed > pop.fa
$qdir/genome_create -P 30 -s 0 -l 10000 -S 0.0005 -C 0.0002 -N 0.02 -n 0.01 \
     -A 0.00005 -L 0.00002 -T 0.0002 -s $seed > pop.fa

# Use the first 20 to build a pangenome
rm seq_0*
perl -lane 'if (/>/) {s/>//;close(FH),open(FH,">$_");print FH ">$_";next} {print FH $_} END {close(FH)}' pop.fa

ls -1 seq_* | head -20 > fofn.train
minigraph -l 1000 -d 1000 -cxggs `cat fofn.train` > pop.gfa
echo "Node count:" `egrep -c '^S' pop.gfa`

# Build the nodeseq index
head -40 pop.fa > train.fa
ssh seq4d "cd `pwd`;/nfs/sam_scratch/jkb/conda22.old/bin/GraphAligner -g pop.gfa -f train.fa -x vg -a pop.gaf --threads 32"
$qdir/gaf2nodeseq.pl.tmp pop.gaf train.fa pop.gfa 75 > pop.gfa.ns75
$qdir/gaf2nodeseq.pl.tmp pop.gaf train.fa pop.gfa 50 > pop.gfa.ns50
$qdir/gaf2nodeseq.pl.tmp pop.gaf train.fa pop.gfa 25 > pop.gfa.ns25

# Use the last 10 as the test set
ls -1 seq_* | tail -10 > fofn.test
echo === Test set are in fofn.test ===

if [ "$test" != "" ]
then
    for i in `cat fofn.test`
    do
	echo === $i ===
	$qdir/align_synth_gfa-k3.sh $i pop.gfa pop.gfa.ns{75,50,25}
    done
fi
