#!/bin/bash

qdir=${QDIR:-.}
seed=$1

test=$2

# Produce sequences

# Challenging: ~400 nodes
#$qdir/genome_create -P 50 -l 20000  -S 0.01 -C 0.001 -N 0.002 -n 0.002 \
#     -A 0.001 -L 0.0002 -T 0.0004 -s $seed > pop.fa

# Complex: 100-200 nodes
#$qdir/genome_create -P 50 -l 10000  -S 0.01 -C 0.0005 -N 0.005 -n 0.005 \
#     -A 0.0005 -L 0.0001 -T 0.0002 -s $seed > pop.fa

# Medium: 40-80 nodes
$qdir/genome_create -P 50 -l 5000  -S 0.001 -C 0.001 -N 0.01 -n 0.01 \
     -A 0.0005 -L 0.0001 -T 0.0002 -s $seed > pop.fa

# Simpler: 10-50 nodes
#$qdir/genome_create -P 50 -l 5000 -S 0.001 -C 0.0005 -N 0.01 -n 0.01 \
#     -A 0.00005 -L 0.00002 -T 0.0002 -s $seed > pop.fa

# Use the first 30 to build a pangenome
rm seq_0*
perl -lane 'if (/>/) {s/>//;close(FH),open(FH,">$_");print FH ">$_";next} {print FH $_} END {close(FH)}' pop.fa

ls -1 seq_* | head -30 > fofn.train
minigraph -l 1000 -d 1000 -cxggs `cat fofn.train` > pop.gfa
echo "Node count:" `egrep -c '^S' pop.gfa`

# Build the nodeseq index
head -40 pop.fa > train.fa
ssh seq4d "cd `pwd`;/nfs/sam_scratch/jkb/conda22.old/bin/GraphAligner -g pop.gfa -f train.fa -x vg -a pop.gaf --threads 32"
$qdir/gaf2nodeseq2.pl pop.gaf train.fa pop.gfa 75 > pop.gfa.ns75
$qdir/gaf2nodeseq2.pl pop.gaf train.fa pop.gfa 50 > pop.gfa.ns50
$qdir/gaf2nodeseq2.pl pop.gaf train.fa pop.gfa 35 > pop.gfa.ns35

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
