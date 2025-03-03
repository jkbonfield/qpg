#!/bin/bash

seed=$1
k1=$2
k2=$3
k3=$4


# Produce sequences

# Challenging: ~400-800 nodes
#genome_create -P 50 -l 20000  -S 0.01 -C 0.001 -N 0.002 -n 0.002 \
#     -A 0.001 -L 0.0002 -T 0.0004 -s $seed > pop.fa
#genome_create -P 50 -l 20000  -S 0.01 -C 0.001 -N 0.002 -n 0.002 \
#     -A 0.001 -L 0.0002 -T 0.0004 -s $seed > pop.fa

# Complex: 100-200 nodes
#genome_create -P 50 -l 10000  -S 0.01 -C 0.0005 -N 0.005 -n 0.005 \
#     -A 0.0005 -L 0.0001 -T 0.0002 -s $seed > pop.fa
#genome_create -P 50 -l 10000  -S 0.01 -C 0.0005 -N 0.005 -n 0.005 \
#    -A 0.0005 -L 0.0001 -T 0.00005 -I 0.0001 -s $seed > pop.fa

# Medium: 40-80 nodes
#genome_create -P 50 -l 5000  -S 0.001 -C 0.001 -N 0.01 -n 0.01 \
#    -A 0.0005 -L 0.0001 -T 0.00005 -I 0.0002 -s $seed > pop.fa
echo === Creating population
genome_create -P 50 -l 5000  -S 0.001 -C 0.0003 -N 0.01 -n 0.01 \
    -A 0.0005 -L 0.0001 -T 0.00005 -I 0.0001 -s $seed > pop.fa 2>pop.create_genome.err

# Simpler: 10-50 nodes
#genome_create -P 50 -l 5000 -S 0.001 -C 0.0005 -N 0.01 -n 0.01 \
#     -A 0.00005 -L 0.00002 -T 0.0002 -s $seed > pop.fa

# Use the first 30 to build a pangenome
echo === Building GFA
perl -lane 'if (/>/) {s/>//;close(FH),open(FH,">$_");print FH ">$_";next} {print FH $_} END {close(FH)}' pop.fa

ls -1 seq_* | head -30 > fofn.train
minigraph -l 1000 -d 1000 -cxggs `cat fofn.train` > pop.gfa 2>pop.minigraph.err
echo "Node count:" `egrep -c '^S' pop.gfa`

# Build the nodeseq index
echo === Creating nodeseq files, kmers $k1, $k2, $k3
head -40 pop.fa > train.fa
/nfs/sam_scratch/jkb/conda22.old/bin/GraphAligner -g pop.gfa -f train.fa -x vg -a pop.gaf >pop.GraphAligner.out
gaf2nodeseq2.pl pop.gaf train.fa pop.gfa $k1 > pop.gfa.ns$k1
gaf2nodeseq2.pl pop.gaf train.fa pop.gfa $k2 > pop.gfa.ns$k2
gaf2nodeseq2.pl pop.gaf train.fa pop.gfa $k3 > pop.gfa.ns$k3

# Use the last 10 as the test set
ls -1 seq_* | tail -10 > fofn.test
