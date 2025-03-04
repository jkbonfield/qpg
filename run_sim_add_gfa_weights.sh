#!/bin/bash

. ${CONFIG:-$QDIR/sim_path_config_hifi.sh}

gfa=$1
query=$2
k1=$3
k2=$4
k3=$5

nodeseq1=$gfa.ns$k1
nodeseq2=$gfa.ns$k2
nodeseq3=$gfa.ns$k3

# Shred the FASTA file into small bits, simulating sequencing
echo === Shredding genome $query
shred_fa=$query.shred.fa
shred.pl -s 1 -l $shred_len -e $shred_err -d $shred_depth $query > $shred_fa

# Map kmers in the shredded fasta to the graph nodeseq
echo === Mapping kmers: k=$k1,$k2,$k3
kmer2node2 -U -k$k1 $nodeseq1 $shred_fa | grep Node > $query.nodes.$k1
kmer2node2 -U -k$k2 $nodeseq2 $shred_fa | grep Node > $query.nodes.$k2
kmer2node2    -k$k3 $nodeseq3 $shred_fa | grep Node > $query.nodes.$k3
merge_kmer2node.pl $query.nodes.$k1 $query.nodes.$k2 $query.nodes.$k3 > $query.nodes

# Annotate the GFA with kmer counts
tag_gfa.pl $gfa < $query.nodes > $query.gfa
