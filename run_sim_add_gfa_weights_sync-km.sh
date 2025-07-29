#!/bin/bash

gfa=$1
query=$2
k1=$3
k2=$4
k3=$5
mode=$6

. ${CONFIG:-$QDIR/config_illumina.sh} $mode


nodeseq1=$gfa.ns$k1
nodeseq2=$gfa.ns$k2
nodeseq3=$gfa.ns$k3

# Shred the FASTA file into small bits, simulating sequencing
echo === Shredding genome $query
shred_fa=$query.shred.fa
shred.pl -s 1 -l $shred_len -e $shred_err -d $shred_depth $query > $shred_fa

# De novo assemble
eval syncasm $syncasm_opts $shred_fa
awk '/^S/ {printf(">%s\n",$2);print $3}' syncasm.asm.utg.gfa > $query.asm.fa

# Map kmers in the shredded fasta to the graph nodeseq
echo === Mapping kmers: k=$k1,$k2,$k3
#kmer2node="kmer2node2"
kmer2node="kmer2node4 -G $gfa";

eval $kmer2node -U -k$k1 $nodeseq1 $query.asm.fa -E /dev/stdout | egrep 'Edge|Node' \
| perl -lane '$"="\t";$F[-1]*=30;print "@F"' > $query.nodes.syncasm.$k1
eval $kmer2node -U -k$k2 $nodeseq1 $query.asm.fa -E /dev/stdout | egrep 'Edge|Node' \
| perl -lane '$"="\t";$F[-1]*=30;print "@F"' > $query.nodes.syncasm.$k2
eval $kmer2node -U -k$k3 $nodeseq1 $query.asm.fa -E /dev/stdout | egrep 'Edge|Node' \
| perl -lane '$"="\t";$F[-1]*=30;print "@F"' > $query.nodes.syncasm.$k3

eval $kmer2node -U -k$k1 $nodeseq1 $shred_fa -E /dev/stdout | egrep 'Edge|Node' > $query.nodes.$k1
eval $kmer2node -U -k$k2 $nodeseq2 $shred_fa -E /dev/stdout | egrep 'Edge|Node' > $query.nodes.$k2
eval $kmer2node -U -k$k3 $nodeseq3 $shred_fa -E /dev/stdout | egrep 'Edge|Node' > $query.nodes.$k3
merge_kmer2node.pl $gfa $query.nodes.* > $query.nodes

# Annotate the GFA with kmer counts
tag_gfa_km.pl $gfa < $query.nodes > $query.gfa

exit 0
