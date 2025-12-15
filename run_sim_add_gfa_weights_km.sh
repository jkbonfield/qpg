#!/bin/bash

gfa=$1
query=$2

. ${CONFIG:-$QDIR/config_illumina.sh}
kmers=${kmers:-75 35 20}

# Build the nodeseq indices
set -- $kmers; k1=$1
if [ ! -e $gfa.ns$k1 ]
then
    echo "=== Creating nodeseq files, kmers $kmers"
    GraphAligner -g $gfa -f train.fa -x vg -a pop.gaf >pop.GraphAligner.out
    gaf2nodeseq2.pl pop.gaf train.fa $gfa $k1 > $gfa.ns$k1
fi

# Shred the FASTA file into small bits, simulating sequencing
echo === Shredding genome $query with depth $shred_depth
shred_fa=$query.shred.fa
shred.pl -s 1 -l $shred_len -e $shred_err -d $shred_depth $query > $shred_fa

# Map kmers in the shredded fasta to the graph nodeseq
echo === Mapping kmers: k=$kmers
for k in $kmers
do
    eval kmer2node4 -G $gfa -U -k$k -K$k1 $gfa.ns$k1 $shred_fa -E /dev/stdout | egrep 'Edge|Node' > $query.nodes.$k
done
merge_kmer2node.pl $gfa $query.nodes.* > $query.nodes

# Annotate the GFA with kmer counts
tag_gfa_km.pl $gfa < $query.nodes > $query.gfa

exit 0
