#!/bin/bash

. ${CONFIG:-$QDIR/config_illumina.sh}

gfa=$1
query=$2

# Shred the FASTA file into small bits, simulating sequencing
echo === Shredding genome $query
shred_fa=$query.shred.fa
shred.pl -s 1 -l $shred_len -e $shred_err -d $shred_depth $query > $shred_fa

# Map shredded data to the GFA with GraphAligner
echo === Mapping with Giraffe

# See https://github.com/vgteam/vg/issues/4160

if [ ! -e $gfa.giraffe.gbz ]
then
    vg autoindex --workflow giraffe --prefix $gfa -g $gfa
fi
vg giraffe -Z $gfa.giraffe.gbz --named-coordinates -f $shred_fa -o gaf > $shred_fa.gaf

# Annotate the GFA with kmer counts
tag_gfa_ga.pl $gfa $shred_fa.gaf > $query.gfa
