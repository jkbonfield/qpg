#!/bin/bash

. ${CONFIG:-$QDIR/config_illumina.sh}

gfa=$1
query=$2

# Shred the FASTA file into small bits, simulating sequencing
echo === Shredding genome $query
shred_fa=$query.shred.fa
shred.pl -s 1 -l $shred_len -e $shred_err -d $shred_depth $query > $shred_fa

# Map shredded data to the GFA with GraphAligner
echo === Mapping with GraphAligner
eval /nfs/sam_scratch/jkb/conda22.old/bin/GraphAligner $ga_opts -g $gfa -f $shred_fa -x vg -a $shred_fa.gaf

# Annotate the GFA with kmer counts
tag_gfa_ga.pl $gfa 10 $shred_fa.gaf > $query.gfa
