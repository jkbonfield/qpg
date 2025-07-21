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
vg giraffe ${GIRAFFE_OPTS:-} -Z $gfa.giraffe.gbz --named-coordinates -f $shred_fa -o gaf > $shred_fa.gaf

# Annotate the GFA with kmer counts
tag_gfa_ga.pl $gfa 10 $shred_fa.gaf > $query.gfa


# --annotate vg with GIRAFFE_OPTS=?, using pathfinder
# -M1  10015.0 10127.5    91.8    91.5     1.6     2.3     1.5     0.1    96.9
# -M2  10015.0 10091.7    91.7    91.6     1.7     2.2     1.5     0.1    97.0
# -M5  10015.0 10091.7    91.7    91.6     1.7     2.2     1.5     0.1    97.0
# trim 10015.0  9891.2    89.8    90.2     1.6     2.1     3.1     0.1    95.4
# (M1)
