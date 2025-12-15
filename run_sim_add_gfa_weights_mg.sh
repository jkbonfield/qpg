#!/bin/bash
. ${CONFIG:-$QDIR/config_illumina.sh}

# Adds weights using minigraph to align and amend the GFA
gfa=$1
query=$2

# Shred the FASTA file into small bits, simulating sequencing
echo === Shredding genome $query with depth $shred_depth
shred_fa=$query.shred.fa
shred.pl -s 1 -l $shred_len -e $shred_err -d $shred_depth $query > $shred_fa

# Length specific minigraph options.  Long-read vs short-read
if [ $shred_len -gt 500 ]
then
    mg_opts="-x lr -w 15 -l 100 -d 100 -m 50,40 -n 5,5"
else
    mg_opts="-x sr"
fi
mg_opts="$mg_opts -j 0.01 --cov ${MINIGRAPH_OPTS}"


# Run minigraph to get a GFA
eval minigraph $mg_opts $gfa $shred_fa 2>mg.err > $query.mg

# Quantise dc:f field to kmer-count by multplying up the length.
# For edges we use it as EC direct, but in both cases it must be integer.
tag_gfa_mg.pl $query.mg > $query.gfa
