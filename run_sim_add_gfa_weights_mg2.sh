#!/bin/bash

### TODO ###
# Use minigraph --vc and tag_gfa_ga.pl style processing to turn
# minigraph >s1>s2... paths into kmer counts, instead of using --cov
############

. ${CONFIG:-$QDIR/config_illumina.sh}

# Adds weights using minigraph to align and amend the GFA

gfa=$1
query=$2

# Shred the FASTA file into small bits, simulating sequencing
echo === Shredding genome $query
shred_fa=$query.shred.fa
shred.pl -s 1 -l $shred_len -e $shred_err -d $shred_depth $query > $shred_fa

# Length specific minigraph options.  Long-read vs short-read
if [ $shred_len -gt 500 ]
then
    mg_opts="-x lr -w 15 -l 100 -d 100 -m 50,40 -n 5,5"
else
    mg_opts="-x sr"
fi
mg_opts="$mg_opts -j 0.01 --vc ${MINIGRAPH_OPTS}"


# Run minigraph to get a GFA
echo === Aligning data
eval minigraph $mg_opts $gfa $shred_fa 2>/dev/null > $query.mg

# Annotate GFA
tag_gfa_ga.pl $gfa $query.mg > $query.gfa
