#!/bin/bash

. ${CONFIG:-$QDIR/config_illumina.sh}

gfa=$1
query=$2

# Shred the FASTA file into small bits, simulating sequencing
echo === Shredding genome $query
shred_fa=$query.shred.fa
shred.pl -s 1 -l $shred_len -e $shred_err -d $shred_depth $query > $shred_fa

# Assemble and add these to the shredded input so we map both
minimap2_opts="-k17 -w7 -X -m50 -f 0.001"
miniasm_opts="-s100 -F0.95 -c3 -m120 -e10 -d100 -h10"
eval minimap2 $minimap2_opts -t8 $shred_fa $shred_fa | gzip -1 > $query.shred.paf.gz
eval miniasm $miniasm_opts -f $shred_fa $query.shred.paf.gz > $query.asm.gfa
awk '/^S/ {printf(">%s\n",$2);print $3}' $query.asm.gfa > $query.asm.fa

(cat $shred_fa; for i in `seq 1 30`; do cat $query.asm.fa | sed 's/>.*/&'_$i'/'; done) > $query.frags.fa
#cat $shred_fa $query.asm.fa > $query.frags.fa
frags_fa=$query.frags.fa

# Map shredded data to the GFA with GraphAligner
echo === Mapping with GraphAligner
eval GraphAligner $ga_opts -g $gfa -f $frags_fa -x vg -a $shred_fa.gaf

# Annotate the GFA with kmer counts
tag_gfa_ga.pl $gfa 10 $shred_fa.gaf > $query.gfa
