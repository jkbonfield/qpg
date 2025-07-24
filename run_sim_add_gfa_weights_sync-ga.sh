#!/bin/bash

. ${CONFIG:-$QDIR/config_illumina.sh}

gfa=$1
query=$2

# Shred the FASTA file into small bits, simulating sequencing
echo === Shredding genome $query
shred_fa=$query.shred.fa
shred.pl -s 1 -l $shred_len -e $shred_err -d $shred_depth $query > $shred_fa

# Assemble and add these to the shredded input so we map both
eval syncasm $syncasm_opts $shred_fa
awk '/^S/ {printf(">%s\n",$2);print $3}' syncasm.asm.utg.gfa > $query.asm.fa

(cat $shred_fa; for i in `seq 1 30`; do cat $query.asm.fa; done) > $query.frags.fa
#cat $shred_fa $query.asm.fa > $query.frags.fa
frags_fa=$query.frags.fa

# Map shredded data to the GFA with GraphAligner
echo === Mapping with GraphAligner
eval GraphAligner $ga_opts -g $gfa -f $frags_fa -x vg -a $shred_fa.gaf

# Annotate the GFA with kmer counts
tag_gfa_ga.pl $gfa 10 $shred_fa.gaf > $query.gfa
