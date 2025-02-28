#!/bin/bash

# New args
query=$1
shred=$2
t=$3
idx=$4


echo === Evaluating result: GFA node sequence
truth=$query

echo "Try bwa in shell"
bwa index $truth 
bwa mem $truth $query.path_seq.$t.$idx
echo "------"

echo $truth $query.path_seq.$t.$idx
echo "start candidate stats.sh" 1>&2
candidate_stats.sh $truth $query.path_seq.$t.$idx | tee $query.eval_seq.$t.$idx

#--- Sample-specific consensus sequence by remapping back to candidate seq
# To get higher QUAL
echo === Evaluating result: Sample consensus
echo "samtools view" 1>&2
samtools view $shred -o $query.shred.fq

#bwa index $query.path_seq
#bwa mem -t8 $query.path_seq _ar.fastq |samtools sort -o _ar.bam
minimap2 -x lr:hq -a $query.path_seq.$t.$idx $query.shred.fq | samtools sort -o $query.bam.$t.$idx

samtools consensus -C0 $query.bam.$t.$idx -o $query.path_cons.$t.$idx
candidate_stats.sh $truth $query.path_cons.$t.$idx | tee $query.eval_cons.$t.$idx

#--- TODO: round 2: remap back to the consensus to improve further?
