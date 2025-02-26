#!/bin/bash

gfa=$1
path=$2
query=$3
shred=$4

#--- Candidate sequence from the GFA
# Turn the path into a candidate sequence
# FIXME: for now we expect pathfinder output format
pathfinder2seq.pl $gfa $path > $query.path_seq

echo === Evaluating result: GFA node sequence
truth=$query

echo $truth $query.path_seq
candidate_stats.sh $truth $query.path_seq | tee $query.eval_seq

#--- Sample-specific consensus sequence by remapping back to candidate seq
# To get higher QUAL
echo === Evaluating result: Sample consensus
samtools view $shred -o $query.shred.fq

#bwa index $query.path_seq
#bwa mem -t8 $query.path_seq _ar.fastq |samtools sort -o _ar.bam
minimap2 -x lr:hq -a $query.path_seq $query.shred.fq | samtools sort -o $query.bam

samtools consensus -C0 $query.bam -o $query.path_cons
candidate_stats.sh $truth $query.path_cons | tee $query.eval_cons

#--- TODO: round 2: remap back to the consensus to improve further?
