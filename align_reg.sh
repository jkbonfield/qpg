#!/bin/sh -x

qdir=$HOME/work/quantum
pathfinder3=/lustre/scratch127/qpg/cz3/QuantumTangle/pathfinder/pathfinder
pathfinder=/nfs/users/nfs_j/jkb/lustre/quantum/pathfinder/pathfinder

file=$1
region=$2
gfa=$3
k=51

/bin/rm -f _ar.*

#gfatools view -R $region $gfa > _ar.gfa
gfatools view -R $region $gfa > _ar.full.gfa
$qdir/split_gfa_nodes.pl _ar.full.gfa 10000 > _ar.gfa
#cp _ar.full.gfa _ar.gfa

samtools view --write-index -o _ar.bam $file $region
samtools fasta _ar.bam > _ar.fasta
samtools consensus -r $region -C0 _ar.bam -o _ar.cons

#$qdir/gaf2nodeseq.pl /dev/null /dev/null _ar.gfa $k > _ar.nodeseq
$qdir/gfa2nodeseq.pl _ar.gfa $k > _ar.nodeseq

$qdir/kmer2node2 -U -k$k _ar.nodeseq _ar.fasta | grep Node > _ar.nodes
awk '/Node/ && $NF>=1' _ar.nodes

$qdir/tag_gfa.pl _ar.gfa < _ar.nodes > _ar.gfa_tagged
$pathfinder -a _ar.gfa_tagged 2>/dev/null | $qdir/pathfinder2seq.pl _ar.gfa > _ar.called_seq
$qdir/candidate_stats.pl _ar.cons _ar.called_seq

ls -latr _ar*
dotter _ar.cons _ar.called_seq

