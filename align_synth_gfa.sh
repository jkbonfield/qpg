#!/bin/sh

qdir=$HOME/work/quantum
pathfinder3=/lustre/scratch127/qpg/cz3/QuantumTangle/pathfinder/pathfinder
pathfinder=/nfs/users/nfs_j/jkb/lustre/quantum/pathfinder/pathfinder

query=$1
gfa=$2
nodeseq=$3

. $qdir/sim_path_config_hifi.sh

/bin/rm -f _ar.*

#$qdir/split_gfa_nodes.pl _ar.full.gfa 10000 > _ar.gfa
#gfa=_ar.gfa

echo Shredding genome $2
$qdir/shred.pl -s 1 -l $shred_len -e $shred_err -d $shred_depth $query > _ar.fasta

echo Mapping kmers
$qdir/kmer2node2 -U -k$kmer2node_kmer $nodeseq _ar.fasta | grep Node > _ar.nodes
#awk '/Node/ && $NF>=1' _ar.nodes
$qdir/tag_gfa.pl $gfa < _ar.nodes > _ar.gfa_tagged

echo Finding path
$pathfinder3 -a _ar.gfa_tagged 2>/dev/null | $qdir/pathfinder2seq.pl $gfa > _ar.called_seq

echo Evaluating result
$qdir/candidate_stats.pl $truth _ar.called_seq

ls -latr $query _ar*
#dotter $query _ar.called_seq

