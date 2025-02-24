#!/bin/sh

qdir=$HOME/work/quantum

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
#$qdir/kmer2node3 -U -k$kmer2node_kmer $nodeseq _ar.fasta | grep Node > _ar.nodes
#$qdir/kmer2node3 -D -k$kmer2node_kmer $nodeseq _ar.fasta | grep Node > _ar.nodes

$qdir/tag_gfa.pl $gfa < _ar.nodes > _ar.gfa_tagged

echo Finding path
$pathfinder_jkb.jkb2 -C20  -a _ar.gfa_tagged 2>/dev/null | tee _ar.gfa_pf | $qdir/pathfinder2seq.pl $gfa > _ar.called_seq
#$pathfinder_jkb.jkb1  -a _ar.gfa_tagged 2>/dev/null | tee _ar.gfa_pf | $qdir/pathfinder2seq.pl $gfa > _ar.called_seq
#$pathfinder_cz3 _ar.gfa_tagged 2>/dev/null | tee _ar.gfa_pf | $qdir/pathfinder2seq.pl $gfa > _ar.called_seq
cat _ar.gfa_pf | sed -n '/PATH/,$p'|awk '{printf("%s ",$3)} END {print "\n"}'

echo Evaluating result
sed 's/^>.*/>truth.fa/' $query > _truth.fa
$qdir/candidate_stats.pl _truth.fa _ar.called_seq

echo Realignment evaluation
# To get higher QUAL
samtools view _ar.fasta -o _ar.fastq
#bwa index _ar.called_seq
#bwa mem -t8 _ar.called_seq _ar.fastq |samtools sort -o _ar.bam
minimap2 -x lr:hq -a _ar.called_seq _ar.fastq | samtools sort -o _ar.bam
samtools consensus -C0 _ar.bam -o _ar.cons
$qdir/candidate_stats.pl _truth.fa _ar.cons

ls -latr $query _ar*
#dotter $query _ar.called_seq

