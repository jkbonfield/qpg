#!/bin/sh
qdir=${QDIR:-.}

query=$1
gfa=$2
nodeseq1=$3
nodeseq2=$4
nodeseq3=$5
k1=75
k2=50
k3=35

. $qdir/sim_path_config_hifi.sh

/bin/rm -f _ar.*

shred_len=150

#$qdir/split_gfa_nodes.pl _ar.full.gfa 10000 > _ar.gfa
#gfa=_ar.gfa

echo Shredding genome $2
$qdir/shred.pl -s 1 -l $shred_len -e $shred_err -d $shred_depth $query > _ar.fasta

echo Mapping kmers: Note kmer2node2 vs kmer2node3 usage
$qdir/kmer2node2 -U -k$k1 $nodeseq1 _ar.fasta | grep Node > _ar.nodes.$k1
$qdir/kmer2node2 -U -k$k2 $nodeseq2 _ar.fasta | grep Node > _ar.nodes.$k2
$qdir/kmer2node2    -k$k3 $nodeseq3 _ar.fasta | grep Node > _ar.nodes.$k3
$qdir/merge_kmer2node.pl _ar.nodes.$k1 _ar.nodes.$k2 _ar.nodes.$k3 > _ar.nodes

#$qdir/kmer2node3 -U -k$kmer2node_kmer $nodeseq _ar.fasta | grep Node > _ar.nodes
#$qdir/kmer2node3 -D -k$kmer2node_kmer $nodeseq _ar.fasta | grep Node > _ar.nodes
#awk '/Node/ && $NF>=1' _ar.nodes
$qdir/tag_gfa.pl $gfa < _ar.nodes > _ar.gfa_tagged

echo Finding path
# -a seems routinely to produce shorter paths
#$pathfinder_cz3 _ar.gfa_tagged 2>/dev/null | tee _ar.gfa_pf | $qdir/pathfinder2seq.pl $gfa > _ar.called_seq
$pathfinder_jkb -C40  _ar.gfa_tagged 2>/dev/null | tee _ar.gfa_pf | $qdir/pathfinder2seq.pl $gfa > _ar.called_seq
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

