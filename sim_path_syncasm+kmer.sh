#!/bin/sh
. ${CONFIG:-./config_illumina.sh}

# Path reconstruction with a simulated sequence and pathfinder

PATH=/software/badger/opt/pangenomes/bin/:$PATH
LD_LIBRARY_PATH=/software/badger/opt/pangenomes/lib
export LD_LIBRARY_PATH

seed=${1:-$$}
tmp=/tmp/sim.syncasm-kmer.$seed
if [ -e $tmp ]
then
  rm -rf $tmp
fi
mkdir $tmp

echo Using random seed $seed, output directory $tmp

# Create a genome.
# TODO: add options here to control genome size and complexity
echo Creating gake genome
./genome_create -s $seed -l $genome_len > $tmp/true.fa 2>>$tmp/stderr

# Shred genome
echo Shredding genome
./shred.pl -s 1 -l $shred_len -e $shred_err -d $shred_depth $tmp/true.fa > $tmp/shred.fa
./shred.pl -s 1 -l $shred_len -e 0 -d $shred_depth $tmp/true.fa > $tmp/shred_perfect.fa

# Assemble.
echo Assembling perfect data with syncasm
syncasm -k ${syncasm_kmer} $tmp/shred_perfect.fa -o $tmp/syncasm 2>>$tmp/stderr
gfa=$tmp/syncasm.utg.gfa

#gfatools asm  -t 10,50000 -b 100000 -t 10,100000  -b 1000000 -t 10,200000 -b 1\
#000000 -u $tmp/syncasm.utg.gfa > $tmp/syncasm.gfatools.gfa
#gfa=$tmp/syncasm.gfatools.gfa

echo Kmer indexing and matching of perfect graph vs erroneous reads
k1=${KMER1:-$kmer2node_kmer}

./gfa2nodeseq.pl $gfa $k1 > $tmp/kmer$k1.nodeseq

eval ${KMER2NODE:-./kmer2node2} -k $k1 -K $k1 $tmp/kmer$k1.nodeseq $tmp/shred.fa |grep Node > $tmp/kmer$k1.nodes
cat $tmp/kmer$k1.nodes

cp $tmp/kmer$k1.nodes $tmp/kmers.nodes
./tag_gfa.pl $gfa < $tmp/kmers.nodes > $gfa.tagged
gfa=$gfa.tagged

# Find a path and generate the candidate sequence
echo Pathfinding
eval $pathfinder $pathfinder_opts $gfa | ./pathfinder2seq.pl $gfa > $tmp/candidate.fa

# Report
./candidate_stats.pl $tmp/true.fa $tmp/candidate.fa
