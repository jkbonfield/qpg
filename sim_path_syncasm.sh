#!/bin/sh

# Path reconstruction with a simulated sequence and pathfinder

pathfinder=/lustre/scratch127/qpg/cz3/QuantumTangle/pathfinder/pathfinder

PATH=/software/badger/opt/pangenomes/bin/:$PATH
LD_LIBRARY_PATH=/software/badger/opt/pangenomes/lib
export LD_LIBRARY_PATH

seed=${1:-$$}
tmp=/tmp/sim.syncasm.$seed
if [ -e $tmp ]
then
  rm -rf $tmp
fi
mkdir $tmp

echo Using random seed $seed, output directory $tmp

# Create a genome.
# TODO: add options here to control genome size and complexity
echo Creating gake genome
./genome_create -s $seed > $tmp/true.fa 2>>$tmp/stderr

# Shred genome
echo Shredding genome
./shred.pl -s 1 -l 2000 -e 1e-5 -d 30 $tmp/true.fa > $tmp/shred.fa

# Assemble.
# TODO: options for kmer size 
echo Assembling with syncasm
syncasm -k 301 $tmp/shred.fa -o $tmp/syncasm 2>>$tmp/stderr
gfa=$tmp/syncasm.utg.final.gfa

# Find a path and generate the candidate sequence
echo Pathfinding
$pathfinder $gfa | ./pathfinder2seq.pl $gfa > $tmp/candidate.fa

# Compare
echo Compare /tmp/true.fa /tmp/candidate.fa
echo dotter /tmp/sim.1/true.fa /tmp/sim.1/candidate.fa

# Count fragments
minimap2  --no-long-join --secondary=no $tmp/true.fa $tmp/candidate.fa \
    > $tmp/minimap2_compare.paf
wc -l $tmp/minimap2_compare.paf
