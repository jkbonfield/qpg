#!/bin/sh
. ${CONFIG:-./sim_path_config.sh}

# Path reconstruction with a simulated sequence and pathfinder

PATH=/software/badger/opt/pangenomes/bin/:$PATH
LD_LIBRARY_PATH=/software/badger/opt/pangenomes/lib
export LD_LIBRARY_PATH

seed=${1:-$$}
tmp=/tmp/sim.mdbg.$seed
if [ -e $tmp ]
then
  rm -rf $tmp
fi
mkdir $tmp

echo Using random seed $seed, output directory $tmp

# Create a genome.
# TODO: add options here to control genome size and complexity
echo Creating fake genome
./genome_create -s $seed -l $genome_len > $tmp/true.fa 2>>$tmp/stderr

# Shred genome
echo Shredding genome
./shred.pl -s 1 -l $shred_len -e $shred_err -d $shred_depth $tmp/true.fa > $tmp/shred.fa

# Assemble.
echo Assembling with mdbg
eval rust-mdbg $mdbg_opts --prefix $tmp/mdbg $tmp/shred.fa 2>>$tmp/stderr
gfa=$tmp/mdbg.gfa

# Populate GFA so sequence isn't "*".
echo /software/badger/opt/pangenomes/share/rust-mdbg/utils/magic_simplify $tmp/mdbg
/software/badger/opt/pangenomes/share/rust-mdbg/utils/magic_simplify $tmp/mdbg
gfa=$tmp/mdbg.msimpl.gfa

# Find a path and generate the candidate sequence
echo Pathfinding
eval $pathfinder $pathfinder_opts $gfa | ./pathfinder2seq.pl $gfa > $tmp/candidate.fa

# # Count fragments.
# # Also consider e.g. "-r 50,20k" (default 500,20k) to spot smaller gaps
# minimap2 --no-long-join --secondary=no $tmp/true.fa $tmp/candidate.fa \
#     > $tmp/minimap2_compare.paf
# wc -l $tmp/minimap2_compare.paf

# Report
echo === mdbg.msimpl
./candidate_stats.pl $tmp/true.fa $tmp/mdbg.msimpl.fa

echo "=== pathfinder (ususally poorer accuracy and low coverage)"
./candidate_stats.pl $tmp/true.fa $tmp/candidate.fa
