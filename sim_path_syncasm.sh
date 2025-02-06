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
# 96282.8 1 7.45455 35 0.454545 98.4218 0
#syncasm -k 301 $tmp/shred.fa -o $tmp/syncasm 2>>$tmp/stderr

# y       = xx      xx      =        y
# 97152.4 1 10.6364 48.3636 0.454545 98.5627 0
#syncasm -k 201 $tmp/shred.fa -o $tmp/syncasm 2>>$tmp/stderr

# xxx     xx      yy      yyy     yyy       yy
# 70675.6 2.90909 4.45455 11.9091 0.0909091 99.3755 0
#syncasm -k 501 $tmp/shred.fa -o $tmp/syncasm 2>>$tmp/stderr

# opt2:
# y       y yy      yy      yy       y
# 96633.3 1 6.81818 26.1818 0.272727 98.9536 0
syncasm -k 401 $tmp/shred.fa -o $tmp/syncasm 2>>$tmp/stderr

# poor coverage
# 82349.5 1.36364 5.09091 19.8182 0.363636 98.96 0
# syncasm -k 451 $tmp/shred.fa -o $tmp/syncasm 2>>$tmp/stderr

# xxx     xxxx    yyyy    xxxx    yy yy
# 59127.9 22.7273 1.36364 5.54545 0 99.8227 0
#syncasm -k 401 -c20 $tmp/shred.fa -o $tmp/syncasm 2>>$tmp/stderr
# x       ~x      =       =       =       ~x
# 95572.8 1.09091 6.81818 26.1818 0.272727 98.9455 0
#syncasm -k 401 -c10 $tmp/shred.fa -o $tmp/syncasm 2>>$tmp/stderr

# xx      = y       y       x        ~x
# 92859.5 1 6.63636 22.1818 0.363636 98.87 0
#syncasm -k 401 -s25 $tmp/shred.fa -o $tmp/syncasm 2>>$tmp/stderr

# Lower/higher arc cov (def 0.35)
# yy      = =       x       x        ~x
# 97777.2 1 6.81818 27.4545 0.363636 98.8964 0
#syncasm -k 401 -a 0.25 $tmp/shred.fa -o $tmp/syncasm 2>>$tmp/stderr
# xxx     x       y yy      x        ~
# 88016.8 1.63636 6 23.3636 0.363636 98.9455 0
# syncasm -k 401 -a 0.45 $tmp/shred.fa -o $tmp/syncasm 2>>$tmp/stderr

# --max-tip 10000: identical
# --max-tip 1000:  identical

#gfa=$tmp/syncasm.utg.final.gfa
gfa=$tmp/syncasm.utg.gfa

#gfatools asm  -t 10,50000 -b 100000 -t 10,100000  -b 1000000 -t 10,200000 -b 1\
#000000 -u $tmp/syncasm.utg.gfa > $tmp/syncasm.gfatools.gfa
#gfa=$tmp/syncasm.gfatools.gfa

# Find a path and generate the candidate sequence
echo Pathfinding
$pathfinder $gfa | ./pathfinder2seq.pl $gfa > $tmp/candidate.fa

# Report
./candidate_stats.pl $tmp/true.fa $tmp/candidate.fa
