#!/bin/sh
. ${CONFIG:-./config_illumina.sh}

# Path reconstruction with a simulated sequence and pathfinder

PATH=/software/badger/opt/pangenomes/bin/:$PATH
LD_LIBRARY_PATH=/software/badger/opt/pangenomes/lib
export LD_LIBRARY_PATH

seed=${1:-$$}
tmp=/tmp/sim.miniasm.$seed
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
# TODO: options for kmer size 
echo Assembling with minimap2 / miniasm
eval minimap2 $minimap2_opts -t8 $tmp/shred.fa $tmp/shred.fa \
| gzip -1 > $tmp/miniasm.paf.gz
#1: POOR
# cov     ncont   nbreak  ind diff   id
# 53209.3 10.1818 12.2727 1 0.454545 99.5718 0
#miniasm -s1000 -f $tmp/shred.fa $tmp/miniasm.paf.gz > $tmp/miniasm.gfa

#2: OK
# good    good    =       poor    ok       ~=
# 57205.9 8.18182 12.5455 1.81818 0.545455 99.5845 0
#miniasm -s1000 -c5 -m250 -e10  -f $tmp/shred.fa $tmp/miniasm.paf.gz > $tmp/miniasm.gfa

# increasing -s dramatically decreases nbreak and increases ncontig
# ~       xx      yyy     y        x       ~
# 57880.5 20.5455 2.72727 0.909091 1.27273 99.5155 0
#miniasm -s1500 -c5 -m250 -e10  -f $tmp/shred.fa $tmp/miniasm.paf.gz > $tmp/miniasm.gfa

# x       =       =       y        y        y
# 54527.2 8.18182 12.1818 0.909091 0.454545 99.62 0
#miniasm -s1000 -c5 -m200 -e10  -f $tmp/shred.fa $tmp/miniasm.paf.gz > $tmp/miniasm.gfa

# increasing -m gives better coverage. Decreasing lower. 250 good middle
# y     y x       x x        x
# 58414 8 13.4545 2 0.636364 99.4945 0
#miniasm -s1000 -c5 -m300 -e10  -f $tmp/shred.fa $tmp/miniasm.paf.gz > $tmp/miniasm.gfa

# increasing -c (min depth) is better cov/break, but more contig
# yy      xxx     yy      y       x        y
# 60116.8 16.3636 8.09091 1.54545 0.727273 99.6136 0
#miniasm -s1000 -c10 -m250 -e10  -f $tmp/shred.fa $tmp/miniasm.paf.gz > $tmp/miniasm.gfa

# Bigger -c as above, but compensate with lower -s; not worth it
# x       x       y       =       x        x
# 55943.8 9.18182 10.0909 1.81818 0.636364 99.5527 0
#miniasm -s750 -c10 -m250 -e10  -f $tmp/shred.fa $tmp/miniasm.paf.gz > $tmp/miniasm.gfa

# 10: Opposite: lower -c, bigger -s.  Good
# yyy     y x       y       x x
# 62415.3 7 13.6364 1.45455 1 99.5264 0
#miniasm -s1250 -c3 -m250 -e10  -f $tmp/shred.fa $tmp/miniasm.paf.gz > $tmp/miniasm.gfa

# vs 10: decrease -m is slightly more accurate but lower coverage.
# x     y       ~       x       y        ~
# 60411 6.45455 13.1818 1.54545 0.909091 99.5218 0
# miniasm -s1250 -c3 -m200 -e10  -f $tmp/shred.fa $tmp/miniasm.paf.gz > $tmp/miniasm.gfa

# vs 10; lower -e is better coverage but a bit less accuracy
# ~y      x       x       x       y        y
# 62986.6 7.36364 14.4545 2.09091 0.818182 99.5891 0
#miniasm -s1250 -c3 -m250 -e7  -f $tmp/shred.fa $tmp/miniasm.paf.gz > $tmp/miniasm.gfa

# vs 10; bigger -e is better accuracy, fewer contigs, but a bit less coverage
# xx      yy      y       =       y        ~y
# 60908.3 6.63636 12.8182 1.45455 0.818182 99.5391 0
#miniasm -s1250 -c3 -m250 -e15  -f $tmp/shred.fa $tmp/miniasm.paf.gz > $tmp/miniasm.gfa

# 11: Bigger -s is again much better coverage and accuracy, more more contig
# Accuracy matters more I think
# yyy     xxx     yyy     y      y         y
# 70130.4 16.5455 5.27273 1.27273 0.727273 99.6691 0
#miniasm -s1500 -c3 -m250 -e15  -f $tmp/shred.fa $tmp/miniasm.paf.gz > $tmp/miniasm.gfa

# 12: vs 11: -e15 vs -e10 is better coverage/accuracy, a few more contigs again
# yy      xx      y       y       yy       y
# 73036.2 18.9091 5.09091 1.18182 0.545455 99.6891 0
# miniasm -s1500 -c3 -m250 -e10  -f $tmp/shred.fa $tmp/miniasm.paf.gz > $tmp/miniasm.gfa

# vs 11:  -e15 vs -e5 is better coverage/accuracy, a few more contigs again
# But it's starting to get more breaks again.
# yyy     xxx     x       x       y       y
# 75395.4 21.7273 5.63636 1.36364 0.636364 99.6873 0
#miniasm -s1500 -c3 -m250 -e5  -f $tmp/shred.fa $tmp/miniasm.paf.gz > $tmp/miniasm.gfa

# yy      xx      =       =       yy       y
# 72791.6 18.6364 5.27273 1.27273 0.545455 99.6836 0
#miniasm -s1500 -c3 -m300 -e10  -f $tmp/shred.fa $tmp/miniasm.paf.gz > $tmp/miniasm.gfa

# vs 12, -g500 (vs1000).  Identical
# vs 12, reducing -o1000 (was as -s). Identical
# vs 12, -h500.  Identical

# end-to-end match ratio (default -I0.8); minimal cov vs ncontig. Not a win
# xx      y       y       x       =        ~
# 68324.1 17.8182 4.90909 1.27273 0.545455 99.6827 0
# miniasm -s1500 -I0.7 -c3 -m250 -e10  -f $tmp/shred.fa $tmp/miniasm.paf.gz > $tmp/miniasm.gfa

# lower -F; aggr drop ratio (default -F0.8)
# ~x      =       x       x       y        x
# 72212.1 18.9091 5.36364 1.27273 0.454545 99.6991 0
#miniasm -s1500 -F0.7 -c3 -m250 -e10  -f $tmp/shred.fa $tmp/miniasm.paf.gz > $tmp/miniasm.gfa

# Higher -F
# y       ~x y       y       =        y
# 73288.4 19 4.81818 1.09091 0.545455 99.7055 0
#miniasm -s1500 -F0.85 -c3 -m250 -e10  -f $tmp/shred.fa $tmp/miniasm.paf.gz > $tmp/miniasm.gfa
# 13:
# y       ~x      y       y x        ~x
# 73301.4 19.0909 4.63636 1 0.636364 99.6809 0
#miniasm -s1500 -F0.9 -c3 -m250 -e10  -f $tmp/shred.fa $tmp/miniasm.paf.gz > $tmp/miniasm.gfa

# 14: Big -m increase
# y       y       y       y        =        ~  vs 12
# ~       y       =       y        y        ~  vs 13
# 73239.7 18.8182 4.63636 0.727273 0.545455 99.6636 0
eval miniasm $miniasm_opts  -f $tmp/shred.fa $tmp/miniasm.paf.gz > $tmp/miniasm.gfa

# Even bigger -m, plus explicit -o (doc as same as -s)
# y       x       yy      y        x        ~
# 74432.1 19.3636 3.63636 0.636364 0.636364 99.62 0
#miniasm -s1500 -o1500 -F0.9 -c3 -m1000 -e10  -f $tmp/shred.fa $tmp/miniasm.paf.gz > $tmp/miniasm.gfa

# -c3 to -c5: POOR (but very accurate)
# 59400.8 21 1.54545 0.545455 1.09091 99.5491 0
# miniasm -s1500 -o1500 -F0.9 -c5 -m1000 -e10  -f $tmp/shred.fa $tmp/miniasm.paf.gz > $tmp/miniasm.gfa

# -c3 to -c2: many fewer contigs, but less coverage and accuracy
# xx      yyy     xx      x        x        ~x
# 70417.5 12.0909 6.72727 0.818182 0.727273 99.6373 0
#miniasm -s1500 -o1500 -F0.9 -c2 -m1000 -e10  -f $tmp/shred.fa $tmp/miniasm.paf.gz > $tmp/miniasm.gfa

# -e10 to -e5; too many contigs.  Try -e15 again?
# y       xx      y       =        x        ~y
# 74692.4 21.8182 4.54545 0.727273 0.636364 99.6664 0   vs 14
#miniasm -s1500 -F0.9 -c3 -m500 -e5  -f $tmp/shred.fa $tmp/miniasm.paf.gz > $tmp/miniasm.gfa
# xx      yy      =       =        x       ~x
# 71126.7 16.8182 4.63636 0.727273 0.727273 99.6236 0
#miniasm -s1500 -F0.9 -c3 -m500 -e15  -f $tmp/shred.fa $tmp/miniasm.paf.gz > $tmp/miniasm.gfa

# -n3 vs -n2.  Identical to 14.

gfa=$tmp/miniasm.gfa


# Find a path and generate the candidate sequence
echo Pathfinding
eval $pathfinder $pathfinder_opts $gfa | ./pathfinder2seq.pl $gfa > $tmp/candidate.fa

# Report
./candidate_stats.pl $tmp/true.fa $tmp/candidate.fa
