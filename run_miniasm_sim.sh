#!/bin/bash

if [ $# -lt 2 ]
then
    echo Usage: run_miniasm_sim.sh seed [out_prefix] 1>&2
    echo Set CONFIG environment variable to change configuration. 1>&2
    exit 1
fi

seed=$1
prefix=${2:-sim_}
export QDIR=${QDIR:-`pwd`}
export PATH=$QDIR:$PATH

. ${CONFIG:-$QDIR/config_illumina.sh}

out_dir=`printf "$prefix%05d" $seed`
rm -rf $out_dir 2>/dev/null
mkdir $out_dir

cd $out_dir

(
# Create fake genomes
# Even though we don't need it, we use run_sim_create_gfa.sh here so
# we get a consistent list of test genomes.  TODO: hive this off to
# its own script for efficiency.
export use_syncasm=1
run_sim_create_gfa.sh $seed 10 10

## Creates seq_* and fofn.test
#eval genome_create $genome_opts -s $seed | tail -20 | \
#perl -lane 'if (/>/) {s/>//;close(FH),open(FH,">$_");print FH ">$_";next} {print FH $_} END {close(FH)}'
#ls -1 seq_* | tail -10 > fofn.test

# Foreach test genome, assemble and eval
for i in `cat fofn.test`
do
    # Shred
    echo === $i ===
    shred_fa=$i.shred.fa
    shred.pl -s 1 -l $shred_len -e $shred_err -d $shred_depth $i > $shred_fa

    # Assemble and turn GFA to contigs fasta
    minimap2_opts=${minimap2_opts:-"-k19 -Xw5 -m100"}
    #miniasm_opts=${miniasm_opts:-"-s100 -F0.9 -c3 -m100 -e10"}
    # Higher -f is fewer contigs and better %used, but worse covered%
    # Higher -k is                                     poor  covered%
    # Higher -k is more  contigs and better %used, and more  covered%
    # Higher -w is fewer contigs and better %used, but worse covered%
    minimap2_opts="-k17 -w7 -X -m50 -f 0.0005"
    # slower, but maybe a bit better
    minimap2_opts="-k17 -w7 -X -m50 -f 0.001"

    # Higher -F is fewer contigs and better %used, but slightly worse covered%
    # Lower  -m is                   poorer %used, poorer covered%?
    # Higher -m is                   BETTER %used, poorer covered%
    miniasm_opts="-s100 -F0.95 -c3 -m120 -e10 -d100 -h10"
    eval minimap2 $minimap2_opts -t8 $shred_fa $shred_fa | gzip -1 > $i.shred.paf.gz
    eval miniasm $miniasm_opts -f $shred_fa $i.shred.paf.gz > $i.asm.gfa
    awk '/^S/ {printf(">%s\n",$2);print $3}' $i.asm.gfa > $i.asm.fa

    # Evaluate
    candidate_stats.sh $i $i.asm.fa > $i.eval_seq
done

# Summary
echo
echo === Summary ===
head -1 `ls -1 *.eval_seq|head -1`
cat *.eval_seq|awk '!/Per/'
) 2>sim.err | tee sim.out
