#!/bin/bash

if [ $# -lt 2 ]
then
    echo Usage: run_syncasm_sim.sh seed [out_prefix] 1>&2
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
run_sim_create_gfa.sh $seed 10

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
    eval syncasm $syncasm_opts $shred_fa
    awk '/^S/ {printf(">%s\n",$2);print $3}' syncasm.asm.utg.gfa > $i.asm.fa

    # Evaluate
    candidate_stats.sh $i $i.asm.fa > $i.eval_seq
done

# Summary
echo
echo === Summary ===
head -1 `ls -1 *.eval_seq|head -1`
cat *.eval_seq|awk '!/Per/'
) 2>sim.err | tee sim.out
