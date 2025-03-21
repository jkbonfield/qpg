#!/bin/bash

if [ $# -lt 2 ]
then
    echo Usage: run_gfa_sim.sh seed solver [out_prefix] 1>&2
    exit 1
fi

seed=$1
solver=$2
prefix=${3:-sim_}
export QDIR=${QDIR:-`pwd`}
export PATH=$QDIR:$PATH

. ${CONFIG:-$QDIR/sim_path_config_hifi.sh}

out_dir=`printf "$prefix%05d" $seed`
rm -rf $out_dir 2>/dev/null
mkdir $out_dir

#k1=75
#k2=50
#k3=35

#k1=75
#k2=35
#k3=15

k1=75
k2=35
k3=20

cd $out_dir

(
# Create fake genomes and use the training set to create a pangenome
# Creates:
#    pop.gfa
#    pop.gfa.ns$k1
#    pop.gfa.ns$k2
#    pop.gfa.ns$k3
#    fofn.test
#    fofn.train
run_sim_create_gfa.sh $seed $k1 $k2 $k3

# Foreach test genome, not used in pangenome creation, find path and eval
for i in `cat fofn.test`
do
    if [ "x$use_mg" = "x1" ]
    then
        # Add weights to the GFA via minigraph
        # Creates:
        #     $i.gfa (primary output; annotated pop.gfa)
        #     $i.shred.fa
        #     $i.mg
        run_sim_add_gfa_weights_mg.sh pop.gfa $i
    else
        # Add weights to the GFA via kmer2node.
        # Creates:
        #     $i.gfa (primary output; annotated pop.gfa)
        #     $i.shred.fa
        #     $i.nodes.$k1
        #     $i.nodes.$k2
        #     $i.nodes.$k3
        #     $i.nodes
        run_sim_add_gfa_weights.sh pop.gfa $i $k1 $k2 $k3
    fi


    # Find a path
    # Creates:
    #     $i.path
    echo === Running solver $solver
    run_sim_solver_$solver.sh $i.gfa > $i.path
    sed -n '/PATH/,$p' $i.path \
	| awk '/^\[/ {printf("%s ",$3)} END {print "\n"}' 1>&2

    # Evaluate the path
    # Creates:
    #     $i.path_seq
    #     $i.bam
    #     $i.path_cons
    #     $i.eval_seq
    #     $i.eval_cons
    run_sim_evaluate_path.sh pop.gfa $i.path $i $i.shred.fa
done

# Summary
echo
echo === Summary ===
head -1 `ls -1 *.eval_cons|head -1`
cat *.eval_cons|awk '!/Per/'
) 2>sim.err | tee sim.out
