#!/bin/bash

help() {
    echo Usage: run_gfa_sim.sh [options] [seed solver [out_prefix]] 1>&2
    echo Options:
    echo "    -c,--config    FILE    Use FILE as configuration"
    echo "    -s,--seed      INT     Specificy random number seed"
    echo "    -p,--prefix    STR     Use STR as the output dir prefix [sim_]"
    echo "       --solver    STR     Specify the solver"
    echo "    -a,--annotate  STR     GFA node weight algorithm"
    echo "       --shred_len INT     Shotgun read length"
    echo "       --shred_err FLOAT   Shotgun read error rate (fraction)"
    echo ""
    echo "The old API is still supported with fixed argument order."
}

config=${CONFIG:-$QDIR/config_hifi_km.sh}
. $config

prefix="sim_"
seed=1
solver=pathfinder

while true
do
    case "$1" in
	'-h')
	    help
	    exit 0
	    ;;
	'-c'|'--config')
	    # Run now so we can override with other arguments
	    config=$2
	    . $config
	    shift 2
	    continue
	    ;;
	'-s'|'--seed')
	    seed=$2
	    shift 2
	    continue
	    ;;
	'--solver')
	    solver=$2
	    shift 2
	    continue
	    ;;
	'-p'|'--prefix')
	    prefix=$2
	    shift 2
	    continue
	    ;;
	'-a'|'--annotate')
	    annotate=$2
	    shift 2
	    continue
	    ;;
	'--shred-len')
	    shred_len=$2
	    shift 2
	    continue
	    ;;
	'--shred-err')
	    shred_err=0.001
	    shift 2
	    continue
	    ;;
	*)
	    break
	    ;;
    esac
done

# Old CLI syntax also used and overrides any --opt setting
if [ $# -gt 0 ]; then
    seed=$1
    shift
fi
if [ $# -gt 0 ]; then
    solver=$1
    shift
fi
if [ $# -gt 0 ]; then
    prefix=$1
    shift
fi

export QDIR=${QDIR:-`pwd`}
export PATH=$QDIR:$PATH
export CONFIG=$config; # still used in some other scripts

echo "Config:   $config"
echo "Solver:   $solver"
echo "Seed:     $seed"
echo "Annotate: $annotate"
echo "prefix:   $prefix"
echo ""

out_dir=`printf "$prefix%05d" $seed`
rm -rf $out_dir 2>/dev/null
mkdir $out_dir

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
    # Add weights to the GFA via minigraph
    # Creates:
    #     $i.gfa (primary output; annotated pop.gfa)
    #     $i.shred.fa
    #     $i.mg
    # For kmer2node also creates:
    #     $i.nodes.$k1
    #     $i.nodes.$k2
    #     $i.nodes.$k3
    #     $i.nodes
    eval run_sim_add_gfa_weights_${annotate}.sh pop.gfa $i $k1 $k2 $k3

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
head -1 `ls -1 *.eval_seq|head -1`
cat *.eval_seq|awk '!/Per/'
) 2>sim.err | tee sim.out
