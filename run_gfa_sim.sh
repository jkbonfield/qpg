#!/bin/bash
help() {
    echo Usage: run_gfa_sim.sh [options] [seed solver [out_prefix]] 1>&2
    echo Options:
    echo "    -c,--config    FILE    Use FILE as configuration"
    echo "    -s,--seed      INT     Specify random number seed [1]"
    echo "    -p,--prefix    STR     Use STR as the output dir prefix [sim_]"
    echo "       --solver    STR     Specify the solver [pathfinder]"
    echo "    -a,--annotate  STR     GFA node weight algorithm [km]"
    echo "       --shred_len INT     Shotgun read length"
    echo "       --shred_err FLOAT   Shotgun read error rate (fraction)"
    echo "    -t,--times     INT_LIST   Time limits provided to QUBO solvers"
    echo "    -j,--jobs      INT     Number of runs of QUBO solvers"
    echo "    -n,--training  INT     Number of strings to use as training set [10]"
    echo ""
    echo "The old API is still supported with fixed argument order."
}

config=${CONFIG:-$QDIR/config_hifi_km.sh}
. $config

prefix="sim_"
seed=1
solver=pathfinder
num_training=10

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
	'-t'|'--times')
	    time_limits=$2
	    shift 2
	    continue
	    ;;
	'-j'|'--jobs')
	    num_jobs=$2
	    shift 2
	    continue
	    ;;
	'-n'|'--training')
	    num_training=$2
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


export QDIR=${QDIR:-$(pwd)}
export PATH=$QDIR:$PATH
export CONFIG=$config; # still used in some other scripts

echo "Config:   $config"
echo "Solver:   $solver"
echo "Seed:     $seed"
echo "Annotate: $annotate"
echo "prefix:   $prefix"
echo ""

out_dir=$(printf "$prefix%05d" "$seed")
echo "$out_dir"
rm -rf "$out_dir" 2>/dev/null
mkdir "$out_dir"

k1=75
k2=35
k3=20

cd "$out_dir" || exit 1

(
# Create fake genomes and use the training set to create a pangenome
# Creates:
#    pop.gfa
#    pop.gfa.ns$k1
#    pop.gfa.ns$k2
#    pop.gfa.ns$k3
#    fofn.test
#    fofn.train
run_sim_create_gfa.sh "$seed" $k1 $k2 $k3 "$num_training"

# Foreach test genome, not used in pangenome creation, find path and eval
for i in $(cat fofn.test)
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
    echo "Annotate: run_sim_add_gfa_weights_${annotate}.sh pop.gfa $i $k1 $k2 $k3"
    eval run_sim_add_gfa_weights_${annotate}.sh pop.gfa $i $k1 $k2 $k3

    # Find a path
    # Creates:
    #     $i.path
    qubo_solvers="mqlib gurobi dwave"
    if [[ " $qubo_solvers " =~  $solver  ]]; then
        if [ -z "${time_limits}" ]; then
            echo "Default time limits: 3,5"
            time_limits="3,5"
        fi
        if [ -z "${num_jobs}" ]; then
            echo "Default jobs: 2"
            num_jobs=2
        fi
        run_sim_solver_qubo.sh "$i".gfa "$solver" "$i" "$time_limits" "$num_jobs"
        echo "Finished sim solver qubo"
        for t in ${time_limits//,/ }; do
            for ((idx=0;idx<num_jobs;idx++)); do
                path2seq.pl "$i".gfa "$i.gaf.$t.$idx" > "$i".path_seq."$t".$idx
            done
        done
        echo "Finished path2seq"
    else
        time_limits=0
        num_jobs=1
        run_sim_solver_"$solver".sh "$i".gfa > "$i".path
        pathfinder2seq.pl pop.gfa "$i".path > "$i".path_seq.0.0
        sed -n '/PATH/,$p' "$i".path \
        | awk '/^\[/ {printf("%s ",$3)} END {print "\n"}' 1>&2
    fi


    # Evaluate the paths
    # Creates:
    #     $i.path_seq
    #     $i.bam
    #     $i.path_cons
    #     $i.eval_seq
    #     $i.eval_cons
    # New args
    echo "Start evaluate"
    for t in ${time_limits//,/ }; do
        for ((idx=0;idx<num_jobs;idx++)); do
            echo
            echo
            echo "=== Evaluate $t $idx ==="
            run_sim_evaluate_path.sh "$i" "$i".shred.fa "$t" $idx
        done
    done
done

# Summary
for t in ${time_limits//,/ }; do
    for ((idx=0;idx<num_jobs;idx++)); do
        echo 
        echo
        echo "=== Summary $t $idx==="
        head -1 $(ls -1 *.eval_cons.$t.$idx | head -1)
        cat ./*.eval_cons."$t".$idx | awk '!/Per/'  
    done
done

) 2>sim.err | tee sim.out

exit 0
