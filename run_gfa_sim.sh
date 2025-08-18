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
    echo "       --edge2node         Use edge2node version"
    echo "       --trim-edges        Use trim_edges.pl"
    echo "       --pathfinder        Use pathfinder to get subgraphs"
    echo ""
    echo "The old API is still supported with fixed argument order."
}

config=${CONFIG:-$QDIR/config_hifi_km.sh}
. $config

prefix="sim_"
seed=1
solver=pathfinder
num_training=10
edge2node=0
trimedges=0
pathfinder_copy_numbers=0

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
	'-c1'|'--const1')
	    const1=$2
	    shift 2
	    continue
	    ;;
	'-c2'|'--const2')
	    const2=$2
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
    '--edge2node')
        edge2node=1
        shift 1
        continue
        ;;
    '--trim-edges')
        trimedges=1
        shift 1
        continue
        ;;
    '--pathfinder')
        pathfinder_copy_numbers=1
        shift 1
        continue
        ;;
    '--pathfinder_graph')
        pathfinder_graph=1
        shift 1
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
echo "Edge2node:$edge2node"
echo "TrimEdges:$trimedges"
echo "PathfinderGraph:$pathfinder_graph"
echo "Pathfinder:$pathfinder_copy_numbers"
echo ""

out_dir=$(printf "$prefix%05d" "$seed")
echo "$out_dir"
rm -rf "$out_dir" 2>/dev/null
mkdir "$out_dir"

cd "$out_dir" || exit 1

(
# Create fake genomes and use the training set to create a pangenome
# Creates:
#    pop.gfa
#    fofn.test
#    fofn.train
run_sim_create_gfa.sh "$seed" "$num_training"

# Foreach test genome, not used in pangenome creation, find path and eval
for i in $(cat fofn.test)
do
    # Add weights to the GFA via minigraph
    # Creates:
    #     $i.gfa (primary output; annotated pop.gfa)
    #     $i.shred.fa
    #     $i.mg
    echo "Annotate: run_sim_add_gfa_weights_${annotate}.sh pop.gfa $i"
    eval run_sim_add_gfa_weights_${annotate}.sh pop.gfa "$i"

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

        gfa_file_name="$i".gfa
        if [ "$edge2node" -eq 1 ]; then
            echo "Using edge2node"
	        gfa_edge_to_node.pl < "$i".gfa > "$i".edge2node.gfa
            gfa_file_name="$i".edge2node.gfa
        fi

        if [ "$trimedges" -eq 1 ]; then
            echo ">>> Using trim_edges.pl"
            trim_edges.pl $gfa_file_name > $i.edited.gfa
            gfa_file_name=$i.edited.gfa
        fi

        run_sim_solver_qubo.sh -f "$gfa_file_name" -s "$solver" -q "$i" \
        -t "$time_limits" -j "$num_jobs" -a "$annotate" \
        -c1 "$const1" -c2 "$const2" \
        --edge2node "$edge2node" --pathfinder "$pathfinder_copy_numbers" --pathfinder_graph "$pathfinder_graph"
        echo "Finished sim solver qubo"
    else
        time_limits=0
        num_jobs=1
        gfa=$i.gfa
        if [ "$trimedges" -eq 1 ]; then
            echo ">>> Using trim_edges.pl"
            trim_edges.pl $i.gfa > $i.edited.gfa
            gfa=$i.edited.gfa
        fi
        echo "Start $solver"
        run_sim_solver_"$solver".sh $gfa > "$i".path
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
        echo "=== Summary $t $idx ==="

        echo "eval seq"
        head -1 $(ls -1 *.eval_seq.$t.$idx | head -1)
        cat ./*.eval_seq."$t".$idx | awk '!/Per/'  

	    echo "eval cons"
        head -1 $(ls -1 *.eval_cons.$t.$idx | head -1)
        cat ./*.eval_cons."$t".$idx | awk '!/Per/'  
    done
done

) 2>sim.err | tee sim.out

exit 0
