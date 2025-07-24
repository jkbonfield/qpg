#!/bin/bash

help() {
    echo Usage: run_sim_solver_qubo.sh [options] 1>&2
    echo Options:
    echo "    -f,--filename  FILE    Use FILE as graph"
    echo "    -s,--solver    STR     Specify the solver"
    echo "    -q,--query     STR     Name of sequence aligned against the graph"
    echo "    -t,--times     INT_LIST   Time limits provided to QUBO solvers"
    echo "    -j,--jobs      INT     Number of runs of QUBO solvers"
    echo "       --pathfinder INT    Use Pathfinder to get subgraphs and copy numbers if eq 1"
    echo "       --edge2node INT    Use edge2node to get copy numbers if eq 1"
    echo ""
}

edge2node=0
pathfinder_copy_numbers=0

while true
do
    case "$1" in
	'-h')
	    help
	    exit 0
	    ;;
	'-f'|'--filename')
	    gfa_filepath=$2
	    shift 2
	    continue
	    ;;
	'-s'|'--solver')
	    solver=$2
	    shift 2
	    continue
	    ;;
	'-q'|'--query')
	    query=$2
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
	'--pathfinder')
	    pathfinder_copy_numbers=$2
	    shift 2
	    continue
	    ;;
    '--edge2node')
        edge2node=$2
        shift 2
        continue
        ;;
	*)
	    break
	    ;;
    esac
done


outdir="."

. ${CONFIG:-$QDIR/config_illumina.sh}

QUBO_DIR=/software/qpg/qubo
PATH=$PATH:$QUBO_DIR
source $QUBO_DIR/qubo_venv/bin/activate

echo "Gfa:         $gfa_filepath"
echo "Query:       $query"
echo "Solver:      $solver"
echo "Time limits: $time_limits"
echo "Num jobs:    $num_jobs"
echo "Edge2node:   $edge2node"
echo "Pathfinder:  $pathfinder_copy_numbers"
echo ""

if [[ " dwave " =~  $solver  ]]; then
    penalties="10,5,1"
else
    penalties="100,50,1"
fi
echo $penalties

if [ "$edge2node" -eq 1 ]; then
    echo "Solve with edge2node"

    copy_numbers=$(perl -e '
    use strict;
    while (<>) {
        next unless /^S/;
        m/sp:f:([0-9.]*).*sm:f:([0-9.]*)/;
        print int($1/'$shred_depth' + .5), ",";
        print int($2/'$shred_depth' + .5), ",";
    }
    ' $gfa_filepath)

    
    echo "copy_numbers" >> sim.err
    echo $copy_numbers >> sim.err

    python3 "$QUBO_DIR/build_edge2node_qubo_matrix.py" -f "$gfa_filepath" -d "$outdir" -c "$copy_numbers" -p "$penalties"
    python3 "$QUBO_DIR/oriented_max_path.py" -s "$solver" -f "$gfa_filepath" -d "$outdir" -j "$num_jobs" -t "$time_limits" -o "$query.gaf" "--edge2node"


    for t in ${time_limits//,/ }; do
        for ((idx=0;idx<num_jobs;idx++)); do
            echo ">contig_1" >> "$query.path_seq.$t.$idx"
            path2seq.pl "$gfa_filepath" "$query.gaf.$t.$idx" >> "$query.path_seq.$t.$idx"
        done
    done


elif [ "$pathfinder_copy_numbers" -eq 1 ]; then
    run_sim_solver_qubo_with_pathfinder.sh $gfa_filepath $query $outdir $penalties $solver $num_jobs $time_limits

    
else
    echo "Default solve"

    copy_numbers=$(perl -e '
    use strict;
    while (<>) {
        next unless /^S/;
        m/SC:f:([0-9.]*)/;
        print int($1/30 + 0.8), ",";
    }
    ' $gfa_filepath)
    
    # print int($1/'$shred_depth' + 0.8), ",";
    python3 "$QUBO_DIR/build_oriented_qubo_matrix.py" -f "$gfa_filepath" -d "$outdir" -c "$copy_numbers" -p "$penalties"
    python3 "$QUBO_DIR/oriented_max_path.py" -s "$solver" -f "$gfa_filepath" -d "$outdir" -j "$num_jobs" -t "$time_limits" -o "$query.gaf"

    for t in ${time_limits//,/ }; do
        for ((idx=0;idx<num_jobs;idx++)); do
            echo ">contig_1" >> "$query.path_seq.$t.$idx"
            path2seq.pl "$gfa_filepath" "$query.gaf.$t.$idx" >> "$query.path_seq.$t.$idx"
        done
    done
fi


exit 0
