#!/bin/bash

help() {
    echo Usage: run_sim_solver_qubo.sh [options] 1>&2
    echo Options:
    echo "    -f,--filename  FILE    Use FILE as graph"
    echo "    -s,--solver    STR     Specify the solver"
    echo "    -q,--query     STR     Name of sequence aligned against the graph"
    echo "    -t,--times     INT_LIST   Time limits provided to QUBO solvers"
    echo "    -j,--jobs      INT     Number of runs of QUBO solvers"
    echo "    -a,--annotator STR     The annotator used to build the gfa"
    echo "       --pathfinder_graph INT    Use Pathfinder to get subgraphs and copy numbers if eq 1"
    echo "       --pathfinder INT    Use Pathfinder to get subgraphs if eq 1"
    echo "       --edge2node INT     Use edge2node to get copy numbers if eq 1 [DEPRECATED]"
    echo "       --subgraph  D W     Split graphs where nodes are more than D edges from"
    echo "                               nodes with weight >= W."
    echo ""
}

edge2node=0
pathfinder_copy_numbers=0
pathfinder_graph=0
subgraph_D=0
subgraph_W=0
const1=0.6
const2=1.0

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
	'-a'|'--annotator')
	    annotator=$2
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
	'--pathfinder_graph')
	    pathfinder_graph=$2
	    shift 2
	    continue
	    ;;
	'--pathfinder')
	    pathfinder_copy_numbers=$2
	    shift 2
	    continue
	    ;;
	'--subgraph')
	    subgraph_D=$2
	    subgraph_W=$3
	    shift 3
	    continue;
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

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
QUBO_DIR=$SCRIPT_DIR/qubo/qubo_solvers/oriented_tangle
echo $QUBO_DIR

echo "Gfa:                    $gfa_filepath"
echo "Query:                  $query"
echo "Solver:                 $solver"
echo "Time limits:            $time_limits"
echo "Num jobs:               $num_jobs"
echo "Annotator:              $annotator"
echo "Edge2node:              $edge2node"
echo "Pathfinder graph only:  $pathfinder_graph"
echo "Pathfinder:             $pathfinder_copy_numbers"
echo "Subgraph:               $subgraph_D $subgraph_W"
echo ""

if [[ " dwave " =~  $solver  ]]; then
    penalties="40,10,1"
else
    penalties="200,50,1"
fi
echo $penalties

if [ "$edge2node" -eq 1 ]; then
    echo "Solve with edge2node is deprecated"
    exit 1

    # copy_numbers=$(perl -e '
    # use strict;
    # while (<>) {
    #     next unless /^S/;
    #     m/sp:f:([0-9.]*).*sm:f:([0-9.]*)/;
    #     print int($1/'$shred_depth' + .5), ",";
    #     print int($2/'$shred_depth' + .5), ",";
    # }
    # ' $gfa_filepath)

    
    # echo "copy_numbers" >> sim.err
    # echo $copy_numbers >> sim.err

    # python3 "$QUBO_DIR/build_edge2node_qubo_matrix.py" -f "$gfa_filepath" -d "$outdir" -c "$copy_numbers" -p "$penalties"
    # python3 "$QUBO_DIR/oriented_max_path.py" -s "$solver" -f "$gfa_filepath" -d "$outdir" -j "$num_jobs" -t "$time_limits" -o "$query.gaf" "--edge2node"


    # for t in ${time_limits//,/ }; do
    #     for ((idx=0;idx<num_jobs;idx++)); do
    #         echo ">contig_1" >> "$query.path_seq.$t.$idx"
    #         path2seq.pl "$gfa_filepath" "$query.gaf.$t.$idx" >> "$query.path_seq.$t.$idx"
    #     done
    # done

elif [ "$pathfinder_graph" -eq 1 ]; then
    echo "Solve with pathfinder graph"
    run_sim_solver_qubo_with_pathfinder_graph.sh "$gfa_filepath" "$query" "$outdir" "$penalties" "$solver" "$num_jobs" "$time_limits" "$annotator" "$const1" "$const2"


elif [ "$pathfinder_copy_numbers" -eq 1 ]; then
    echo "Solve with pathfinder copy numbers"
    run_sim_solver_qubo_with_pathfinder.sh "$gfa_filepath" "$query" "$outdir" "$penalties" "$solver" "$num_jobs" "$time_limits"

    
else
    echo "Default solve"
    if [ "$subgraph_D" -ne 0 ]; then
	partition_graph.pl $gfa_filepath $subgraph_D $subgraph_W
	gfa_list=`echo $gfa_filepath.sub_graph.*`
    else
	gfa_list=$gfa_filepath
    fi

    counter=0
    for gfa_filepath in $gfa_list
    do
        counter=$((counter+1))
        echo gfa_filepath=$gfa_filepath
        copy_numbers=$(perl -e '
    use strict;
    open(my $gfa, "<", shift(@ARGV)) || die;
    while (<$gfa>) {
        next unless /^S/;
        m/SC:f:([0-9.]*)/;
        #print int($1/'$shred_depth' + $ARGV[0]), ",";
        print $1/'$shred_depth'*$ARGV[1] + $ARGV[0], ",";
    }
    ' "$gfa_filepath" "$const1" "$const2")
    # print int($1/30 + $ARGV[0]), ",";
    
        # print int($1/'$shred_depth' + 0.8), ",";
        python3 "$QUBO_DIR/build_oriented_qubo_matrix.py" -f "$gfa_filepath" -d "$outdir" -c "$copy_numbers" -p "$penalties"
        python3 "$QUBO_DIR/oriented_max_path.py" -s "$solver" -f "$gfa_filepath" -d "$outdir" -j "$num_jobs" -t "$time_limits" -o "$query.gaf"


        for t in ${time_limits//,/ }; do
            for ((idx=0;idx<num_jobs;idx++)); do
                fragment_content=""
                in_fragment=false
                while IFS= read -r line || [[ -n "$line" ]]; do
                    if [[ "$line" == "Begin fragment" ]]; then
                        if [ "$in_fragment" = true ] && [ -n "$fragment_content" ]; then
                            tmp_file=$(mktemp)
                            echo -e "$fragment_content" > "$tmp_file"
                            echo ">contig_$counter" >> "$query.path_seq.$t.$idx"
                            path2seq.pl "$gfa_filepath" "$tmp_file" >> "$query.path_seq.$t.$idx"
                            counter=$((counter+1))
                            rm "$tmp_file"
                        fi

                        in_fragment=true
                        fragment_content=""
                        continue 
                    fi

                    if [ "$in_fragment" = true ]; then
                        fragment_content+="$line"$'\n'
                    fi
                done  < "$query.gaf.$t.$idx"
                
                if [ "$in_fragment" = true ] && [ -n "$fragment_content" ]; then
                    tmp_file=$(mktemp)
                    echo -e "$fragment_content" > "$tmp_file"
                    echo ">contig_$counter" >> "$query.path_seq.$t.$idx"
                    path2seq.pl "$gfa_filepath" "$tmp_file" >> "$query.path_seq.$t.$idx"
                    rm "$tmp_file"
                fi
            done
        done
    done; #gfa_list
fi

# Tidy up intermediary files
rm *.pkl mqlib_input*.txt

exit 0
