#!/bin/bash

gfa_filepath=$1
query=$2 
outdir=$3 
penalties=$4
solver=$5 
num_jobs=$6 
time_limits=$7
annotator=$8
const1=$9
const2=${10}

. ${CONFIG:-$QDIR/config_illumina.sh}

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
QUBO_DIR=$SCRIPT_DIR/qubo/qubo_solvers/oriented_tangle
# source $QUBO_DIR/qubo_venv/bin/activate

echo "Solve with pathfinder copy numbers"
eval $pathfinder $pathfinder_opts "$gfa_filepath" 2>"$gfa_filepath.pf.err" > "$query.path"

awk '
    BEGIN { 
        in_subgraph_table = 0;
        subgraph_idx = 1; 
    }
    /^PATH/ {
        in_subgraph_table = 0;
    }
    /^SUBGRAPH/ {
        in_subgraph_table = 1;
        print "Subgraph", subgraph_idx
        subgraph_idx = subgraph_idx + 1
        next; # Skip the header line
    }

    in_subgraph_table == 1 {
        if (NF >= 5) {
            print $3, $6;
        }
    }
    ' "$query.path" > "$query.subgraph"

    awk '
    BEGIN {
        in_subgraph_table = 0;           
        current_subgraph_name = "";      
        node_list = "";                  
    }

    /^Subgraph/ {
        if (in_subgraph_table == 1 && current_subgraph_name != "") {
            print "SUBGRAPH_START:" current_subgraph_name;
            print "NODES_LIST:" node_list;
            print "SUBGRAPH_END";
        }

        in_subgraph_table = 1;           
        current_subgraph_name = $2;      # Assume the subgraph name is the second field on the line
        node_list = "";    
        next;                            
    }

    in_subgraph_table == 1 {
        if (NF >= 2 && $2 > 0) {
            if (node_list == "") {
                node_list = $1;
            } else {
                node_list = node_list " " $1; 
            }
        } else if (NF == 0) { 
            if (current_subgraph_name != "") {
                print "SUBGRAPH_START:" current_subgraph_name;
                print "NODES_LIST:" node_list;
                print "SUBGRAPH_END";
            }
            in_subgraph_table = 0;       
            current_subgraph_name = "";  
            node_list = "";    
        }
    }

    END {
        if (in_subgraph_table == 1 && current_subgraph_name != "") {
            print "SUBGRAPH_START:" current_subgraph_name;
            print "NODES_LIST:" node_list;
            print "SUBGRAPH_END";
        }
    }
    ' "$query".subgraph | while IFS= read -r line; do

        if [[ "$line" == "SUBGRAPH_START:"* ]]; then
            current_subgraph_name="${line#SUBGRAPH_START:}" 
            current_nodes_str=""
        elif [[ "$line" == "NODES_LIST:"* ]]; then
            current_nodes_str="${line#NODES_LIST:}"
        elif [[ "$line" == "SUBGRAPH_END" ]]; then
            output_gfa_file="${query}.subgraph.${current_subgraph_name}.gfa"

            touch "$output_gfa_file" 


            if [[ -n "$current_nodes_str" ]]; then 

                awk -v nodes_to_include="$current_nodes_str" '
                    BEGIN {
                        split(nodes_to_include, node_set_arr, " ");
                        for (i in node_set_arr) {
                            node_set[node_set_arr[i]] = 1;
                        }
                    }

                    $1 == "S" {
                        if (node_set[$2]) {
                            print $0;
                        }
                    }

                    $1 == "L" {
                        if (node_set[$2] && node_set[$4]) {
                            print $0;
                        }
                    }
                ' "$gfa_filepath" >> "$output_gfa_file"

                # if [[ " km " =~  $annotator  ]]; then
                #     const1=0.6; const2=0.8;
                # elif [[ " mg " =~  $annotator  ]]; then
                #     const1=0.8; const2=0.8;
                # else
                #     const1=0.6; const2=0.8;
                # fi

                copy_numbers=$(perl -e '
                use strict;
                open(my $gfa, "<", shift(@ARGV)) || die;
                while (<$gfa>) {
                    next unless /^S/;
                    m/SC:f:([0-9.]*)/;
                    #print ((($1/30+.1)**0.8 + $ARGV[0] + int($1/30+$ARGV[1]))/2);
                    #print ((3*(($1/30+.1)**0.8 + $ARGV[0]) + 7*(int($1/30+$ARGV[1])))/10);
                    #print (($1/30)**$ARGV[0]+$ARGV[1]);
                    print (($1/30)+$ARGV[0]); #P1
                    print ",";
                }
                ' "$output_gfa_file" "$const1" "$const2")

                python3 "$QUBO_DIR/build_oriented_qubo_matrix.py" -f "$output_gfa_file" -d "$outdir" -c "$copy_numbers" -p "$penalties"
                python3 "$QUBO_DIR/oriented_max_path.py" -s "$solver" -f "$output_gfa_file" -d "$outdir" -j "$num_jobs" -t "$time_limits" -o "$query.subgraph.$current_subgraph_name.gaf"

                for t in ${time_limits//,/ }; do
                    for ((idx=0;idx<num_jobs;idx++)); do
                        fragment_content=""
                        in_fragment=false
                        counter=0
                        while IFS= read -r line || [[ -n "$line" ]]; do
                            if [[ "$line" == "Begin fragment" ]]; then
                                if [ "$in_fragment" = true ] && [ -n "$fragment_content" ]; then
                                    tmp_file=$(mktemp)
                                    echo -e "$fragment_content" > "$tmp_file"
                                    echo ">contig_$current_subgraph_name.$counter" >> "$query.path_seq.$t.$idx"
                                    path2seq.pl "$output_gfa_file" "$tmp_file" >> "$query.path_seq.$t.$idx"
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
                        done  < "$query.subgraph.$current_subgraph_name.gaf.$t.$idx"
                        
                        if [ "$in_fragment" = true ] && [ -n "$fragment_content" ]; then
                            tmp_file=$(mktemp)
                            echo -e "$fragment_content" > "$tmp_file"
                            echo ">contig_$current_subgraph_name.$counter" >> "$query.path_seq.$t.$idx"
                            path2seq.pl "$output_gfa_file" "$tmp_file" >> "$query.path_seq.$t.$idx"
                            rm "$tmp_file"
                        fi
                    done
                done

                current_subgraph_name=""
                current_nodes_str=""
            fi
        fi
    done
