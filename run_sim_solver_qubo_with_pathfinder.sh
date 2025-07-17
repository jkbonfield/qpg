#!/bin/bash

gfa_filepath=$1
query=$2 
outdir=$3 
penalties=$4
solver=$5 
num_jobs=$6 
time_limits=$7

. ${CONFIG:-$QDIR/config_illumina.sh}

QUBO_DIR=/software/qpg/qubo
PATH=$PATH:$QUBO_DIR
source $QUBO_DIR/qubo_venv/bin/activate

echo "Solve with pathfinder copy numbers"
echo $pathfinder $pathfinder_opts $gfa_filepath
eval $pathfinder $pathfinder_opts $gfa_filepath 2>$gfa_filepath.pf.err > "$query".path

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
    ' "$query".path > "$query".subgraph

    awk '
    BEGIN {
        in_subgraph_table = 0;           # Flag to indicate if we are currently inside a subgraph table
        current_subgraph_name = "";      # Stores the name of the current subgraph being processed
        node_list = "";                  # Accumulates space-separated node names for the current subgraph
        data_list = "";                  # Accumulates space-separated second column data for the current subgraph
    }

    /^Subgraph/ {
        if (in_subgraph_table == 1 && current_subgraph_name != "") {
            print "SUBGRAPH_START:" current_subgraph_name;
            print "NODES_LIST:" node_list;
            print "DATA_LIST:" data_list;
            print "SUBGRAPH_END";
        }

        in_subgraph_table = 1;           
        current_subgraph_name = $2;      # Assume the subgraph name is the second field on the line
        node_list = "";    
        data_list = "";              
        next;                            
    }

    in_subgraph_table == 1 {
        if (NF >= 1) {
            if (node_list == "") {
                node_list = $1;
            } else {
                node_list = node_list " " $1; 
            }
            if (NF >= 2) {
                if (data_list == "") {
                    data_list = $2;
                } else {
                    data_list = data_list " " $2; 
                }
            } else {
                if (data_list == "") {
                    data_list = ""; 
                } else {
                    data_list = data_list " ";
                }
            }
        } else if (NF == 0) { 
            if (current_subgraph_name != "") {
                print "SUBGRAPH_START:" current_subgraph_name;
                print "NODES_LIST:" node_list;
                print "DATA_LIST:" data_list;
                print "SUBGRAPH_END";
            }
            in_subgraph_table = 0;       
            current_subgraph_name = "";  
            node_list = "";    
            data_list = "";          
        }
    }

    END {
        if (in_subgraph_table == 1 && current_subgraph_name != "") {
            print "SUBGRAPH_START:" current_subgraph_name;
            print "NODES_LIST:" node_list;
            print "DATA_LIST:" data_list;
            print "SUBGRAPH_END";
        }
    }
    ' "$query".subgraph | while IFS= read -r line; do

        if [[ "$line" == "SUBGRAPH_START:"* ]]; then
            current_subgraph_name="${line#SUBGRAPH_START:}" # Extract the subgraph name
            current_nodes_str=""
            current_data_str=""
        elif [[ "$line" == "NODES_LIST:"* ]]; then
            current_nodes_str="${line#NODES_LIST:}"
        elif [[ "$line" == "DATA_LIST:"* ]]; then
            current_data_str="${line#DATA_LIST:}"
        elif [[ "$line" == "SUBGRAPH_END" ]]; then
            output_gfa_file="${query}.subgraph.${current_subgraph_name}.gfa"
            output_data_file="${query}.copy_numbers.${current_subgraph_name}.txt"

            echo "Creating $output_gfa_file and $output_data_file..."
            > "$output_gfa_file" 
            > "$output_data_file" 

            if [ -n "$current_data_str" ]; then
                echo "$current_data_str" | tr ' ' ',' >> "$output_data_file"
            fi

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

                python3 "$QUBO_DIR/build_oriented_qubo_matrix.py" -f "$output_gfa_file" -d "$outdir" -c "$(cat $output_data_file)" -p "$penalties"
                python3 "$QUBO_DIR/oriented_max_path.py" -s "$solver" -f "$output_gfa_file" -d "$outdir" -j "$num_jobs" -t "$time_limits" -o "$query.subgraph.$current_subgraph_name.gaf"

                for t in ${time_limits//,/ }; do
                    for ((idx=0;idx<num_jobs;idx++)); do
                        echo ">contig_$current_subgraph_name" >> "$query.path_seq.$t.$idx"
                        path2seq.pl "$output_gfa_file" "$query.subgraph.$current_subgraph_name.gaf.$t.$idx" >> "$query.path_seq.$t.$idx"
                    done
                done

                current_subgraph_name=""
                current_nodes_str=""
                current_data_str=""
            fi
        fi
    done
