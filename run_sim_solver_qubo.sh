#!/bin/bash
gfa_filepath="$1"
solver="$2"
query="$3"
time_limits="$4"
num_jobs="$5"
edge2node="$6"

outdir="."

. ${CONFIG:-$QDIR/config_hifi_km.sh} # $mode

QUBO_DIR=/software/qpg/qubo
PATH=$PATH:$QUBO_DIR
source $QUBO_DIR/qubo_venv/bin/activate

echo $gfa_filepath




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
else

    copy_numbers=$(perl -e '
    use strict;
    while (<>) {
        next unless /^S/;
        m/SC:f:([0-9.]*)/;
        print int($1/50 + 0.8), ",";
    }
    ' $gfa_filepath)
    # print int($1/'$shred_depth' + 0.8), ",";
    python3 "$QUBO_DIR/build_oriented_qubo_matrix.py" -f "$gfa_filepath" -d "$outdir" -c "$copy_numbers" -p "$penalties"
    python3 "$QUBO_DIR/oriented_max_path.py" -s "$solver" -f "$gfa_filepath" -d "$outdir" -j "$num_jobs" -t "$time_limits" -o "$query.gaf"
fi


exit 0
