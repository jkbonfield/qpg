#!/bin/bash
gfa_filepath="$1"
solver="$2"
query="$3"
time_limits="$4"
num_jobs="$5"
# mode="$6"

outdir="."

. ${CONFIG:-$QDIR/config_hifi_km.sh} # $mode

QUBO_DIR=/software/qpg/qubo
PATH=$PATH:$QUBO_DIR
source $QUBO_DIR/qubo_venv/bin/activate

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

if [[ " dwave " =~  $solver  ]]; then
    penalties="10,5,1"
else
    penalties="100,50,1"
fi
echo $penalties

python3 $QUBO_DIR/build_edge2node_qubo_matrix.py -f $gfa_filepath -d $outdir -c $copy_numbers -p $penalties
python3 $QUBO_DIR/oriented_max_path.py -s $solver -f "$gfa_filepath" -d "$outdir" -j "$num_jobs" -t $time_limits -o "$query.gaf" --edge2node

exit 0
