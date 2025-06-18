#!/bin/bash
gfa_filepath="$1"
solver="$2"
query="$3"
time_limits="$4"
num_jobs="$5"
mode="$6"

outdir="."

. ${CONFIG:-$QDIR/sim_path_config_hifi.sh} $mode

QUBO_DIR=/nfs/users/nfs_j/jc59/quantumwork/pangenome/qubo_solvers
source $QUBO_DIR/.venv/bin/activate

copy_numbers=$(perl -e '
use strict;
while (<>) {
    next unless /^S/;
    m/dc:f:([0-9.]*)/;
    print int($1/'$shred_depth' + .8), ",";
}
' $gfa_filepath)

echo $copy_numbers >> sim.err
if [[ " dwave " =~ " $solver " ]]; then
    penalties="10,5,1"
else
    penalties="100,50,1"
fi
echo $penalties

python3 $QUBO_DIR/qubo_solvers/oriented_tangle/build_oriented_qubo_matrix.py -f $gfa_filepath -d $outdir -c $copy_numbers -p $penalties
python3 $QUBO_DIR/qubo_solvers/oriented_tangle/oriented_max_path.py -s $solver -f "$gfa_filepath" -d "$outdir" -j "$num_jobs" -t $time_limits -o "$query.gaf"

exit 0