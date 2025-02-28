#!/bin/bash
gfa_filepath="$1"
solver="$2"
query="$3"
time_limits="$4"
num_jobs="$5"

outdir="."

QUBO_DIR=/nfs/users/nfs_j/jc59/quantumwork/pangenome/qubo_solvers
source $QUBO_DIR/.venv/bin/activate

python3 $QUBO_DIR/qubo_solvers/oriented_tangle/build_oriented_qubo_matrix.py -f $gfa_filepath -d $outdir
python3 $QUBO_DIR/qubo_solvers/oriented_tangle/oriented_max_path.py -s $solver -f "$gfa_filepath" -d "$outdir" -j "$num_jobs" -t $time_limits -o "$query.gaf"


# TODO: read file output and write to path string
