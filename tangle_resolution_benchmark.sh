#!/bin/bash

qpg_dir="$HOME/quantumwork/pangenome/modules/qpg"
max_seed="$1"
time_limits="$2"
num_jobs="$3"
num_training="$4"
solvers="$5"


# parallel --header  : "$qpg_dir/run_gfa_sim.sh" -s {seed} -c "$qpg_dir"/config_illumina_{annotate}.sh -a {annotate} --solver {solver} -p {solver}.{annotate}. -t $time_limits -j $num_jobs -n $num_training --pathfinder ::: seed $(seq 1 $max_seed) ::: annotate mg km ga ::: solver pathfinder gurobi mqlib

# parallel --header  : "$qpg_dir/run_gfa_sim.sh" -s {seed} -c "$qpg_dir"/config_illumina_{annotate}.sh -a {annotate} --solver {solver} -p {solver}.{annotate}. -t $time_limits -j $num_jobs -n $num_training --pathfinder ::: seed $(seq 1 $max_seed) ::: annotate mg km ga ::: solver pathfinder dwave

parallel --header  : "$qpg_dir/run_gfa_sim.sh" -s {seed} -c "$qpg_dir"/config_illumina_{annotate}.sh -a {annotate} --solver {solver} -p {solver}.{annotate}. -t $time_limits -j $num_jobs -n $num_training --pathfinder ::: seed $(seq 1 $max_seed) ::: annotate mg km ga ::: solver ${solvers}