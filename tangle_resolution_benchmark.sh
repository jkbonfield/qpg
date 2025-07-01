#!/bin/bash

qpg_dir="$HOME/quantumwork/pangenome/modules/qpg"

parallel -j 128 --header  : "$qpg_dir/run_gfa_sim.sh" -s {seed} -c "$qpg_dir"/config_illumina_{annotate}.sh -a {annotate} --solver {solver} -p {solver}.{annotate}. -t 1,3 -j 2 -n 2 ::: seed $(seq 1 5) ::: annotate ga km mg ::: solver pathfinder gurobi mqlib