#!/bin/bash
max_seed="$1"
time_limits="$2"
num_jobs="$3"
num_training="$4"
solvers="$5"

export QDIR=${QDIR:-$(pwd)}

for const1 in 0.6 0.7 0.8 0.9; do
for const2 in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9; do
    dirname="power.const.int_$const1.$const2"
    rm -rf "$dirname" 2>/dev/null
    mkdir $dirname
    cd $dirname || exit 1
    parallel --header  : "$QDIR/run_gfa_sim.sh" -s {seed} -c "$QDIR"/config_illumina_{annotate}.sh \
    -c1 $const1 -c2 $const2 \
    -a {annotate} --solver {solver} -p {solver}.{annotate}. -t "$time_limits" -j "$num_jobs" -n "$num_training" --pathfinder_graph \
    ::: seed $(seq 1 "$max_seed") ::: annotate mg km ga ::: solver ${solvers}

    "$QDIR/tangle_resolution_benchmark_stats.sh" "$max_seed" "$time_limits" "$num_jobs" "$num_training" "$solvers" "cons"
    "$QDIR/tangle_resolution_benchmark_parse_stats.sh" "cons"
    "$QDIR/tangle_resolution_benchmark_violin.sh" "$time_limits" "$solvers" "cons"
    rm ./*summary.txt
    rm ./*avg.txt
    rm ./*max.txt
    cd ".."
done
done

