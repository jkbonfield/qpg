#!/bin/bash
max_seed="$1"
time_limits="$2"
num_jobs="$3"
num_training="$4"
solvers="$5"
pathfinder="$6"

export QDIR=${QDIR:-$(pwd)}

if [ "$pathfinder" -eq 1 ]; then

    parallel --header  : "$QDIR/run_gfa_sim.sh" -s {seed} -c "$QDIR"/config_illumina_{annotate}.sh \
    -a {annotate} --solver {solver} -p {solver}.{annotate}. -t "$time_limits" -j "$num_jobs" -n "$num_training" --pathfinder_graph \
    ::: seed $(seq 1 "$max_seed") ::: annotate mg km ga ::: solver ${solvers}

else
    parallel --header  : "$QDIR/run_gfa_sim.sh" -s {seed} -c "$QDIR"/config_illumina_{annotate}.sh \
    -a {annotate} --solver {solver} -p {solver}.{annotate}. -t "$time_limits" -j "$num_jobs" -n "$num_training" \
    ::: seed $(seq 1 "$max_seed") ::: annotate mg km ga ::: solver ${solvers}
fi

"$QDIR/tangle_resolution_benchmark_stats.sh" "$max_seed" "$time_limits" "$num_jobs" "$num_training" "$solvers" "cons"
"$QDIR/tangle_resolution_benchmark_parse_stats.sh" "cons"
"$QDIR/tangle_resolution_benchmark_violin.sh" "$time_limits" "$solvers" "cons"