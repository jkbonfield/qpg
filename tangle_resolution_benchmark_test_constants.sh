#!/bin/bash
max_seed="$1"
time_limits="$2"
num_jobs="$3"
num_training="$4"
solvers="$5"

export QDIR=${QDIR:-$(pwd)}

for const1 in 0.5; do
for const2 in 0.8 0.9; do
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
    find . -type f ! -name '*.gfa' -a  ! -name '*subgraph' -a ! -name 'sim.*' -a ! -name '*cons.avg.parsed.txt' -a ! -name '*.cons.violin.txt' -delete
    cd ".."
done
done

# mqlib km:

# power.const.int_0.9.0.9/mqlib.km.cons.avg.parsed.txt
# 15 90.1973 94.6763 2.9333 1.6000 0.5000 0.3000 98.8383

# power.const.int_0.8.0.9/mqlib.km.cons.avg.parsed.txt
# 15 86.6480 94.6137 2.8000 1.6000 0.6333 0.3667 98.8097



# mqlib ga:

# power.const.int_0.9.0.6/mqlib.ga.cons.avg.parsed.txt
# 15 87.6187 94.9703 1.7000 1.7333 0.7000 0.2333 98.8377

# power.const.int_0.6.0.7/mqlib.ga.cons.avg.parsed.txt
# 15 87.7097 94.4943 1.7000 1.9000 0.8000 0.4333 98.3590

#power.const.int_0.4.0.7/mqlib.ga.cons.avg.parsed.txt
# 15 87.1830 94.4887 1.4667 1.8667 0.7000 0.4000 98.4787


# mqlib mg

# power.const.int_0.9.0.9/mqlib.mg.cons.avg.parsed.txt
# 15 81.5530 95.4750 3.9000 1.8333 0.3000 0.1333 99.3090

# power.const.int_0.7.0.9/mqlib.mg.cons.avg.parsed.txt
# 15 81.0590 95.8533 3.9000 1.9333 0.2667 0.1000 99.3483

# power.const.int_0.6.0.9/mqlib.mg.cons.avg.parsed.txt
# 15 81.3320 96.0403 3.9000 1.9667 0.3000 0.1000 99.3540

# power.const.int_0.5.0.9/mqlib.mg.cons.avg.parsed.txt
# 15 80.6057 96.6263 3.9000 1.6667 0.3000 0.1333 99.2727

# power.const.int_0.3.0.9/mqlib.mg.cons.avg.parsed.txt
# 15 80.2380 96.9163 3.9333 1.7667 0.3000 0.1000 99.3300