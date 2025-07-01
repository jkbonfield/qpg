#!/bin/bash
export PATH="$HOME/quantumwork/pangenome/bin:$PATH"

max_seed="$1"
time_limits="$2"
num_jobs="$3"
num_training="$4"

for annotate in ga km mg; do

solver="pathfinder"
max_file_name="$solver.$annotate.max.txt"
avg_file_name="$solver.$annotate.avg.txt"
rm "$max_file_name" 2> /dev/null
rm "$avg_file_name" 2> /dev/null

summary_file_name="$solver.$annotate.summary.txt"
rm "$summary_file_name" 2> /dev/null

for seed in $(seq 1 $max_seed); do
    grep -A $(( 4 + 2 * num_training )) "Summary " --no-group-separator < "$(printf "$solver.$annotate.%05d" "$seed")/sim.out" 2>/dev/null >> "$summary_file_name"
done


tangle_resolution_benchmark_max.sh "$summary_file_name" "$num_jobs" "$num_training" "seq"  >> "$max_file_name"

{
    echo "Average stats"
    average_gfa_sim.sh "$max_file_name"
} >> "$avg_file_name"

echo "===============" >> "$avg_file_name"



for solver in mqlib gurobi; do
    max_file_name="$solver.$annotate.max.txt"
    avg_file_name="$solver.$annotate.avg.txt"
    rm "$max_file_name" 2> /dev/null
    rm "$avg_file_name" 2> /dev/null


    for t in ${time_limits//,/ }; do
        summary_file_name="$solver.$annotate.$t.summary.txt"
        rm "$summary_file_name" 2> /dev/null

        for seed in $(seq 1 $max_seed); do
            grep -A $(( 4 + 2 * num_training )) "Summary $t " --no-group-separator < "$(printf "$solver.$annotate.%05d" "$seed")/sim.out" 2>/dev/null >> "$summary_file_name"
        done


        tangle_resolution_benchmark_max.sh "$summary_file_name" "$num_jobs" "$num_training" "seq"  >> "$max_file_name"

        {
            echo "Average stats for best runs with time limit $t"
            average_gfa_sim.sh "$max_file_name"
        } >> "$avg_file_name"

        echo "===============" >> "$avg_file_name"
    done
done


done

exit 0