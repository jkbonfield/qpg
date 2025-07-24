#!/bin/bash
export PATH="$HOME/quantumwork/pangenome/modules/qpg:$PATH"

max_seed="$1"
time_limits="$2"
num_jobs="$3"
num_training="$4"
solvers="$5"
data_type="$6"

for annotate in ga km mg; do

solver="pathfinder"
max_file_name="$solver.$annotate.$data_type.max.txt"
avg_file_name="$solver.$annotate.$data_type.avg.txt"
rm "$max_file_name" 2> /dev/null
rm "$avg_file_name" 2> /dev/null

summary_file_name="$solver.$annotate.$data_type.summary.txt"
rm "$summary_file_name" 2> /dev/null

for seed in $(seq 1 $max_seed); do
    grep -A $(( 4 + 2 * num_training )) "Summary " --no-group-separator < "$(printf "$solver.$annotate.%05d" "$seed")/sim.out" 2>/dev/null >> "$summary_file_name"
done


tangle_resolution_benchmark_max.sh "$summary_file_name" "$num_training" "1" "$data_type"  >> "$max_file_name"

{
    echo "Average stats"
    tangle_resolution_benchmark_avg.sh "$max_file_name"
} >> "$avg_file_name"

echo "===============" >> "$avg_file_name"



for solver in $solvers; do    
    avg_file_name="$solver.$annotate.$data_type.avg.txt"
    rm "$avg_file_name" 2> /dev/null


    for t in ${time_limits//,/ }; do
        summary_file_name="$solver.$annotate.$data_type.$t.summary.txt"
        max_file_name="$solver.$annotate.$data_type.$t.max.txt"
        rm "$summary_file_name" 2> /dev/null
        rm "$max_file_name" 2> /dev/null

        for seed in $(seq 1 $max_seed); do
            grep -A $(( 4 + 2 * num_training )) "Summary $t " --no-group-separator < "$(printf "$solver.$annotate.%05d" "$seed")/sim.out" 2>/dev/null >> "$summary_file_name"
        done

        tangle_resolution_benchmark_max.sh "$summary_file_name" "$num_training" "$num_jobs" "$data_type"  >> "$max_file_name"

        {
            echo "Average stats for best runs with time limit $t"
            tangle_resolution_benchmark_avg.sh "$max_file_name"
        } >> "$avg_file_name"

        echo "===============" >> "$avg_file_name"
    done
done


done

exit 0