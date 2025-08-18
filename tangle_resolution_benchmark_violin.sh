#!/bin/bash
time_limits=$1
solvers=$2
data_type="$3"

for annotate in km mg ga; do

    for solver in $solvers; do
    [[ " pathfinder " =~ $solver ]] && continue;   
    rm "$solver.$annotate.$data_type.violin.txt" 2> /dev/null
    for t in ${time_limits//,/ }; do
        {
            echo "$t"
            awk '/seq/ && NF {print ($4 + $5)/2}' < $solver.$annotate.$data_type.$t.max.txt
            echo '-------'
        } >> "$solver.$annotate.$data_type.violin.txt"
    done
    done

    solver=pathfinder
    rm "$solver.$annotate.violin.txt" 2> /dev/null
    {
        echo 0
        awk '/seq/ && NF {print ($4 + $5)/2}' < $solver.$annotate.$data_type.max.txt
        echo '-------'
    } >> "$solver.$annotate.$data_type.violin.txt"


done