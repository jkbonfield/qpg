#!/bin/bash
data_type="$1"
for annotate in ga km mg; do
for solver in pathfinder mqlib gurobi dwave; do
# for solver in pathfinder dwave; do


INPUT_FILE="./$solver.$annotate.$data_type.avg.txt"


if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found or is not readable." >&2
else
    awk '
    BEGIN {
        patterns[0] = "Average stats"
        patterns[1] = "covered"
        patterns[2] = "used"
        patterns[3] = "ncontig"
        patterns[4] = "nbreaks"
        patterns[5] = "nindel"
        patterns[6] = "ndiff"
        patterns[7] = "identity"
        num_patterns = 8 

        current_pattern_idx = 0

        delete extracted_values
    }

    {
        if (index($0, patterns[current_pattern_idx]) == 1) {
            # If a match is found:

            val = $NF

            gsub(/[^0-9.]/, "", val)

            extracted_values[current_pattern_idx] = val

            current_pattern_idx++

            if (current_pattern_idx == num_patterns) {

                for (i = 0; i < num_patterns; i++) {
                    printf "%s%s", extracted_values[i], (i == num_patterns - 1 ? "" : " ")
                }
                printf "\n" 

                current_pattern_idx = 0
                delete extracted_values
            }
        }
    }
    ' "$INPUT_FILE" > $solver.$annotate.$data_type.avg.parsed.txt
fi


done
done