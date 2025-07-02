#!/bin/bash

for annotate in ga km mg; do
for solver in pathfinder mqlib gurobi; do

INPUT_FILE="./$solver.$annotate.avg.txt"


if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found or is not readable." >&2
    exit 1
fi

awk '
BEGIN {
    patterns[0] = "Average stats"
    patterns[1] = "covered"
    patterns[2] = "used"
    patterns[3] = "nbreaks"
    patterns[4] = "nindel"
    patterns[5] = "ndiff"
    patterns[6] = "identity"
    num_patterns = 7 

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
' "$INPUT_FILE" > $solver.$annotate.avg.parsed.txt


done
done