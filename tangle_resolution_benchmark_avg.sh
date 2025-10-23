#!/bin/bash

input_file="$1"
temp_file="merged_data.tmp"

column_names=("-" "len-t" "len-q" "covered" "used" "ncontig" "nbreaks" "nindel" "ndiff" "identity")

# Remove the temp file if it exists
rm -f "$temp_file"

# Extract numeric data from all tables and merge them
awk '!/^=/ && !/^-/ && !/^#/ && !/^Per/ && NF {print}' "$input_file" >> "$temp_file"

# Calculate column means and standard deviations
awk -v col_names="${column_names[*]}" '
BEGIN {
    split(col_names, names, " ");  # Split column names into an array
}
{
    for (i = 2; i <= NF; i++) {
        sum[i]  += $i
        sumsq[i] += ($i * $i)
        count[i]++
    }
}
END {
    print "Column\tMean\tStdDev"
    for (i = 2; i in sum; i++) {
        mean = sum[i] / count[i]
        variance = (sumsq[i] / count[i]) - (mean * mean)
        stddev = (variance > 0 ? sqrt(variance) : 0)
        col_name = (i <= length(names)) ? names[i] : "Column" i
        printf "%s\t%.2f\t%.2f\n", col_name, mean, stddev
    }
}' "$temp_file"

# Clean up temporary file
rm -f "$temp_file"
