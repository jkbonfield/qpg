#!/bin/bash

input_file="$1" 
table_size_new="$2"  # Number of rows per table
batch_size="$3"  # Number of rows per table
data_version="$4"
temp_file="merged_data.tmp"
temp_file2="merged_data2.tmp"

# Remove the temp file if it exists
rm -f "$temp_file"
rm -f "$temp_file2"

# Extract numeric data from all tables and merge them

grep -A $(( 1 + table_size_new )) "eval $data_version" < "$input_file" 2>/dev/null >> "$temp_file"
awk '!/^=/ && !/^-/ && !/^#/ && !/^Per/ && !/^eval/ && NF {print}' "$temp_file" >> "$temp_file2"


awk -v table_size_new=$table_size_new -v batch_size=$batch_size '
{
    if (NF < 4 || $3 == 0) next;  # Skip lines with insufficient columns

    row_idx = (NR - 1) % table_size_new # Row index within each table (0: num_training - 1)
    table_idx = int(((NR-1) % (table_size_new * batch_size)) / table_size_new) # Table index within each batch
    batch_idx = int( (NR-1) / (table_size_new * batch_size))
    
    avg = ($4 + $5) / 2  # Average of 3rd and 4th columns

    # Check if this row is the best so far for this row index in the batch
    if (!(batch_idx, row_idx) in max_avg || avg > max_avg[batch_idx, row_idx]) {
        max_avg[batch_idx, row_idx] = avg
        best_row[batch_idx, row_idx] = $0
    }
}
END {
    for (batch = 0; batch * table_size_new * batch_size < NR; batch++) {
        print "==============="
        print batch
        for (row = 0; row < table_size_new; row++) {
            if ((batch, row) in best_row) {
                print best_row[batch, row]
            }
        }
        print "----------------"
    }
}' "$temp_file2"


# Clean up temporary file
rm -f "$temp_file"
rm -f "$temp_file2"