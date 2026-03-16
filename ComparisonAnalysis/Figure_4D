######################################
#Figure 4 D
######################################
#!/bin/bash
# process_all_matrices.sh

# Directory containing matrix files
MATRIX_DIR="."
OUTPUT_DIR="middle_values"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Find all gzipped matrix files
matrix_files=$(ls matrix_*.gz)

echo "Found $(echo "$matrix_files" | wc -w) matrix files to process"

# Process each file
for file in $matrix_files; do
    # Extract base name without .gz extension
    base_name="${file%.gz}"
    output_file="${OUTPUT_DIR}/${base_name}_middle_8_values.txt"
    
    echo "Processing: $file -> $output_file"
    
    # Extract middle 8 values (columns 37-44 as in your example)
    zcat "$file" | 
    awk '!/^@/ {
        # Print metadata: chrom, start, end, ID
        printf("%s\t%s\t%s\t%s", $1, $2, $3, $4);
        
        # Print middle 10 values (columns 37-46)
        for(i=37;i<=44;i++) printf("\t%f", $i);
        printf("\n")
    }' > "$output_file"
    
    # Count lines processed
    line_count=$(wc -l < "$output_file")
    echo "  Extracted $line_count regions"
done

echo "Processing complete! Files saved in: $OUTPUT_DIR"
