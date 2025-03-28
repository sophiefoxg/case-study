#!/bin/bash

# Input files
input_csv="Infected_vs_Control.csv"
mapping="mart_export_ENS.txt"
output_csv="Infected_vs_Control_ENS.csv"

# Step 1: Create gene → ENSMUSG mapping (first match only)
tail -n +2 "$mapping" | awk -F'\t' '!seen[$1]++ && $2 != "" {print $1, $2}' > ens_mapped.txt

# Step 2: Copy header, changing first column name
head -n 1 "$input_csv" | sed 's/^"[^"]*"/"Mouse_gene_ID"/' > "$output_csv"

# Step 3: Process data rows
tail -n +2 "$input_csv" | while IFS=',' read -r gene rest; do
    # Remove quotes from the gene ID
    clean_gene=$(echo "$gene" | tr -d '"')

    # Find mapped Ensembl ID
    ens=$(grep -m1 "^$clean_gene " ens_mapped.txt | awk '{print $2}')

    # Add quotes back to mapped ID (or NA)
    echo "\"${ens:-NA}\",$rest" >> "$output_csv"
done

