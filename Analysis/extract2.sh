#!/bin/bash

input_csv="all_significant.csv"
mapping_file="matched_genes.txt"
output_csv="all_significant_ensembl.csv"

awk -F',' '
    BEGIN {
        OFS = ",";
        # Build map from full gene name (prefix + "_" + gene ID) → ENSMUSG
        while ((getline < "'"$mapping_file"'") > 0) {
            key = $1 "_" $2;  # e.g., MGP_BALBcJ_G0000004
            map[key] = $3;
        }
        close("'"$mapping_file"'");

        # Print header
        print "Gene,log2FoldChange,padj";
    }

    NR > 1 {
        ensm = (map[$1] ? map[$1] : "NA");
        print ensm, $2, $3;
    }
' "$input_csv" > "$output_csv"

