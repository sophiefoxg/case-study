#!/bin/bash

# Define variables
workdir=$(pwd)

# Process each .featurecount file to extract gene counts
for f in ${workdir}/*.featurecount
do
    base=$(basename "$f" | sed 's/.featurecount//')

    # Extract GeneID and count columns, remove headers
    output_file="${workdir}/${base}_counts.tsv"

    awk -F '\t' 'NR>2 {print $1 "\t" $7}' "$f" > "$output_file"
done

# Define output CSV file
csv_file="${workdir}/count_matrix.csv"

# Write CSV header manually
echo "GeneID,Hpb_1,Hpb_2,Hpb_3,Hpb_4,Hpb_5,Naive_1,Naive_2,Naive_3,Naive_4,Naive_5" > "$csv_file"

# Merge data and convert tabs to commas
paste <(awk -F '\t' '{print $1}' Hpb_1_counts.tsv) \
      <(awk -F '\t' '{print $2}' Hpb_1_counts.tsv) \
      <(awk -F '\t' '{print $2}' Hpb_2_counts.tsv) \
      <(awk -F '\t' '{print $2}' Hpb_3_counts.tsv) \
      <(awk -F '\t' '{print $2}' Hpb_4_counts.tsv) \
      <(awk -F '\t' '{print $2}' Hpb_5_counts.tsv) \
      <(awk -F '\t' '{print $2}' Naive_1_counts.tsv) \
      <(awk -F '\t' '{print $2}' Naive_2_counts.tsv) \
      <(awk -F '\t' '{print $2}' Naive_3_counts.tsv) \
      <(awk -F '\t' '{print $2}' Naive_4_counts.tsv) \
      <(awk -F '\t' '{print $2}' Naive_5_counts.tsv) | tr '\t' ',' >> "$csv_file"

