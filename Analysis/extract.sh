#!/bin/bash

workdir=$(pwd)


#upregulated matching


awk '/gene_id/ && /projection_parent_gene/ {match($0, /gene_id "([^"]+)"/, g); match($0, /projection_parent_gene "([^"]+)\./, p); print g[1], p[1]}' ~/mydata/case-study/Pre-processing/Genome/Mus_musculus_balbcj.BALB_cJ_v1.113.gtf > ~/mydata/case-study/Pre-processing/Genome/extracted_genes.gtf

# Loop through each Gene ID in the first column of all_DEGs_upregulated.csv
tail -n +2 upregulated.csv | cut -d ',' -f1 | while read gene_id; do
    # Use grep to find the matching row in extracted_genes.gtf
    grep -m1 "$gene_id" ~/mydata/case-study/Pre-processing/Genome/extracted_genes.gtf >> matched_genes.txt

done



#Making new csv file

# Create header for the new CSV file
echo "projection_parent_gene,Log2FoldChange,p_value" > upregulated_ensembl.csv

# Loop through each line in all_DEGs_upregulated.csv (excluding header)
tail -n +2 upregulated.csv | while IFS=',' read -r gene_id log2fc pval; do
    # Find corresponding ENSMUSG ID using grep
    ensmusg=$(grep -m1 "^$gene_id " matched_genes.txt | awk '{print $2}')
    
    # If a match is found, print to the new file; otherwise, assign "NA"
    echo "${ensmusg:-NA},$log2fc,$pval" >> upregulated_ensembl.csv
done



#downregulated matching

awk '/gene_id/ && /projection_parent_gene/ {match($0, /gene_id "([^"]+)"/, g); match($0, /projection_parent_gene "([^"]+)\./, p); print g[1], p[1]}' ~/mydata/case-study/Pre-processing/Genome/Mus_musculus_balbcj.BALB_cJ_v1.113.gtf > ~/mydata/case-study/Pre-processing/Genome/extracted_genes.gtf

# Loop through each Gene ID in the first column of all_DEGs_downregulated.csv
tail -n +2 downregulated.csv | cut -d ',' -f1 | while read gene_id; do
    # Use grep to find the matching row in extracted_genes.gtf
    grep -m1 "$gene_id" ~/mydata/case-study/Pre-processing/Genome/extracted_genes.gtf >> matched_genes.txt

done

#making new csv file

# Create header for the new CSV file
echo "projection_parent_gene,Log2FoldChange,p_value" > downregulated_ensembl.csv

# Loop through each line in all_DEGs_upregulated.csv (excluding header)
tail -n +2 downregulated.csv | while IFS=',' read -r gene_id log2fc pval; do
    # Find corresponding ENSMUSG ID using grep
    ensmusg=$(grep -m1 "^$gene_id " matched_genes.txt | awk '{print $2}')

    # If a match is found, print to the new file; otherwise, assign "NA"
    echo "${ensmusg:-NA},$log2fc,$pval" >> downregulated_ensembl.csv
done


