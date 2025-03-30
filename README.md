# Differential Gene Expression Analysis #

## Case Study Project : data was supplied by a client, in this case the investigation was to see how a combination of a poor diet and infection with a parasite exacerbates inflammation in the colon ##
1. downloaded all raw reads and input into a folder named "rawdata"
2. downloaded Genome from: [https://www.ensembl.org/Mus_musculus/Info/Index](https://www.ensembl.org/Mus_musculus_BALB_cJ/Info/Index)
   Mus_musculus_balbcj.BALB_cJ_v1.113.gtf
    Mus_musculus_balbcj.BALB_cJ_v1.dna_sm.toplevel.fa

input files:
Genome/
  Mus_musculus_balbcj.BALB_cJ_v1.113.gtf
    Mus_musculus_balbcj.BALB_cJ_v1.dna_sm.toplevel.fa

rawdata/
    Hpb_1_1.fq.gz  Hpb_3_2.fq.gz  Naive_1_1.fq.gz  Naive_3_2.fq.gz
    Hpb_1_2.fq.gz  Hpb_4_1.fq.gz  Naive_1_2.fq.gz  Naive_4_1.fq.gz
    Hpb_2_1.fq.gz  Hpb_4_2.fq.gz  Naive_2_1.fq.gz  Naive_4_2.fq.gz
    Hpb_2_2.fq.gz  Hpb_5_1.fq.gz  Naive_2_2.fq.gz  Naive_5_1.fq.gz
    Hpb_3_1.fq.gz  Hpb_5_2.fq.gz  Naive_3_1.fq.gz  Naive_5_2.fq.gz

# PRE-PROCESSING   #
In UNIX system:
Execute all slurm scripts in order 1-4

#### 1-fastQCP.slurm, 
Quality checks, trimming (adapter removal), trimmed quality check
  -> OUTPUT: raw data quality check multiqc_report_raw.html, trimmed reads e.g."Hpb_1_trim_R1.fq.gz" , trimmed data quality check multiqc_report_trim.html

#### 2-star_index_genome.slurm 
Indexes genome using STAR from: 
      Mus_musculus_balbcj.BALB_cJ_v1.113.gtf
      Mus_musculus_balbcj.BALB_cJ_v1.dna_sm.toplevel.fa
->  OUTPUT: SA files, Genome genomeparameters.txt etc.

#### 3-star.slurm 
Aligning trimmed data to genome (STAR), QC alignment
->  OUTPUT: STAR alignment files for each sample "Aligned.out.bam", multiqc_report_star.html

#### 4-featurecounts.slurm 
Creates feature count files, quantifying gene expression at the gene level for each sample
-> OUTPUT: Naive_1.featureCount 


#### 5. extraction.sh
Processes of featurecount information into "sample".tsv files and a count_matrix.csv file
metadata files:
Targets.txt file created manually
for the tidied data (samples "4" removed) use Targets2.txt




# RNA-seq ANALYSIS 
#### Deseq.R:
Packages: devtools v2.4.5, SARTools v1.8.1, DESeq2 1.42.1

Processes feature count files, filters poor quality samples and unaligned reads.It normalises counts and uses parametric dispersion estimation for statistical modelling. generates diagnostic plots and summary reports of differentially expressed genes that can be seen in the results file. 
"Control" set as the reference condition and differentially expressed genes are identified using an adjusted p-value threshold of 0.05 (Benjamini-Hochberg method).
Using .html files for easy comprehension. DietvsWorm2.html is the cleaned dataset, using the file Targets2.txt.
 -> Outputs: DietvsWorm.html DietvsWorm.RData OR DietvsWorm2.html
*Note input files are located in Pre-Processing directory


#### DEGs.R
using the results from Deseq, DEGs are indicated using log2foldchange values of >1 or <-1, with 0.05 cutoff for p.adjusted values to prevent false discovery
these were extracted into CSV files for downstream analysis.
-> Outputs: upregulated.csv downregulated.csv and all_significant.csv

#### extract.sh
The strain of mouse used in this analysis was very specific, IDs had to be coverted for any biological interpretation.
they were in the format MGP_BALBcJ_GOOOOOOOO4
using the online biomart tool, these IDs from upregulated and downregulated files were inputted and a list
of ENSMUG IDs were Provided. Note some of the terms were NA and some were matched to the same ENSMUG
##### this script replaced "MGP" terms with the correseponding "ENSMUG" IDs using a file called matched_genes.txt 
-> Output:upregulated_geneID.csv + downregulated_geneID.csv


#### Volcano.R
packages required: ggplot2 v3.5.1 and ggrepel v0.9.6
Using the Infected_vs_Control.csv file outputted by DESEQ script, it plots upregulated and downregulated genes using the same log fold change and P.adjusted cut-offs. it labels the top 10 upregulated and downregulated genes using the files outputted from the extract.sh script
-> Output: Volcano.png (Results folder)

This script has the option to look at specific gene changes in the form of a barplot. As the client was interested in ALOX5 and ALOX15 these were visualised using their log2foldchange and the significance (P value). 
-> Output: ALOX5 and ALOX15.png (Results folder)


#### Bio_processes.R
Packages:BiocManager_1.30.25, clusterProfiler_4.12.6 , org.Mm.eg.db_3.19.1, AnnotationDbi_1.66.0  , enrichplot_1.24.4 , ggplot2_3.5.1, DOSE_3.30.5, ggrepel_0.9.6 , biomaRt_2.60.1 , GOplot_1.0.2 , pathview_1.44.0 
performs:
###### Gene Ontology Enrichment (Biological Processes)
-> Output: Dotplot, Barplot, .csv, Cnet plot

###### KEGG Enrichment (Biological Processes)
-> Output: Dotplot, Pathview analysis of specific pathways of interest for this analysis I looked at coagulation and wound healing by request of client. Script can be edited to look at different pathways


###### ENS_rename.sh
Using the biomart convert tool I was able to download the strain version gene IDs to the corresponding IDs on the GRCm39 genome.
this file is called :
mart_export_ENS.txt

AKA using a script that converted them into a new .csv. with each “MGP_BALBcJ” being linked to an “ENSMUG” identifier

It reads in the Infected_vs_Control.csv from the Deseq and DEG script. and replaces the IDs to allow for easier GSEA.
->Output: Infected_vs_Control_ENS.csv

###### GSEA.R
#### KEGG and GO gene set enrichment analyses of entire dataset
packages: "clusterProfiler", "org.Mm.eg.db", "enrichplot", "DOSE", "biomaRt"
Using the Infected_vs_Control_ENS.csv, this script ranks the gene list, removes duplicate IDs and adds in adding in jitter to ensure each term has a distinct ranking value (noise to break ties between terms), which preserves the statistical assumptions of GSEA and enables accurate calculation of enrichment scores and p-values. 
it allows for multiple visualisations for both GO and KEGG: inclusing enrichment map plots, cnet plots etc.



#### GSEAapp.R
Shiny app that allows client to set her own parameters for Log2Fold Change. For each GO biological Process term that was deemed significant, you can see what genes in that set are being up or downregulated. results are plotted using a graph.


For full explanation of how to use this app, please look at the app explanation document 









