# Differential Gene Expression Analysis #

## Case Study Project where data was supplied by a client, in this case the invistigation was to see how a combination of a poor diet and infection with a parasite exacerbates inflammation in the colon ##
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
Execute all slurm scripts 

#### 1-fastQCP.slurm, Quality checks, trimming (adapter removal), quality check
  -> OUTPUT: raw data quality check multiqc_report_raw.html, trimmed reads e.g."Hpb_1_trim_R1.fq.gz" , trimmed data quality check multiqc_report_trim.html

#### 2-star_index_genome.slurm , indexes genome using STAR from: 
      Mus_musculus_balbcj.BALB_cJ_v1.113.gtf
      Mus_musculus_balbcj.BALB_cJ_v1.dna_sm.toplevel.fa
->  OUTPUT: SA files, genome parameters etc.

#### 3-star.slurm , aligning trimmed data to genome (STAR), QC alignment
->  OUTPUT: STAR alignment files for each sample "Aligned.out.bam", multiqc_report_star.html

#### 4-featurecounts.slurm , creates feature count files for each sample

#### 5. extraction.sh
extraction of featurecount information into "sample".tsv files and a count_matrix.csv file
Targets.txt file created manually
for the tidied data (samples "4" removed) use Targets2.txt




# RNA-seq ANALYSIS 
#### Deseq.R:
Packages: devtools v2.4.5, SARTools v1.8.1, DESeq2 1.42.1

Processes feature count files, filters poor quality, unaligned samples. normalises counts
identifies differentially expressed genes using parametric dispersion estimation. generates diagnostic plots and summary reports
 ->outputs: DietvsWorm.html DietvsWorm.RData



#### DEGs.R
using the results from seseq, DEGs are indicated using log2foldchange values of >1 or <-1, with 0.05 cutoff for padjusted values to prevent false dscovery
these were extracted into CSV files for downstream analysis.
-> Outputs: upregulated.csv downregulated.csv and all_significant.csv

#### extract.sh
The strain of mouse used in this analysis was very specific, IDs had to be coverted for any biological interpretation.
they were in the format MGP_BALBcJ_GOOOOOOOO4
using the online biomart tool, these IDs from upregulated and downregulated files were inputted and a list
of ENSMUG IDs were Provided. Note some of the terms were NA and some were matched to the same ENSMUG
##### this script replaced "MGP" terms with the correseponding "ENSMUG" IDs using a file called matched_genes.txt 


Volcano.R

Bio_processes.R



Packages and depdendencies:

