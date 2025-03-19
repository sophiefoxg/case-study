README


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

PRE-PROCESSING   
In command line:
Execute all slurm scripts 
1-fastQCP.slurm, Quality checks, trimming (adapter removal), quality check
  -> OUTPUT: raw data quality check multiqc_report_raw.html, trimmed reads e.g."Hpb_1_trim_R1.fq.gz" , trimmed data quality check multiqc_report_trim.html
2-star_index_genome.slurm , indexes genome using STAR from: 
      Mus_musculus_balbcj.BALB_cJ_v1.113.gtf
      Mus_musculus_balbcj.BALB_cJ_v1.dna_sm.toplevel.fa
->  OUTPUT: SA files, genome parameters etc
3-star.slurm , aligninging trimmed data to genome (STAR), QC alignment
->  OUTPUT: STAR alignment files for each sample, multiqc_report_star.html
4-featurecounts.slurm , creates feature count files for each sample

5. extraction.sh
extraction of featurecount information into "sample".tsv files and a count_matrix.csv file
Targets.txt file created manually

ALL THAT'S REQUIRED FOR ANALYSIS IS:

Before data tidying:
Targets.txt
featureCount (folder containing .tsv files)


After data tidying
Targets2.txt
featureCount2


RNA-seq ANALYSIS
R
Deseq.R:

DEGs.R

Bio_processes.R

Volcano.R


