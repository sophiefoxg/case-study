#!/bin/bash
#SBATCH --job-name=star_genome_generate
#SBATCH --output=star_genome_generate_%j.out
#SBATCH --error=star_genome_generate_%j.err
#SBATCH --partition=defq   # Using 'defq' based on your HPC partition
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8      # Number of CPU cores
#SBATCH --mem=64G              # Memory allocation
#SBATCH --time=3-00:00:00      # Time limit (3 days)
#SBATCH --mail-user=Fox-GmuerSE@cardiff.ac.uk
#SBATCH --mail-type=END,FAIL

# Load STAR module
module load star/2.7.11b-conda

workdir=$(pwd)

# Generate Genome Index
STAR    --runThreadN 8 \
        --limitGenomeGenerateRAM 9594470667 \
        --runMode genomeGenerate \
        --genomeDir ~/mydata/case-study/Pre-processing/Genome \
        --genomeFastaFiles ~/mydata/case-study/Pre-processing/Genome/Mus_musculus_balbcj.BALB_cJ_v1.fa \
        --sjdbGTFfile ~/mydata/case-study/Pre-processing/Genome/Mus_musculus_balbcj.BALB_cJ_v1.gtf \
        --sjdbOverhang 74

# Unload STAR module
module unload star/2.7.11b-conda

