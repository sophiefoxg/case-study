#!/bin/bash
#SBATCH --job-name=RNAseq_pipeline
#SBATCH --output=pipeline.out
#SBATCH --error=pipeline.err
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=Fox-GmuerSE@cardiff.ac.uk
#SBATCH --mail-type=END,FAIL

# Set working directory
workdir=$(pwd)

# Step 1: Quality Control and trimming (FastQCP)
echo "Running FastQC on raw reads..."
sbatch 1-fastQCP.slurm
wait

# Step 2: Genome Indexing (STAR)
echo "Generating STAR genome index..."
sbatch 2-star_index_genome.slurm
wait

# Step 3: Read Alignment (STAR)
echo "Aligning reads to reference genome..."
sbatch 3-star.slurm
wait

# Step 4: Gene Quantification (FeatureCounts)
echo "Counting reads per gene..."
sbatch 4-featurecounts.slurm
wait

#Step 5: Reformatting
echo "FeatureCounts2 to .CSV"
sbatch extraction.sh
wait

echo "Converting to GO terms"
sbatch extract.sh
wait


echo "Pre-Processing complete"
