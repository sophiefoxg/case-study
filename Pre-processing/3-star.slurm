#!/bin/bash
#SBATCH --job-name=star
#SBATCH --output=star.out
#SBATCH --error=star.err
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --mail-user=Fox-GmuerSE@cardiff.ac.uk
#SBATCH --mail-type=END,FAIL

## Load some Modules
module load star/2.7.11b-conda

## Useful shortcuts
#varibles
workdir=$(pwd)

rawdir=${workdir}/rawdata
trimdir=${workdir}/trimdata
mkdir ${workdir}/star
stardir=${workdir}/star
genomedir=~/mydata/case-study/Pre-processing/Genome/

for f in ${workdir}/rawdata/*_1.fq.gz
do
R1=$(basename $f | cut -f1 -d.)
base=$(echo $R1 | sed 's/_1//')

# map forward and reverse reads to genome
STAR   --outMultimapperOrder Random \
       --outSAMmultNmax 1 \
       --runThreadN 8 \
       --runMode alignReads \
       --outSAMtype BAM Unsorted \
       --quantMode GeneCounts \
       --readFilesCommand zcat \
       --outFileNamePrefix ${stardir}/${base}-unsort. \
       --genomeDir ${genomedir} \
       --readFilesIn ${trimdir}/${base}_trim_R1.fq.gz ${trimdir}/${base}_trim_R2.fq.gz

done

module unload star/2.7.11b-conda

