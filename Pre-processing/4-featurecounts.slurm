#!/bin/bash
#SBATCH --job-name=fc
#SBATCH --output=fc.out
#SBATCH --error=fc.err
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --mail-user=Fox-GmuerSE@cardiff.ac.uk
#SBATCH --mail-type=END,FAIL

# Write jobscript to output file (good for reproducibility)
cat $0

# Load some modules
module load subread/2.0.6-abbqxcc

## Useful shortcuts
workdir=$(pwd)

rawdir=${workdir}/rawdata
trimdir=${workdir}/trimdata
stardir=${workdir}/star
genomedir=~/mydata/case-study/Pre-processing/Genome
markdir=${workdir}/markdup
mkdir ${workdir}/featureCounts
fcdir=${workdir}/featureCounts

for f in ${workdir}/rawdata/*_1.fq.gz
do
R1=$(basename $f | cut -f1 -d.)
base=$(echo $R1 | sed 's/_1//')

featureCounts \
        -T 4 -p -F GTF -t exon -g gene_id \
	-a ${genomedir}/Mus_musculus_balbcj.BALB_cJ_v1.113.gtf \
	-o ${fcdir}/${base}.featurecount \
	${stardir}/${base}*.bam


done
