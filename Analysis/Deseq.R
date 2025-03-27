################################################################################
## RNA-seq analysis script using SARTools and DESeq2
### Author: Sophie Fox-Gmuer
### Date: Adapted from November 28th, 2019 (Hugo Varet) 
################################################################################

# Install and load  packages devtools 2.4.5 SARTools 1.8.1
install.packages("devtools")
devtools::install_github("PF2-pasteur-fr/SARTools", build_opts="--no-resave-data")
install.packages("SARTools")

library(devtools)
library(SARTools)


#optional
rm(list=ls())  # Clear workspace 

#####set working directory (this is an example)
workDir<-("~/mydata/case-study/Analysis")

#####CHANGE THE project TITLE ####
projectName <- "DietvsWorm"

author <- "Sophie Fox-Gmuer"     # Author name for documentation purposes


# Path to the experimental design file
targetFile <- "~/mydata/case-study/Pre-processing/featureCounts/Targets.txt"  

##### Change path to feature count files 
rawDir <- "~/mydata/case-study/Pre-processing/featureCounts" 



# Analysis settings
varInt <- "Infection"       # Variable of interest in the dataset
condRef <- "Control"        # Reference condition for comparisons

batch <- NULL               # Optional batch effect variable (set to NULL if not applicable)

# Filter unwanted features
featuresToRemove <- c("alignment_not_unique", "not_aligned", "too_low_aQual")

# DESeq2 statistical options
fitType <- "parametric"      # Method for dispersion estimation
cooksCutoff <- TRUE          # Enable/disable outlier detection
independentFiltering <- TRUE # Enable/disable independent filtering
alpha <- 0.05                # Significance threshold for adjusted p-values
pAdjustMethod <- "BH"        # Method for multiple testing correction

# Transformation and visualisation options
typeTrans <- "VST"           # Data transformation for visualisation (VST or rlog)
locfunc <- "median"          # Method for estimating size factors
colors <- c("#f3c300", "#875692", "#f38400", "#a1caf1", "#be0032", "#c2b280",
            "#848482", "#008856", "#e68fac", "#0067a5", "#f99379", "#604e97")
forceCairoGraph <- FALSE     # Set to TRUE for Cairo graphics

################################################################################
######
################################################################################

setwd(workDir)  # Set the working directory
if (forceCairoGraph) options(bitmapType="cairo")  # Set Cairo for graphics rendering if needed

# Validate parameters and settings
checkParameters.DESeq2(
  projectName=projectName, author=author, targetFile=targetFile,
  rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
  condRef=condRef, batch=batch, fitType=fitType, cooksCutoff=cooksCutoff,
  independentFiltering=independentFiltering, alpha=alpha, pAdjustMethod=pAdjustMethod,
  typeTrans=typeTrans, locfunc=locfunc, colors=colors
)

# Load target file with experimental design
target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)

# Import raw count data
counts <- loadCountData(target=target, rawDir=rawDir)

# Generate initial diagnostic plots
majSequences <- descriptionPlots(counts=counts, group=target[, varInt], col=colors)

# Perform DESeq2 analysis
out.DESeq2 <- run.DESeq2(
  counts=counts, target=target, varInt=varInt, batch=batch, locfunc=locfunc,
  fitType=fitType, pAdjustMethod=pAdjustMethod, cooksCutoff=cooksCutoff,
  independentFiltering=independentFiltering, alpha=alpha
)

# PCA and clustering for exploratory analysis
exploreCounts(object=out.DESeq2$dds, group=target[, varInt], typeTrans=typeTrans, col=colors)

# Summarise and visualise the analysis
summaryResults <- summarizeResults.DESeq2(
  out.DESeq2, group=target[, varInt], col=colors,
  independentFiltering=independentFiltering, cooksCutoff=cooksCutoff, alpha=alpha
)

# Save R session to allow further downstream analysis without re-running the script
save.image(file=paste0(gsub(" ", "", projectName), ".RData"))

# Generate a detailed HTML report summarising info 

writeReport.DESeq2(
  target=target, counts=counts, out.DESeq2=out.DESeq2, summaryResults=summaryResults,
  majSequences=majSequences, workDir=workDir,projectName=projectName, author=author,
  targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
  condRef=condRef, batch=batch, fitType=fitType, cooksCutoff=cooksCutoff,
  independentFiltering=independentFiltering, alpha=alpha, pAdjustMethod=pAdjustMethod,
  typeTrans=typeTrans, locfunc=locfunc, colors=colors
)


