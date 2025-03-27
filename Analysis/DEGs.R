################################################################################
### Identify and Filter Differentially Expressed (DE) Genes
### Author: Sophie Fox-Gmuer (edited from Sarah Christofedes)
### This script loads saved DESeq2 results (.Rdata), identifies DE genes, 
### and organises them for analysis. 

setwd("~/mydata/case-study/Analysis")

##Mark genes that are differentially expressed and filter 


load("~/mydata/case-study/Analysis/DietvsWorm2.RData")

#Put the DESeq2 results objects into a list called deData
deData<-out.DESeq2$results

#Create a factor indicating which genes are DE
#For a gene to be DE in a pairwise comparison it needs to be a) significant overall and b) have a log2FC >1
indicateDEGS<-function(deData){
  DE<-deData$padj <= 0.05 
  DE[is.na(DE)]<-FALSE
  DE[deData$log2FoldChange < 1 & deData$log2FoldChange > -1]<-FALSE 
  return(DE)
}

for (i in names(deData)){
  deData[[i]]$DE<-indicateDEGS(deData[[i]])
}

####Filter the Infected_vs_Control to see DE genes 
View(as.data.frame(deData$Infected_vs_Control[deData$Infected_vs_Control$DE==T,]))

write.csv(as.data.frame(deData$Infected_vs_Control), "Infected_vs_Control.csv")

#looking at data
data<-read.csv("Infected_vs_Control.csv")


print(head(data))


# Extract relevant columns and clean gene names
data_all <- read.csv("Infected_vs_Control.csv")
colnames(data_all)[1] <- "Gene"
data_all <- data_all[, c("Gene", "log2FoldChange", "padj")]


head(data_all)
significant_all <- data_all[!is.na(data_all$padj) & data_all$padj < 0.05, ]  
upregulated_all <- significant_all[significant_all$log2FoldChange > 1, ]  
downregulated_all <- significant_all[significant_all$log2FoldChange < -1, ] 

write.csv(significant_all, "all_significant.csv", row.names = FALSE, quote=FALSE)
write.csv(upregulated_all, "upregulated.csv", row.names = FALSE, quote = FALSE)
write.csv(downregulated_all, "downregulated.csv", row.names = FALSE, quote = FALSE)
