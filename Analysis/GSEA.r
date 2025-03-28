# GSEA 

# Install necessary packages (only runs if missing)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

packages <- c("clusterProfiler", "org.Mm.eg.db", "enrichplot", "DOSE", "biomaRt")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
}

# Load libraries
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(DOSE)
library(biomaRt)


ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")


Genes_set_all <- read.csv("~/mydata/case-study/Analysis/Infected_vs_Control_ENS.csv")


head(Genes_set_all)

# Convert Ensembl IDs to Entrez IDs
gene_ids <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  filters = "ensembl_gene_id",
  values = Genes_set_all$Mouse_gene_ID,
  mart = ensembl
)

# Merge Entrez IDs with expression data
Genes_set_all_merged <- merge(Genes_set_all, gene_ids, by.x = "Mouse_gene_ID", by.y = "ensembl_gene_id")
Genes_set_all_merged <- Genes_set_all_merged[!is.na(Genes_set_all_merged$entrezgene_id), ]

# Prepare ranked gene list for GSEA
gene_list <- setNames(Genes_set_all_merged$log2FoldChange, Genes_set_all_merged$entrezgene_id)

# Fix duplicate IDs and ties in stats
gene_list <- gene_list[!duplicated(names(gene_list))]
gene_list <- gene_list + rnorm(length(gene_list), mean = 0, sd = 1e-6)
gene_list <- sort(gene_list, decreasing = TRUE)

# GSEA GO
gsea_go <- gseGO(
  geneList = gene_list,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  keyType = "ENTREZID",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = TRUE,
  eps = 0
)


# GSEA KEGG 
gsea_kegg <- gseKEGG(
  geneList = gene_list,
  organism = "mmu",
  keyType = "ncbi-geneid",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = TRUE,
  eps = 0
)



# Save 
write.csv(as.data.frame(gsea_go), file = "GSEA_GO_results.csv", row.names = FALSE)
write.csv(as.data.frame(gsea_kegg), file = "GSEA_KEGG_results.csv", row.names = FALSE)



# Visualisations

# Dotplots (showing both up- and downregulated pathways)
dotplot(gsea_go, showCategory = 10, split = ".sign") + ggtitle("GO GSEA - Dotplot") + facet_grid(.~.sign)
dotplot(gsea_kegg, showCategory = 15) + ggtitle("KEGG GSEA - Dotplot")


# Compute term similarity matrix
gsea_go <- pairwise_termsim(gsea_go)
gsea_kegg <- pairwise_termsim(gsea_kegg)

# Enrichment map plots
emapplot(gsea_go, showCategory = 10) + ggtitle("GO Enrichment Map")
emapplot(gsea_kegg, showCategory = 10) + ggtitle("KEGG Enrichment Map")

# Category network plots
cnetplot(gsea_go, categorySize = "pvalue", foldChange = gene_list, showCategory = 5) + ggtitle("GO Cnetplot")
cnetplot(gsea_kegg, categorySize = "pvalue", foldChange = gene_list, showCategory = 5) + ggtitle("KEGG Cnetplot")
