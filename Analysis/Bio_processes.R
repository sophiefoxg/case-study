################################################################################
### Gene Set Enrichment Visualization for GO, KEGG, and GSEA (Mouse Data)
### Author: Sophie Fox-Gmuer
################################################################################

# 0️⃣ Install required packages if not already installed
packages <- c("BiocManager", "clusterProfiler", "org.Mm.eg.db", "AnnotationDbi", "enrichplot", "ggplot2", "DOSE", "ggrepel", "biomaRt", "GOplot")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "GOplot") install.packages(pkg) 
    else BiocManager::install(pkg)
  }
}
lapply(packages, library, character.only = TRUE)


# Create results folder if it doesn't exist
results_dir <- "~/mydata/case-study/Results"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

#  Load and prepare data
upregulated <- read.csv("~/mydata/case-study/Analysis/upregulated_ensembl.csv", stringsAsFactors = FALSE)
downregulated <- read.csv("~/mydata/case-study/Analysis/downregulated_ensembl.csv", stringsAsFactors = FALSE)

convert_to_gene_name <- function(df) {
  gene_names <- mapIds(org.Mm.eg.db,  
                       keys = df$projection_parent_gene,
                       column = "SYMBOL",  
                       keytype = "ENSEMBL",
                       multiVals = "first")  
  df$GeneName <- gene_names
  df <- df[!is.na(df$GeneName), ]  
  df <- df[, c("GeneName", "Log2FoldChange", "p_value")]  
  return(df)
}

convert_to_entrez <- function(df) {
  entrez_ids <- mapIds(org.Mm.eg.db,  
                       keys = df$GeneName,
                       column = "ENTREZID",
                       keytype = "SYMBOL",
                       multiVals = "first")  
  df$EntrezID <- entrez_ids
  df <- df[!is.na(df$EntrezID), ]  
  df <- df[, c("GeneName", "EntrezID", "Log2FoldChange", "p_value")]  
  return(df)
}

# Apply conversions
upregulated_geneID <- convert_to_entrez(convert_to_gene_name(upregulated))
downregulated_geneID <- convert_to_entrez(convert_to_gene_name(downregulated))

#  GO Enrichment (BP)
go_up <- enrichGO(gene = upregulated_geneID$GeneName, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
go_down <- enrichGO(gene = downregulated_geneID$GeneName, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)

# GO Upregulated Dotplot
go_up_dotplot <- dotplot(go_up, showCategory = 20, title = "GO Biological Processes Upregulated") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90")
  ) +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40))

ggsave(filename = file.path(results_dir, "GO_Upregulated_Dotplot.png"), plot = go_up_dotplot, width = 10, height = 8, dpi = 300)

# GO Upregulated Barplot
go_up_barplot <- barplot(go_up, showCategory = 10, title = "Barplot: GO Upregulated", font.size = 10)

ggsave(filename = file.path(results_dir, "GO_Upregulated_Barplot.png"), plot = go_up_barplot, width = 10, height = 8, dpi = 300)

# Save GO Biological Process results
go_up_df <- as.data.frame(go_up)


write.csv(go_up_df, "~/mydata/case-study/Results/go_upregulated_BP.csv", row.names = FALSE)


#concept
library(enrichplot)


cnet_up <- cnetplot(
  go_up,
  showCategory = 5,
  circular = TRUE,
  color.params = list(
    foldChange = NULL,     
    edge = TRUE            
  )
)
#  named vector of log2FC for your upregulated genes
fc_vector <- upregulated_geneID$Log2FoldChange
names(fc_vector) <- upregulated_geneID$GeneName


cnet_up <- cnetplot(
  go_up,
  showCategory = 5,
  circular = TRUE,
  color.params = list(
    foldChange = fc_vector,
    edge = TRUE
  )
)

ggsave("~/mydata/case-study/Results/GO_Upregulated_Cnetplot.png", cnet_up, width = 10, height = 8, dpi = 300)




#  KEGG Enrichment
kegg_up <- enrichKEGG(gene = upregulated_geneID$EntrezID, organism = "mmu", pvalueCutoff = 0.05)
kegg_down <- enrichKEGG(gene = downregulated_geneID$EntrezID, organism = "mmu", pvalueCutoff = 0.05)

#  KEGG Visualizations


# KEGG Upregulated Dotplot
kegg_up_dotplot <- dotplot(kegg_up, showCategory = 15, title = "KEGG Upregulated")

ggsave(filename = file.path(results_dir, "KEGG_Upregulated_Dotplot.png"), plot = kegg_up_dotplot, width = 10, height = 8, dpi = 300)

kegg_up_df <- as.data.frame(kegg_up)

write.csv(kegg_up_df, "~/mydata/case-study/Results/kegg_upregulated.csv", row.names = FALSE)
head(kegg_up_df)

if (!requireNamespace("pathview", quietly = TRUE)) BiocManager::install("pathview")
library(pathview)

# Prepare logFC named vector for pathview
logFC_vector <- upregulated_geneID$Log2FoldChange
names(logFC_vector) <- upregulated_geneID$EntrezID

head(logFC_vector)



pathview(
  gene.data    = logFC_vector,
  pathway.id   = "mmu04610",  
  species      = "mmu",
  out.suffix   = "Upregulated",
  kegg.native  = TRUE,
  same.layer   = FALSE,
  output.path  = results_dir
)
# Run pathview for the selected KEGG pathway
pathview(
  gene.data    = logFC_vector,
  pathway.id   = "mmu04610",  # Complement and coagulation cascades
  species      = "mmu",
  out.suffix   = "Upregulated",
  kegg.native  = TRUE,
  same.layer   = FALSE,
  output.path  = results_dir
)



gene_list <- setNames(upregulated_geneID$Log2FoldChange, upregulated_geneID$GeneName)
gene_list_up <- sort(gene_list_up, decreasing = TRUE)

gsea_results_up <- gseGO(
  geneList = gene_list_up,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  scoreType = "pos"
)
head (gsea_results_up)
write.csv(as.data.frame(gsea_results_up), "GSEA_GO_Results.csv", row.names = FALSE)


# Combine Up and Downregulated Genes
combined_gene_data <- rbind(upregulated_geneID, downregulated_geneID)

combined_gene_data <- combined_gene_data[!duplicated(combined_gene_data$GeneName), ]

# ranked list
gene_list_combined <- setNames(combined_gene_data$Log2FoldChange, combined_gene_data$GeneName)
gene_list_combined <- sort(gene_list_combined, decreasing = TRUE)

# Run GSEA on Combined Gene List
gsea_results_combined <- gseGO(
  geneList = gene_list_combined,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  scoreType = "std"  # accepts both pos and neg scores
)

head(gsea_results_combined)




#Save GSEA results
write.csv(as.data.frame(gsea_results_combined), "GSEA.csv", row.names = FALSE)


