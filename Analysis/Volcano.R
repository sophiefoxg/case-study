################################################################################
### Volcano plots for differentially expressed genes 
### Author: Sophie Fox-Gmuer
### Loads the results from the DEG filtering script to create a volcano plot and labels the up and downregulated genes. 
# Note: Portions of this script were generated or optimized with ChatGPT.
#
################################################################################

# Loading libraries
install.packages("ggplot2")
library(ggplot2)
install.packages("ggrepel")
library(ggrepel)

# Load differential expression data
deg <- read.csv("~/mydata/case-study/Analysis/Infected_vs_Control.csv")
deg <- deg[!is.na(deg$padj) & !is.na(deg$log2FoldChange), ]

# Log2fold change and Padj cutoffs
deg$expression <- ifelse(deg$log2FoldChange > 1 & deg$padj < 0.05, "Upregulated",
                         ifelse(deg$log2FoldChange < -1 & deg$padj < 0.05, "Downregulated", "Neutral"))

#  gene labels 
up_labels <- read.csv("~/mydata/case-study/Analysis/upregulated_geneID.csv")
down_labels <- read.csv("~/mydata/case-study/Analysis/downregulated_geneID.csv")

# Merge labels with DEG data based on Log2FoldChange and p_value
deg <- merge(deg, up_labels, by.x = c("log2FoldChange", "padj"), by.y = c("Log2FoldChange", "p_value"), all.x = TRUE)
deg <- merge(deg, down_labels, by.x = c("log2FoldChange", "padj"), by.y = c("Log2FoldChange", "p_value"), all.x = TRUE, suffixes = c("_up", "_down"))


deg$label <- ifelse(!is.na(deg$GeneName_up), deg$GeneName_up, 
                    ifelse(!is.na(deg$GeneName_down), deg$GeneName_down, ""))

# Select top upregulated and downregulated genes for labeling
ten_up <- deg[deg$expression == "Upregulated", ]
ten_up <- ten_up[order(-ten_up$log2FoldChange), ][1:10, ]

ten_down <- deg[deg$expression == "Downregulated", ]
ten_down <- ten_down[order(ten_down$log2FoldChange), ][1:10, ]

labels <- rbind(ten_down, ten_up)

#  volcano plot 
volcano <- ggplot(deg, aes(x = log2FoldChange, y = -log10(padj), color = expression)) +
  geom_point(alpha = 0.8, size = 2) +
  geom_label_repel(
    data = labels, aes(label = label), 
    size = 3, color = "black", fill = "white", 
    box.padding = 0.5, point.padding = 0.3, 
    segment.color = "black", max.overlaps = Inf
  ) +
 
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Neutral" = "grey50")) +
  labs(title = "Differentially Expressed Genes",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.title = element_blank(),
    legend.position = "top",
    panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.25)
  )

# Plot +save
plot(volcano)

ggsave(filename = file.path(results_dir, "Volcano_plot.png"), plot = volcano, width = 10, height =8, dpi =300 )

#ALOX15 and ALOX5

View(data.frame(upregulated_geneID))


library(ggplot2)
library(dplyr)

# Read the CSV file
up_genes <- read.csv("~/mydata/case-study/Analysis/upregulated_geneID.csv")

# Filter for Gene of choice
alox_genes <- up_genes %>%
  filter(GeneName %in% c("Alox5", "Alox15"))

# Add significance stars
alox_genes$signif <- cut(alox_genes$p_value,
                         breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                         labels = c("***", "**", "*", "ns"))

# Bar plot
barplot<- ggplot(alox_genes, aes(x = GeneName, y = Log2FoldChange, fill = GeneName)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = round(Log2FoldChange, 2)), vjust = 1.5, color = "white", size = 4.5, fontface = "bold") +
  geom_text(aes(label = signif), vjust = -0.5, size = 5) +
  labs(title = "ALOX Genes Differential Expression",
       x = "Gene", y = "Log2 Fold Change") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

#save
ggsave(filename = file.path(results_dir, "ALOX5 + ALOX15.png"), plot= barplot, width = 10, height = 8, dpi = 300)



