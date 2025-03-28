#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

#note portions of the script were optimised using chatgpt

# GSEA Explorer Shiny App with Gene Names 

library(shiny)
library(clusterProfiler)
library(dplyr)
library(DT)
library(ggplot2)
library(org.Mm.eg.db)

# Load GSEA results (with ENTREZ IDs in core_enrichment)
gsea_data <- read.csv("~/mydata/case-study/Analysis/GSEA_GO_results.csv")
expression_data <- read.csv("~/mydata/case-study/Analysis/Infected_vs_Control_ENS.csv")

# Map Ensembl IDs in expression data to gene symbols
ensembl_to_symbol <- AnnotationDbi::select(org.Mm.eg.db, 
                                           keys = expression_data$Mouse_gene_ID, 
                                           columns = c("SYMBOL"), 
                                           keytype = "ENSEMBL")
expression_data <- merge(expression_data, ensembl_to_symbol, by.x = "Mouse_gene_ID", by.y = "ENSEMBL")

# Helper function to convert Entrez IDs to gene symbols
decode_entrez_to_symbol <- function(entrez_ids) {
  symbols <- AnnotationDbi::mapIds(org.Mm.eg.db, 
                                   keys = entrez_ids, 
                                   column = "SYMBOL", 
                                   keytype = "ENTREZID", 
                                   multiVals = "first")
  return(symbols)
}

extract_genes <- function(core_enrichment_str) {
  unlist(strsplit(as.character(core_enrichment_str), "/"))
}

# UI
ui <- fluidPage(
  titlePanel("GSEA Pathway Explorer"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("pathway", "Choose a GO Pathway:", 
                  choices = setNames(gsea_data$ID, gsea_data$Description)),
      numericInput("padj_filter", "Max adjusted p-value:", value = 0.05, step = 0.01),
      numericInput("logfc_up", "Min log2FC (Upregulated):", value = 0.5, step = 0.1),
      numericInput("logfc_down", "Max log2FC (Downregulated):", value = -0.5, step = 0.1),
      helpText("Select a GO term and apply filters to view the key genes driving enrichment.")
    ),
    
    mainPanel(
      h4("Leading-edge Genes in Pathway"),
      DTOutput("gene_table"),
      h4("Expression of Leading-edge Genes"),
      plotOutput("expression_plot")
    )
  )
)

# Server
server <- function(input, output, session) {
  
  selected_pathway_data <- reactive({
    gsea_data[gsea_data$ID == input$pathway, ]
  })
  
  leading_edge_genes <- reactive({
    entrez_ids <- extract_genes(selected_pathway_data()$core_enrichment)
    gene_symbols <- decode_entrez_to_symbol(entrez_ids)
    expression_data_filtered <- expression_data %>% 
      filter(SYMBOL %in% gene_symbols) %>%
      filter(padj <= input$padj_filter) %>%
      filter(log2FoldChange >= input$logfc_up | log2FoldChange <= input$logfc_down)
    expression_data_filtered
  })
  
  output$gene_table <- renderDT({
    leading_edge_genes() %>% 
      select(SYMBOL, log2FoldChange, padj) %>% 
      rename(Gene = SYMBOL)
  })
  
  output$expression_plot <- renderPlot({
    df <- leading_edge_genes()
    ggplot(df, aes(x = reorder(SYMBOL, log2FoldChange), y = log2FoldChange)) +
      geom_col(fill = "steelblue") +
      coord_flip() +
      labs(x = "Gene", y = "log2 Fold Change", title = "Gene Expression") +
      theme_minimal()
  })
}

# Run the app
shinyApp(ui, server)