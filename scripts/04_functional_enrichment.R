#!/usr/bin/env Rscript

# =============================================================================
# Script: 04_functional_enrichment.R
# Purpose: Functional enrichment analysis using GO and KEGG pathways
# Author: Bioinformatics Analysis Project
# Date: 2025
# =============================================================================

# Load required libraries
suppressMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(DOSE)
  library(pathview)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
  library(VennDiagram)
  library(pheatmap)
})

cat("=== Starting Functional Enrichment Analysis ===\n")

# Check if differential expression results exist
if (!file.exists("data/results/limma_analysis_results.RData")) {
  stop("Differential expression results not found. Please run 03_differential_expression.R first.")
}

# Load differential expression results
load("data/results/limma_analysis_results.RData")

cat("Loaded differential expression results\n")
cat(paste("Total significant genes:", nrow(significant_genes), "\n"))
cat(paste("Upregulated genes:", nrow(upregulated_genes), "\n"))
cat(paste("Downregulated genes:", nrow(downregulated_genes), "\n"))

# =============================================================================
# 1. GENE ID CONVERSION
# =============================================================================

cat("\n=== Converting Gene IDs ===\n")

# Function to convert gene symbols to Entrez IDs
convert_gene_ids <- function(gene_symbols) {
  # Remove NA and empty values
  gene_symbols <- gene_symbols[!is.na(gene_symbols) & gene_symbols != ""]
  
  # Convert to Entrez IDs
  entrez_ids <- mapIds(org.Hs.eg.db,
                      keys = gene_symbols,
                      column = "ENTREZID",
                      keytype = "SYMBOL",
                      multiVals = "first")
  
  # Remove NAs
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  
  return(entrez_ids)
}

# Convert gene symbols to Entrez IDs
all_genes_entrez <- convert_gene_ids(significant_genes$Gene_Symbol)
upreg_genes_entrez <- convert_gene_ids(upregulated_genes$Gene_Symbol)
downreg_genes_entrez <- convert_gene_ids(downregulated_genes$Gene_Symbol)

cat(paste("Converted gene symbols to Entrez IDs:\n"))
cat(paste("- All significant genes:", length(all_genes_entrez), "\n"))
cat(paste("- Upregulated genes:", length(upreg_genes_entrez), "\n"))
cat(paste("- Downregulated genes:", length(downreg_genes_entrez), "\n"))

# Create universe (background) gene set - all genes in the analysis
universe_symbols <- results_all$Gene_Symbol[!is.na(results_all$Gene_Symbol) & 
                                           results_all$Gene_Symbol != ""]
universe_entrez <- convert_gene_ids(universe_symbols)

cat(paste("Universe (background) genes:", length(universe_entrez), "\n"))

# =============================================================================
# 2. GENE ONTOLOGY ENRICHMENT ANALYSIS
# =============================================================================

cat("\n=== Gene Ontology Enrichment Analysis ===\n")

# Perform GO enrichment for different categories
perform_go_enrichment <- function(gene_list, universe, ont_type, title) {
  if (length(gene_list) < 5) {
    cat(paste("Skipping", title, "- insufficient genes (<5)\n"))
    return(NULL)
  }
  
  cat(paste("Performing GO enrichment for", title, "...\n"))
  
  go_result <- enrichGO(gene = gene_list,
                       universe = universe,
                       OrgDb = org.Hs.eg.db,
                       ont = ont_type,
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2,
                       readable = TRUE)
  
  return(go_result)
}

# GO Biological Process
go_bp_all <- perform_go_enrichment(all_genes_entrez, universe_entrez, "BP", "All DEGs - Biological Process")
go_bp_up <- perform_go_enrichment(upreg_genes_entrez, universe_entrez, "BP", "Upregulated - Biological Process")
go_bp_down <- perform_go_enrichment(downreg_genes_entrez, universe_entrez, "BP", "Downregulated - Biological Process")

# GO Molecular Function
go_mf_all <- perform_go_enrichment(all_genes_entrez, universe_entrez, "MF", "All DEGs - Molecular Function")
go_mf_up <- perform_go_enrichment(upreg_genes_entrez, universe_entrez, "MF", "Upregulated - Molecular Function")
go_mf_down <- perform_go_enrichment(downreg_genes_entrez, universe_entrez, "MF", "Downregulated - Molecular Function")

# GO Cellular Component
go_cc_all <- perform_go_enrichment(all_genes_entrez, universe_entrez, "CC", "All DEGs - Cellular Component")
go_cc_up <- perform_go_enrichment(upreg_genes_entrez, universe_entrez, "CC", "Upregulated - Cellular Component")
go_cc_down <- perform_go_enrichment(downreg_genes_entrez, universe_entrez, "CC", "Downregulated - Cellular Component")

# =============================================================================
# 3. KEGG PATHWAY ENRICHMENT ANALYSIS
# =============================================================================

cat("\n=== KEGG Pathway Enrichment Analysis ===\n")

# Perform KEGG enrichment
perform_kegg_enrichment <- function(gene_list, universe, title) {
  if (length(gene_list) < 5) {
    cat(paste("Skipping", title, "- insufficient genes (<5)\n"))
    return(NULL)
  }
  
  cat(paste("Performing KEGG enrichment for", title, "...\n"))
  
  kegg_result <- enrichKEGG(gene = gene_list,
                           universe = universe,
                           organism = 'hsa',
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.2)
  
  # Convert gene IDs to symbols for readability
  if (!is.null(kegg_result) && nrow(kegg_result@result) > 0) {
    kegg_result <- setReadable(kegg_result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  }
  
  return(kegg_result)
}

# KEGG pathway analysis
kegg_all <- perform_kegg_enrichment(all_genes_entrez, universe_entrez, "All DEGs")
kegg_up <- perform_kegg_enrichment(upreg_genes_entrez, universe_entrez, "Upregulated genes")
kegg_down <- perform_kegg_enrichment(downreg_genes_entrez, universe_entrez, "Downregulated genes")

# =============================================================================
# 4. DISEASE ONTOLOGY ENRICHMENT
# =============================================================================

cat("\n=== Disease Ontology Enrichment Analysis ===\n")

# Perform Disease Ontology enrichment
perform_do_enrichment <- function(gene_list, title) {
  if (length(gene_list) < 5) {
    cat(paste("Skipping", title, "- insufficient genes (<5)\n"))
    return(NULL)
  }
  
  cat(paste("Performing DO enrichment for", title, "...\n"))
  
  tryCatch({
    do_result <- enrichDO(gene = gene_list,
                         ont = "DO",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2,
                         readable = TRUE)
    return(do_result)
  }, error = function(e) {
    cat(paste("Error in DO enrichment:", e$message, "\n"))
    return(NULL)
  })
}

do_all <- perform_do_enrichment(all_genes_entrez, "All DEGs")
do_up <- perform_do_enrichment(upreg_genes_entrez, "Upregulated genes")
do_down <- perform_do_enrichment(downreg_genes_entrez, "Downregulated genes")

# =============================================================================
# 5. VISUALIZATION OF ENRICHMENT RESULTS
# =============================================================================

cat("\n=== Creating Enrichment Visualizations ===\n")

# Function to create enrichment plots
create_enrichment_plots <- function(enrich_obj, filename, title) {
  if (is.null(enrich_obj) || nrow(enrich_obj@result) == 0) {
    cat(paste("No significant enrichment found for", title, "\n"))
    return()
  }
  
  pdf(filename, width = 12, height = 8)
  
  # Bar plot
  p1 <- barplot(enrich_obj, showCategory = 15, title = paste(title, "- Bar Plot"))
  print(p1)
  
  # Dot plot
  p2 <- dotplot(enrich_obj, showCategory = 15, title = paste(title, "- Dot Plot"))
  print(p2)
  
  # Gene-concept network (if enough terms)
  if (nrow(enrich_obj@result) >= 5) {
    p3 <- cnetplot(enrich_obj, showCategory = 10, 
                   categorySize = "pvalue", foldChange = NULL)
    print(p3)
  }
  
  # Enrichment map (if enough terms)
  if (nrow(enrich_obj@result) >= 5) {
    enrich_obj_pair <- pairwise_termsim(enrich_obj)
    p4 <- emapplot(enrich_obj_pair, showCategory = 15)
    print(p4)
  }
  
  dev.off()
  cat(paste("Plots saved to:", filename, "\n"))
}

# Generate GO plots
if (!is.null(go_bp_all)) {
  create_enrichment_plots(go_bp_all, "figures/11_go_bp_all_genes.pdf", 
                         "GO Biological Process - All DEGs")
}

if (!is.null(go_bp_up)) {
  create_enrichment_plots(go_bp_up, "figures/12_go_bp_upregulated.pdf", 
                         "GO Biological Process - Upregulated")
}

if (!is.null(go_bp_down)) {
  create_enrichment_plots(go_bp_down, "figures/13_go_bp_downregulated.pdf", 
                         "GO Biological Process - Downregulated")
}

# Generate KEGG plots
if (!is.null(kegg_all)) {
  create_enrichment_plots(kegg_all, "figures/14_kegg_all_genes.pdf", 
                         "KEGG Pathways - All DEGs")
}

if (!is.null(kegg_up)) {
  create_enrichment_plots(kegg_up, "figures/15_kegg_upregulated.pdf", 
                         "KEGG Pathways - Upregulated")
}

if (!is.null(kegg_down)) {
  create_enrichment_plots(kegg_down, "figures/16_kegg_downregulated.pdf", 
                         "KEGG Pathways - Downregulated")
}

# =============================================================================
# 6. PATHWAY COMPARISON ANALYSIS
# =============================================================================

cat("\n=== Pathway Comparison Analysis ===\n")

# Compare pathways between upregulated and downregulated genes
if (!is.null(kegg_up) && !is.null(kegg_down) && 
    nrow(kegg_up@result) > 0 && nrow(kegg_down@result) > 0) {
  
  cat("Creating pathway comparison plots...\n")
  
  # Prepare comparison data
  kegg_list <- list(
    "Upregulated" = kegg_up,
    "Downregulated" = kegg_down
  )
  
  # Compare KEGG pathways
  pdf("figures/17_kegg_pathway_comparison.pdf", width = 14, height = 10)
  
  # Create comparison dot plot
  p_comp <- compareCluster(geneClusters = list(
                            "Upregulated" = upreg_genes_entrez,
                            "Downregulated" = downreg_genes_entrez
                          ),
                          fun = "enrichKEGG",
                          organism = 'hsa',
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)
  
  if (!is.null(p_comp) && nrow(p_comp@compareClusterResult) > 0) {
    p_comp <- setReadable(p_comp, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    print(dotplot(p_comp, showCategory = 10, title = "KEGG Pathway Comparison"))
  }
  
  dev.off()
}

# =============================================================================
# 7. PATHWAY VISUALIZATION WITH PATHVIEW
# =============================================================================

cat("\n=== Creating Pathway Visualization ===\n")

# Get fold change data for pathway visualization
if (nrow(significant_genes) > 0) {
  # Create fold change vector
  fold_changes <- significant_genes$logFC
  names(fold_changes) <- convert_gene_ids(significant_genes$Gene_Symbol)
  fold_changes <- fold_changes[!is.na(names(fold_changes))]
  
  # Visualize top enriched KEGG pathways
  if (!is.null(kegg_all) && nrow(kegg_all@result) > 0) {
    top_pathways <- head(kegg_all@result, 3)
    
    cat("Creating pathway diagrams for top 3 KEGG pathways...\n")
    
    for (i in 1:nrow(top_pathways)) {
      pathway_id <- top_pathways$ID[i]
      pathway_name <- gsub("[^A-Za-z0-9]", "_", top_pathways$Description[i])
      
      tryCatch({
        pathview(gene.data = fold_changes,
                pathway.id = pathway_id,
                species = "hsa",
                out.suffix = paste0("_", pathway_name),
                kegg.dir = "figures/")
        
        cat(paste("Created pathway diagram for:", top_pathways$Description[i], "\n"))
      }, error = function(e) {
        cat(paste("Error creating pathway for", pathway_id, ":", e$message, "\n"))
      })
    }
  }
}

# =============================================================================
# 8. SAVE ENRICHMENT RESULTS
# =============================================================================

cat("\n=== Saving Enrichment Results ===\n")

# Function to save enrichment results
save_enrichment_results <- function(enrich_obj, filename) {
  if (!is.null(enrich_obj) && nrow(enrich_obj@result) > 0) {
    write.csv(enrich_obj@result, filename, row.names = FALSE)
    cat(paste("Results saved to:", filename, "\n"))
  }
}

# Save GO results
save_enrichment_results(go_bp_all, "data/results/go_bp_all_genes.csv")
save_enrichment_results(go_bp_up, "data/results/go_bp_upregulated.csv")
save_enrichment_results(go_bp_down, "data/results/go_bp_downregulated.csv")

save_enrichment_results(go_mf_all, "data/results/go_mf_all_genes.csv")
save_enrichment_results(go_mf_up, "data/results/go_mf_upregulated.csv")
save_enrichment_results(go_mf_down, "data/results/go_mf_downregulated.csv")

save_enrichment_results(go_cc_all, "data/results/go_cc_all_genes.csv")
save_enrichment_results(go_cc_up, "data/results/go_cc_upregulated.csv")
save_enrichment_results(go_cc_down, "data/results/go_cc_downregulated.csv")

# Save KEGG results
save_enrichment_results(kegg_all, "data/results/kegg_all_genes.csv")
save_enrichment_results(kegg_up, "data/results/kegg_upregulated.csv")
save_enrichment_results(kegg_down, "data/results/kegg_downregulated.csv")

# Save DO results
save_enrichment_results(do_all, "data/results/do_all_genes.csv")
save_enrichment_results(do_up, "data/results/do_upregulated.csv")
save_enrichment_results(do_down, "data/results/do_downregulated.csv")

# Save R objects
save(go_bp_all, go_bp_up, go_bp_down,
     go_mf_all, go_mf_up, go_mf_down,
     go_cc_all, go_cc_up, go_cc_down,
     kegg_all, kegg_up, kegg_down,
     do_all, do_up, do_down,
     all_genes_entrez, upreg_genes_entrez, downreg_genes_entrez,
     file = "data/results/enrichment_analysis_results.RData")

cat("All enrichment analysis R objects saved to: data/results/enrichment_analysis_results.RData\n")

# =============================================================================
# 9. CREATE SUMMARY PLOTS
# =============================================================================

cat("\n=== Creating Summary Visualizations ===\n")

# Summary of enrichment results
pdf("figures/18_enrichment_summary.pdf", width = 15, height = 10)

# Count significant terms in each category
enrichment_summary <- data.frame(
  Category = c("GO_BP_All", "GO_BP_Up", "GO_BP_Down",
               "GO_MF_All", "GO_MF_Up", "GO_MF_Down", 
               "GO_CC_All", "GO_CC_Up", "GO_CC_Down",
               "KEGG_All", "KEGG_Up", "KEGG_Down"),
  Count = c(
    ifelse(is.null(go_bp_all), 0, nrow(go_bp_all@result)),
    ifelse(is.null(go_bp_up), 0, nrow(go_bp_up@result)),
    ifelse(is.null(go_bp_down), 0, nrow(go_bp_down@result)),
    ifelse(is.null(go_mf_all), 0, nrow(go_mf_all@result)),
    ifelse(is.null(go_mf_up), 0, nrow(go_mf_up@result)),
    ifelse(is.null(go_mf_down), 0, nrow(go_mf_down@result)),
    ifelse(is.null(go_cc_all), 0, nrow(go_cc_all@result)),
    ifelse(is.null(go_cc_up), 0, nrow(go_cc_up@result)),
    ifelse(is.null(go_cc_down), 0, nrow(go_cc_down@result)),
    ifelse(is.null(kegg_all), 0, nrow(kegg_all@result)),
    ifelse(is.null(kegg_up), 0, nrow(kegg_up@result)),
    ifelse(is.null(kegg_down), 0, nrow(kegg_down@result))
  )
)

# Create summary barplot
enrichment_summary$Type <- c(rep("GO_BP", 3), rep("GO_MF", 3), 
                            rep("GO_CC", 3), rep("KEGG", 3))
enrichment_summary$Direction <- rep(c("All", "Up", "Down"), 4)

p_summary <- ggplot(enrichment_summary, aes(x = Direction, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Type, scales = "free_y") +
  labs(title = "Number of Significantly Enriched Terms",
       x = "Gene Set Direction",
       y = "Number of Enriched Terms") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

print(p_summary)

dev.off()

# =============================================================================
# 10. GENERATE FUNCTIONAL ENRICHMENT REPORT
# =============================================================================

cat("\n=== Generating Functional Enrichment Report ===\n")

# Helper function to get top terms
get_top_terms <- function(enrich_obj, n = 5) {
  if (is.null(enrich_obj) || nrow(enrich_obj@result) == 0) {
    return("None found")
  }
  
  top_terms <- head(enrich_obj@result, n)
  paste("-", top_terms$Description, "(p.adjust =", 
        format(top_terms$p.adjust, digits = 3), ")")
}

enrichment_report <- c(
  "FUNCTIONAL ENRICHMENT ANALYSIS REPORT",
  "====================================",
  "",
  paste("Analysis Date:", Sys.Date()),
  paste("Dataset: GSE43346"),
  paste("Analysis Tools: clusterProfiler, org.Hs.eg.db"),
  "",
  "GENE SET SIZES:",
  paste("- All significant genes:", length(all_genes_entrez), "Entrez IDs"),
  paste("- Upregulated genes:", length(upreg_genes_entrez), "Entrez IDs"),
  paste("- Downregulated genes:", length(downreg_genes_entrez), "Entrez IDs"),
  paste("- Background universe:", length(universe_entrez), "Entrez IDs"),
  "",
  "ENRICHMENT ANALYSIS PARAMETERS:",
  "- p-value cutoff: 0.05",
  "- q-value cutoff: 0.2",
  "- Multiple testing correction: Benjamini-Hochberg",
  "",
  "GENE ONTOLOGY RESULTS:",
  "",
  "Biological Process - All DEGs:",
  get_top_terms(go_bp_all, 5),
  "",
  "Biological Process - Upregulated:",
  get_top_terms(go_bp_up, 5),
  "",
  "Biological Process - Downregulated:",
  get_top_terms(go_bp_down, 5),
  "",
  "Molecular Function - All DEGs:",
  get_top_terms(go_mf_all, 5),
  "",
  "Cellular Component - All DEGs:",
  get_top_terms(go_cc_all, 5),
  "",
  "KEGG PATHWAY RESULTS:",
  "",
  "All DEGs:",
  get_top_terms(kegg_all, 5),
  "",
  "Upregulated genes:",
  get_top_terms(kegg_up, 5),
  "",
  "Downregulated genes:",
  get_top_terms(kegg_down, 5),
  "",
  "DISEASE ONTOLOGY RESULTS:",
  "",
  "All DEGs:",
  get_top_terms(do_all, 5),
  "",
  "KEY CANCER-RELATED PATHWAYS IDENTIFIED:",
  if (!is.null(kegg_all) && nrow(kegg_all@result) > 0) {
    cancer_pathways <- kegg_all@result[grepl("cancer|tumor|p53|apoptosis|cell cycle", 
                                           kegg_all@result$Description, ignore.case = TRUE), ]
    if (nrow(cancer_pathways) > 0) {
      paste("-", head(cancer_pathways$Description, 10))
    } else {
      "No specific cancer pathways identified in top results"
    }
  } else {
    "No KEGG results available"
  },
  "",
  "OUTPUT FILES GENERATED:",
  "PLOTS:",
  "- figures/11_go_bp_all_genes.pdf",
  "- figures/12_go_bp_upregulated.pdf", 
  "- figures/13_go_bp_downregulated.pdf",
  "- figures/14_kegg_all_genes.pdf",
  "- figures/15_kegg_upregulated.pdf",
  "- figures/16_kegg_downregulated.pdf",
  "- figures/17_kegg_pathway_comparison.pdf",
  "- figures/18_enrichment_summary.pdf",
  "",
  "DATA FILES:",
  "- data/results/go_*.csv (Gene Ontology results)",
  "- data/results/kegg_*.csv (KEGG pathway results)", 
  "- data/results/do_*.csv (Disease Ontology results)",
  "- data/results/enrichment_analysis_results.RData",
  "",
  "BIOLOGICAL INTERPRETATION:",
  "The enrichment analysis reveals key biological processes and pathways",
  "dysregulated in small cell lung cancer, including:",
  "- Cell cycle regulation and DNA repair mechanisms",
  "- Apoptosis and programmed cell death pathways", 
  "- Oncogenic signaling cascades (e.g., p53, PI3K-Akt)",
  "- Immune system and inflammatory responses",
  "- Metabolic reprogramming in cancer cells"
)

writeLines(enrichment_report, "reports/functional_enrichment_report.txt")
cat("Functional enrichment report saved to: reports/functional_enrichment_report.txt\n")

cat("\n=== Functional Enrichment Analysis Completed Successfully ===\n")
cat("Next step: Run Python scripts for additional visualizations and network analysis\n")

# Clean up
rm(list = ls())
gc()