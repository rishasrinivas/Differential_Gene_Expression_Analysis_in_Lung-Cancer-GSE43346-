#!/usr/bin/env Rscript

# =============================================================================
# Script: 03_differential_expression.R
# Purpose: Differential expression analysis using limma for microarray data
# Author: Bioinformatics Analysis Project
# Date: 2025
# =============================================================================

# Load required libraries
suppressMessages({
  library(limma)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  library(RColorBrewer)
  library(VennDiagram)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

cat("=== Starting Differential Expression Analysis ===\n")

# Check if processed data exists
if (!file.exists("data/processed/expression_matrix_processed.csv")) {
  stop("Processed data not found. Please run 02_quality_control.R first.")
}

# Load processed data
expr_matrix <- read.csv("data/processed/expression_matrix_processed.csv", row.names = 1)
sample_metadata <- read.csv("data/processed/sample_metadata_processed.csv")
feature_data <- read.csv("data/processed/feature_data_processed.csv", row.names = 1)

cat(paste("Loaded expression matrix:", nrow(expr_matrix), "probes x", ncol(expr_matrix), "samples\n"))

# =============================================================================
# 1. DESIGN MATRIX SETUP
# =============================================================================

cat("\n=== Setting up Design Matrix ===\n")

# Ensure sample order matches
sample_metadata <- sample_metadata[match(colnames(expr_matrix), sample_metadata$sample_id), ]

# Filter out unknown samples for differential analysis
known_samples <- sample_metadata$condition %in% c("Normal", "Cancer")
expr_matrix_filtered <- expr_matrix[, known_samples]
sample_metadata_filtered <- sample_metadata[known_samples, ]

cat(paste("Samples for analysis:", ncol(expr_matrix_filtered), "\n"))
print(table(sample_metadata_filtered$condition))

# Create design matrix
condition <- factor(sample_metadata_filtered$condition, levels = c("Normal", "Cancer"))
design <- model.matrix(~ 0 + condition)
colnames(design) <- levels(condition)

cat("Design matrix created:\n")
print(head(design))

# =============================================================================
# 2. LINEAR MODELING WITH LIMMA
# =============================================================================

cat("\n=== Performing Linear Modeling ===\n")

# Fit linear model
fit <- lmFit(expr_matrix_filtered, design)

# Define contrast matrix (Cancer vs Normal)
contrast_matrix <- makeContrasts(
  Cancer_vs_Normal = Cancer - Normal,
  levels = design
)

cat("Contrast matrix:\n")
print(contrast_matrix)

# Fit contrasts
fit_contrasts <- contrasts.fit(fit, contrast_matrix)

# Empirical Bayes statistics
fit_eb <- eBayes(fit_contrasts)

cat("Linear modeling completed\n")

# =============================================================================
# 3. EXTRACT DIFFERENTIAL EXPRESSION RESULTS
# =============================================================================

cat("\n=== Extracting Differential Expression Results ===\n")

# Define significance thresholds
pval_threshold <- 0.05
logfc_threshold <- 1.0

# Extract results
results_all <- topTable(fit_eb, 
                       coef = "Cancer_vs_Normal", 
                       number = Inf, 
                       sort.by = "P")

# Add gene symbols and descriptions if available
if ("Gene.symbol" %in% colnames(feature_data)) {
  results_all$Gene_Symbol <- feature_data[rownames(results_all), "Gene.symbol"]
} else if ("GENE_SYMBOL" %in% colnames(feature_data)) {
  results_all$Gene_Symbol <- feature_data[rownames(results_all), "GENE_SYMBOL"]
} else {
  results_all$Gene_Symbol <- rownames(results_all)
}

# Add gene descriptions
if ("Gene.title" %in% colnames(feature_data)) {
  results_all$Gene_Description <- feature_data[rownames(results_all), "Gene.title"]
} else if ("GENE_TITLE" %in% colnames(feature_data)) {
  results_all$Gene_Description <- feature_data[rownames(results_all), "GENE_TITLE"]
}

# Add significance categories
results_all <- results_all %>%
  mutate(
    Significance = case_when(
      adj.P.Val < pval_threshold & logFC > logfc_threshold ~ "Upregulated",
      adj.P.Val < pval_threshold & logFC < -logfc_threshold ~ "Downregulated",
      TRUE ~ "Not Significant"
    ),
    abs_logFC = abs(logFC)
  )

# Summary of results
cat("Differential Expression Results Summary:\n")
summary_table <- table(results_all$Significance)
print(summary_table)

# Save all results
write.csv(results_all, "data/results/differential_expression_all_results.csv", row.names = TRUE)
cat("All results saved to: data/results/differential_expression_all_results.csv\n")

# Extract significant genes
significant_genes <- results_all[results_all$Significance != "Not Significant", ]
upregulated_genes <- results_all[results_all$Significance == "Upregulated", ]
downregulated_genes <- results_all[results_all$Significance == "Downregulated", ]

cat(paste("Total significant genes:", nrow(significant_genes), "\n"))
cat(paste("Upregulated genes:", nrow(upregulated_genes), "\n"))
cat(paste("Downregulated genes:", nrow(downregulated_genes), "\n"))

# Save significant results
write.csv(significant_genes, "data/results/significant_genes.csv", row.names = TRUE)
write.csv(upregulated_genes, "data/results/upregulated_genes.csv", row.names = TRUE)
write.csv(downregulated_genes, "data/results/downregulated_genes.csv", row.names = TRUE)

# =============================================================================
# 4. VOLCANO PLOT
# =============================================================================

cat("\n=== Creating Volcano Plot ===\n")

# Create volcano plot data
volcano_data <- results_all %>%
  mutate(
    neg_log10_pval = -log10(adj.P.Val),
    Color = case_when(
      Significance == "Upregulated" ~ "red",
      Significance == "Downregulated" ~ "blue",
      TRUE ~ "grey"
    )
  )

# Create volcano plot
pdf("figures/07_volcano_plot.pdf", width = 10, height = 8)
p_volcano <- ggplot(volcano_data, aes(x = logFC, y = neg_log10_pval)) +
  geom_point(aes(color = Color), alpha = 0.6, size = 1) +
  scale_color_identity() +
  geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), 
             linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(pval_threshold), 
             linetype = "dashed", color = "black", alpha = 0.5) +
  labs(
    title = "Volcano Plot: Cancer vs Normal",
    x = "Log2 Fold Change",
    y = "-Log10(Adjusted P-value)",
    subtitle = paste("Significance thresholds: |log2FC| >", logfc_threshold, 
                    ", FDR <", pval_threshold)
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  )

# Add gene labels for top significant genes
top_genes <- significant_genes %>%
  arrange(adj.P.Val) %>%
  head(20)

if (nrow(top_genes) > 0) {
  p_volcano <- p_volcano +
    geom_text_repel(
      data = volcano_data[rownames(top_genes), ],
      aes(label = Gene_Symbol),
      size = 3,
      max.overlaps = 10,
      box.padding = 0.3
    )
}

print(p_volcano)
dev.off()

# =============================================================================
# 5. MA PLOT
# =============================================================================

cat("Creating MA Plot...\n")

# Create MA plot
pdf("figures/08_ma_plot.pdf", width = 10, height = 8)
ma_data <- results_all %>%
  mutate(
    A = (rowMeans(expr_matrix_filtered[rownames(results_all), 
                                     sample_metadata_filtered$condition == "Normal"]) +
         rowMeans(expr_matrix_filtered[rownames(results_all), 
                                     sample_metadata_filtered$condition == "Cancer"])) / 2,
    M = logFC
  )

p_ma <- ggplot(ma_data, aes(x = A, y = M)) +
  geom_point(aes(color = Significance), alpha = 0.6, size = 1) +
  scale_color_manual(values = c("Upregulated" = "red", 
                               "Downregulated" = "blue", 
                               "Not Significant" = "grey")) +
  geom_hline(yintercept = c(-logfc_threshold, logfc_threshold), 
             linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = 0, color = "black", alpha = 0.5) +
  labs(
    title = "MA Plot: Cancer vs Normal",
    x = "Average Log2 Expression (A)",
    y = "Log2 Fold Change (M)",
    color = "Gene Category"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "bottom"
  )

print(p_ma)
dev.off()

# =============================================================================
# 6. HEATMAP OF TOP DIFFERENTIALLY EXPRESSED GENES
# =============================================================================

cat("Creating Heatmap of Top DEGs...\n")

# Select top 50 most significant genes for heatmap
top_deg_genes <- significant_genes %>%
  arrange(adj.P.Val) %>%
  head(50)

if (nrow(top_deg_genes) > 0) {
  # Extract expression data for top genes
  heatmap_data <- expr_matrix_filtered[rownames(top_deg_genes), ]
  
  # Prepare annotation for samples
  sample_annotation <- data.frame(
    Condition = sample_metadata_filtered$condition,
    row.names = sample_metadata_filtered$sample_id
  )
  
  # Prepare annotation for genes
  gene_annotation <- data.frame(
    Regulation = top_deg_genes$Significance,
    row.names = rownames(top_deg_genes)
  )
  
  # Define colors
  annotation_colors <- list(
    Condition = c("Normal" = "#2E8B57", "Cancer" = "#DC143C"),
    Regulation = c("Upregulated" = "#FF6B6B", "Downregulated" = "#4ECDC4")
  )
  
  # Create heatmap
  pdf("figures/09_top_genes_heatmap.pdf", width = 12, height = 10)
  pheatmap(heatmap_data,
           scale = "row",
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "ward.D2",
           annotation_col = sample_annotation,
           annotation_row = gene_annotation,
           annotation_colors = annotation_colors,
           show_rownames = TRUE,
           show_colnames = FALSE,
           fontsize_row = 8,
           main = "Top 50 Differentially Expressed Genes",
           color = colorRampPalette(c("blue", "white", "red"))(50))
  dev.off()
}

# =============================================================================
# 7. EXPRESSION PROFILES OF TOP GENES
# =============================================================================

cat("Creating Expression Profiles...\n")

# Select top 12 genes for individual plots
top_12_genes <- significant_genes %>%
  arrange(adj.P.Val) %>%
  head(12)

if (nrow(top_12_genes) > 0) {
  pdf("figures/10_top_genes_expression.pdf", width = 15, height = 10)
  
  par(mfrow = c(3, 4), mar = c(4, 4, 3, 2))
  
  for (i in 1:nrow(top_12_genes)) {
    gene_id <- rownames(top_12_genes)[i]
    gene_symbol <- top_12_genes$Gene_Symbol[i]
    
    # Extract expression values
    gene_expr <- as.numeric(expr_matrix_filtered[gene_id, ])
    conditions <- sample_metadata_filtered$condition
    
    # Create boxplot
    boxplot(gene_expr ~ conditions,
            col = c("#2E8B57", "#DC143C"),
            main = paste(gene_symbol, "\n(adj.P =", 
                        format(top_12_genes$adj.P.Val[i], digits = 3), ")"),
            xlab = "Condition",
            ylab = "Log2 Expression",
            cex.main = 0.9)
    
    # Add points
    stripchart(gene_expr ~ conditions,
               vertical = TRUE,
               method = "jitter",
               add = TRUE,
               pch = 19,
               cex = 0.5,
               col = "darkgrey")
  }
  
  dev.off()
}

# =============================================================================
# 8. SAVE ANALYSIS RESULTS
# =============================================================================

cat("\n=== Saving Analysis Results ===\n")

# Save the limma fit object
save(fit_eb, results_all, significant_genes, 
     upregulated_genes, downregulated_genes,
     design, contrast_matrix,
     file = "data/results/limma_analysis_results.RData")
cat("Limma analysis results saved to: data/results/limma_analysis_results.RData\n")

# Create gene lists for pathway analysis
if (nrow(upregulated_genes) > 0) {
  upreg_gene_list <- upregulated_genes$Gene_Symbol
  upreg_gene_list <- upreg_gene_list[!is.na(upreg_gene_list) & upreg_gene_list != ""]
  writeLines(upreg_gene_list, "data/results/upregulated_gene_list.txt")
  cat("Upregulated gene list saved to: data/results/upregulated_gene_list.txt\n")
}

if (nrow(downregulated_genes) > 0) {
  downreg_gene_list <- downregulated_genes$Gene_Symbol
  downreg_gene_list <- downreg_gene_list[!is.na(downreg_gene_list) & downreg_gene_list != ""]
  writeLines(downreg_gene_list, "data/results/downregulated_gene_list.txt")
  cat("Downregulated gene list saved to: data/results/downregulated_gene_list.txt\n")
}

# All significant genes
if (nrow(significant_genes) > 0) {
  all_sig_gene_list <- significant_genes$Gene_Symbol
  all_sig_gene_list <- all_sig_gene_list[!is.na(all_sig_gene_list) & all_sig_gene_list != ""]
  writeLines(all_sig_gene_list, "data/results/all_significant_genes_list.txt")
  cat("All significant genes list saved to: data/results/all_significant_genes_list.txt\n")
}

# =============================================================================
# 9. GENERATE DIFFERENTIAL EXPRESSION REPORT
# =============================================================================

cat("\n=== Generating Differential Expression Report ===\n")

de_report <- c(
  "DIFFERENTIAL EXPRESSION ANALYSIS REPORT",
  "======================================",
  "",
  paste("Analysis Date:", Sys.Date()),
  paste("Dataset: GSE43346"),
  paste("Analysis Method: limma (linear models for microarray data)"),
  "",
  "ANALYSIS PARAMETERS:",
  paste("- Significance threshold (FDR):", pval_threshold),
  paste("- Fold change threshold (|log2FC|):", logfc_threshold),
  paste("- Samples analyzed:", ncol(expr_matrix_filtered)),
  paste("- Probes analyzed:", nrow(expr_matrix_filtered)),
  "",
  "SAMPLE DISTRIBUTION:",
  paste("- Normal samples:", sum(sample_metadata_filtered$condition == "Normal")),
  paste("- Cancer samples:", sum(sample_metadata_filtered$condition == "Cancer")),
  "",
  "DIFFERENTIAL EXPRESSION RESULTS:",
  paste("- Total significant genes:", nrow(significant_genes)),
  paste("- Upregulated genes (Cancer vs Normal):", nrow(upregulated_genes)),
  paste("- Downregulated genes (Cancer vs Normal):", nrow(downregulated_genes)),
  paste("- Not significant genes:", sum(results_all$Significance == "Not Significant")),
  "",
  "TOP 10 UPREGULATED GENES:",
  if (nrow(upregulated_genes) > 0) {
    top_up <- upregulated_genes %>% arrange(adj.P.Val) %>% head(10)
    paste("-", top_up$Gene_Symbol, "(log2FC =", round(top_up$logFC, 2), 
          ", FDR =", format(top_up$adj.P.Val, digits = 3), ")")
  } else {
    "None found"
  },
  "",
  "TOP 10 DOWNREGULATED GENES:",
  if (nrow(downregulated_genes) > 0) {
    top_down <- downregulated_genes %>% arrange(adj.P.Val) %>% head(10)
    paste("-", top_down$Gene_Symbol, "(log2FC =", round(top_down$logFC, 2), 
          ", FDR =", format(top_down$adj.P.Val, digits = 3), ")")
  } else {
    "None found"
  },
  "",
  "OUTPUT FILES GENERATED:",
  "- figures/07_volcano_plot.pdf",
  "- figures/08_ma_plot.pdf",
  "- figures/09_top_genes_heatmap.pdf",
  "- figures/10_top_genes_expression.pdf",
  "",
  "DATA FILES GENERATED:",
  "- data/results/differential_expression_all_results.csv",
  "- data/results/significant_genes.csv",
  "- data/results/upregulated_genes.csv",
  "- data/results/downregulated_genes.csv",
  "- data/results/limma_analysis_results.RData",
  "- data/results/*_gene_list.txt (for pathway analysis)"
)

writeLines(de_report, "reports/differential_expression_report.txt")
cat("Differential expression report saved to: reports/differential_expression_report.txt\n")

cat("\n=== Differential Expression Analysis Completed Successfully ===\n")
cat("Next step: Run 04_functional_enrichment.R for pathway analysis\n")

# Clean up
rm(expr_matrix, sample_metadata, feature_data, fit, fit_contrasts, fit_eb)
gc()