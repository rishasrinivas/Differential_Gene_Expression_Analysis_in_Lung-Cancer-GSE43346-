#!/usr/bin/env Rscript

# =============================================================================
# Script: 02_quality_control.R
# Purpose: Quality control and preprocessing of GSE43346 dataset
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
  library(corrplot)
  library(VennDiagram)
  library(gridExtra)
})

cat("=== Starting Quality Control Analysis ===\n")

# Load the data
if (!file.exists("data/raw/expression_matrix.csv")) {
  stop("Expression matrix not found. Please run 01_data_download.R first.")
}

# Load expression matrix
expr_matrix <- read.csv("data/raw/expression_matrix.csv", row.names = 1)
sample_metadata <- read.csv("data/raw/sample_metadata.csv")
feature_data <- read.csv("data/raw/feature_data.csv", row.names = 1)

cat(paste("Loaded expression matrix:", nrow(expr_matrix), "probes x", ncol(expr_matrix), "samples\n"))

# =============================================================================
# 1. DATA PREPROCESSING
# =============================================================================

cat("\n=== Data Preprocessing ===\n")

# Ensure sample order matches between expression matrix and metadata
sample_metadata <- sample_metadata[match(colnames(expr_matrix), sample_metadata$sample_id), ]

# Check for missing values
missing_values <- sum(is.na(expr_matrix))
cat(paste("Missing values in expression matrix:", missing_values, "\n"))

if (missing_values > 0) {
  cat("Handling missing values by removing probes with >10% missing values...\n")
  # Remove probes with more than 10% missing values
  missing_percent <- rowSums(is.na(expr_matrix)) / ncol(expr_matrix)
  expr_matrix <- expr_matrix[missing_percent <= 0.1, ]
  
  # Impute remaining missing values with row median
  for (i in 1:nrow(expr_matrix)) {
    expr_matrix[i, is.na(expr_matrix[i, ])] <- median(as.numeric(expr_matrix[i, ]), na.rm = TRUE)
  }
  
  cat(paste("After filtering: ", nrow(expr_matrix), "probes retained\n"))
}

# Log2 transformation (if data is not already log-transformed)
# Check if data needs log transformation
if (max(expr_matrix, na.rm = TRUE) > 100) {
  cat("Applying log2 transformation...\n")
  expr_matrix <- log2(expr_matrix + 1)
} else {
  cat("Data appears to be already log-transformed\n")
}

# =============================================================================
# 2. QUALITY CONTROL VISUALIZATIONS
# =============================================================================

cat("\n=== Generating Quality Control Plots ===\n")

# Create a color palette for conditions
condition_colors <- c("Normal" = "#2E8B57", "Cancer" = "#DC143C", "Unknown" = "#808080")
sample_colors <- condition_colors[sample_metadata$condition]

# 2.1 Sample distribution (boxplots)
pdf("figures/01_sample_distribution.pdf", width = 12, height = 8)
par(mar = c(8, 4, 4, 2))
boxplot(expr_matrix, 
        las = 2, 
        col = sample_colors,
        main = "Sample Expression Distribution",
        ylab = "Log2 Expression",
        cex.axis = 0.7)
legend("topright", 
       legend = names(condition_colors), 
       fill = condition_colors, 
       cex = 0.8)
dev.off()

# 2.2 Density plots
pdf("figures/02_density_plots.pdf", width = 10, height = 6)
plot(density(expr_matrix[, 1]), 
     main = "Sample Expression Density", 
     xlab = "Log2 Expression", 
     ylim = c(0, max(density(expr_matrix[, 1])$y) * 1.2))

for (i in 1:ncol(expr_matrix)) {
  lines(density(expr_matrix[, i]), col = sample_colors[i])
}

legend("topright", 
       legend = names(condition_colors), 
       col = condition_colors, 
       lty = 1, 
       cex = 0.8)
dev.off()

# 2.3 Principal Component Analysis (PCA)
cat("Performing PCA analysis...\n")
pca_result <- prcomp(t(expr_matrix), scale. = TRUE)
pca_data <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  PC3 = pca_result$x[, 3],
  Sample = colnames(expr_matrix),
  Condition = sample_metadata$condition
)

# Calculate variance explained
var_explained <- round((pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100, 1)

# PCA plot
pdf("figures/03_pca_analysis.pdf", width = 12, height = 5)
p1 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = condition_colors) +
  labs(
    title = "Principal Component Analysis",
    x = paste0("PC1 (", var_explained[1], "% variance)"),
    y = paste0("PC2 (", var_explained[2], "% variance)")
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- ggplot(pca_data, aes(x = PC1, y = PC3, color = Condition)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = condition_colors) +
  labs(
    title = "PC1 vs PC3",
    x = paste0("PC1 (", var_explained[1], "% variance)"),
    y = paste0("PC3 (", var_explained[3], "% variance)")
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(p1, p2, ncol = 2)
dev.off()

# 2.4 Sample correlation heatmap
cat("Computing sample correlations...\n")
sample_cor <- cor(expr_matrix, method = "spearman")

pdf("figures/04_sample_correlation.pdf", width = 10, height = 10)
pheatmap(sample_cor,
         annotation_col = data.frame(Condition = sample_metadata$condition,
                                    row.names = sample_metadata$sample_id),
         annotation_colors = list(Condition = condition_colors),
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Sample Correlation Heatmap",
         color = colorRampPalette(c("blue", "white", "red"))(50))
dev.off()

# 2.5 Hierarchical clustering
pdf("figures/05_hierarchical_clustering.pdf", width = 12, height = 8)
sample_dist <- dist(t(expr_matrix))
sample_hc <- hclust(sample_dist, method = "ward.D2")
plot(sample_hc, 
     labels = paste(sample_metadata$condition, 1:nrow(sample_metadata), sep = "_"),
     main = "Hierarchical Clustering of Samples",
     xlab = "",
     cex = 0.6)
dev.off()

# =============================================================================
# 3. OUTLIER DETECTION
# =============================================================================

cat("\n=== Outlier Detection ===\n")

# Detect outliers using PCA
# Calculate distances from center
center_distances <- sqrt(pca_data$PC1^2 + pca_data$PC2^2)
outlier_threshold <- quantile(center_distances, 0.95)
outliers <- which(center_distances > outlier_threshold)

if (length(outliers) > 0) {
  cat(paste("Potential outliers detected:", length(outliers), "samples\n"))
  cat("Outlier samples:\n")
  print(pca_data[outliers, c("Sample", "Condition")])
  
  # Create outlier visualization
  pdf("figures/06_outlier_detection.pdf", width = 8, height = 6)
  plot(pca_data$PC1, pca_data$PC2, 
       col = sample_colors, 
       pch = 19,
       main = "Outlier Detection in PCA Space",
       xlab = paste0("PC1 (", var_explained[1], "% variance)"),
       ylab = paste0("PC2 (", var_explained[2], "% variance)"))
  
  # Highlight outliers
  points(pca_data$PC1[outliers], pca_data$PC2[outliers], 
         col = "black", pch = 1, cex = 2, lwd = 2)
  
  legend("topright", 
         legend = c(names(condition_colors), "Outliers"), 
         col = c(condition_colors, "black"), 
         pch = c(19, 19, 19, 1),
         cex = 0.8)
  dev.off()
} else {
  cat("No significant outliers detected\n")
}

# =============================================================================
# 4. BATCH EFFECT ASSESSMENT
# =============================================================================

cat("\n=== Batch Effect Assessment ===\n")

# If batch information is available in metadata
if ("batch" %in% colnames(sample_metadata) || any(grepl("batch", colnames(sample_metadata), ignore.case = TRUE))) {
  cat("Batch information found. Analyzing potential batch effects...\n")
  # Add batch effect analysis here if needed
} else {
  cat("No explicit batch information found in metadata\n")
}

# =============================================================================
# 5. PROBE FILTERING AND ANNOTATION
# =============================================================================

cat("\n=== Probe Filtering and Annotation ===\n")

# Remove probes with low variance (bottom 10%)
probe_var <- apply(expr_matrix, 1, var, na.rm = TRUE)
var_threshold <- quantile(probe_var, 0.1, na.rm = TRUE)
high_var_probes <- probe_var > var_threshold

cat(paste("Removing", sum(!high_var_probes), "low-variance probes\n"))
expr_matrix_filtered <- expr_matrix[high_var_probes, ]
feature_data_filtered <- feature_data[high_var_probes, ]

cat(paste("After filtering:", nrow(expr_matrix_filtered), "probes retained\n"))

# =============================================================================
# 6. SAVE PROCESSED DATA
# =============================================================================

cat("\n=== Saving Processed Data ===\n")

# Save processed expression matrix
write.csv(expr_matrix_filtered, "data/processed/expression_matrix_processed.csv", row.names = TRUE)
cat("Processed expression matrix saved to: data/processed/expression_matrix_processed.csv\n")

# Save updated sample metadata
write.csv(sample_metadata, "data/processed/sample_metadata_processed.csv", row.names = FALSE)
cat("Sample metadata saved to: data/processed/sample_metadata_processed.csv\n")

# Save filtered feature data
write.csv(feature_data_filtered, "data/processed/feature_data_processed.csv", row.names = TRUE)
cat("Feature data saved to: data/processed/feature_data_processed.csv\n")

# Save PCA results
save(pca_result, pca_data, var_explained, file = "data/processed/pca_results.RData")
cat("PCA results saved to: data/processed/pca_results.RData\n")

# =============================================================================
# 7. GENERATE QC REPORT
# =============================================================================

cat("\n=== Generating Quality Control Report ===\n")

qc_summary <- list(
  total_samples_initial = ncol(expr_matrix),
  total_probes_initial = nrow(expr_matrix),
  total_samples_final = ncol(expr_matrix_filtered),
  total_probes_final = nrow(expr_matrix_filtered),
  normal_samples = sum(sample_metadata$condition == "Normal"),
  cancer_samples = sum(sample_metadata$condition == "Cancer"),
  unknown_samples = sum(sample_metadata$condition == "Unknown"),
  outliers_detected = length(outliers),
  variance_explained_pc1 = var_explained[1],
  variance_explained_pc2 = var_explained[2],
  data_transformation = ifelse(max(expr_matrix, na.rm = TRUE) > 100, "log2", "none")
)

# Create QC report
qc_report_text <- c(
  "QUALITY CONTROL SUMMARY REPORT",
  "================================",
  "",
  paste("Date:", Sys.Date()),
  paste("Dataset: GSE43346"),
  "",
  "SAMPLE INFORMATION:",
  paste("- Total samples:", qc_summary$total_samples_final),
  paste("- Normal samples:", qc_summary$normal_samples),
  paste("- Cancer samples:", qc_summary$cancer_samples),
  paste("- Unknown samples:", qc_summary$unknown_samples),
  "",
  "PROBE INFORMATION:",
  paste("- Initial probes:", qc_summary$total_probes_initial),
  paste("- Final probes (after filtering):", qc_summary$total_probes_final),
  paste("- Probes removed:", qc_summary$total_probes_initial - qc_summary$total_probes_final),
  "",
  "QUALITY METRICS:",
  paste("- Outliers detected:", qc_summary$outliers_detected),
  paste("- PC1 variance explained:", paste0(qc_summary$variance_explained_pc1, "%")),
  paste("- PC2 variance explained:", paste0(qc_summary$variance_explained_pc2, "%")),
  paste("- Data transformation:", qc_summary$data_transformation),
  "",
  "FILES GENERATED:",
  "- figures/01_sample_distribution.pdf",
  "- figures/02_density_plots.pdf", 
  "- figures/03_pca_analysis.pdf",
  "- figures/04_sample_correlation.pdf",
  "- figures/05_hierarchical_clustering.pdf",
  if(length(outliers) > 0) "- figures/06_outlier_detection.pdf" else NULL,
  "",
  "PROCESSED DATA FILES:",
  "- data/processed/expression_matrix_processed.csv",
  "- data/processed/sample_metadata_processed.csv",
  "- data/processed/feature_data_processed.csv",
  "- data/processed/pca_results.RData"
)

writeLines(qc_report_text, "reports/quality_control_report.txt")
cat("Quality control report saved to: reports/quality_control_report.txt\n")

cat("\n=== Quality Control Analysis Completed Successfully ===\n")
cat("Next step: Run 03_differential_expression.R for differential expression analysis\n")

# Clean up
rm(expr_matrix, sample_metadata, feature_data, expr_matrix_filtered, 
   feature_data_filtered, pca_result, pca_data, sample_cor, qc_summary)
gc()