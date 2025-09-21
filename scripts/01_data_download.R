#!/usr/bin/env Rscript

# =============================================================================
# Script: 01_data_download.R
# Purpose: Download and initial processing of GSE43346 dataset from GEO
# Author: Bioinformatics Analysis Project
# Date: 2025
# =============================================================================

# Load required libraries
suppressMessages({
  library(GEOquery)
  library(Biobase)
  library(dplyr)
  library(readr)
})

# Print session information for reproducibility
cat("=== R Session Information ===\n")
sessionInfo()

# Create directory structure
dir.create("data", showWarnings = FALSE)
dir.create("data/raw", showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
dir.create("data/results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE)
dir.create("reports", showWarnings = FALSE)

cat("\n=== Starting GSE43346 Data Download ===\n")

# Download GSE43346 dataset
gse_id <- "GSE43346"
cat(paste("Downloading", gse_id, "from GEO database...\n"))

tryCatch({
  # Download the GEO dataset
  gset <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE)
  
  # Extract expression set
  if (length(gset) > 1) {
    # If multiple platforms, select the first one
    eset <- gset[[1]]
    cat("Multiple platforms detected. Using the first platform.\n")
  } else {
    eset <- gset[[1]]
  }
  
  cat(paste("Successfully downloaded", gse_id, "\n"))
  
  # Display basic information about the dataset
  cat("\n=== Dataset Information ===\n")
  cat(paste("Number of samples:", ncol(eset), "\n"))
  cat(paste("Number of probes:", nrow(eset), "\n"))
  cat(paste("Platform:", annotation(eset), "\n"))
  
  # Extract expression matrix
  expr_matrix <- exprs(eset)
  
  # Extract phenotype data
  pheno_data <- pData(eset)
  
  # Extract feature data (gene annotations)
  feature_data <- fData(eset)
  
  # Display sample information
  cat("\n=== Sample Information ===\n")
  print(head(pheno_data[, c("title", "source_name_ch1", "characteristics_ch1")]))
  
  # Create sample metadata
  sample_metadata <- pheno_data %>%
    select(geo_accession, title, source_name_ch1) %>%
    mutate(
      sample_id = geo_accession,
      sample_name = title,
      tissue_type = source_name_ch1
    )
  
  # Classify samples based on tissue type
  # Based on GSE43346: normal lung tissue vs SCLC samples
  sample_metadata <- sample_metadata %>%
    mutate(
      condition = case_when(
        grepl("normal|Normal|NORMAL", tissue_type, ignore.case = TRUE) ~ "Normal",
        grepl("cancer|tumor|Cancer|Tumor|SCLC|sclc", tissue_type, ignore.case = TRUE) ~ "Cancer",
        grepl("lung", tissue_type, ignore.case = TRUE) & !grepl("cancer|tumor", tissue_type, ignore.case = TRUE) ~ "Normal",
        TRUE ~ "Unknown"
      )
    )
  
  # Display condition summary
  cat("\n=== Sample Condition Summary ===\n")
  print(table(sample_metadata$condition))
  
  # Quality check: ensure we have both conditions
  if (!"Normal" %in% sample_metadata$condition || !"Cancer" %in% sample_metadata$condition) {
    cat("WARNING: Could not properly classify samples into Normal and Cancer groups.\n")
    cat("Please check the sample metadata manually.\n")
    
    # Display unique tissue types for manual inspection
    cat("\nUnique tissue types found:\n")
    print(unique(sample_metadata$tissue_type))
  }
  
  # Save raw data
  cat("\n=== Saving Raw Data ===\n")
  
  # Save expression matrix
  write.csv(expr_matrix, "data/raw/expression_matrix.csv", row.names = TRUE)
  cat("Expression matrix saved to: data/raw/expression_matrix.csv\n")
  
  # Save sample metadata
  write.csv(sample_metadata, "data/raw/sample_metadata.csv", row.names = FALSE)
  cat("Sample metadata saved to: data/raw/sample_metadata.csv\n")
  
  # Save feature data (gene annotations)
  write.csv(feature_data, "data/raw/feature_data.csv", row.names = TRUE)
  cat("Feature data saved to: data/raw/feature_data.csv\n")
  
  # Save the entire ExpressionSet object for later use
  saveRDS(eset, "data/raw/gse43346_expressionset.rds")
  cat("ExpressionSet object saved to: data/raw/gse43346_expressionset.rds\n")
  
  # Create a summary report
  cat("\n=== Creating Summary Report ===\n")
  
  summary_report <- list(
    dataset_id = gse_id,
    download_date = Sys.Date(),
    platform = annotation(eset),
    total_samples = ncol(eset),
    total_probes = nrow(eset),
    normal_samples = sum(sample_metadata$condition == "Normal"),
    cancer_samples = sum(sample_metadata$condition == "Cancer"),
    unknown_samples = sum(sample_metadata$condition == "Unknown")
  )
  
  # Save summary report
  writeLines(
    c(
      paste("Dataset:", summary_report$dataset_id),
      paste("Download Date:", summary_report$download_date),
      paste("Platform:", summary_report$platform),
      paste("Total Samples:", summary_report$total_samples),
      paste("Total Probes:", summary_report$total_probes),
      paste("Normal Samples:", summary_report$normal_samples),
      paste("Cancer Samples:", summary_report$cancer_samples),
      paste("Unknown Samples:", summary_report$unknown_samples)
    ),
    "data/raw/download_summary.txt"
  )
  
  cat("Download summary saved to: data/raw/download_summary.txt\n")
  
}, error = function(e) {
  cat("ERROR: Failed to download or process the dataset.\n")
  cat("Error message:", e$message, "\n")
  cat("\nPlease check your internet connection and try again.\n")
  cat("You may also need to install required Bioconductor packages:\n")
  cat("BiocManager::install(c('GEOquery', 'Biobase'))\n")
  
  quit(status = 1)
})

cat("\n=== Data Download Completed Successfully ===\n")
cat("Next step: Run 02_quality_control.R for data preprocessing and quality assessment\n")

# Clean up
rm(gset, eset, expr_matrix, pheno_data, feature_data, sample_metadata, summary_report)
gc()
