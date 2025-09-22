# Differential Gene Expression Analysis in Lung Cancer (GSE43346)

## Project Overview

This project investigates differential gene expression in small cell lung cancer (SCLC) using the publicly available microarray dataset GSE43346 from Gene Expression Omnibus (GEO). The study applies statistical methods to identify significantly upregulated and downregulated genes between SCLC and normal lung tissue samples, followed by functional enrichment analysis to reveal key biological pathways involved in tumor progression.

## Dataset Information

- **GEO Accession**: GSE43346
- **Platform**: Affymetrix Human Gene Expression Array
- **Samples**: 
  - 43 normal lung tissue samples
  - 23 small cell lung cancer (SCLC) samples
- **Total**: 66 samples

## Research Objectives

1. **Primary Analysis**:
   - Identify significantly differentially expressed genes (DEGs) between SCLC and normal tissue
   - Apply statistical significance thresholds (adjusted p-value < 0.05, |log2FC| > 1)

2. **Functional Analysis**:
   - Perform Gene Ontology (GO) enrichment analysis
   - Conduct KEGG pathway analysis
   - Identify key biological processes and molecular functions

3. **Visualization**:
   - Generate volcano plots for DEG visualization
   - Create heatmaps for expression patterns
   - Produce pathway network visualizations

4. **Biomarker Discovery**:
   - Identify potential therapeutic targets
   - Discover candidate biomarkers for SCLC diagnosis and prognosis

## Repository Structure

```
lung-cancer-deg-analysis/
├── README.md                          # This file
├── requirements.txt                   # Python dependencies
├── data/
│   ├── raw/                          # Raw data (not included - download from GEO)
│   ├── processed/                    # Processed data files
│   └── results/                      # Analysis results
├── scripts/
│   ├── 01_data_download.R            # GEO data download
│   ├── 02_quality_control.R          # QC and preprocessing
│   ├── 03_differential_expression.R  # DESeq2 analysis
│   ├── 04_functional_enrichment.R    # GO/KEGG analysis
│   ├── 05_visualization.py           # Python visualizations
│   └── 06_network_analysis.py        # Network analysis
├── figures/                          # Generated plots and figures
├── reports/                          # Analysis reports
```

## Key Features

### Statistical Analysis
- **Differential Expression**: limma package for microarray analysis
- **Multiple Testing Correction**: Benjamini-Hochberg FDR correction
- **Effect Size Filtering**: Log2 fold change threshold (|log2FC| > 1)
- **Quality Control**: PCA, sample clustering, and batch effect assessment

### Functional Enrichment
- **Gene Ontology (GO)**: Biological processes, molecular functions, cellular components
- **KEGG Pathways**: Disease-relevant pathway analysis
- **Enrichment Tools**: DAVID, Enrichr, g:Profiler integration

### Visualization Tools
- **Volcano Plots**: DEG visualization with significance thresholds
- **Heatmaps**: Expression pattern clustering
- **PCA Plots**: Sample grouping and outlier detection
- **Pathway Networks**: Cytoscape-compatible network files

## Getting Started

### Prerequisites

- R (>= 4.0.0)
- Python (>= 3.8)
- Required R packages: GEOquery, limma, ggplot2, dplyr, pheatmap
- Required Python packages: pandas, numpy, matplotlib, seaborn, networkx

### Installation

1. **Clone the repository**:
```bash
git clone https://github.com/yourusername/lung-cancer-deg-analysis.git
cd lung-cancer-deg-analysis
```

2. **Set up R environment**:
```r
install.packages(c("BiocManager", "devtools"))
BiocManager::install(c("GEOquery", "limma", "org.Hs.eg.db", "clusterProfiler"))
install.packages(c("ggplot2", "dplyr", "pheatmap", "VennDiagram"))
```

3. **Set up Python environment**:
```bash
pip install -r requirements.txt
```

### Usage

1. **Download and preprocess data**:
```bash
Rscript scripts/01_data_download.R
Rscript scripts/02_quality_control.R
```

2. **Perform differential expression analysis**:
```bash
Rscript scripts/03_differential_expression.R
```

3. **Functional enrichment analysis**:
```bash
Rscript scripts/04_functional_enrichment.R
```

4. **Generate visualizations**:
```bash
python scripts/05_visualization.py
python scripts/06_network_analysis.py
```

## Expected Results

### Statistical Outcomes
- **Total significant genes**: 9901
- **Upregulated genes (Cancer vs Normal)**: 4133
- **Downregulated genes (Cancer vs Normal)**: 5768
- **Not significant genes**: 39306

### Key Pathways Expected
KEY CANCER-RELATED PATHWAYS IDENTIFIED:
- Cell cycle
- MicroRNAs in cancer
- Prostate cancer
- p53 signaling pathway
- Bladder cancer
- Pathways in cancer
- Colorectal cancer
- Proteoglycans in cancer
- Small cell lung cancer
- Transcriptional misregulation in cancer
- 
NETWORK CONSTRUCTION:
- **Input genes**: 200
- **Correlation threshold**: 0.7
- **Network nodes**: 200
- **Network edges**: 17707
- **Network density**: 0.8898

## Industry Applications

This analysis workflow is directly applicable to:

- **Pharmaceutical Industry**: Drug target identification and validation
- **Precision Medicine**: Personalized treatment strategy development  
- **Clinical Diagnostics**: Biomarker discovery for early detection
- **Biotechnology**: Therapeutic development and companion diagnostics

## Scientific Impact

The project contributes to:
- Understanding molecular mechanisms of SCLC progression
- Identification of novel therapeutic targets
- Development of prognostic biomarkers
- Advancement of personalized cancer treatment approaches


## Contact

- **Author**: Risha Srinivas
  
## Acknowledgments

- Gene Expression Omnibus (GEO) for providing public datasets
- The original data contributors of GSE43346
- Bioconductor and R community for statistical tools
- Python scientific computing community

---

*This project is part of a comprehensive bioinformatics workflow for cancer transcriptomics analysis and serves as a template for similar studies in other cancer types.*
