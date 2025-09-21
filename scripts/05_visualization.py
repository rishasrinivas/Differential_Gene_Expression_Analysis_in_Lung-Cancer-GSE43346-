#!/usr/bin/env python3

"""
Script: 05_visualization.py
Purpose: Advanced visualizations for lung cancer differential expression analysis
Author: Bioinformatics Analysis Project
Date: 2025
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import warnings
warnings.filterwarnings('ignore')

# Set style parameters
plt.style.use('default')
sns.set_palette("husl")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 12

print("=== Starting Advanced Visualization Analysis ===")

# =============================================================================
# 1. LOAD DATA
# =============================================================================

print("\n=== Loading Data ===")

try:
    # Load differential expression results
    de_results = pd.read_csv('data/results/differential_expression_all_results.csv', index_col=0)
    significant_genes = pd.read_csv('data/results/significant_genes.csv', index_col=0)
    upregulated_genes = pd.read_csv('data/results/upregulated_genes.csv', index_col=0)
    downregulated_genes = pd.read_csv('data/results/downregulated_genes.csv', index_col=0)
    
    # Load expression data
    expr_data = pd.read_csv('data/processed/expression_matrix_processed.csv', index_col=0)
    sample_metadata = pd.read_csv('data/processed/sample_metadata_processed.csv')
    
    print(f"Loaded data:")
    print(f"- DE results: {de_results.shape[0]} genes")
    print(f"- Significant genes: {significant_genes.shape[0]} genes")
    print(f"- Expression matrix: {expr_data.shape[0]} genes x {expr_data.shape[1]} samples")
    print(f"- Sample metadata: {sample_metadata.shape[0]} samples")
    
except FileNotFoundError as e:
    print(f"Error loading data: {e}")
    print("Please run the R scripts first (01-04) to generate the required data files.")
    exit(1)

# =============================================================================
# 2. ENHANCED VOLCANO PLOT
# =============================================================================

print("\n=== Creating Enhanced Volcano Plot ===")

def create_enhanced_volcano_plot():
    """Create an enhanced volcano plot with annotations and statistics."""
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
    
    # Prepare data
    volcano_data = de_results.copy()
    volcano_data['neg_log10_pval'] = -np.log10(volcano_data['adj.P.Val'].clip(lower=1e-300))
    volcano_data['abs_logFC'] = abs(volcano_data['logFC'])
    
    # Define significance thresholds
    pval_threshold = 0.05
    logfc_threshold = 1.0
    
    # Color mapping
    color_map = {'Normal': '#2E8B57', 'Cancer': '#DC143C', 'Unknown': '#808080'}
    colors = [color_map.get(condition, '#808080') for condition in conditions]
    
    # PCA plots
    # PC1 vs PC2
    axes[0, 0].scatter(pca_result[:, 0], pca_result[:, 1], c=colors, alpha=0.7, s=60)
    axes[0, 0].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
    axes[0, 0].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
    axes[0, 0].set_title('PCA: PC1 vs PC2')
    axes[0, 0].grid(True, alpha=0.3)
    
    # PC2 vs PC3
    axes[0, 1].scatter(pca_result[:, 1], pca_result[:, 2], c=colors, alpha=0.7, s=60)
    axes[0, 1].set_xlabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
    axes[0, 1].set_ylabel(f'PC3 ({pca.explained_variance_ratio_[2]:.1%} variance)')
    axes[0, 1].set_title('PCA: PC2 vs PC3')
    axes[0, 1].grid(True, alpha=0.3)
    
    # Variance explained plot
    cumsum_var = np.cumsum(pca.explained_variance_ratio_)
    axes[0, 2].plot(range(1, 11), pca.explained_variance_ratio_[:10], 'bo-', linewidth=2, markersize=8)
    axes[0, 2].plot(range(1, 11), cumsum_var[:10], 'ro-', linewidth=2, markersize=8)
    axes[0, 2].set_xlabel('Principal Component')
    axes[0, 2].set_ylabel('Variance Explained')
    axes[0, 2].set_title('PCA Variance Explained')
    axes[0, 2].legend(['Individual', 'Cumulative'])
    axes[0, 2].grid(True, alpha=0.3)
    
    # t-SNE plot
    axes[1, 0].scatter(tsne_result[:, 0], tsne_result[:, 1], c=colors, alpha=0.7, s=60)
    axes[1, 0].set_xlabel('t-SNE 1')
    axes[1, 0].set_ylabel('t-SNE 2')
    axes[1, 0].set_title('t-SNE Visualization')
    axes[1, 0].grid(True, alpha=0.3)
    
    # UMAP plot
    axes[1, 1].scatter(umap_result[:, 0], umap_result[:, 1], c=colors, alpha=0.7, s=60)
    axes[1, 1].set_xlabel('UMAP 1')
    axes[1, 1].set_ylabel('UMAP 2')
    axes[1, 1].set_title('UMAP Visualization')
    axes[1, 1].grid(True, alpha=0.3)
    
    # Sample distances heatmap
    from scipy.spatial.distance import pdist, squareform
    distances = squareform(pdist(expr_matrix, metric='euclidean'))
    im = axes[1, 2].imshow(distances, cmap='viridis', aspect='auto')
    axes[1, 2].set_title('Sample Distance Matrix')
    axes[1, 2].set_xlabel('Samples')
    axes[1, 2].set_ylabel('Samples')
    
    # Add colorbar
    plt.colorbar(im, ax=axes[1, 2])
    
    # Create legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=color_map[condition], label=condition) 
                      for condition in color_map.keys() if condition in conditions]
    fig.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(0.98, 0.98))
    
    plt.tight_layout()
    plt.savefig('figures/22_dimensionality_reduction.png', dpi=300, bbox_inches='tight')
    plt.savefig('figures/22_dimensionality_reduction.pdf', bbox_inches='tight')
    plt.close()
    
    print("Dimensionality reduction plots saved to figures/22_dimensionality_reduction.*")

# Check if required libraries are available
try:
    create_pca_plots()
except ImportError as e:
    print(f"Skipping PCA plots due to missing library: {e}")
    print("Install missing libraries with: pip install scikit-learn umap-learn")

# =============================================================================
# 6. GENE SET OVERLAP ANALYSIS
# =============================================================================

print("\n=== Creating Gene Set Overlap Analysis ===")

def create_overlap_analysis():
    """Create Venn diagrams and upset plots for gene set overlaps."""
    
    from matplotlib_venn import venn2, venn3
    
    # Create sets
    up_genes = set(upregulated_genes['Gene_Symbol'].dropna())
    down_genes = set(downregulated_genes['Gene_Symbol'].dropna())
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Venn diagram for up/down regulation
    if len(up_genes) > 0 and len(down_genes) > 0:
        venn2([up_genes, down_genes], 
              set_labels=('Upregulated', 'Downregulated'),
              ax=axes[0, 0])
        axes[0, 0].set_title('Gene Set Overlap: Up vs Down Regulation')
    
    # Expression distribution comparison
    if len(significant_genes) > 0:
        up_fc = upregulated_genes['logFC'] if len(upregulated_genes) > 0 else []
        down_fc = downregulated_genes['logFC'] if len(downregulated_genes) > 0 else []
        
        axes[0, 1].hist(up_fc, bins=30, alpha=0.7, label='Upregulated', color='red')
        axes[0, 1].hist(down_fc, bins=30, alpha=0.7, label='Downregulated', color='blue')
        axes[0, 1].set_xlabel('Log2 Fold Change')
        axes[0, 1].set_ylabel('Frequency')
        axes[0, 1].set_title('Distribution of Fold Changes')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)
    
    # P-value distribution
    axes[1, 0].hist(de_results['adj.P.Val'], bins=50, alpha=0.7, color='purple')
    axes[1, 0].set_xlabel('Adjusted P-value')
    axes[1, 0].set_ylabel('Frequency')
    axes[1, 0].set_title('Distribution of Adjusted P-values')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Significance categories pie chart
    sig_counts = de_results['Significance'].value_counts()
    axes[1, 1].pie(sig_counts.values, labels=sig_counts.index, autopct='%1.1f%%',
                   colors=['#FF6B6B', '#4ECDC4', '#95A5A6'])
    axes[1, 1].set_title('Gene Categories Distribution')
    
    plt.tight_layout()
    plt.savefig('figures/23_gene_set_analysis.png', dpi=300, bbox_inches='tight')
    plt.savefig('figures/23_gene_set_analysis.pdf', bbox_inches='tight')
    plt.close()
    
    print("Gene set overlap analysis saved to figures/23_gene_set_analysis.*")

try:
    create_overlap_analysis()
except ImportError as e:
    print(f"Skipping Venn diagrams due to missing library: {e}")
    print("Install with: pip install matplotlib-venn")

# =============================================================================
# 7. PATHWAY ENRICHMENT VISUALIZATION
# =============================================================================

print("\n=== Creating Pathway Enrichment Visualization ===")

def create_pathway_visualization():
    """Create pathway enrichment visualization from saved results."""
    
    # Try to load enrichment results
    try:
        kegg_results = pd.read_csv('data/results/kegg_all_genes.csv')
        go_results = pd.read_csv('data/results/go_bp_all_genes.csv')
    except FileNotFoundError:
        print("Enrichment results not found. Skipping pathway visualization.")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(20, 12))
    
    # KEGG pathway enrichment
    if not kegg_results.empty:
        top_kegg = kegg_results.head(15)
        y_pos = np.arange(len(top_kegg))
        
        axes[0, 0].barh(y_pos, -np.log10(top_kegg['p.adjust']), 
                       color='steelblue', alpha=0.7)
        axes[0, 0].set_yticks(y_pos)
        axes[0, 0].set_yticklabels(top_kegg['Description'], fontsize=10)
        axes[0, 0].set_xlabel('-Log10(Adjusted P-value)')
        axes[0, 0].set_title('KEGG Pathway Enrichment')
        axes[0, 0].grid(True, alpha=0.3, axis='x')
        
        # KEGG pathway gene ratio plot
        gene_ratios = top_kegg['GeneRatio'].apply(lambda x: eval(x) if isinstance(x, str) else x)
        axes[0, 1].barh(y_pos, gene_ratios, color='lightcoral', alpha=0.7)
        axes[0, 1].set_yticks(y_pos)
        axes[0, 1].set_yticklabels(top_kegg['Description'], fontsize=10)
        axes[0, 1].set_xlabel('Gene Ratio')
        axes[0, 1].set_title('KEGG Pathway Gene Ratio')
        axes[0, 1].grid(True, alpha=0.3, axis='x')
    
    # GO biological process enrichment
    if not go_results.empty:
        top_go = go_results.head(15)
        y_pos_go = np.arange(len(top_go))
        
        axes[1, 0].barh(y_pos_go, -np.log10(top_go['p.adjust']), 
                       color='forestgreen', alpha=0.7)
        axes[1, 0].set_yticks(y_pos_go)
        axes[1, 0].set_yticklabels(top_go['Description'], fontsize=9)
        axes[1, 0].set_xlabel('-Log10(Adjusted P-value)')
        axes[1, 0].set_title('GO Biological Process Enrichment')
        axes[1, 0].grid(True, alpha=0.3, axis='x')
        
        # GO pathway gene count plot
        axes[1, 1].barh(y_pos_go, top_go['Count'], color='gold', alpha=0.7)
        axes[1, 1].set_yticks(y_pos_go)
        axes[1, 1].set_yticklabels(top_go['Description'], fontsize=9)
        axes[1, 1].set_xlabel('Gene Count')
        axes[1, 1].set_title('GO Biological Process Gene Count')
        axes[1, 1].grid(True, alpha=0.3, axis='x')
    
    plt.tight_layout()
    plt.savefig('figures/24_pathway_enrichment.png', dpi=300, bbox_inches='tight')
    plt.savefig('figures/24_pathway_enrichment.pdf', bbox_inches='tight')
    plt.close()
    
    print("Pathway enrichment visualization saved to figures/24_pathway_enrichment.*")

create_pathway_visualization()

# =============================================================================
# 8. SUMMARY DASHBOARD
# =============================================================================

print("\n=== Creating Summary Dashboard ===")

def create_summary_dashboard():
    """Create a comprehensive summary dashboard."""
    
    fig = plt.figure(figsize=(20, 16))
    gs = fig.add_gridspec(4, 4, hspace=0.3, wspace=0.3)
    
    # Title
    fig.suptitle('Lung Cancer Differential Expression Analysis - Summary Dashboard', 
                 fontsize=20, fontweight='bold', y=0.98)
    
    # 1. Gene categories pie chart (top-left)
    ax1 = fig.add_subplot(gs[0, 0])
    sig_counts = de_results['Significance'].value_counts()
    ax1.pie(sig_counts.values, labels=sig_counts.index, autopct='%1.1f%%',
            colors=['#FF6B6B', '#4ECDC4', '#95A5A6'])
    ax1.set_title('Gene Categories', fontweight='bold')
    
    # 2. Top upregulated genes (top-middle)
    ax2 = fig.add_subplot(gs[0, 1])
    if len(upregulated_genes) > 0:
        top_up = upregulated_genes.head(10)
        y_pos = np.arange(len(top_up))
        ax2.barh(y_pos, top_up['logFC'], color='red', alpha=0.7)
        ax2.set_yticks(y_pos)
        ax2.set_yticklabels(top_up['Gene_Symbol'], fontsize=8)
        ax2.set_xlabel('Log2 Fold Change')
        ax2.set_title('Top Upregulated Genes', fontweight='bold')
        ax2.grid(True, alpha=0.3, axis='x')
    
    # 3. Top downregulated genes (top-right)
    ax3 = fig.add_subplot(gs[0, 2:])
    if len(downregulated_genes) > 0:
        top_down = downregulated_genes.head(10)
        y_pos = np.arange(len(top_down))
        ax3.barh(y_pos, top_down['logFC'], color='blue', alpha=0.7)
        ax3.set_yticks(y_pos)
        ax3.set_yticklabels(top_down['Gene_Symbol'], fontsize=8)
        ax3.set_xlabel('Log2 Fold Change')
        ax3.set_title('Top Downregulated Genes', fontweight='bold')
        ax3.grid(True, alpha=0.3, axis='x')
    
    # 4. Volcano plot (middle-left)
    ax4 = fig.add_subplot(gs[1, :2])
    volcano_data = de_results.copy()
    volcano_data['neg_log10_pval'] = -np.log10(volcano_data['adj.P.Val'].clip(lower=1e-300))
    
    colors = {'Upregulated': '#FF6B6B', 'Downregulated': '#4ECDC4', 'Not Significant': '#95A5A6'}
    
    for category in colors.keys():
        data_subset = volcano_data[volcano_data['Significance'] == category]
        ax4.scatter(data_subset['logFC'], data_subset['neg_log10_pval'],
                   c=colors[category], alpha=0.6, s=20, label=category)
    
    ax4.axvline(x=1, color='black', linestyle='--', alpha=0.5)
    ax4.axvline(x=-1, color='black', linestyle='--', alpha=0.5)
    ax4.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
    ax4.set_xlabel('Log2 Fold Change')
    ax4.set_ylabel('-Log10(Adjusted P-value)')
    ax4.set_title('Volcano Plot', fontweight='bold')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # 5. Expression heatmap (middle-right)
    ax5 = fig.add_subplot(gs[1, 2:])
    if len(significant_genes) > 0:
        top_genes = significant_genes.head(20)
        heatmap_data = expr_data.loc[top_genes.index].values
        
        im = ax5.imshow(heatmap_data, cmap='RdBu_r', aspect='auto')
        ax5.set_title('Top 20 DEGs Expression', fontweight='bold')
        ax5.set_xlabel('Samples')
        ax5.set_ylabel('Genes')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax5, fraction=0.046, pad=0.04)
        cbar.set_label('Expression Level')
    
    # 6. Statistics summary (bottom-left)
    ax6 = fig.add_subplot(gs[2, :2])
    ax6.axis('off')
    
    stats_text = f"""
    ANALYSIS STATISTICS
    
    Total Genes Analyzed: {len(de_results):,}
    Significant Genes: {len(significant_genes):,}
    Upregulated: {len(upregulated_genes):,}
    Downregulated: {len(downregulated_genes):,}
    
    Samples:
    • Normal: {sum(sample_metadata['condition'] == 'Normal')}
    • Cancer: {sum(sample_metadata['condition'] == 'Cancer')}
    
    Significance Thresholds:
    • FDR < 0.05
    • |Log2FC| > 1.0
    """
    
    ax6.text(0.1, 0.9, stats_text, transform=ax6.transAxes, fontsize=12,
            verticalalignment='top', bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue"))
    
    # 7. Sample correlation heatmap (bottom-right)
    ax7 = fig.add_subplot(gs[2, 2:])
    
    # Calculate sample correlations (using subset for speed)
    subset_expr = expr_data.iloc[:1000].T  # Use top 1000 genes for correlation
    correlation_matrix = subset_expr.corr()
    
    im = ax7.imshow(correlation_matrix.values, cmap='coolwarm', vmin=-1, vmax=1)
    ax7.set_title('Sample Correlations', fontweight='bold')
    ax7.set_xlabel('Samples')
    ax7.set_ylabel('Samples')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax7, fraction=0.046, pad=0.04)
    cbar.set_label('Correlation')
    
    # 8. Analysis workflow (bottom)
    ax8 = fig.add_subplot(gs[3, :])
    ax8.axis('off')
    
    workflow_text = """
    ANALYSIS WORKFLOW: Data Download → Quality Control → Differential Expression → Functional Enrichment → Visualization → Network Analysis
    """
    
    ax8.text(0.5, 0.5, workflow_text, transform=ax8.transAxes, fontsize=14,
            horizontalalignment='center', verticalalignment='center',
            bbox=dict(boxstyle="round,pad=0.5", facecolor="lightyellow"))
    
    plt.savefig('figures/25_summary_dashboard.png', dpi=300, bbox_inches='tight')
    plt.savefig('figures/25_summary_dashboard.pdf', bbox_inches='tight')
    plt.close()
    
    print("Summary dashboard saved to figures/25_summary_dashboard.*")

create_summary_dashboard()

# =============================================================================
# 9. SAVE VISUALIZATION SUMMARY
# =============================================================================

print("\n=== Saving Visualization Summary ===")

visualization_files = [
    "figures/19_enhanced_volcano_plot.*",
    "figures/20_interactive_volcano.html",
    "figures/21_expression_heatmap_clustered.*",
    "figures/22_dimensionality_reduction.*",
    "figures/23_gene_set_analysis.*",
    "figures/24_pathway_enrichment.*",
    "figures/25_summary_dashboard.*"
]

summary_text = f"""
ADVANCED VISUALIZATION ANALYSIS SUMMARY
=====================================

Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}
Dataset: GSE43346 (Small Cell Lung Cancer)

VISUALIZATIONS CREATED:

1. Enhanced Volcano Plot (figures/19_enhanced_volcano_plot.*)
   - Standard and density-enhanced volcano plots
   - Annotated with top significant genes
   - Color-coded by significance categories

2. Interactive Volcano Plot (figures/20_interactive_volcano.html)
   - Plotly-based interactive visualization
   - Hover information for each gene
   - Zoomable and exportable

3. Expression Heatmap (figures/21_expression_heatmap_clustered.*)
   - Top 100 most significant genes
   - Hierarchical clustering of samples and genes
   - Z-score normalized expression values

4. Dimensionality Reduction (figures/22_dimensionality_reduction.*)
   - PCA analysis with variance explained
   - t-SNE visualization
   - UMAP clustering
   - Sample distance matrix

5. Gene Set Analysis (figures/23_gene_set_analysis.*)
   - Venn diagrams for gene overlaps
   - Fold change distributions
   - P-value distributions
   - Category proportions

6. Pathway Enrichment (figures/24_pathway_enrichment.*)
   - KEGG pathway enrichment bars
   - GO biological process enrichment
   - Gene ratios and counts

7. Summary Dashboard (figures/25_summary_dashboard.*)
   - Comprehensive overview of all results
   - Key statistics and visualizations
   - Analysis workflow summary

DATA PROCESSED:
- Total genes: {len(de_results):,}
- Significant genes: {len(significant_genes):,}
- Upregulated: {len(upregulated_genes):,}
- Downregulated: {len(downregulated_genes):,}
- Samples: {len(sample_metadata)}

All visualizations are publication-ready and saved in both PNG (300 DPI) and PDF formats.
"""

with open('reports/visualization_summary.txt', 'w') as f:
    f.write(summary_text)

print("Visualization summary saved to reports/visualization_summary.txt")
print("\n=== Advanced Visualization Analysis Completed Successfully ===")
print("Next step: Run 06_network_analysis.py for gene network analysis")
    colors = {
        'Upregulated': '#FF6B6B',
        'Downregulated': '#4ECDC4', 
        'Not Significant': '#95A5A6'
    }
    
    # Plot 1: Standard volcano plot
    for category in colors.keys():
        data_subset = volcano_data[volcano_data['Significance'] == category]
        ax1.scatter(data_subset['logFC'], 
                   data_subset['neg_log10_pval'],
                   c=colors[category],
                   alpha=0.6,
                   s=30,
                   label=f"{category} ({len(data_subset)})")
    
    # Add threshold lines
    ax1.axvline(x=logfc_threshold, color='black', linestyle='--', alpha=0.5)
    ax1.axvline(x=-logfc_threshold, color='black', linestyle='--', alpha=0.5)
    ax1.axhline(y=-np.log10(pval_threshold), color='black', linestyle='--', alpha=0.5)
    
    # Annotate top genes
    top_genes = significant_genes.head(20)
    for idx, row in top_genes.iterrows():
        if idx in volcano_data.index:
            ax1.annotate(row['Gene_Symbol'], 
                        xy=(volcano_data.loc[idx, 'logFC'], 
                           volcano_data.loc[idx, 'neg_log10_pval']),
                        xytext=(5, 5), 
                        textcoords='offset points',
                        fontsize=8,
                        alpha=0.8)
    
    ax1.set_xlabel('Log2 Fold Change')
    ax1.set_ylabel('-Log10(Adjusted P-value)')
    ax1.set_title('Enhanced Volcano Plot: Cancer vs Normal')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Density-enhanced volcano plot
    scatter = ax2.scatter(volcano_data['logFC'], 
                         volcano_data['neg_log10_pval'],
                         c=volcano_data['abs_logFC'],
                         cmap='viridis',
                         alpha=0.6,
                         s=30)
    
    ax2.axvline(x=logfc_threshold, color='white', linestyle='--', alpha=0.8)
    ax2.axvline(x=-logfc_threshold, color='white', linestyle='--', alpha=0.8)
    ax2.axhline(y=-np.log10(pval_threshold), color='white', linestyle='--', alpha=0.8)
    
    ax2.set_xlabel('Log2 Fold Change')
    ax2.set_ylabel('-Log10(Adjusted P-value)')
    ax2.set_title('Volcano Plot: Colored by Absolute Fold Change')
    ax2.grid(True, alpha=0.3)
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax2)
    cbar.set_label('Absolute Log2 Fold Change')
    
    plt.tight_layout()
    plt.savefig('figures/19_enhanced_volcano_plot.png', dpi=300, bbox_inches='tight')
    plt.savefig('figures/19_enhanced_volcano_plot.pdf', bbox_inches='tight')
    plt.close()
    
    print("Enhanced volcano plot saved to figures/19_enhanced_volcano_plot.*")

create_enhanced_volcano_plot()

# =============================================================================
# 3. INTERACTIVE VOLCANO PLOT
# =============================================================================

print("\n=== Creating Interactive Volcano Plot ===")

def create_interactive_volcano():
    """Create an interactive volcano plot using Plotly."""
    
    volcano_data = de_results.copy()
    volcano_data['neg_log10_pval'] = -np.log10(volcano_data['adj.P.Val'].clip(lower=1e-300))
    
    # Create hover text
    volcano_data['hover_text'] = (
        volcano_data['Gene_Symbol'] + '<br>' +
        'Log2FC: ' + volcano_data['logFC'].round(3).astype(str) + '<br>' +
        'Adj P-val: ' + volcano_data['adj.P.Val'].apply(lambda x: f'{x:.2e}') + '<br>' +
        'Significance: ' + volcano_data['Significance']
    )
    
    # Color mapping
    color_map = {
        'Upregulated': '#FF6B6B',
        'Downregulated': '#4ECDC4',
        'Not Significant': '#95A5A6'
    }
    
    fig = go.Figure()
    
    for category in color_map.keys():
        data_subset = volcano_data[volcano_data['Significance'] == category]
        
        fig.add_trace(go.Scatter(
            x=data_subset['logFC'],
            y=data_subset['neg_log10_pval'],
            mode='markers',
            name=f"{category} ({len(data_subset)})",
            marker=dict(
                color=color_map[category],
                size=6,
                opacity=0.7
            ),
            text=data_subset['hover_text'],
            hovertemplate='%{text}<extra></extra>'
        ))
    
    # Add threshold lines
    fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="black", opacity=0.5)
    fig.add_vline(x=1, line_dash="dash", line_color="black", opacity=0.5)
    fig.add_vline(x=-1, line_dash="dash", line_color="black", opacity=0.5)
    
    fig.update_layout(
        title="Interactive Volcano Plot: Cancer vs Normal",
        xaxis_title="Log2 Fold Change",
        yaxis_title="-Log10(Adjusted P-value)",
        width=800,
        height=600,
        showlegend=True
    )
    
    fig.write_html("figures/20_interactive_volcano.html")
    print("Interactive volcano plot saved to figures/20_interactive_volcano.html")

create_interactive_volcano()

# =============================================================================
# 4. EXPRESSION HEATMAP WITH CLUSTERING
# =============================================================================

print("\n=== Creating Expression Heatmap with Clustering ===")

def create_expression_heatmap():
    """Create a comprehensive expression heatmap."""
    
    # Select top 100 most significant genes for heatmap
    top_genes = significant_genes.head(100)
    
    if len(top_genes) == 0:
        print("No significant genes found for heatmap")
        return
    
    # Extract expression data for top genes
    heatmap_data = expr_data.loc[top_genes.index].T
    
    # Prepare sample annotation
    sample_colors = sample_metadata.set_index('sample_id')['condition'].map({
        'Normal': '#2E8B57',
        'Cancer': '#DC143C'
    })
    
    # Create the heatmap
    fig, ax = plt.subplots(figsize=(15, 12))
    
    # Z-score normalization
    heatmap_data_zscore = heatmap_data.apply(lambda x: (x - x.mean()) / x.std(), axis=0)
    
    # Create clustermap
    g = sns.clustermap(heatmap_data_zscore.T, 
                      cmap='RdBu_r',
                      center=0,
                      figsize=(15, 12),
                      row_cluster=True,
                      col_cluster=True,
                      col_colors=sample_colors,
                      cbar_kws={'label': 'Z-score'},
                      xticklabels=False,
                      yticklabels=True)
    
    # Adjust layout
    g.ax_row_dendrogram.set_visible(True)
    g.ax_col_dendrogram.set_visible(True)
    
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), fontsize=6)
    
    g.fig.suptitle('Expression Heatmap: Top 100 Significant Genes', 
                   fontsize=16, y=0.95)
    
    plt.savefig('figures/21_expression_heatmap_clustered.png', 
                dpi=300, bbox_inches='tight')
    plt.savefig('figures/21_expression_heatmap_clustered.pdf', 
                bbox_inches='tight')
    plt.close()
    
    print("Expression heatmap saved to figures/21_expression_heatmap_clustered.*")

create_expression_heatmap()

# =============================================================================
# 5. PCA AND DIMENSIONALITY REDUCTION PLOTS
# =============================================================================

print("\n=== Creating PCA and Dimensionality Reduction Plots ===")

def create_pca_plots():
    """Create comprehensive PCA analysis plots."""
    
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE
    import umap
    
    # Prepare data
    expr_matrix = expr_data.values.T  # Transpose for PCA (samples as rows)
    conditions = sample_metadata['condition'].values
    
    # PCA
    pca = PCA(n_components=10)
    pca_result = pca.fit_transform(expr_matrix)
    
    # t-SNE
    tsne = TSNE(n_components=2, random_state=42, perplexity=min(30, len(expr_matrix)-1))
    tsne_result = tsne.fit_transform(expr_matrix)
    
    # UMAP
    umap_reducer = umap.UMAP(random_state=42)
    umap_result = umap_reducer.fit_transform(expr_matrix)
    
    # Create subplots
    fig, axes = plt.subplots(2, 3, figsize=(20, 12))
    
    # Color mapping
