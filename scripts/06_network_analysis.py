#!/usr/bin/env python3

"""
Script: 06_network_analysis.py
Purpose: Gene network analysis and protein-protein interaction analysis
Author: Bioinformatics Analysis Project
Date: 2025
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from scipy.stats import pearsonr
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import warnings
warnings.filterwarnings('ignore')

print("=== Starting Gene Network Analysis ===")

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
    
    print(f"Loaded data successfully:")
    print(f"- Significant genes: {len(significant_genes)}")
    print(f"- Expression matrix: {expr_data.shape}")
    
except FileNotFoundError as e:
    print(f"Error loading data: {e}")
    print("Please run the previous scripts (01-05) to generate the required data files.")
    exit(1)

# =============================================================================
# 2. GENE CO-EXPRESSION NETWORK ANALYSIS
# =============================================================================

print("\n=== Building Gene Co-expression Network ===")

def build_coexpression_network(expression_data, gene_list, correlation_threshold=0.8):
    """Build co-expression network from expression data."""
    
    # Select genes of interest
    if isinstance(gene_list, pd.DataFrame):
        genes_of_interest = list(gene_list.index)
    else:
        genes_of_interest = gene_list
    
    # Subset expression data
    expr_subset = expression_data.loc[genes_of_interest]
    
    # Calculate correlation matrix
    correlation_matrix = expr_subset.T.corr()
    
    # Create adjacency matrix based on correlation threshold
    adjacency_matrix = (abs(correlation_matrix) >= correlation_threshold).astype(int)
    np.fill_diagonal(adjacency_matrix.values, 0)  # Remove self-loops
    
    # Create network graph
    G = nx.from_pandas_adjacency(adjacency_matrix)
    
    # Add node attributes (fold change, p-value)
    for node in G.nodes():
        if node in de_results.index:
            G.nodes[node]['logFC'] = de_results.loc[node, 'logFC']
            G.nodes[node]['adj_pval'] = de_results.loc[node, 'adj.P.Val']
            G.nodes[node]['gene_symbol'] = de_results.loc[node, 'Gene_Symbol']
        else:
            G.nodes[node]['logFC'] = 0
            G.nodes[node]['adj_pval'] = 1
            G.nodes[node]['gene_symbol'] = node
    
    # Add edge attributes (correlation values)
    for edge in G.edges():
        gene1, gene2 = edge
        correlation_value = correlation_matrix.loc[gene1, gene2]
        G.edges[edge]['correlation'] = correlation_value
        G.edges[edge]['abs_correlation'] = abs(correlation_value)
    
    return G, correlation_matrix

# Build network for significant genes (top 200 for computational efficiency)
top_sig_genes = significant_genes.head(200) if len(significant_genes) > 200 else significant_genes

network, corr_matrix = build_coexpression_network(
    expr_data, 
    top_sig_genes, 
    correlation_threshold=0.7
)

print(f"Built co-expression network:")
print(f"- Nodes: {network.number_of_nodes()}")
print(f"- Edges: {network.number_of_edges()}")
print(f"- Density: {nx.density(network):.3f}")

# =============================================================================
# 3. NETWORK ANALYSIS AND METRICS
# =============================================================================

print("\n=== Analyzing Network Properties ===")

def analyze_network_properties(G):
    """Calculate network properties and centrality measures."""
    
    if G.number_of_nodes() == 0:
        print("Network is empty. Skipping network analysis.")
        return None
    
    # Basic network properties
    properties = {
        'nodes': G.number_of_nodes(),
        'edges': G.number_of_edges(),
        'density': nx.density(G),
        'is_connected': nx.is_connected(G)
    }
    
    if not nx.is_connected(G):
        properties['connected_components'] = nx.number_connected_components(G)
        properties['largest_cc_size'] = len(max(nx.connected_components(G), key=len))
        
        # Use largest connected component for further analysis
        largest_cc = max(nx.connected_components(G), key=len)
        G_cc = G.subgraph(largest_cc).copy()
    else:
        G_cc = G.copy()
        properties['connected_components'] = 1
        properties['largest_cc_size'] = G.number_of_nodes()
    
    # Calculate centrality measures
    if G_cc.number_of_nodes() > 0:
        degree_centrality = nx.degree_centrality(G_cc)
        betweenness_centrality = nx.betweenness_centrality(G_cc)
        closeness_centrality = nx.closeness_centrality(G_cc)
        
        # Eigenvector centrality (may fail for some graphs)
        try:
            eigenvector_centrality = nx.eigenvector_centrality(G_cc, max_iter=1000)
        except:
            eigenvector_centrality = {node: 0 for node in G_cc.nodes()}
        
        # Add centrality as node attributes
        nx.set_node_attributes(G, degree_centrality, 'degree_centrality')
        nx.set_node_attributes(G, betweenness_centrality, 'betweenness_centrality')
        nx.set_node_attributes(G, closeness_centrality, 'closeness_centrality')
        nx.set_node_attributes(G, eigenvector_centrality, 'eigenvector_centrality')
    
    return properties, G_cc

network_props, largest_component = analyze_network_properties(network)

if network_props:
    print("Network Properties:")
    for key, value in network_props.items():
        print(f"- {key}: {value}")

# =============================================================================
# 4. IDENTIFY HUB GENES
# =============================================================================

print("\n=== Identifying Hub Genes ===")

def identify_hub_genes(G, top_n=20):
    """Identify hub genes based on centrality measures."""
    
    if G.number_of_nodes() == 0:
        return pd.DataFrame()
    
    # Create dataframe with centrality measures
    hub_data = []
    
    for node in G.nodes():
        node_data = {
            'gene_id': node,
            'gene_symbol': G.nodes[node].get('gene_symbol', node),
            'degree': G.degree(node),
            'degree_centrality': G.nodes[node].get('degree_centrality', 0),
            'betweenness_centrality': G.nodes[node].get('betweenness_centrality', 0),
            'closeness_centrality': G.nodes[node].get('closeness_centrality', 0),
            'eigenvector_centrality': G.nodes[node].get('eigenvector_centrality', 0),
            'logFC': G.nodes[node].get('logFC', 0),
            'adj_pval': G.nodes[node].get('adj_pval', 1)
        }
        hub_data.append(node_data)
    
    hub_df = pd.DataFrame(hub_data)
    
    # Calculate composite hub score
    centrality_cols = ['degree_centrality', 'betweenness_centrality', 
                      'closeness_centrality', 'eigenvector_centrality']
    
    # Normalize centrality measures
    scaler = StandardScaler()
    hub_df[centrality_cols] = scaler.fit_transform(hub_df[centrality_cols])
    
    # Calculate composite score
    hub_df['hub_score'] = hub_df[centrality_cols].mean(axis=1)
    
    # Sort by hub score
    hub_df = hub_df.sort_values('hub_score', ascending=False)
    
    return hub_df.head(top_n)

hub_genes = identify_hub_genes(network)

if not hub_genes.empty:
    print(f"Top 10 Hub Genes:")
    for _, row in hub_genes.head(10).iterrows():
        print(f"- {row['gene_symbol']} (score: {row['hub_score']:.3f}, degree: {row['degree']})")
    
    # Save hub genes
    hub_genes.to_csv('data/results/hub_genes.csv', index=False)
    print("Hub genes saved to: data/results/hub_genes.csv")

# =============================================================================
# 5. NETWORK CLUSTERING
# =============================================================================

print("\n=== Performing Network Clustering ===")

def perform_network_clustering(G):
    """Perform community detection on the network."""
    
    if G.number_of_nodes() < 3:
        print("Network too small for clustering analysis")
        return None, None
    
    # Community detection using different algorithms
    try:
        # Louvain algorithm
        communities_louvain = nx.community.louvain_communities(G, seed=42)
        
        # Greedy modularity communities
        communities_greedy = nx.community.greedy_modularity_communities(G)
        
        print(f"Community Detection Results:")
        print(f"- Louvain communities: {len(communities_louvain)}")
        print(f"- Greedy modularity communities: {len(communities_greedy)}")
        
        # Add community assignments as node attributes
        for i, community in enumerate(communities_louvain):
            for node in community:
                G.nodes[node]['louvain_community'] = i
        
        for i, community in enumerate(communities_greedy):
            for node in community:
                G.nodes[node]['greedy_community'] = i
        
        return communities_louvain, communities_greedy
        
    except Exception as e:
        print(f"Error in community detection: {e}")
        return None, None

communities_louvain, communities_greedy = perform_network_clustering(network)

# =============================================================================
# 6. NETWORK VISUALIZATION
# =============================================================================

print("\n=== Creating Network Visualizations ===")

def visualize_network(G, communities=None, save_prefix="network"):
    """Create network visualization."""
    
    if G.number_of_nodes() == 0:
        print("Network is empty. Skipping visualization.")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(20, 16))
    fig.suptitle('Gene Co-expression Network Analysis', fontsize=16, fontweight='bold')
    
    # Use largest connected component for visualization if network is disconnected
    if not nx.is_connected(G):
        largest_cc = max(nx.connected_components(G), key=len)
        G_vis = G.subgraph(largest_cc)
    else:
        G_vis = G
    
    # Layout for visualization
    if G_vis.number_of_nodes() > 100:
        pos = nx.spring_layout(G_vis, k=1, iterations=50, seed=42)
    else:
        pos = nx.spring_layout(G_vis, k=2, iterations=100, seed=42)
    
    # 1. Network colored by fold change
    ax1 = axes[0, 0]
    node_colors = [G_vis.nodes[node].get('logFC', 0) for node in G_vis.nodes()]
    node_sizes = [G_vis.degree(node) * 20 + 50 for node in G_vis.nodes()]
    
    nx.draw_networkx_nodes(G_vis, pos, node_color=node_colors, 
                          node_size=node_sizes, cmap='RdBu_r', 
                          alpha=0.8, ax=ax1)
    nx.draw_networkx_edges(G_vis, pos, alpha=0.3, width=0.5, ax=ax1)
    
    ax1.set_title('Network colored by Log2 Fold Change')
    ax1.axis('off')
    
    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap='RdBu_r', 
                              norm=plt.Normalize(vmin=min(node_colors), 
                                               vmax=max(node_colors)))
    sm.set_array([])
    cbar1 = plt.colorbar(sm, ax=ax1, fraction=0.046, pad=0.04)
    cbar1.set_label('Log2 Fold Change')
    
    # 2. Network colored by degree centrality
    ax2 = axes[0, 1]
    centrality_colors = [G_vis.nodes[node].get('degree_centrality', 0) for node in G_vis.nodes()]
    
    nx.draw_networkx_nodes(G_vis, pos, node_color=centrality_colors, 
                          node_size=node_sizes, cmap='viridis', 
                          alpha=0.8, ax=ax2)
    nx.draw_networkx_edges(G_vis, pos, alpha=0.3, width=0.5, ax=ax2)
    
    ax2.set_title('Network colored by Degree Centrality')
    ax2.axis('off')
    
    # Add colorbar
    sm2 = plt.cm.ScalarMappable(cmap='viridis', 
                               norm=plt.Normalize(vmin=min(centrality_colors), 
                                                vmax=max(centrality_colors)))
    sm2.set_array([])
    cbar2 = plt.colorbar(sm2, ax=ax2, fraction=0.046, pad=0.04)
    cbar2.set_label('Degree Centrality')
    
    # 3. Community structure (if available)
    ax3 = axes[1, 0]
    if communities:
        community_colors = plt.cm.Set3(np.linspace(0, 1, len(communities)))
        node_colors_comm = []
        
        for node in G_vis.nodes():
            for i, community in enumerate(communities):
                if node in community:
                    node_colors_comm.append(community_colors[i])
                    break
            else:
                node_colors_comm.append('gray')
        
        nx.draw_networkx_nodes(G_vis, pos, node_color=node_colors_comm, 
                              node_size=node_sizes, alpha=0.8, ax=ax3)
        nx.draw_networkx_edges(G_vis, pos, alpha=0.3, width=0.5, ax=ax3)
        
        ax3.set_title(f'Network Communities ({len(communities)} communities)')
    else:
        ax3.text(0.5, 0.5, 'Community detection\nnot available', 
                transform=ax3.transAxes, ha='center', va='center')
        ax3.set_title('Network Communities')
    
    ax3.axis('off')
    
    # 4. Degree distribution
    ax4 = axes[1, 1]
    degrees = [G.degree(node) for node in G.nodes()]
    ax4.hist(degrees, bins=20, alpha=0.7, color='skyblue', edgecolor='black')
    ax4.set_xlabel('Degree')
    ax4.set_ylabel('Frequency')
    ax4.set_title('Degree Distribution')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'figures/26_{save_prefix}_visualization.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'figures/26_{save_prefix}_visualization.pdf', bbox_inches='tight')
    plt.close()
    
    print(f"Network visualization saved to figures/26_{save_prefix}_visualization.*")

if network.number_of_nodes() > 0:
    visualize_network(network, communities_louvain, "coexpression")

# =============================================================================
# 7. CORRELATION HEATMAP OF TOP GENES
# =============================================================================

print("\n=== Creating Correlation Heatmap ===")

def create_correlation_heatmap(correlation_matrix, top_genes_df, save_name="correlation_heatmap"):
    """Create correlation heatmap for top genes."""
    
    if correlation_matrix.empty or top_genes_df.empty:
        print("No data available for correlation heatmap")
        return
    
    # Select top genes for heatmap (limit to 50 for readability)
    n_genes = min(50, len(top_genes_df))
    top_genes_list = list(top_genes_df.index[:n_genes])
    
    # Get gene symbols for labeling
    gene_symbols = []
    for gene_id in top_genes_list:
        if gene_id in de_results.index:
            symbol = de_results.loc[gene_id, 'Gene_Symbol']
            gene_symbols.append(symbol if pd.notna(symbol) else gene_id)
        else:
            gene_symbols.append(gene_id)
    
    # Create correlation matrix subset
    corr_subset = correlation_matrix.loc[top_genes_list, top_genes_list]
    
    # Create the heatmap
    fig, ax = plt.subplots(figsize=(15, 12))
    
    # Create annotation matrix for significant correlations
    mask = np.abs(corr_subset) < 0.5  # Mask weak correlations
    
    sns.heatmap(corr_subset, 
                annot=False, 
                cmap='RdBu_r', 
                center=0,
                vmin=-1, vmax=1,
                square=True,
                mask=mask,
                xticklabels=gene_symbols,
                yticklabels=gene_symbols,
                cbar_kws={'label': 'Correlation Coefficient'},
                ax=ax)
    
    plt.xticks(rotation=45, ha='right', fontsize=8)
    plt.yticks(rotation=0, fontsize=8)
    plt.title(f'Gene Co-expression Correlation Heatmap (Top {n_genes} Genes)', 
              fontsize=14, fontweight='bold', pad=20)
    
    plt.tight_layout()
    plt.savefig(f'figures/27_{save_name}.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'figures/27_{save_name}.pdf', bbox_inches='tight')
    plt.close()
    
    print(f"Correlation heatmap saved to figures/27_{save_name}.*")

if not corr_matrix.empty and not significant_genes.empty:
    create_correlation_heatmap(corr_matrix, significant_genes)

# =============================================================================
# 8. PATHWAY-BASED NETWORK ANALYSIS
# =============================================================================

print("\n=== Creating Pathway-based Network ===")

def create_pathway_network(enrichment_file, gene_set, network_name="pathway"):
    """Create network based on pathway enrichment results."""
    
    try:
        # Load pathway enrichment results
        pathway_results = pd.read_csv(enrichment_file)
        
        if pathway_results.empty:
            print(f"No pathway results found in {enrichment_file}")
            return
        
        # Create pathway-gene bipartite network
        G_pathway = nx.Graph()
        
        # Add pathway and gene nodes
        for _, row in pathway_results.head(20).iterrows():  # Top 20 pathways
            pathway_id = row['ID'] if 'ID' in row else row['Description']
            pathway_name = row['Description']
            
            # Add pathway node
            G_pathway.add_node(f"pathway_{pathway_id}", 
                             node_type='pathway', 
                             name=pathway_name,
                             pvalue=row['p.adjust'])
            
            # Add gene nodes and edges
            if 'geneID' in row and pd.notna(row['geneID']):
                genes = str(row['geneID']).split('/')
                
                for gene in genes[:10]:  # Limit to first 10 genes per pathway
                    if gene in gene_set.index:
                        G_pathway.add_node(gene, 
                                         node_type='gene',
                                         logFC=gene_set.loc[gene, 'logFC'] if 'logFC' in gene_set.columns else 0,
                                         symbol=gene_set.loc[gene, 'Gene_Symbol'] if 'Gene_Symbol' in gene_set.columns else gene)
                        
                        # Add edge between pathway and gene
                        G_pathway.add_edge(f"pathway_{pathway_id}", gene)
        
        # Visualize pathway network
        if G_pathway.number_of_nodes() > 0:
            fig, ax = plt.subplots(figsize=(16, 12))
            
            pos = nx.spring_layout(G_pathway, k=3, iterations=50, seed=42)
            
            # Separate pathway and gene nodes
            pathway_nodes = [n for n in G_pathway.nodes() if G_pathway.nodes[n]['node_type'] == 'pathway']
            gene_nodes = [n for n in G_pathway.nodes() if G_pathway.nodes[n]['node_type'] == 'gene']
            
            # Draw pathway nodes
            nx.draw_networkx_nodes(G_pathway, pos, nodelist=pathway_nodes,
                                 node_color='lightcoral', node_size=800, 
                                 alpha=0.8, ax=ax, label='Pathways')
            
            # Draw gene nodes
            gene_colors = []
            for gene in gene_nodes:
                logfc = G_pathway.nodes[gene].get('logFC', 0)
                gene_colors.append(logfc)
            
            nx.draw_networkx_nodes(G_pathway, pos, nodelist=gene_nodes,
                                 node_color=gene_colors, node_size=300,
                                 cmap='RdBu_r', alpha=0.8, ax=ax)
            
            # Draw edges
            nx.draw_networkx_edges(G_pathway, pos, alpha=0.3, ax=ax)
            
            # Add labels for pathways
            pathway_labels = {node: G_pathway.nodes[node]['name'][:30] + '...' 
                            if len(G_pathway.nodes[node]['name']) > 30 
                            else G_pathway.nodes[node]['name']
                            for node in pathway_nodes}
            
            nx.draw_networkx_labels(G_pathway, pos, labels=pathway_labels,
                                  font_size=8, ax=ax)
            
            ax.set_title(f'Pathway-Gene Network ({network_name.title()})', 
                        fontsize=14, fontweight='bold')
            ax.axis('off')
            
            plt.tight_layout()
            plt.savefig(f'figures/28_{network_name}_pathway_network.png', 
                       dpi=300, bbox_inches='tight')
            plt.savefig(f'figures/28_{network_name}_pathway_network.pdf', 
                       bbox_inches='tight')
            plt.close()
            
            print(f"Pathway network saved to figures/28_{network_name}_pathway_network.*")
            
    except FileNotFoundError:
        print(f"Pathway file {enrichment_file} not found. Skipping pathway network.")
    except Exception as e:
        print(f"Error creating pathway network: {e}")

# Create pathway networks for KEGG and GO results
create_pathway_network('data/results/kegg_all_genes.csv', significant_genes, "kegg")
create_pathway_network('data/results/go_bp_all_genes.csv', significant_genes, "go_bp")

# =============================================================================
# 9. EXPORT NETWORK FOR CYTOSCAPE
# =============================================================================

print("\n=== Exporting Network for Cytoscape ===")

def export_for_cytoscape(G, filename_prefix="network"):
    """Export network files for Cytoscape visualization."""
    
    if G.number_of_nodes() == 0:
        print("Network is empty. Skipping Cytoscape export.")
        return
    
    # Export node attributes
    node_data = []
    for node in G.nodes():
        node_attrs = G.nodes[node].copy()
        node_attrs['node_id'] = node
        node_data.append(node_attrs)
    
    node_df = pd.DataFrame(node_data)
    node_df.to_csv(f'data/results/{filename_prefix}_nodes.csv', index=False)
    
    # Export edge attributes
    edge_data = []
    for edge in G.edges():
        edge_attrs = G.edges[edge].copy()
        edge_attrs['source'] = edge[0]
        edge_attrs['target'] = edge[1]
        edge_data.append(edge_attrs)
    
    edge_df = pd.DataFrame(edge_data)
    edge_df.to_csv(f'data/results/{filename_prefix}_edges.csv', index=False)
    
    # Export as GraphML (preferred format for Cytoscape)
    nx.write_graphml(G, f'data/results/{filename_prefix}.graphml')
    
    print(f"Network exported for Cytoscape:")
    print(f"- Nodes: data/results/{filename_prefix}_nodes.csv")
    print(f"- Edges: data/results/{filename_prefix}_edges.csv") 
    print(f"- GraphML: data/results/{filename_prefix}.graphml")

if network.number_of_nodes() > 0:
    export_for_cytoscape(network, "coexpression_network")

# =============================================================================
# 10. NETWORK ANALYSIS SUMMARY REPORT
# =============================================================================

print("\n=== Generating Network Analysis Report ===")

def generate_network_report():
    """Generate comprehensive network analysis report."""
    
    report_lines = [
        "GENE NETWORK ANALYSIS REPORT",
        "===========================",
        "",
        f"Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"Dataset: GSE43346 (Small Cell Lung Cancer)",
        "",
        "NETWORK CONSTRUCTION:",
        f"- Input genes: {len(top_sig_genes)}",
        f"- Correlation threshold: 0.7",
        f"- Network nodes: {network.number_of_nodes()}",
        f"- Network edges: {network.number_of_edges()}",
        f"- Network density: {nx.density(network):.4f}",
        ""
    ]
    
    if network_props:
        report_lines.extend([
            "NETWORK PROPERTIES:",
            f"- Connected: {network_props.get('is_connected', 'Unknown')}",
            f"- Connected components: {network_props.get('connected_components', 'Unknown')}",
            f"- Largest component size: {network_props.get('largest_cc_size', 'Unknown')}",
            ""
        ])
    
    if not hub_genes.empty:
        report_lines.extend([
            "TOP HUB GENES:",
            *[f"- {row['gene_symbol']} (score: {row['hub_score']:.3f}, degree: {row['degree']})" 
              for _, row in hub_genes.head(10).iterrows()],
            ""
        ])
    
    if communities_louvain:
        report_lines.extend([
            "COMMUNITY DETECTION:",
            f"- Louvain communities: {len(communities_louvain)}",
            f"- Average community size: {np.mean([len(c) for c in communities_louvain]):.1f}",
            ""
        ])
    
    report_lines.extend([
        "FILES GENERATED:",
        "VISUALIZATIONS:",
        "- figures/26_coexpression_network_visualization.*",
        "- figures/27_correlation_heatmap.*",
        "- figures/28_*_pathway_network.* (if pathway data available)",
        "",
        "DATA FILES:",
        "- data/results/hub_genes.csv",
        "- data/results/coexpression_network_nodes.csv",
        "- data/results/coexpression_network_edges.csv", 
        "- data/results/coexpression_network.graphml",
        "",
        "BIOLOGICAL INTERPRETATION:",
        "The network analysis reveals:",
        "- Hub genes that may play central roles in lung cancer",
        "- Gene modules (communities) with similar expression patterns",
        "- Co-expression relationships that suggest functional interactions",
        "- Potential therapeutic targets based on network centrality",
        "",
        "NEXT STEPS:",
        "1. Import GraphML file into Cytoscape for advanced visualization",
        "2. Validate hub genes through literature review",
        "3. Perform experimental validation of key network relationships",
        "4. Integrate with protein-protein interaction data"
    ])
    
    return "\n".join(report_lines)

report_text = generate_network_report()

with open('reports/network_analysis_report.txt', 'w') as f:
    f.write(report_text)

print("Network analysis report saved to: reports/network_analysis_report.txt")

print("\n=== Gene Network Analysis Completed Successfully ===")
print("All analysis files have been generated and saved.")
print("\nFor advanced network visualization:")
print("1. Open Cytoscape")
print("2. Import data/results/coexpression_network.graphml") 
print("3. Use the node and edge attribute files for styling")
