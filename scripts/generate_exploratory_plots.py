import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA

# 1. Load the data
df = pd.read_csv("vst_normalized_counts.csv", index_col=0)

# Set global seaborn style
sns.set_theme(style="white", context="talk")

# 2. PCA Plot
pca = PCA(n_components=2)
# Data must be samples as rows, genes as columns for PCA
pca_result = pca.fit_transform(df.T)

pca_df = pd.DataFrame({
    'PC1': pca_result[:, 0],
    'PC2': pca_result[:, 1],
    'Sample': df.columns
})

# Variance explained
var_exp = pca.explained_variance_ratio_ * 100

plt.figure(figsize=(10, 8))
# Modern distinct colors for the 4 samples
colors = ['#E63946', '#F4A261', '#2A9D8F', '#264653']

sns.scatterplot(
    x='PC1', y='PC2', 
    hue='Sample', 
    data=pca_df,
    palette=colors,
    s=400, # Large dots
    edgecolor='black',
    linewidth=1.5,
    alpha=0.9
)

plt.title('Principal Component Analysis (PCA)', pad=20, fontweight='bold')
plt.xlabel(f'PC1 ({var_exp[0]:.1f}% Variance)', fontweight='bold')
plt.ylabel(f'PC2 ({var_exp[1]:.1f}% Variance)', fontweight='bold')
plt.legend(title='', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
plt.grid(True, linestyle='--', alpha=0.3)
plt.tight_layout()
plt.savefig("7_PCA_plot_Modern.png", dpi=300, bbox_inches='tight')
plt.close()
print("Saved 7_PCA_plot_Modern.png")

# 3. Top 50 Variable Genes Heatmap
# Calculate variance across the 4 samples for each gene
gene_vars = df.var(axis=1)
# Get top 50
top50_genes = gene_vars.nlargest(50).index
df_top50 = df.loc[top50_genes]

# Normalize rows (Z-score) for heatmap visualization
df_top50_z = df_top50.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

# Plot modern clustered heatmap
# We don't want cluster labels to clutter the plot, just the heatmap
plt.figure(figsize=(10, 12))
g = sns.clustermap(
    df_top50_z,
    cmap="mako", # Modern cool colormap
    figsize=(8, 10),
    cbar_kws={'label': 'Z-Score (Expression)'},
    dendrogram_ratio=(0.1, 0.2),
    linewidths=0.5,
    linecolor='black',
    tree_kws={'linewidth': 1.5}
)
g.ax_heatmap.set_yticklabels([]) # Hide gene names as there are 50
g.ax_heatmap.set_ylabel("Top 50 Most Variable Genes", fontweight='bold')
g.ax_heatmap.set_xlabel("Samples", fontweight='bold')
g.fig.suptitle("Expression Heatmap (4 Samples)", y=1.02, fontweight='bold', fontsize=18)
plt.savefig("11_Heatmap_Top50_Variable_Modern.png", dpi=300, bbox_inches='tight')
plt.close()
print("Saved 11_Heatmap_Top50_Variable_Modern.png")

# 4. Sample Distance Matrix
from scipy.spatial.distance import pdist, squareform

dist_matrix = pdist(df.T, metric='euclidean')
sq_dist = squareform(dist_matrix)
dist_df = pd.DataFrame(sq_dist, index=df.columns, columns=df.columns)

plt.figure(figsize=(8, 6))
sns.heatmap(
    dist_df, 
    cmap="rocket_r", # Modern colormap
    annot=True, 
    fmt=".1f",
    linewidths=1,
    linecolor='black',
    cbar_kws={'label': 'Euclidean Distance'}
)
plt.title("Sample Distance Matrix", pad=20, fontweight='bold')
plt.tight_layout()
plt.savefig("8_Sample_Distance_Matrix_Modern.png", dpi=300, bbox_inches='tight')
plt.close()
print("Saved 8_Sample_Distance_Matrix_Modern.png")
