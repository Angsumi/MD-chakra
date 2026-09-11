import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import os

print("Generating exact-style KEGG Pathway Analysis dotplot...")

os.makedirs("visualizations", exist_ok=True)
os.makedirs("new plots/pdf", exist_ok=True)
os.makedirs("new plots/png", exist_ok=True)
os.makedirs("docs/assets", exist_ok=True)

# Select high-confidence representative pathways matching the exact order & style
pathway_data = [
    {'pathway': 'Phospholipase D signaling pathway', 'count': 12, 'ratio': 0.0692, 'pvalue': 0.00035},
    {'pathway': 'Thyroid hormone signaling pathway', 'count': 11, 'ratio': 0.0635, 'pvalue': 0.00042},
    {'pathway': 'Oxytocin signaling pathway', 'count': 11, 'ratio': 0.0635, 'pvalue': 0.00078},
    {'pathway': 'Choline metabolism in cancer', 'count': 11, 'ratio': 0.0635, 'pvalue': 0.00045},
    {'pathway': 'Platelet activation', 'count': 10, 'ratio': 0.0577, 'pvalue': 0.00062},
    {'pathway': 'Cell cycle', 'count': 9, 'ratio': 0.0519, 'pvalue': 0.00165},
    {'pathway': 'Phosphatidylinositol signaling pathway', 'count': 8, 'ratio': 0.0461, 'pvalue': 0.00115},
    {'pathway': 'B cell receptor signaling pathway', 'count': 7, 'ratio': 0.0404, 'pvalue': 0.00105},
    {'pathway': 'Lysine degradation', 'count': 6, 'ratio': 0.0346, 'pvalue': 0.00095},
    {'pathway': 'African trypanosomiasis', 'count': 5, 'ratio': 0.0288, 'pvalue': 0.00098}
]

df = pd.DataFrame(pathway_data)

# Reverse so highest gene ratio is at the top (matching reference image)
df = df.iloc[::-1].reset_index(drop=True)

# -------------------------------------------------------------------------
# Matplotlib Exact Academic Styling (ggplot2 style)
# -------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8.8, 7.2), dpi=300)

# Background and Grid styling matching reference
ax.set_facecolor('#f0f0f0') # Classic ggplot2 light grey
fig.patch.set_facecolor('white')
ax.grid(True, color='white', linewidth=1.2, linestyle='-')
ax.set_axisbelow(True)

# Spines removal/styling
for spine in ax.spines.values():
    spine.set_visible(False)

# Custom Colormap: Bright Red -> Magenta/Pink -> Dark Purple/Navy (matching reference image gradient)
colors = ['#ff0011', '#e6005c', '#c71585', '#6b177a', '#2c105e']
cmap = LinearSegmentedColormap.from_list('custom_kegg', colors, N=256)

# Normalization for p-value colorbar
norm = matplotlib.colors.Normalize(vmin=0.0002, vmax=0.0018)

# Point sizes scaling (Selection counts)
sizes = (df['count'] ** 2) * 2.8

# Scatter plot
scatter = ax.scatter(
    df['ratio'], 
    range(len(df)), 
    s=sizes, 
    c=df['pvalue'], 
    cmap=cmap, 
    norm=norm, 
    edgecolors='none',
    zorder=5
)

# Y-ticks with pathway names
ax.set_yticks(range(len(df)))
ax.set_yticklabels(df['pathway'], fontsize=11, fontfamily='DejaVu Sans', color='#111111')
ax.tick_params(axis='y', which='both', length=0, pad=12)

# X-ticks & label
ax.set_xlabel('Gene ratio (selection counts/selection size)', fontsize=11.5, fontfamily='DejaVu Sans', labelpad=10, color='#222222')
ax.tick_params(axis='x', which='both', color='#999999', labelsize=10.5)
ax.set_xlim(0.026, 0.074)
ax.set_xticks([0.03, 0.04, 0.05, 0.06, 0.07])

# Title
ax.set_title('Pathway analysis', fontsize=13.5, fontfamily='DejaVu Sans', pad=14, color='#222222')

# -------------------------------------------------------------------------
# Right Side Legends (Continuous Colorbar for p-value + Discrete Size Circles)
# -------------------------------------------------------------------------
# Position the main plot to leave space on right for legends
plt.subplots_adjust(left=0.44, right=0.80, top=0.92, bottom=0.10)

# 1. Colorbar (p-value)
cax = fig.add_axes([0.83, 0.52, 0.038, 0.36])
cb = plt.colorbar(scatter, cax=cax, orientation='vertical')
cb.ax.set_title('p value', fontsize=11, fontstyle='italic', fontfamily='DejaVu Sans', pad=10, color='#111111', ha='center')
cb.ax.tick_params(labelsize=9.5, color='#444444')
cb.outline.set_visible(False)
cb.set_ticks([0.0005, 0.0010, 0.0015])
cb.set_ticklabels(['0.0005', '0.0010', '0.0015'])

# 2. Size Legend (Selection counts)
# Create legend circles matching reference (6, 8, 10, 12)
counts_legend = [6, 8, 10, 12]
fig.text(0.83, 0.44, 'Selection counts', fontsize=10.5, fontfamily='DejaVu Sans', color='#111111', ha='left')

legend_y_positions = [0.38, 0.30, 0.21, 0.11]
for cnt, y_pos in zip(counts_legend, legend_y_positions):
    sz = (cnt ** 2) * 2.8
    leg_ax = fig.add_axes([0.83, y_pos - 0.025, 0.08, 0.06])
    leg_ax.axis('off')
    leg_ax.scatter([0.25], [0.5], s=sz, c='#111111', edgecolors='none')
    leg_ax.text(0.65, 0.45, str(cnt), fontsize=9.5, fontfamily='DejaVu Sans', color='#111111', va='center')
    leg_ax.set_xlim(0, 1)
    leg_ax.set_ylim(0, 1)

# Save Outputs
output_png = "visualizations/16_KEGG_Pathway_Dotplot.png"
output_pdf = "visualizations/16_KEGG_Pathway_Dotplot.pdf"
docs_png = "docs/assets/16_KEGG_Pathway_Dotplot.png"
docs_pdf = "docs/assets/16_KEGG_Pathway_Dotplot.pdf"

plt.savefig(output_png, dpi=300, bbox_inches='tight')
plt.savefig(output_pdf, bbox_inches='tight')
plt.savefig(docs_png, dpi=300, bbox_inches='tight')
plt.savefig(docs_pdf, bbox_inches='tight')
plt.close()

print("Exact ggplot2-style KEGG Pathway Analysis Dotplot generated and saved successfully!")
