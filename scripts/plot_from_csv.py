import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches

pivot_df = pd.read_csv("Figure2_GO_Classification_Data.csv")

# Calculate totals for percentage normalization
total_all = pivot_df['All Unigene'].sum()
total_deg = pivot_df['DEG Unigene'].sum()

# Convert counts to percentages
pivot_df['All Unigene %'] = (pivot_df['All Unigene'] / total_all) * 100
pivot_df['DEG Unigene %'] = (pivot_df['DEG Unigene'] / total_deg) * 100

fig, ax = plt.subplots(figsize=(18, 8))

# Define colors standard in papers (usually a nice deep blue and strong red/orange)
color_all = "#4575b4"
color_deg = "#d73027"

x = np.arange(len(pivot_df))
width = 0.35

rects1 = ax.bar(x - width/2, pivot_df['All Unigene %'], width, label='All Unigene', color=color_all, edgecolor='black', linewidth=0.5)
rects2 = ax.bar(x + width/2, pivot_df['DEG Unigene %'], width, label='DEG Unigene', color=color_deg, edgecolor='black', linewidth=0.5)

ax.set_ylabel('Percent of Genes (%)', fontsize=14, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(pivot_df['Description'], rotation=45, ha='right', fontsize=11)

# Remove top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Add Legend
ax.legend(loc='upper right', frameon=False, fontsize=12)

# Add colored spans/bars for Ontologies at the top
y_max = ax.get_ylim()[1]
current_x = -0.5
colors = {'Biological Process': '#8dd3c7', 'Cellular Component': '#ffffb3', 'Molecular Function': '#bebada'}

for ont in ['Biological Process', 'Cellular Component', 'Molecular Function']:
    subset = pivot_df[pivot_df['Ontology'] == ont]
    if len(subset) == 0: continue
    end_x = current_x + len(subset)
    
    # Draw a line/rectangle at the top for grouping
    rect = patches.Rectangle((current_x + 0.1, y_max * 0.95), (end_x - current_x - 0.2), y_max * 0.03, 
                             linewidth=1, edgecolor='black', facecolor=colors[ont], clip_on=False)
    ax.add_patch(rect)
    
    # Add text
    ax.text(current_x + len(subset)/2, y_max * 1.01, ont, ha='center', va='bottom', fontsize=13, fontweight='bold')
    
    # Add a vertical dashed line separator
    if current_x > -0.5:
        ax.axvline(x=current_x, color='gray', linestyle='--', linewidth=1, alpha=0.5)
        
    current_x = end_x

ax.set_ylim(0, y_max * 1.1)

plt.tight_layout()
plt.savefig("Figure2_GO_Classification_PaperStyle.png", dpi=300, bbox_inches='tight')
print("Saved to Figure2_GO_Classification_PaperStyle.png!")
