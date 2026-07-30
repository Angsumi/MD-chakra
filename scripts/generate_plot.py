import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Read the CSV
df = pd.read_csv('Table 1.csv')

# Drop the Reference Genome column as it doesn't have per-sample metrics for most rows
df = df.drop(columns=['Reference Genome (All Genes)'])

# Clean up percentage signs
for col in ['NPFM1', 'NPFM2', 'NPFM3', 'NPFM4']:
    df[col] = df[col].astype(str).str.replace('%', '').astype(float)

# Setup the figure
fig, axes = plt.subplots(3, 3, figsize=(15, 15))
fig.suptitle('Statistics of RNA-seq data and Expressed Genes', fontsize=20, y=0.95)
axes = axes.flatten()

metrics = df['Metric'].values
samples = ['NPFM1', 'NPFM2', 'NPFM3', 'NPFM4']
colors = ['#FF9999', '#66B2FF', '#99FF99', '#FFCC99']

for i, metric in enumerate(metrics):
    ax = axes[i]
    values = df.iloc[i][samples].values
    
    # Create colorful bar plot
    bars = ax.bar(samples, values, color=colors, edgecolor='black')
    
    ax.set_title(metric, fontsize=14, pad=10)
    ax.set_ylabel('Value')
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Add values on top of bars
    for bar in bars:
        height = bar.get_height()
        # Format the text depending on the value magnitude
        if height > 1e9:
            text = f"{height/1e9:.1f}B"
        elif height > 1e6:
            text = f"{height/1e6:.1f}M"
        elif height > 1e3:
            text = f"{height/1e3:.1f}K"
        else:
            text = f"{height:.1f}"
            
        ax.text(bar.get_x() + bar.get_width()/2., height,
                text,
                ha='center', va='bottom', fontsize=10, rotation=0)

plt.tight_layout(rect=[0, 0.03, 1, 0.93])
plt.savefig('table 1 plot.png', dpi=300, bbox_inches='tight')
plt.savefig('table 1 plot.pdf', bbox_inches='tight')
print("Plots generated successfully!")
