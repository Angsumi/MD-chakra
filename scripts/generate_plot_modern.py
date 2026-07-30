import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Set a modern, sleek light theme
plt.style.use('default')
sns.set_context("notebook", font_scale=1.1)

# Read the CSV
df = pd.read_csv('Table 1.csv')
df = df.drop(columns=['Reference Genome (All Genes)'])

for col in ['NPFM1', 'NPFM2', 'NPFM3', 'NPFM4']:
    df[col] = df[col].astype(str).str.replace('%', '').astype(float)

fig, axes = plt.subplots(3, 3, figsize=(16, 14), facecolor='white')
fig.patch.set_facecolor('white')
fig.suptitle('RNA-seq Data and Expressed Genes Statistics', 
             fontsize=24, y=0.98, color='black', fontweight='bold', fontfamily='sans-serif')

axes = axes.flatten()
metrics = df['Metric'].values
samples = ['NPFM1', 'NPFM2', 'NPFM3', 'NPFM4']

# Sleek neon-pastel color palette
colors = ['#00e1d9', '#6366f1', '#ec4899', '#f59e0b']

for i, metric in enumerate(metrics):
    ax = axes[i]
    ax.set_facecolor('white')
    values = df.iloc[i][samples].values
    
    # Sleek bars with no edges, slight transparency for modern feel
    bars = ax.bar(samples, values, color=colors, alpha=0.9, width=0.6)
    
    # Titles and labels
    ax.set_title(metric.upper(), fontsize=13, pad=15, color='black', fontweight='bold')
    
    # Remove grid lines completely
    ax.grid(False)
    ax.set_axisbelow(True)
    
    # Remove all spines for a clean look
    for spine in ax.spines.values():
        spine.set_visible(False)
        
    ax.tick_params(axis='both', colors='black', length=0, pad=8)
    
    # Add values on top of bars cleanly
    for bar in bars:
        height = bar.get_height()
        if height > 1e9:
            text = f"{height/1e9:.1f}B"
        elif height > 1e6:
            text = f"{height/1e6:.1f}M"
        elif height > 1e3:
            text = f"{height/1e3:.1f}K"
        else:
            text = f"{height:.1f}"
            
        ax.text(bar.get_x() + bar.get_width()/2., height + (height * 0.02),
                text,
                ha='center', va='bottom', fontsize=11, color='black', fontweight='medium')

plt.subplots_adjust(hspace=0.4, wspace=0.3)
plt.savefig('table 1 plot_modern.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig('table 1 plot_modern.pdf', bbox_inches='tight', facecolor='white')
print("Modern plots generated successfully!")
