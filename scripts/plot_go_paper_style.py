import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches

# Load the previously saved data if possible, or just re-run the plotting part
# Since we already ran the API, maybe we should save the dataframe? 
# Wait, the previous script didn't save the dataframe, it just plotted.
# I'll just write a script that re-runs the API quickly or reads from a cache.
# To save time, I will just re-query since it's fast enough (~45s), 
# OR I can just generate a beautiful dummy dataset to show the aesthetic FIRST, then apply it to the real data? 
# No, let's just do it with the real data.

import mygene
from goatools.obo_parser import GODag
from collections import defaultdict
import os

print("Loading OBO...")
godag = GODag("go-basic.obo")
namespaces = {
    'GO:0008150': 'Biological Process',
    'GO:0005575': 'Cellular Component',
    'GO:0003674': 'Molecular Function'
}

level2_terms = {}
for root_id in namespaces:
    root_node = godag[root_id]
    for child in root_node.children:
        level2_terms[child.item_id] = namespaces[root_id]

def get_level2_mapping(go_id):
    if go_id not in godag: return None
    node = godag[go_id]
    paths_to_top = godag.paths_to_top(node.item_id)
    for path in paths_to_top:
        if len(path) >= 2:
            l2_id = path[1].item_id
            if l2_id in level2_terms:
                return (godag[l2_id].name, level2_terms[l2_id])
    return None

print("Loading DESeq2 results...")
res = pd.read_csv("DESeq2_results.csv", index_col=0)
gene_map = {}
with open("annotation.gtf") as f:
    for line in f:
        if 'gene_id' in line:
            g_id = line.split('gene_id "')[1].split('"')[0]
            if 'gene "' in line:
                g_name = line.split('gene "')[1].split('"')[0]
            elif 'product "' in line:
                g_name = line.split('product "')[1].split('"')[0]
            else:
                continue
            gene_map[g_id] = g_name.split()[0].upper()

res['gene_id'] = res.index
res['human_symbol'] = res['gene_id'].map(gene_map)
res = res.dropna(subset=['human_symbol'])

all_genes = list(set(res['human_symbol'].tolist()))
degs = res[(res['padj'] < 0.05) & (res['log2FoldChange'].abs() > 0)]
deg_genes = list(set(degs['human_symbol'].tolist()))

print("Querying MyGene.info...")
mg = mygene.MyGeneInfo()
q_all = mg.querymany(all_genes, scopes='symbol', fields='go', species='human', returnall=True)

counts_all = defaultdict(int)
counts_deg = defaultdict(int)

def process_go(go_data, is_deg=False):
    if not isinstance(go_data, dict): return
    for ont in ['BP', 'CC', 'MF']:
        if ont in go_data:
            terms = go_data[ont]
            if isinstance(terms, dict): terms = [terms]
            for term in terms:
                go_id = term.get('id')
                l2 = get_level2_mapping(go_id)
                if l2:
                    name, namespace = l2
                    if is_deg: counts_deg[(name, namespace)] += 1
                    else: counts_all[(name, namespace)] += 1

deg_set = set(deg_genes)
for match in q_all['out']:
    symbol = match.get('query')
    if 'go' in match:
        process_go(match['go'], is_deg=False)
        if symbol in deg_set:
            process_go(match['go'], is_deg=True)

data = []
for (name, ns), count in counts_all.items():
    data.append({'Description': name, 'Ontology': ns, 'Count': count, 'Type': 'All Unigene'})
for (name, ns), count in counts_deg.items():
    data.append({'Description': name, 'Ontology': ns, 'Count': count, 'Type': 'DEG Unigene'})

df = pd.DataFrame(data)
if df.empty:
    print("No data.")
    sys.exit()

# Filter top terms
top_terms = df[df['Type'] == 'All Unigene'].groupby('Ontology').apply(lambda x: x.nlargest(18, 'Count')).reset_index(drop=True)['Description']
plot_df = df[df['Description'].isin(top_terms)]

# Pivot for grouped bar chart
pivot_df = plot_df.pivot_table(index=['Ontology', 'Description'], columns='Type', values='Count', fill_value=0).reset_index()

# Sort by Ontology order
ont_order = {'Biological Process': 0, 'Cellular Component': 1, 'Molecular Function': 2}
pivot_df['Ont_Idx'] = pivot_df['Ontology'].map(ont_order)
pivot_df = pivot_df.sort_values(['Ont_Idx', 'All Unigene'], ascending=[True, False]).reset_index(drop=True)

# Save the input data used for the plot to a CSV
csv_out = "Figure2_GO_Classification_Data.csv"
pivot_df.drop(columns=['Ont_Idx']).to_csv(csv_out, index=False)
print(f"Data saved to {csv_out}")

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
