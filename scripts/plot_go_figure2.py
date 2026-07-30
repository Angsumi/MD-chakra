import sys
import pandas as pd
import mygene
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import requests
import os

# 1. Get the latest go-basic.obo
if not os.path.exists("go-basic.obo"):
    print("Downloading go-basic.obo...")
    r = requests.get("http://purl.obolibrary.org/obo/go/go-basic.obo")
    with open("go-basic.obo", "w") as f:
        f.write(r.text)

from goatools.obo_parser import GODag
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

# 2. Extract genes
print("Loading DESeq2 results...")
res = pd.read_csv("DESeq2_results.csv", index_col=0)

print("Parsing GTF...")
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

print(f"Total Unique Genes: {len(all_genes)}, True DEGs: {len(deg_genes)}")

# 3. Query MyGene.info
mg = mygene.MyGeneInfo()
print("Querying MyGene.info for GO terms...")
q_all = mg.querymany(all_genes, scopes='symbol', fields='go', species='human', returnall=True)

# 4. Count
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

# 5. Plot
data = []
for (name, ns), count in counts_all.items():
    data.append({'Description': name, 'Ontology': ns, 'Count': count, 'Type': 'All Unigene'})
for (name, ns), count in counts_deg.items():
    data.append({'Description': name, 'Ontology': ns, 'Count': count, 'Type': 'DEG Unigene'})

df = pd.DataFrame(data)
if df.empty:
    print("No GO terms found.")
    sys.exit()

top_terms = df[df['Type'] == 'All Unigene'].groupby('Ontology').apply(lambda x: x.nlargest(15, 'Count')).reset_index(drop=True)['Description']
plot_df = df[df['Description'].isin(top_terms)]

plt.figure(figsize=(15, 7))
sns.set_theme(style="whitegrid")
g = sns.catplot(
    data=plot_df, kind="bar",
    x="Description", y="Count", hue="Type", col="Ontology",
    height=6, aspect=1.2, palette=["#4B8BBE", "#E06666"], sharex=False, sharey=False
)
g.set_titles("{col_name}", size=14)
g.set_axis_labels("", "Number of Genes", size=12)
g.set_xticklabels(rotation=45, ha='right', size=10)
plt.tight_layout()
plt.savefig("Figure2_GO_Classification.png", dpi=300, bbox_inches='tight')
print("Saved to Figure2_GO_Classification.png!")
