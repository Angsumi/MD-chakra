import pandas as pd
import re

gtf_file = 'reference_genome/annotation.gtf'

data = []
with open(gtf_file, 'r') as f:
    for line in f:
        if line.startswith('#'): continue
        parts = line.strip().split('\t')
        if len(parts) < 9: continue
        
        chrom = parts[0]
        feature_type = parts[2]
        
        info = parts[8]
        
        # Extract attributes
        gene_id_match = re.search(r'gene_id "([^"]+)"', info)
        gene_name_match = re.search(r'gene "([^"]+)"', info)
        product_match = re.search(r'product "([^"]+)"', info)
        
        gene_id = gene_id_match.group(1) if gene_id_match else None
        gene_name = gene_name_match.group(1) if gene_name_match else None
        description = product_match.group(1) if product_match else None
        
        if gene_id:
            data.append({
                'Gene ID': gene_id,
                'Chromosome': chrom,
                'Gene Name': gene_name,
                'Description': description,
                'feature_type': feature_type
            })

df = pd.DataFrame(data)

# We might have multiple entries per gene_id (e.g., gene, mRNA, CDS)
# Let's aggregate to get the best gene name and description
def agg_func(x):
    res = {}
    res['Chromosome'] = x['Chromosome'].iloc[0]
    
    names = x['Gene Name'].dropna().unique()
    res['Gene Name'] = names[0] if len(names) > 0 else 'Unknown'
    
    desc = x['Description'].dropna().unique()
    res['Description'] = desc[0] if len(desc) > 0 else 'Unknown'
    
    return pd.Series(res)

gene_anno = df.groupby('Gene ID').apply(agg_func).reset_index()
gene_anno.to_csv('processed_data/gene_annotation.csv', index=False)
print("Saved gene annotation to processed_data/gene_annotation.csv")
