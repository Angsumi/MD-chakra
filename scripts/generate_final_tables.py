import pandas as pd
import numpy as np
import os

outdir = "downstream_results"
counts_file = "counts/gene_counts.txt"
anno_file = "processed_data/gene_annotation.csv"

# 1. Load absolute counts
counts_df = pd.read_csv(counts_file, sep='\t', comment='#')
# Rename columns
counts_df.columns = [col.replace('alignments/', '').replace('_sorted.bam', '') for col in counts_df.columns]
counts_df.rename(columns={'Geneid': 'Gene ID'}, inplace=True)

# 2. Load Gene Annotations
anno_df = pd.read_csv(anno_file) if os.path.exists(anno_file) else None

# Merge Annotations with Counts
if anno_df is not None:
    merged = pd.merge(counts_df, anno_df[['Gene ID', 'Gene Name', 'Description']], on='Gene ID', how='left')
else:
    merged = counts_df.copy()
    merged['Gene Name'] = 'Unknown'
    merged['Description'] = 'Unknown'

# Fill NaNs
merged['Gene Name'] = merged['Gene Name'].fillna('Unknown')
merged['Description'] = merged['Description'].fillna('Unknown')

# Define Samples
samples = ['NPFM1', 'NPFM2', 'NPFM3', 'NPFM4']

# Filter only those that exist
samples = [s for s in samples if s in merged.columns]

# Calculate Avg Expression and Std Dev
merged['Avg Expression'] = merged[samples].mean(axis=1)
merged['Standard Deviation'] = merged[samples].std(axis=1)

# Calculate Z-score for heatmap input (though normally calculated on normalized counts)
# We will use normalized FPKM for Z-score if available, otherwise just counts.
fpkm_file = f"{outdir}/6_FPKM_normalized_counts_individual.csv"
if os.path.exists(fpkm_file):
    fpkm = pd.read_csv(fpkm_file)
    fpkm.rename(columns={'Unnamed: 0': 'Gene ID'}, inplace=True)
    
    # Calculate Z-score on FPKM
    fpkm_mean = fpkm[samples].mean(axis=1)
    fpkm_std = fpkm[samples].std(axis=1)
    # Avoid div by zero
    fpkm_std[fpkm_std == 0] = 1e-6
    for s in samples:
        fpkm[f'Z-score_{s}'] = (fpkm[s] - fpkm_mean) / fpkm_std
        merged[f'Z-score_{s}'] = fpkm[f'Z-score_{s}']
else:
    # Use raw counts if FPKM is missing
    mean = merged['Avg Expression']
    std = merged['Standard Deviation']
    std_adj = std.replace(0, 1e-6)
    for s in samples:
        merged[f'Z-score_{s}'] = (merged[s] - mean) / std_adj

# Organize columns for Item 9
# Gene ID | Chromosome | Gene Name | Description | S1 | S2 | S3 | S4 | Avg Expression
cols = ['Gene ID', 'Chr', 'Gene Name', 'Description'] + samples + ['Avg Expression', 'Standard Deviation'] + [f'Z-score_{s}' for s in samples]
item9_table = merged[cols]

item9_table.to_excel(f"{outdir}/9_Annotated_Expression_Table_Individual.xlsx", index=False)
item9_table.to_csv(f"{outdir}/9_Annotated_Expression_Table_Individual.csv", index=False)

# Pairwise fold changes based on normalized counts
if os.path.exists(fpkm_file):
    pairwise_df = fpkm[['Gene ID', 'NPFM1', 'NPFM2', 'NPFM3', 'NPFM4']].copy()
else:
    pairwise_df = merged[['Gene ID', 'NPFM1', 'NPFM2', 'NPFM3', 'NPFM4']].copy()

import itertools
pairs = list(itertools.combinations(samples, 2))
for s1, s2 in pairs:
    # Fold change = s2 / s1
    # Adding a pseudocount to avoid division by zero
    pairwise_df[f'FC_{s2}_vs_{s1}'] = (pairwise_df[s2] + 1) / (pairwise_df[s1] + 1)
    pairwise_df[f'log2FC_{s2}_vs_{s1}'] = np.log2(pairwise_df[f'FC_{s2}_vs_{s1}'])

# Merge with annotations
pairwise_final = pd.merge(pairwise_df, item9_table[['Gene ID', 'Chr', 'Gene Name', 'Description']], on='Gene ID', how='left')

pairwise_final.to_excel(f"{outdir}/10_Pairwise_Comparisons.xlsx", index=False)
pairwise_final.to_csv(f"{outdir}/10_Pairwise_Comparisons.csv", index=False)

print("Tables generation complete.")
