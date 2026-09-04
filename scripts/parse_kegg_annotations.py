import sqlite3
import pandas as pd
import numpy as np
import os
import plotly.express as px
import plotly.graph_objects as go
import zipfile, xml.sax.saxutils

print("Starting KEGG Annotation Parser & Visualization Generator...")

os.makedirs("visualizations", exist_ok=True)
os.makedirs("new plots/html", exist_ok=True)
os.makedirs("new plots/pdf", exist_ok=True)
os.makedirs("new plots/png", exist_ok=True)
os.makedirs("docs/assets", exist_ok=True)
os.makedirs("downstream_results", exist_ok=True)

# 1. Parse Diamond hits
hits_file = "kegg_analysis/diamond_hits.seed_orthologs"
if not os.path.exists(hits_file):
    print(f"Error: {hits_file} not found.")
    exit(1)

cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
hits_df = pd.read_csv(hits_file, sep='\t', names=cols)
print(f"Loaded {len(hits_df)} Diamond hits.")

# Extract protein ID from qseqid (e.g., GFS28097.1)
hits_df['protein_id'] = hits_df['qseqid'].str.split().str[0]

# Build map of protein_id -> Gene ID / Locus tag from proteins.fasta
prot_to_gene = {}
prot_to_desc = {}
with open("kegg_analysis/proteins.fasta") as f:
    for line in f:
        if line.startswith(">"):
            parts = line[1:].strip().split()
            pid = parts[0]
            # Search for NPIL_
            npil_matches = [p for p in parts if 'NPIL_' in p]
            if npil_matches:
                prot_to_gene[pid] = npil_matches[0]
            desc = " ".join(parts[1:])
            prot_to_desc[pid] = desc

hits_df['Gene ID'] = hits_df['protein_id'].map(prot_to_gene)
hits_df['Protein Description'] = hits_df['protein_id'].map(prot_to_desc)

# Connect to eggnog.db to extract KEGG KOs, pathways, and descriptions
db_path = "kegg_analysis/eggnog.db"
print("Querying eggnog.db SQL database...")
conn = sqlite3.connect(db_path)
cursor = conn.cursor()

# Query tables in eggnog.db
cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
tables = [row[0] for row in cursor.fetchall()]
print("Tables in eggnog.db:", tables)

# Query og / annotations table
# In eggnog.db, ortholog annotations contain KEGG_ko, KEGG_pathway, Description, Preferred_name
target_seeds = tuple(hits_df['sseqid'].unique())
print(f"Querying {len(target_seeds)} unique target seed orthologs in SQL...")

# Query in chunks to prevent SQL variable limit
chunk_size = 5000
annotations_list = []

# Check columns of relevant table
for tbl in ['og', 'proteins', 'annotations', 'members']:
    if tbl in tables:
        cursor.execute(f"PRAGMA table_info({tbl})")
        print(f"Columns in {tbl}:", [c[1] for c in cursor.fetchall()])

cursor.execute("SELECT name, sql FROM sqlite_master WHERE type='table'")
for tname, sql in cursor.fetchall():
    print(f"Schema for {tname}:\n{sql}\n")

conn.close()
print("SQL exploration complete.")
