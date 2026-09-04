import sqlite3
import pandas as pd
import numpy as np
import os
import plotly.express as px
import plotly.graph_objects as go
import zipfile, xml.sax.saxutils

print("Building Master KEGG Annotation Tables and Publication Visualizations...")

os.makedirs("visualizations", exist_ok=True)
os.makedirs("new plots/html", exist_ok=True)
os.makedirs("new plots/pdf", exist_ok=True)
os.makedirs("new plots/png", exist_ok=True)
os.makedirs("docs/assets", exist_ok=True)
os.makedirs("downstream_results", exist_ok=True)

# 1. Load Diamond hits
hits_file = "kegg_analysis/diamond_hits.seed_orthologs"
cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
hits_df = pd.read_csv(hits_file, sep='\t', names=cols)
hits_df['protein_id'] = hits_df['qseqid'].str.split().str[0]

# Build protein_id -> Gene ID map from proteins.fasta
prot_to_gene = {}
prot_to_desc = {}
with open("kegg_analysis/proteins.fasta") as f:
    for line in f:
        if line.startswith(">"):
            parts = line[1:].strip().split()
            pid = parts[0]
            npil_matches = [p for p in parts if 'NPIL_' in p]
            if npil_matches:
                prot_to_gene[pid] = npil_matches[0]
            prot_to_desc[pid] = " ".join(parts[1:])

hits_df['Gene ID'] = hits_df['protein_id'].map(prot_to_gene)
hits_df['Protein Description'] = hits_df['protein_id'].map(prot_to_desc)

# 2. Query eggnog.db prots table for exact KEGG KO, Pathway, Enzyme EC, Modules, COG, GO
db_path = "kegg_analysis/eggnog.db"
conn = sqlite3.connect(db_path)
cursor = conn.cursor()

unique_seeds = hits_df['sseqid'].unique().tolist()
print(f"Total unique target proteins to map: {len(unique_seeds)}")

# Query in chunks of 500
chunk_size = 500
kegg_rows = []
for i in range(0, len(unique_seeds), chunk_size):
    chunk = unique_seeds[i:i+chunk_size]
    placeholders = ','.join('?' for _ in chunk)
    query = f"SELECT name, pname, kegg_ko, kegg_pathway, kegg_module, kegg_ec, kegg_brite, gos, ogs FROM prots WHERE name IN ({placeholders})"
    cursor.execute(query, chunk)
    kegg_rows.extend(cursor.fetchall())

conn.close()

kegg_df = pd.DataFrame(kegg_rows, columns=['sseqid', 'Preferred_Name', 'KEGG_KO', 'KEGG_Pathway', 'KEGG_Module', 'KEGG_EC', 'KEGG_BRITE', 'GO_Terms', 'eggNOG_OGs'])
print(f"Retrieved {len(kegg_df)} functional annotation records from SQL.")

# Merge with Diamond hits
merged_df = pd.merge(hits_df, kegg_df, on='sseqid', how='left')

# Load FPKM expression data to merge expression values
fpkm_df = pd.read_csv('downstream_results/6_FPKM_normalized_counts_individual.csv')
fpkm_df.rename(columns={'Unnamed: 0': 'Gene ID'}, inplace=True)
merged_df = pd.merge(merged_df, fpkm_df, on='Gene ID', how='left')

# Format Master KEGG Table
final_cols = [
    'Gene ID', 'protein_id', 'Preferred_Name', 'KEGG_KO', 'KEGG_Pathway', 
    'KEGG_Module', 'KEGG_EC', 'KEGG_BRITE', 'sseqid', 'pident', 'evalue', 'bitscore',
    'NPFM1', 'NPFM2', 'NPFM3', 'NPFM4', 'Protein Description', 'GO_Terms', 'eggNOG_OGs'
]
# Select available columns
avail_cols = [c for c in final_cols if c in merged_df.columns]
kegg_master = merged_df[avail_cols].drop_duplicates(subset=['Gene ID', 'protein_id']).copy()

# Save master CSV
kegg_master.to_csv('downstream_results/KEGG_Pathway_Annotation_Master.csv', index=False)
kegg_master.to_csv('docs/assets/KEGG_Pathway_Annotation_Master.csv', index=False)
print("Saved KEGG_Pathway_Annotation_Master.csv with shape:", kegg_master.shape)

# Helper for XLSX
def csv_to_xlsx(df, xlsx_path, sheet_name='KEGG_Annotation'):
    content_types = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
    <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
    <Default Extension="xml" ContentType="application/xml"/>
    <Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>
    <Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>
    <Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>
</Types>'''
    root_rels = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
    <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
</Relationships>'''
    wb = f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
    <sheets><sheet name="{sheet_name}" sheetId="1" r:id="rId1"/></sheets>
</workbook>'''
    wb_rels = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
    <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>
    <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>
</Relationships>'''
    styles = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
    <fonts count="1"><font><name val="Calibri"/><sz val="11"/></font></fonts>
    <fills count="2"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill></fills>
    <borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>
    <cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>
    <cellXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/></cellXfs>
</styleSheet>'''
    sheet_rows = []
    cols = df.columns.tolist()
    header_cells = [f'<c t="inlineStr"><is><t>{xml.sax.saxutils.escape(str(col))}</t></is></c>' for col in cols]
    sheet_rows.append('<row r="1">' + ''.join(header_cells) + '</row>')
    for r_idx, row in df.iterrows():
        row_num = r_idx + 2
        cells = []
        for val in row:
            if pd.isna(val): cells.append('<c/>')
            elif isinstance(val, (int, float, np.number)): cells.append(f'<c><v>{val}</v></c>')
            else: cells.append(f'<c t="inlineStr"><is><t>{xml.sax.saxutils.escape(str(val))}</t></is></c>')
        sheet_rows.append(f'<row r="{row_num}">' + ''.join(cells) + '</row>')
    sheet_xml = f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
    <sheetData>{''.join(sheet_rows)}</sheetData>
</worksheet>'''
    with zipfile.ZipFile(xlsx_path, 'w', compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr('[Content_Types].xml', content_types)
        z.writestr('_rels/.rels', root_rels)
        z.writestr('xl/workbook.xml', wb)
        z.writestr('xl/_rels/workbook.xml.rels', wb_rels)
        z.writestr('xl/styles.xml', styles)
        z.writestr('xl/worksheets/sheet1.xml', sheet_xml)

csv_to_xlsx(kegg_master, 'downstream_results/KEGG_Pathway_Annotation_Master.xlsx')
csv_to_xlsx(kegg_master, 'docs/assets/KEGG_Pathway_Annotation_Master.xlsx')
print("Generated KEGG_Pathway_Annotation_Master.xlsx")

# -------------------------------------------------------------------------
# 3. Parse KEGG Pathways & Compute Pathway Counts & Statistics
# -------------------------------------------------------------------------
pathway_counts = {}
for pways in kegg_master['KEGG_Pathway'].dropna():
    # Format can be comma-separated mapXXXXX or koXXXXX
    for pw in str(pways).split(','):
        pw = pw.strip()
        if pw and pw != '-':
            pathway_counts[pw] = pathway_counts.get(pw, 0) + 1

# Standard KEGG Pathway Names Mapping
kegg_name_map = {
    'ko01100': 'Metabolic pathways',
    'map01100': 'Metabolic pathways',
    'ko01110': 'Biosynthesis of secondary metabolites',
    'map01110': 'Biosynthesis of secondary metabolites',
    'ko04144': 'Endocytosis',
    'map04144': 'Endocytosis',
    'ko04142': 'Lysosome',
    'map04142': 'Lysosome',
    'ko04141': 'Protein processing in endoplasmic reticulum',
    'map04141': 'Protein processing in endoplasmic reticulum',
    'ko03010': 'Ribosome',
    'map03010': 'Ribosome',
    'ko03040': 'Spliceosome',
    'map03040': 'Spliceosome',
    'ko04120': 'Ubiquitin mediated proteolysis',
    'map04120': 'Ubiquitin mediated proteolysis',
    'ko04010': 'MAPK signaling pathway',
    'map04010': 'MAPK signaling pathway',
    'ko04064': 'NF-kappa B signaling pathway',
    'map04064': 'NF-kappa B signaling pathway',
    'ko04310': 'Wnt signaling pathway',
    'map04310': 'Wnt signaling pathway',
    'ko04151': 'PI3K-Akt signaling pathway',
    'map04151': 'PI3K-Akt signaling pathway',
    'ko00190': 'Oxidative phosphorylation',
    'map00190': 'Oxidative phosphorylation',
    'ko00010': 'Glycolysis / Gluconeogenesis',
    'map00010': 'Glycolysis / Gluconeogenesis',
    'ko00020': 'Citrate cycle (TCA cycle)',
    'map00020': 'Citrate cycle (TCA cycle)',
    'ko00230': 'Purine metabolism',
    'map00230': 'Purine metabolism',
    'ko00240': 'Pyrimidine metabolism',
    'map00240': 'Pyrimidine metabolism',
    'ko04810': 'Regulation of actin cytoskeleton',
    'map04810': 'Regulation of actin cytoskeleton',
    'ko04510': 'Focal adhesion',
    'map04510': 'Focal adhesion',
    'ko04210': 'Apoptosis',
    'map04210': 'Apoptosis',
    'ko04145': 'Phagosome',
    'map04145': 'Phagosome',
    'ko04710': 'Circadian rhythm',
    'map04710': 'Circadian rhythm',
    'ko04910': 'Insulin signaling pathway',
    'map04910': 'Insulin signaling pathway',
    'ko04070': 'Phosphatidylinositol signaling system',
    'map04070': 'Phosphatidylinositol signaling system',
    'ko04071': 'Sphingolipid signaling pathway',
    'map04071': 'Sphingolipid signaling pathway',
    'ko04020': 'Calcium signaling pathway',
    'map04020': 'Calcium signaling pathway',
    'ko04150': 'mTOR signaling pathway',
    'map04150': 'mTOR signaling pathway',
    'ko04068': 'FoxO signaling pathway',
    'map04068': 'FoxO signaling pathway',
    'ko04140': 'Autophagy - animal',
    'map04140': 'Autophagy - animal',
    'ko04390': 'Hippo signaling pathway',
    'map04390': 'Hippo signaling pathway',
    'ko04350': 'TGF-beta signaling pathway',
    'map04350': 'TGF-beta signaling pathway'
}

pw_rows = []
for pw, cnt in sorted(pathway_counts.items(), key=lambda x: x[1], reverse=True):
    pname = kegg_name_map.get(pw, f"Pathway ({pw})")
    pw_rows.append({'Pathway_ID': pw, 'Pathway_Name': pname, 'Gene_Count': cnt})

pw_df = pd.DataFrame(pw_rows)
pw_df.to_csv('downstream_results/KEGG_Pathway_Classification_Stats.csv', index=False)
pw_df.to_csv('docs/assets/KEGG_Pathway_Classification_Stats.csv', index=False)
csv_to_xlsx(pw_df, 'downstream_results/KEGG_Pathway_Classification_Stats.xlsx', sheet_name='KEGG_Stats')
csv_to_xlsx(pw_df, 'docs/assets/KEGG_Pathway_Classification_Stats.xlsx', sheet_name='KEGG_Stats')
print(f"Saved KEGG Pathway Stats with {len(pw_df)} pathways.")

# -------------------------------------------------------------------------
# 4. Generate Publication-Quality KEGG Visualizations
# -------------------------------------------------------------------------
# Top 20 Most Represented Specific Pathways (excluding generic metabolic pathways)
filtered_pw = pw_df[~pw_df['Pathway_Name'].isin(['Metabolic pathways', 'Biosynthesis of secondary metabolites'])].head(20)

fig_kegg_bar = go.Figure()
fig_kegg_bar.add_trace(go.Bar(
    y=filtered_pw['Pathway_Name'][::-1],
    x=filtered_pw['Gene_Count'][::-1],
    orientation='h',
    marker=dict(
        color=filtered_pw['Gene_Count'][::-1],
        colorscale='Viridis',
        line=dict(color='#0f172a', width=1)
    ),
    text=filtered_pw['Gene_Count'][::-1],
    textposition='outside'
))

fig_kegg_bar.update_layout(
    title=dict(
        text="<b>KEGG Functional Pathway Classification (Top 20 Pathways)</b><br><sup>Number of mapped Nephila pilipes unigenes per biological pathway</sup>",
        font=dict(size=18, color="#0f172a", family="Plus Jakarta Sans, sans-serif"),
        x=0.5
    ),
    xaxis=dict(title="<b>Number of Mapped Unigenes</b>", showline=True, linecolor="#cbd5e1", linewidth=1.5, showgrid=True, gridcolor="#f1f5f9"),
    yaxis=dict(title="<b>KEGG Pathway</b>", showline=True, linecolor="#cbd5e1", linewidth=1.5),
    paper_bgcolor="white",
    plot_bgcolor="white",
    width=1000,
    height=700,
    margin=dict(l=280, r=50, t=90, b=60)
)

fig_kegg_bar.write_html("new plots/html/plotly_kegg_pathways.html")
fig_kegg_bar.write_image("visualizations/15_KEGG_Pathway_Classification.png", scale=3)
fig_kegg_bar.write_image("visualizations/15_KEGG_Pathway_Classification.pdf")
fig_kegg_bar.write_image("docs/assets/15_KEGG_Pathway_Classification.png", scale=3)
fig_kegg_bar.write_image("docs/assets/15_KEGG_Pathway_Classification.pdf")
print("Saved KEGG Pathway Bar Chart (HTML, PNG, PDF)")

# Bubble Chart: Pathway representation vs Significance/Counts
fig_bubble = go.Figure()
fig_bubble.add_trace(go.Scatter(
    x=filtered_pw['Gene_Count'],
    y=filtered_pw['Pathway_Name'],
    mode='markers',
    marker=dict(
        size=filtered_pw['Gene_Count'] / 15 + 10,
        color=filtered_pw['Gene_Count'],
        colorscale='Plasma',
        showscale=True,
        colorbar=dict(title="Gene Count")
    )
))

fig_bubble.update_layout(
    title=dict(
        text="<b>KEGG Pathway Representation Dotplot</b><br><sup>Visualizing pathway complexity across spider transcriptome</sup>",
        font=dict(size=18, color="#0f172a", family="Plus Jakarta Sans, sans-serif"),
        x=0.5
    ),
    xaxis=dict(title="<b>Gene Count</b>", showline=True, linecolor="#cbd5e1", linewidth=1.5, showgrid=True, gridcolor="#f1f5f9"),
    yaxis=dict(title="<b>KEGG Pathway</b>", showline=True, linecolor="#cbd5e1", linewidth=1.5),
    paper_bgcolor="white",
    plot_bgcolor="white",
    width=1000,
    height=700,
    margin=dict(l=280, r=50, t=90, b=60)
)

fig_bubble.write_html("new plots/html/plotly_kegg_dotplot.html")
fig_bubble.write_image("visualizations/16_KEGG_Pathway_Dotplot.png", scale=3)
fig_bubble.write_image("visualizations/16_KEGG_Pathway_Dotplot.pdf")
fig_bubble.write_image("docs/assets/16_KEGG_Pathway_Dotplot.png", scale=3)
fig_bubble.write_image("docs/assets/16_KEGG_Pathway_Dotplot.pdf")
print("Saved KEGG Pathway Dotplot (HTML, PNG, PDF)")

print("KEGG Master Deliverables and Visualizations built successfully!")
