import sqlite3
import pandas as pd
import numpy as np
import os
import plotly.express as px
import plotly.graph_objects as go
import zipfile, xml.sax.saxutils

print("Starting Protein Domain Analysis & Visualization Generator...")

os.makedirs("visualizations", exist_ok=True)
os.makedirs("new plots/html", exist_ok=True)
os.makedirs("new plots/pdf", exist_ok=True)
os.makedirs("new plots/png", exist_ok=True)
os.makedirs("docs/assets", exist_ok=True)
os.makedirs("downstream_results", exist_ok=True)

# 1. Load Diamond hits & build mappings
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

# 2. Query Pfam domains and functional annotations from eggnog.db
db_path = "kegg_analysis/eggnog.db"
conn = sqlite3.connect(db_path)
cursor = conn.cursor()

seeds = hits_df['sseqid'].unique().tolist()
pfam_map = {}
pname_map = {}
for i in range(0, len(seeds), 500):
    chunk = seeds[i:i+500]
    placeholders = ','.join('?' for _ in chunk)
    cursor.execute(f'SELECT name, pname, pfam FROM prots WHERE name IN ({placeholders}) AND pfam IS NOT NULL AND pfam != ""', chunk)
    for r in cursor.fetchall():
        pname_map[r[0]] = r[1]
        pfam_map[r[0]] = r[2]

conn.close()

hits_df['Preferred_Name'] = hits_df['sseqid'].map(pname_map)
hits_df['Pfam_Domains'] = hits_df['sseqid'].map(pfam_map)

# Merge with FPKM normalized counts
fpkm_df = pd.read_csv('downstream_results/6_FPKM_normalized_counts_individual.csv')
fpkm_df.rename(columns={'Unnamed: 0': 'Gene ID'}, inplace=True)
merged_df = pd.merge(hits_df, fpkm_df, on='Gene ID', how='left')

# Drop duplicates
domain_master = merged_df[merged_df['Pfam_Domains'].notna()].drop_duplicates(subset=['Gene ID', 'protein_id']).copy()

# Save Master Table
final_cols = ['Gene ID', 'protein_id', 'Preferred_Name', 'Pfam_Domains', 'sseqid', 'pident', 'evalue', 'bitscore', 'NPFM1', 'NPFM2', 'NPFM3', 'NPFM4', 'Protein Description']
domain_master = domain_master[[c for c in final_cols if c in domain_master.columns]]

domain_master.to_csv('downstream_results/Protein_Domain_Annotation_Master.csv', index=False)
domain_master.to_csv('docs/assets/Protein_Domain_Annotation_Master.csv', index=False)
print(f"Saved Protein_Domain_Annotation_Master.csv with {len(domain_master)} rows.")

# Helper for XLSX
def csv_to_xlsx(df, xlsx_path, sheet_name='Pfam_Domains'):
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

csv_to_xlsx(domain_master, 'downstream_results/Protein_Domain_Annotation_Master.xlsx')
csv_to_xlsx(domain_master, 'docs/assets/Protein_Domain_Annotation_Master.xlsx')
print("Generated Protein_Domain_Annotation_Master.xlsx")

# -------------------------------------------------------------------------
# 3. Domain Categorization & Frequency Analysis
# -------------------------------------------------------------------------
# Define functional categories and representative domains
domain_category_map = {
    # Proteases
    'Peptidase_A17': 'Proteases & Peptidases',
    'Trypsin': 'Proteases & Peptidases',
    'Peptidase_C1': 'Proteases & Peptidases',
    'Peptidase_M10': 'Proteases & Peptidases',
    'Peptidase_M12A': 'Proteases & Peptidases',
    'Peptidase_C13': 'Proteases & Peptidases',
    'Peptidase_S9': 'Proteases & Peptidases',
    'Metallophos': 'Proteases & Peptidases',
    'Astacin': 'Proteases & Peptidases',
    'CLP_protease': 'Proteases & Peptidases',
    
    # Ion-channel / Neuroactive / Transport
    'Lig_chan': 'Ion-Channel & Transport',
    'Lig_chan-Glu_bd': 'Ion-Channel & Transport',
    'Ion_trans': 'Ion-Channel & Transport',
    '7tm_1': 'Ion-Channel & Transport',
    '7tm_2': 'Ion-Channel & Transport',
    'MFS_1': 'Ion-Channel & Transport',
    'ABC_tran': 'Ion-Channel & Transport',
    'Neur_chan_LBD': 'Ion-Channel & Transport',
    
    # Toxins & Venom-Associated
    'ShKT': 'Toxins & Venom Components',
    'Kunitz_BPTI': 'Toxins & Venom Components',
    'CBP': 'Toxins & Venom Components',
    'CAP': 'Toxins & Venom Components',
    'Lectin_C': 'Toxins & Venom Components',
    'Scorpion_toxin': 'Toxins & Venom Components',
    'Neurotoxin_1': 'Toxins & Venom Components',
    'Phospholip_A2_1': 'Toxins & Venom Components',
    'Kazal_1': 'Toxins & Venom Components',
    
    # Protein Binding & Interactions
    'fn3': 'Protein Binding & Recognition',
    'Ig_3': 'Protein Binding & Recognition',
    'I-set': 'Protein Binding & Recognition',
    'BTB': 'Protein Binding & Recognition',
    'WD40': 'Protein Binding & Recognition',
    'Ank_2': 'Protein Binding & Recognition',
    'Ank': 'Protein Binding & Recognition',
    'SH3_1': 'Protein Binding & Recognition',
    'PDZ': 'Protein Binding & Recognition',
    'TPR_1': 'Protein Binding & Recognition',
    'LRR_8': 'Protein Binding & Recognition',
    'RRM_1': 'Protein Binding & Recognition',
    
    # Enzymatic & Catalytic Activity
    'Pkinase': 'Enzymatic & Catalytic Activity',
    'Pkinase_Tyr': 'Enzymatic & Catalytic Activity',
    'p450': 'Enzymatic & Catalytic Activity',
    'Helicase_C': 'Enzymatic & Catalytic Activity',
    'DEAD': 'Enzymatic & Catalytic Activity',
    'PTPc_motif': 'Enzymatic & Catalytic Activity',
    'Methyltransf_25': 'Enzymatic & Catalytic Activity',
    'GST_N': 'Enzymatic & Catalytic Activity',
    'GST_C': 'Enzymatic & Catalytic Activity',
    'Aldedh': 'Enzymatic & Catalytic Activity',
    
    # Extracellular & Structural Matrix
    'Chitin_bind_4': 'Extracellular & Structural Matrix',
    'EGF': 'Extracellular & Structural Matrix',
    'EGF_CA': 'Extracellular & Structural Matrix',
    'TSP_1': 'Extracellular & Structural Matrix',
    'Collagen': 'Extracellular & Structural Matrix',
    'Zona_pellucida': 'Extracellular & Structural Matrix',
    'Laminin_G_1': 'Extracellular & Structural Matrix',
    'Fibrinogen_C': 'Extracellular & Structural Matrix'
}

domain_desc_map = {
    'Peptidase_A17': 'Peptidase A17 retrotransposon-like protease',
    'Trypsin': 'Trypsin-like serine protease',
    'Peptidase_C1': 'Papain family cysteine protease',
    'Peptidase_M10': 'Matrix metalloproteinase',
    'Astacin': 'Astacin family metallopeptidase',
    'Metallophos': 'Calcineurin-like phosphoesterase',
    'Lig_chan': 'Ligand-gated ion channel',
    'Lig_chan-Glu_bd': 'Ligand-gated ion channel glutamate-binding',
    'Ion_trans': 'Voltage-gated ion channel domain',
    '7tm_1': '7-transmembrane GPCR rhodopsin family',
    'MFS_1': 'Major facilitator superfamily transporter',
    'ShKT': 'ShK toxin domain (potassium channel blocker)',
    'Kunitz_BPTI': 'Kunitz/Bovine pancreatic trypsin inhibitor',
    'CBP': 'Calcium-binding / Lectin venom domain',
    'CAP': 'CRISP/Allergen/PR-1 venom protein',
    'Lectin_C': 'C-type lectin domain',
    'Phospholip_A2_1': 'Phospholipase A2 venom enzyme',
    'Kazal_1': 'Kazal-type serine protease inhibitor',
    'fn3': 'Fibronectin type III domain',
    'Ig_3': 'Immunoglobulin domain',
    'I-set': 'Immunoglobulin I-set domain',
    'BTB': 'BTB/POZ domain (protein interaction)',
    'WD40': 'WD40 repeat protein-binding scaffold',
    'Ank_2': 'Ankyrin repeat',
    'RRM_1': 'RNA recognition motif',
    'Pkinase': 'Protein kinase domain',
    'p450': 'Cytochrome P450 monooxygenase',
    'Helicase_C': 'Helicase conserved C-terminal domain',
    'DEAD': 'DEAD-box RNA/DNA helicase',
    'GST_N': 'Glutathione S-transferase N-terminal',
    'Chitin_bind_4': 'Chitin-binding cuticle/peritrophin domain',
    'EGF': 'Epidermal growth factor-like domain',
    'TSP_1': 'Thrombospondin type 1 repeat',
    'Collagen': 'Collagen triple helix repeat',
    'Zona_pellucida': 'Zona pellucida sperm-binding domain',
    'Fibrinogen_C': 'Fibrinogen C-terminal domain'
}

# Count domain occurrences
domain_counts = {}
for p in domain_master['Pfam_Domains']:
    for d in str(p).split(','):
        d = d.strip()
        if d:
            domain_counts[d] = domain_counts.get(d, 0) + 1

# Build summary dataframe
domain_summary = []
for d, cnt in domain_counts.items():
    cat = domain_category_map.get(d, 'Other Functional Domains')
    desc = domain_desc_map.get(d, f"{d} domain")
    domain_summary.append({'Domain_ID': d, 'Description': desc, 'Functional_Category': cat, 'Transcript_Count': cnt})

domain_sum_df = pd.DataFrame(domain_summary).sort_values(by='Transcript_Count', ascending=False)
domain_sum_df.to_csv('downstream_results/Protein_Domain_Classification_Stats.csv', index=False)
domain_sum_df.to_csv('docs/assets/Protein_Domain_Classification_Stats.csv', index=False)
csv_to_xlsx(domain_sum_df, 'downstream_results/Protein_Domain_Classification_Stats.xlsx', sheet_name='Domain_Stats')
csv_to_xlsx(domain_sum_df, 'docs/assets/Protein_Domain_Classification_Stats.xlsx', sheet_name='Domain_Stats')
print("Saved Protein_Domain_Classification_Stats.xlsx")

# -------------------------------------------------------------------------
# 4. Generate Publication-Quality Protein Domain Visualization
# -------------------------------------------------------------------------
# Select top representative domains across the 6 key categories
selected_cats = [
    'Proteases & Peptidases',
    'Ion-Channel & Transport',
    'Toxins & Venom Components',
    'Protein Binding & Recognition',
    'Enzymatic & Catalytic Activity',
    'Extracellular & Structural Matrix'
]

cat_colors = {
    'Proteases & Peptidases': '#ef4444',        # Coral Red
    'Ion-Channel & Transport': '#06b6d4',      # Cyan
    'Toxins & Venom Components': '#f59e0b',    # Amber Gold
    'Protein Binding & Recognition': '#6366f1',# Electric Indigo
    'Enzymatic & Catalytic Activity': '#10b981',# Emerald
    'Extracellular & Structural Matrix': '#8b5cf6' # Violet
}

# Pick top 4-5 domains per category
plot_domains = []
for c in selected_cats:
    sub = domain_sum_df[domain_sum_df['Functional_Category'] == c].head(4)
    plot_domains.append(sub)

plot_df = pd.concat(plot_domains).sort_values(by='Transcript_Count', ascending=True)

fig_domain = go.Figure()

for cat in selected_cats:
    sub_df = plot_df[plot_df['Functional_Category'] == cat]
    if len(sub_df) > 0:
        fig_domain.add_trace(go.Bar(
            y=sub_df['Domain_ID'] + " (" + sub_df['Description'] + ")",
            x=sub_df['Transcript_Count'],
            name=cat,
            orientation='h',
            marker=dict(color=cat_colors[cat], line=dict(color='#0f172a', width=1)),
            text=sub_df['Transcript_Count'],
            textposition='outside'
        ))

fig_domain.update_layout(
    title=dict(
        text="<b>Top Protein Domains in the <i>Nephila pilipes</i> Control Transcriptome</b><br><sup>Categorized by Functional Classes: Toxins, Proteases, Ion Channels, Binding, Enzymes & Extracellular Matrix</sup>",
        font=dict(size=18, color="#0f172a", family="Plus Jakarta Sans, sans-serif"),
        x=0.5
    ),
    xaxis=dict(title="<b>Number of Transcripts with Identified Domain</b>", showline=True, linecolor="#cbd5e1", linewidth=1.5, showgrid=True, gridcolor="#f1f5f9"),
    yaxis=dict(title="<b>Pfam Protein Domain & Functional Description</b>", showline=True, linecolor="#cbd5e1", linewidth=1.5),
    paper_bgcolor="white",
    plot_bgcolor="white",
    width=1100,
    height=850,
    legend=dict(
        title="<b>Functional Category</b>",
        orientation="h",
        yanchor="bottom",
        y=1.02,
        xanchor="center",
        x=0.5,
        bgcolor="rgba(255,255,255,0.9)",
        bordercolor="#cbd5e1",
        borderwidth=1
    ),
    margin=dict(l=380, r=60, t=140, b=60)
)

fig_domain.write_html("new plots/html/plotly_protein_domains.html")
fig_domain.write_image("visualizations/17_Protein_Domain_Classification.png", scale=3)
fig_domain.write_image("visualizations/17_Protein_Domain_Classification.pdf")
fig_domain.write_image("docs/assets/17_Protein_Domain_Classification.png", scale=3)
fig_domain.write_image("docs/assets/17_Protein_Domain_Classification.pdf")
print("Saved Protein Domain Classification Plot (HTML, PNG, PDF)")

print("Protein Domain Analysis and Visualizations completed successfully!")
