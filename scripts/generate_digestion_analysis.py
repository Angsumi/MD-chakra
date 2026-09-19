import os
import sys
import pandas as pd
import numpy as np
import zipfile
import xml.sax.saxutils
import matplotlib.pyplot as plt
import seaborn as sns

# Set publication style
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Ensure directories
os.makedirs("downstream_results", exist_ok=True)
os.makedirs("visualizations", exist_ok=True)
os.makedirs("docs/assets", exist_ok=True)
os.makedirs("new plots/pdf", exist_ok=True)
os.makedirs("new plots/png", exist_ok=True)

print("1. Loading Master Annotation & Expression Data...")
kegg_master = pd.read_csv('downstream_results/KEGG_Pathway_Annotation_Master.csv')
domain_master = pd.read_csv('downstream_results/Protein_Domain_Annotation_Master.csv')

merged = pd.merge(
    kegg_master, 
    domain_master[['Gene ID', 'protein_id', 'Pfam_Domains']], 
    on=['Gene ID', 'protein_id'], 
    how='outer'
)

if 'Pfam_Domains_x' in merged.columns and 'Pfam_Domains_y' in merged.columns:
    merged['Pfam_Domains'] = merged['Pfam_Domains_x'].combine_first(merged['Pfam_Domains_y'])
    merged.drop(columns=['Pfam_Domains_x', 'Pfam_Domains_y'], inplace=True, errors='ignore')

fpkm_cols = ['NPFM1', 'NPFM2', 'NPFM3', 'NPFM4']
for col in fpkm_cols:
    merged[col] = pd.to_numeric(merged[col], errors='coerce').fillna(0.0)

merged['Mean_FPKM'] = merged[fpkm_cols].mean(axis=1)
merged['Log2_Mean_FPKM'] = np.log2(merged['Mean_FPKM'] + 1.0)

# Define Digestive Enzyme Classification
def classify_digestion(row):
    desc = str(row.get('Protein Description', '')).lower()
    pref = str(row.get('Preferred_Name', '')).lower()
    pf = str(row.get('Pfam_Domains', '')).lower()
    kegg_path = str(row.get('KEGG_Pathway', '')).lower()
    ko = str(row.get('KEGG_KO', '')).lower()
    go = str(row.get('GO_Terms', '')).lower()
    
    # Exclude non-digestive mobile elements / retrotransposons / integrases
    if any(t in desc for t in ['integrase', 'reverse transcriptase', 'transposon', 'gag-pol', 'retrotransposon', 'polyprotein', 'retrovirus']):
        return None

    # --- 1. Protein Digestion & Proteolysis ---
    # Serine Proteases
    if 'trypsin' in desc or 'trypsin' in pf or 'peptidase_s1' in pf:
        return ('Protein Digestion & Proteolysis', 'Serine Proteases', 'Trypsin / Peptidase S1')
    if 'chymotrypsin' in desc:
        return ('Protein Digestion & Proteolysis', 'Serine Proteases', 'Chymotrypsin')
    if 'subtilisin' in desc or 'peptidase_s8' in pf:
        return ('Protein Digestion & Proteolysis', 'Serine Proteases', 'Subtilisin-like (Peptidase S8)')
    if 'peptidase_s9' in pf or 'prolyl oligopeptidase' in desc:
        return ('Protein Digestion & Proteolysis', 'Serine Proteases', 'Prolyl Oligopeptidase (Peptidase S9)')
    if 'peptidase_s10' in pf or 'carboxypeptidase y' in desc:
        return ('Protein Digestion & Proteolysis', 'Serine Proteases', 'Serine Carboxypeptidase (Peptidase S10)')
    if 'elastase' in desc:
        return ('Protein Digestion & Proteolysis', 'Serine Proteases', 'Elastase')

    # Cysteine Proteases
    if 'cathepsin b' in desc or 'cathepsin l' in desc or 'peptidase_c1' in pf or 'papain' in desc:
        return ('Protein Digestion & Proteolysis', 'Cysteine Proteases', 'Cathepsin B/L / Papain (C1A)')
    if 'cathepsin c' in desc or 'dipeptidyl-peptidase i' in desc:
        return ('Protein Digestion & Proteolysis', 'Cysteine Proteases', 'Cathepsin C / DPP-I')
    if 'legumain' in desc or 'peptidase_c13' in pf:
        return ('Protein Digestion & Proteolysis', 'Cysteine Proteases', 'Legumain / Asparaginyl Endopeptidase (C13)')
    if 'calpain' in desc or 'peptidase_c2' in pf:
        return ('Protein Digestion & Proteolysis', 'Cysteine Proteases', 'Calpain (Peptidase C2)')

    # Metalloproteases
    if 'astacin' in desc or 'astacin' in pf or 'peptidase_m12a' in pf:
        return ('Protein Digestion & Proteolysis', 'Metalloproteases', 'Astacin Metalloprotease (M12A)')
    if 'neprilysin' in desc or 'peptidase_m13' in pf:
        return ('Protein Digestion & Proteolysis', 'Metalloproteases', 'Neprilysin (Peptidase M13)')
    if 'peptidase_m14' in pf or 'carboxypeptidase a' in desc or 'carboxypeptidase b' in desc:
        return ('Protein Digestion & Proteolysis', 'Metalloproteases', 'Carboxypeptidase A/B (Peptidase M14)')
    if 'peptidase_m10' in pf or 'matrix metalloproteinase' in desc or 'mmp' in pref:
        return ('Protein Digestion & Proteolysis', 'Metalloproteases', 'Matrix Metalloproteinase (M10)')
    if 'aminopeptidase n' in desc or 'peptidase_m1' in pf:
        return ('Protein Digestion & Proteolysis', 'Metalloproteases', 'Aminopeptidase N (Peptidase M1)')
    if 'leucine aminopeptidase' in desc or 'peptidase_m17' in pf:
        return ('Protein Digestion & Proteolysis', 'Metalloproteases', 'Leucine Aminopeptidase (M17)')
    if 'dipeptidyl peptidase' in desc:
        return ('Protein Digestion & Proteolysis', 'Metalloproteases', 'Dipeptidyl Peptidase')

    # Aspartic Proteases
    if ('cathepsin d' in desc or 'cathepsin e' in desc or 'pepsin' in desc or 'peptidase_a1' in pf or 'gastricsin' in desc):
        return ('Protein Digestion & Proteolysis', 'Aspartic Proteases', 'Cathepsin D/E & Pepsin (Peptidase A1)')

    # --- 2. Carbohydrate & Chitin Digestion ---
    if 'chitinase' in desc or 'glyco_hydro_18' in pf:
        return ('Carbohydrate & Chitin Digestion', 'Chitinolytic Enzymes', 'Chitinase (GH18)')
    if 'hexosaminidase' in desc or 'glyco_hydro_20' in pf:
        return ('Carbohydrate & Chitin Digestion', 'Chitinolytic Enzymes', 'Beta-Hexosaminidase (GH20)')
    if 'amylase' in desc or 'glyco_hydro_13' in pf:
        return ('Carbohydrate & Chitin Digestion', 'Glycoside Hydrolases', 'Alpha-Amylase (GH13)')
    if 'glucosidase' in desc or 'glyco_hydro_1' in pf or 'glyco_hydro_31' in pf:
        return ('Carbohydrate & Chitin Digestion', 'Glycoside Hydrolases', 'Alpha/Beta-Glucosidase (GH1/GH31)')
    if 'trehalase' in desc or 'glyco_hydro_37' in pf:
        return ('Carbohydrate & Chitin Digestion', 'Glycoside Hydrolases', 'Trehalase (GH37)')
    if 'galactosidase' in desc or 'glyco_hydro_2' in pf or 'glyco_hydro_35' in pf:
        return ('Carbohydrate & Chitin Digestion', 'Glycoside Hydrolases', 'Beta-Galactosidase (GH2/GH35)')
    if 'mannosidase' in desc or 'glyco_hydro_38' in pf or 'glyco_hydro_47' in pf:
        return ('Carbohydrate & Chitin Digestion', 'Glycoside Hydrolases', 'Mannosidase (GH38/GH47)')

    # --- 3. Lipid Digestion & Hydrolysis ---
    if 'lipase_3' in pf or 'pancreatic lipase' in desc or 'gastric lipase' in desc or 'triacylglycerol lipase' in desc:
        return ('Lipid Digestion & Hydrolysis', 'Lipases & Esterases', 'Pancreatic/Triacylglycerol Lipase')
    if 'phospholipase a2' in desc or 'pla2' in desc or 'pa2b' in pf:
        return ('Lipid Digestion & Hydrolysis', 'Phospholipases', 'Phospholipase A2 (PLA2)')
    if 'phospholipase c' in desc or 'plc' in pf:
        return ('Lipid Digestion & Hydrolysis', 'Phospholipases', 'Phospholipase C (PLC)')
    if 'phospholipase d' in desc or 'pld' in pf:
        return ('Lipid Digestion & Hydrolysis', 'Phospholipases', 'Phospholipase D (PLD)')
    if 'carboxylesterase' in desc or 'abhydrolase' in pf:
        return ('Lipid Digestion & Hydrolysis', 'Lipases & Esterases', 'Carboxylesterase / Neutral Lipase')

    # --- 4. Nucleic Acid & Phosphate Digestion ---
    if 'deoxyribonuclease' in desc or 'dnase' in desc:
        return ('Nucleic Acid & Phosphate Digestion', 'Nucleases', 'Deoxyribonuclease (DNase)')
    if 'ribonuclease' in desc or 'rnase' in desc:
        return ('Nucleic Acid & Phosphate Digestion', 'Nucleases', 'Ribonuclease (RNase)')
    if 'alkaline phosphatase' in desc or 'alk_phosphatase' in pf:
        return ('Nucleic Acid & Phosphate Digestion', 'Phosphatases', 'Alkaline Phosphatase')
    if 'acid phosphatase' in desc or 'acid_phosphatase' in pf:
        return ('Nucleic Acid & Phosphate Digestion', 'Phosphatases', 'Acid Phosphatase')

    return None

classified_records = []
for idx, row in merged.iterrows():
    c = classify_digestion(row)
    if c is not None:
        r = row.to_dict()
        r['Functional_Process'] = c[0]
        r['Enzyme_Class'] = c[1]
        r['Enzyme_Family'] = c[2]
        r['Expressed_in_Transcriptome'] = 'Yes' if r['Mean_FPKM'] > 0.1 else 'No'
        classified_records.append(r)

digestion_df = pd.DataFrame(classified_records)

# Reorder columns logically
key_cols = [
    'Gene ID', 'protein_id', 'Preferred_Name', 'Functional_Process', 'Enzyme_Class', 'Enzyme_Family',
    'NPFM1', 'NPFM2', 'NPFM3', 'NPFM4', 'Mean_FPKM', 'Log2_Mean_FPKM', 'Expressed_in_Transcriptome',
    'Pfam_Domains', 'KEGG_KO', 'KEGG_Pathway', 'GO_Terms', 'Protein Description'
]
cols_to_use = [c for c in key_cols if c in digestion_df.columns]
digestion_df = digestion_df[cols_to_use].sort_values(by=['Functional_Process', 'Enzyme_Class', 'Mean_FPKM'], ascending=[True, True, False])

print(f"Total classified digestion enzymes: {len(digestion_df)}")

# Helper function to write multi-sheet Excel without external dependencies
def write_multisheet_xlsx(df_dict, xlsx_path):
    sheets_xml = []
    workbook_sheets = []
    
    for idx, (name, df) in enumerate(df_dict.items(), 1):
        s_id = f"rId{idx}"
        workbook_sheets.append(f'<sheet name="{xml.sax.saxutils.escape(name[:31])}" sheetId="{idx}" r:id="{s_id}"/>')
        
        # Build worksheet xml
        cols = list(df.columns)
        rows_xml = []
        
        # Header row
        header_cells = []
        for c_idx, c_name in enumerate(cols):
            c_letter = chr(65 + c_idx) if c_idx < 26 else chr(64 + c_idx // 26) + chr(65 + c_idx % 26)
            header_cells.append(f'<c r="{c_letter}1" t="inlineStr"><is><t>{xml.sax.saxutils.escape(str(c_name))}</t></is></c>')
        rows_xml.append(f'<row r="1">{"".join(header_cells)}</row>')
        
        # Data rows
        for r_idx, row in enumerate(df.itertuples(index=False), 2):
            cells = []
            for c_idx, val in enumerate(row):
                c_letter = chr(65 + c_idx) if c_idx < 26 else chr(64 + c_idx // 26) + chr(65 + c_idx % 26)
                if pd.isna(val):
                    continue
                elif isinstance(val, (int, np.integer)):
                    cells.append(f'<c r="{c_letter}{r_idx}"><v>{val}</v></c>')
                elif isinstance(val, (float, np.floating)):
                    cells.append(f'<c r="{c_letter}{r_idx}"><v>{val:.6f}</v></c>')
                else:
                    cells.append(f'<c r="{c_letter}{r_idx}" t="inlineStr"><is><t>{xml.sax.saxutils.escape(str(val))}</t></is></c>')
            rows_xml.append(f'<row r="{r_idx}">{"".join(cells)}</row>')
            
        ws_xml = f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
    <sheetData>
        {"".join(rows_xml)}
    </sheetData>
</worksheet>'''
        sheets_xml.append((f"xl/worksheets/sheet{idx}.xml", ws_xml))

    content_types = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
                     '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">',
                     '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>',
                     '<Default Extension="xml" ContentType="application/xml"/>',
                     '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>']
    for idx in range(1, len(df_dict) + 1):
        content_types.append(f'<Override PartName="/xl/worksheets/sheet{idx}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>')
    content_types.append('</Types>')

    rels = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">',
            '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>',
            '</Relationships>']

    wb_rels = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
               '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">']
    for idx in range(1, len(df_dict) + 1):
        wb_rels.append(f'<Relationship Id="rId{idx}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{idx}.xml"/>')
    wb_rels.append('</Relationships>')

    workbook_xml = f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
    <sheets>
        {"".join(workbook_sheets)}
    </sheets>
</workbook>'''

    with zipfile.ZipFile(xlsx_path, 'w', compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr('[Content_Types].xml', "".join(content_types))
        z.writestr('_rels/.rels', "".join(rels))
        z.writestr('xl/_rels/workbook.xml.rels', "".join(wb_rels))
        z.writestr('xl/workbook.xml', workbook_xml)
        for s_path, s_content in sheets_xml:
            z.writestr(s_path, s_content)

# 2. Generate Summary Matrices
family_stats = digestion_df.groupby(['Functional_Process', 'Enzyme_Class', 'Enzyme_Family']).agg(
    Transcript_Count=('protein_id', 'nunique'),
    Expressed_Transcripts=('Expressed_in_Transcriptome', lambda x: (x == 'Yes').sum()),
    Mean_FPKM=('Mean_FPKM', 'mean'),
    Total_FPKM=('Mean_FPKM', 'sum'),
    NPFM1_Sum=('NPFM1', 'sum'),
    NPFM2_Sum=('NPFM2', 'sum'),
    NPFM3_Sum=('NPFM3', 'sum'),
    NPFM4_Sum=('NPFM4', 'sum')
).reset_index().sort_values(by=['Functional_Process', 'Total_FPKM'], ascending=[True, False])

# Pivot matrices for Heatmaps
count_pivot = digestion_df.pivot_table(
    index='Enzyme_Family', 
    columns='Functional_Process', 
    values='protein_id', 
    aggfunc='count', 
    fill_value=0
)

abundance_pivot = digestion_df.pivot_table(
    index='Enzyme_Family', 
    columns='Functional_Process', 
    values='Mean_FPKM', 
    aggfunc='sum', 
    fill_value=0.0
)

sample_pivot = digestion_df.pivot_table(
    index='Enzyme_Family', 
    values=['NPFM1', 'NPFM2', 'NPFM3', 'NPFM4'], 
    aggfunc='sum', 
    fill_value=0.0
)[['NPFM1', 'NPFM2', 'NPFM3', 'NPFM4']]

# Save CSV and XLSX datasets
digestion_df.to_csv('downstream_results/01_Digestion_Enzyme_Master_Annotation.csv', index=False)
digestion_df.to_csv('docs/assets/01_Digestion_Enzyme_Master_Annotation.csv', index=False)

family_stats.to_csv('downstream_results/01_Digestion_Family_Summary_Stats.csv', index=False)
family_stats.to_csv('docs/assets/01_Digestion_Family_Summary_Stats.csv', index=False)

count_pivot.to_csv('downstream_results/01_Digestion_Family_by_Process_Matrix.csv')
count_pivot.to_csv('docs/assets/01_Digestion_Family_by_Process_Matrix.csv')

# Write comprehensive Multi-Sheet Excel
write_multisheet_xlsx({
    'Master_Catalogue': digestion_df,
    'Family_Summary_Stats': family_stats,
    'Family_x_Process_Counts': count_pivot.reset_index(),
    'Family_x_Process_Abundance': abundance_pivot.reset_index(),
    'Sample_Expression_Matrix': sample_pivot.reset_index()
}, 'downstream_results/01_Digestion_Enzyme_Master_Annotation.xlsx')

write_multisheet_xlsx({
    'Master_Catalogue': digestion_df,
    'Family_Summary_Stats': family_stats,
    'Family_x_Process_Counts': count_pivot.reset_index(),
    'Family_x_Process_Abundance': abundance_pivot.reset_index(),
    'Sample_Expression_Matrix': sample_pivot.reset_index()
}, 'docs/assets/01_Digestion_Enzyme_Master_Annotation.xlsx')

print("Saved CSV and Excel tables successfully!")

# -------------------------------------------------------------------------
# VISUALIZATION 1: Clustered Heatmap (Enzyme Families × Functional Processes)
# -------------------------------------------------------------------------
print("Generating Visualization 1: Functional Enzyme Clustered Heatmap...")

# Create dual-layer data: Log2(Abundance + 1) for color, formatted string count for annotation
log_abundance = np.log2(abundance_pivot + 1.0)

plt.figure(figsize=(12, 11))
g = sns.clustermap(
    log_abundance,
    cmap="YlGnBu",
    annot=count_pivot, # Show transcript count in each cell
    fmt="d",
    linewidths=1.2,
    linecolor="#E2E8F0",
    figsize=(12, 11),
    cbar_kws={'label': 'Log2(Cumulative FPKM + 1)'},
    tree_kws={'linewidths': 1.2},
    dendrogram_ratio=(0.18, 0.12),
    cbar_pos=(0.02, 0.82, 0.03, 0.15)
)

g.ax_heatmap.set_title("Nephila pilipes Digestive Enzyme Repertoire\nEnzyme Families × Functional Processes", fontsize=14, fontweight='bold', pad=25)
g.ax_heatmap.set_xlabel("Functional Digestive Processes", fontsize=12, fontweight='bold', labelpad=10)
g.ax_heatmap.set_ylabel("Enzyme Families & Classes", fontsize=12, fontweight='bold', labelpad=10)
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=30, ha='right', fontsize=11, fontweight='medium')
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=10)

plt.tight_layout()
fig1_pdf = "visualizations/01_Digestion_Functional_Enzyme_Heatmap.pdf"
fig1_png = "visualizations/01_Digestion_Functional_Enzyme_Heatmap.png"
g.savefig(fig1_pdf, bbox_inches='tight', dpi=300)
g.savefig(fig1_png, bbox_inches='tight', dpi=300)
g.savefig("docs/assets/01_Digestion_Functional_Enzyme_Heatmap.pdf", bbox_inches='tight', dpi=300)
g.savefig("docs/assets/01_Digestion_Functional_Enzyme_Heatmap.png", bbox_inches='tight', dpi=300)
g.savefig("new plots/pdf/01_Digestion_Functional_Enzyme_Heatmap.pdf", bbox_inches='tight', dpi=300)
g.savefig("new plots/png/01_Digestion_Functional_Enzyme_Heatmap.png", bbox_inches='tight', dpi=300)
plt.close()

# -------------------------------------------------------------------------
# VISUALIZATION 2: High-Resolution Clustered Heatmap Across Replicates (NPFM1-NPFM4)
# -------------------------------------------------------------------------
print("Generating Visualization 2: Top Expressed Digestive Enzymes Clustered Heatmap...")

top_genes = digestion_df[digestion_df['Mean_FPKM'] > 0.05].copy()
if len(top_genes) > 50:
    top_genes = top_genes.head(50)

expr_mat = top_genes[fpkm_cols]
log_expr_mat = np.log2(expr_mat + 1.0)
log_expr_mat.index = top_genes['protein_id'] + " (" + top_genes['Enzyme_Family'].apply(lambda x: x.split('/')[0].strip()[:18]) + ")"

# Class colors for row sidebar
unique_classes = top_genes['Enzyme_Class'].unique()
palette = sns.color_palette("tab10", len(unique_classes))
class_to_color = dict(zip(unique_classes, palette))
row_colors = top_genes['Enzyme_Class'].map(class_to_color)
row_colors.index = log_expr_mat.index

g2 = sns.clustermap(
    log_expr_mat,
    cmap="mako",
    row_colors=row_colors,
    linewidths=0.6,
    linecolor="#CBD5E1",
    figsize=(11, 14),
    cbar_kws={'label': 'Log2(FPKM + 1)'},
    tree_kws={'linewidths': 1.0},
    dendrogram_ratio=(0.15, 0.05),
    cbar_pos=(0.02, 0.85, 0.03, 0.12)
)

g2.ax_heatmap.set_title("Expression Profile of Candidate Digestive Enzymes\n(Individual Biological Replicates NPFM1–NPFM4)", fontsize=13, fontweight='bold', pad=20)
g2.ax_heatmap.set_xlabel("Replicate Samples", fontsize=11, fontweight='bold', labelpad=10)
plt.setp(g2.ax_heatmap.get_xticklabels(), rotation=0, fontsize=11, fontweight='bold')
plt.setp(g2.ax_heatmap.get_yticklabels(), rotation=0, fontsize=8.5)

legend_handles = [plt.Line2D([0], [0], marker='s', color='w', label=cls_name,
                             markerfacecolor=color, markersize=10)
                  for cls_name, color in class_to_color.items()]
g2.ax_col_dendrogram.legend(
    handles=legend_handles, 
    title="Enzyme Class", 
    loc="center", 
    bbox_to_anchor=(0.5, 1.8),
    ncol=min(len(legend_handles), 3),
    frameon=True,
    fontsize=9,
    title_fontsize=10
)

plt.tight_layout()
fig2_pdf = "visualizations/01b_Digestion_Gene_Expression_Heatmap.pdf"
fig2_png = "visualizations/01b_Digestion_Gene_Expression_Heatmap.png"
g2.savefig(fig2_pdf, bbox_inches='tight', dpi=300)
g2.savefig(fig2_png, bbox_inches='tight', dpi=300)
g2.savefig("docs/assets/01b_Digestion_Gene_Expression_Heatmap.pdf", bbox_inches='tight', dpi=300)
g2.savefig("docs/assets/01b_Digestion_Gene_Expression_Heatmap.png", bbox_inches='tight', dpi=300)
g2.savefig("new plots/pdf/01b_Digestion_Gene_Expression_Heatmap.pdf", bbox_inches='tight', dpi=300)
g2.savefig("new plots/png/01b_Digestion_Gene_Expression_Heatmap.png", bbox_inches='tight', dpi=300)
plt.close()

# -------------------------------------------------------------------------
# VISUALIZATION 3: Summary Bar Breakdown (Diversity & Cumulative Abundance)
# -------------------------------------------------------------------------
print("Generating Visualization 3: Enzyme Family Diversity & Abundance Profile...")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8), sharey=True)

plot_stats = family_stats.sort_values(by='Total_FPKM', ascending=True)

proc_colors = {
    'Protein Digestion & Proteolysis': '#3B82F6',
    'Lipid Digestion & Hydrolysis': '#10B981',
    'Carbohydrate & Chitin Digestion': '#F59E0B',
    'Nucleic Acid & Phosphate Digestion': '#8B5CF6'
}
bar_colors = [proc_colors.get(p, '#6B7280') for p in plot_stats['Functional_Process']]

bars1 = ax1.barh(plot_stats['Enzyme_Family'], plot_stats['Transcript_Count'], color=bar_colors, alpha=0.85, edgecolor='black', linewidth=0.5)
ax1.set_xlabel("Number of Identified Transcripts / Genes", fontsize=11, fontweight='bold')
ax1.set_title("Transcript Repertoire Diversity", fontsize=12, fontweight='bold', pad=10)
ax1.grid(axis='x', linestyle='--', alpha=0.5)
for bar in bars1:
    width = bar.get_width()
    ax1.text(width + max(plot_stats['Transcript_Count']) * 0.01, bar.get_y() + bar.get_height()/2,
             f"{int(width)}", ha='left', va='center', fontsize=9, fontweight='bold')

bars2 = ax2.barh(plot_stats['Enzyme_Family'], plot_stats['Total_FPKM'], color=bar_colors, alpha=0.85, edgecolor='black', linewidth=0.5)
ax2.set_xlabel("Cumulative Expression Level (Sum FPKM)", fontsize=11, fontweight='bold')
ax2.set_title("Total Transcriptional Output", fontsize=12, fontweight='bold', pad=10)
ax2.grid(axis='x', linestyle='--', alpha=0.5)
for bar in bars2:
    width = bar.get_width()
    if width > 0.05:
        ax2.text(width + max(plot_stats['Total_FPKM']) * 0.01, bar.get_y() + bar.get_height()/2,
                 f"{width:.1f}", ha='left', va='center', fontsize=9, fontweight='bold')

legend_patches = [plt.Rectangle((0,0),1,1, color=color, label=label) for label, color in proc_colors.items()]
fig.legend(handles=legend_patches, title="Functional Digestive Process", loc='upper center', bbox_to_anchor=(0.5, 0.99), ncol=4, frameon=True, fontsize=10, title_fontsize=11)

plt.subplots_adjust(top=0.88, bottom=0.1, left=0.3, right=0.95, wspace=0.1)

fig3_pdf = "visualizations/01c_Digestion_Family_Repertoire_Breakdown.pdf"
fig3_png = "visualizations/01c_Digestion_Family_Repertoire_Breakdown.png"
fig.savefig(fig3_pdf, bbox_inches='tight', dpi=300)
fig.savefig(fig3_png, bbox_inches='tight', dpi=300)
fig.savefig("docs/assets/01c_Digestion_Family_Repertoire_Breakdown.pdf", bbox_inches='tight', dpi=300)
fig.savefig("docs/assets/01c_Digestion_Family_Repertoire_Breakdown.png", bbox_inches='tight', dpi=300)
fig.savefig("new plots/pdf/01c_Digestion_Family_Repertoire_Breakdown.pdf", bbox_inches='tight', dpi=300)
fig.savefig("new plots/png/01c_Digestion_Family_Repertoire_Breakdown.png", bbox_inches='tight', dpi=300)
plt.close()

print("All digestion analyses, tables, and clustered heatmaps successfully generated!")
