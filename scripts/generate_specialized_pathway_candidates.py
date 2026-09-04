import sqlite3
import pandas as pd
import numpy as np
import os
import plotly.express as px
import plotly.graph_objects as go
import zipfile, xml.sax.saxutils

print("Starting Targeted Functional Candidate Gene & Pathway Analysis...")

os.makedirs("visualizations", exist_ok=True)
os.makedirs("new plots/html", exist_ok=True)
os.makedirs("new plots/pdf", exist_ok=True)
os.makedirs("new plots/png", exist_ok=True)
os.makedirs("docs/assets", exist_ok=True)
os.makedirs("downstream_results", exist_ok=True)

# 1. Load Master KEGG & Domain Annotation tables
kegg_master = pd.read_csv('downstream_results/KEGG_Pathway_Annotation_Master.csv')
domain_master = pd.read_csv('downstream_results/Protein_Domain_Annotation_Master.csv')

# Merge to get unified dataset
merged = pd.merge(
    kegg_master, 
    domain_master[['Gene ID', 'protein_id', 'Pfam_Domains']], 
    on=['Gene ID', 'protein_id'], 
    how='outer'
)

# Calculate Mean FPKM
fpkm_cols = ['NPFM1', 'NPFM2', 'NPFM3', 'NPFM4']
merged['Mean_FPKM'] = merged[fpkm_cols].mean(axis=1)

# Define Process Classification Rules with precise keywords, Pfam domains, GO terms, and KEGG pathways
# The 6 Target Biological Processes:
# 1. Digestion
# 2. Detoxification
# 3. Sex determination
# 4. Immunity
# 5. Silk production
# 6. Toxin biogenesis

process_rules = {
    'Digestion': {
        'keywords': ['trypsin', 'peptidase', 'chymotrypsin', 'lipase', 'amylase', 'pepsin', 'carboxypeptidase', 'aminopeptidase', 'glucosidase', 'protease', 'cathepsin', 'chitinase', 'astacin'],
        'pfams': ['Trypsin', 'Peptidase_C1', 'Peptidase_S9', 'Peptidase_M10', 'Peptidase_M12A', 'Peptidase_M14', 'Peptidase_C13', 'Astacin', 'Lipase_3', 'Glyco_hydro'],
        'kegg_pathways': ['ko04974', 'map04974', 'ko04970', 'map04970', 'ko04972', 'map04972', 'ko04971', 'map04971', 'ko00010', 'map00010'] # Protein/Fat/Carbohydrate digestion & absorption
    },
    'Detoxification': {
        'keywords': ['cytochrome p450', 'cyp4', 'cyp6', 'cyp9', 'glutathione s-transferase', 'gst', 'glutathione peroxidase', 'superoxide dismutase', 'sod', 'catalase', 'abc transporter', 'udp-glucuronosyltransferase', 'ugt', 'carboxylesterase', 'aldehyde oxidase'],
        'pfams': ['p450', 'GST_N', 'GST_C', 'ABC_tran', 'UDPGT', 'Sod_Cu', 'Sod_Fe_N', 'Aldedh'],
        'kegg_pathways': ['ko00980', 'map00980', 'ko00982', 'map00982', 'ko00480', 'map00480', 'ko02010', 'map02010'] # Metabolism of xenobiotics by cytochrome P450, Drug metabolism, Glutathione metabolism, ABC transporters
    },
    'Sex determination': {
        'keywords': ['doublesex', 'dsx', 'transformer', 'tra', 'fem-1', 'fem-2', 'fem-3', 'sex lethal', 'sxl', 'fruitless', 'fru', 'intersex', 'ix', 'mab-3', 'dmrt', 'zinc finger protein mab', 'sox', 'tra-2', 'dead-box helicase', 'transformer-2'],
        'pfams': ['DM', 'SR_protein', 'RRM_1', 'zf-C2H2', 'HMG_box', 'DEAD'],
        'kegg_pathways': ['ko04310', 'map04310', 'ko03040', 'map03040', 'ko04914', 'map04914'] # Wnt signaling, Spliceosome, Progesterone-mediated oocyte maturation
    },
    'Immunity': {
        'keywords': ['toll', 'imd', 'spatzle', 'cactus', 'dorsal', 'relish', 'peptidoglycan recognition', 'pgrp', 'gram-negative', 'gnbp', 'scavenger receptor', 'defensin', 'cecropin', 'lectin', 'phenoloxidase', 'prophenoloxidase', 'serpin', 'caspase', 'autophagy'],
        'pfams': ['TIR', 'LRR_8', 'LRR_1', 'C1q', 'PGRP', 'Lectin_C', 'Serpin', 'Caspase', 'DEFT', 'Scavenger'],
        'kegg_pathways': ['ko04624', 'map04624', 'ko04064', 'map04064', 'ko04620', 'map04620', 'ko04145', 'map04145', 'ko04210', 'map04210'] # Toll and Imd signaling pathway, NF-kappa B, Toll-like receptor, Phagosome, Apoptosis
    },
    'Silk production': {
        'keywords': ['spidroin', 'masp', 'masp1', 'masp2', 'misp', 'flag', 'flagelliform', 'pyriform', 'pysp', 'acsp', 'aciniform', 'tubuliform', 'tusp', 'major ampullate', 'minor ampullate', 'dragline', 'silk', 'fibroin', 'collagen', 'zona pellucida', 'cuticle'],
        'pfams': ['Spidroin_N', 'Spidroin_C', 'Collagen', 'Zona_pellucida', 'Chitin_bind_4', 'TSP_1', 'fn3', 'EGF'],
        'kegg_pathways': ['ko04512', 'map04512', 'ko04510', 'map04510', 'ko04141', 'map04141'] # ECM-receptor interaction, Focal adhesion, Protein processing in ER
    },
    'Toxin biogenesis': {
        'keywords': ['toxin', 'venom', 'neurotoxin', 'shkt', 'kunitz', 'crisp', 'latrotoxin', 'latrodectus', 'phospholipase a2', 'hyaluronidase', 'icck', 'knottin', 'calcium channel toxin', 'spider toxin', 'necrosis', 'sphingomyelinase'],
        'pfams': ['ShKT', 'Kunitz_BPTI', 'CAP', 'Phospholip_A2_1', 'CBP', 'Latrotoxin', 'Kazal_1', 'Scorpion_toxin', 'Neurotoxin_1', 'Ank_2'],
        'kegg_pathways': ['ko04020', 'map04020', 'ko04724', 'map04724', 'ko04726', 'map04726', 'ko04210', 'map04210'] # Calcium signaling pathway, Glutamatergic synapse, Serotonergic synapse, Apoptosis
    }
}

# Tag genes to processes
assigned_rows = []

for idx, row in merged.iterrows():
    pdesc = str(row['Protein Description']).lower()
    pname = str(row['Preferred_Name']).lower()
    pfam_str = str(row['Pfam_Domains'])
    kegg_pw_str = str(row['KEGG_Pathway'])
    kegg_ko_str = str(row['KEGG_KO']).lower()
    
    assigned_procs = set()
    for proc, rules in process_rules.items():
        # Match keywords in description or preferred name
        if any(kw in pdesc or kw in pname for kw in rules['keywords']):
            assigned_procs.add(proc)
        # Match Pfam domains
        elif any(pf in pfam_str for pf in rules['pfams']):
            assigned_procs.add(proc)
        # Match KEGG pathways
        elif any(kp in kegg_pw_str for kp in rules['kegg_pathways']):
            assigned_procs.add(proc)
            
    for p in assigned_procs:
        r_copy = row.to_dict()
        r_copy['Biological_Process'] = p
        assigned_rows.append(r_copy)

candidates_df = pd.DataFrame(assigned_rows)
print(f"Total candidate gene instances identified across the 6 processes: {len(candidates_df)}")

# Summary by process
proc_counts = candidates_df['Biological_Process'].value_counts()
print("Candidate breakdown by Biological Process:\n", proc_counts)

# Clean and sort candidates
cols_order = [
    'Biological_Process', 'Gene ID', 'protein_id', 'Preferred_Name', 
    'KEGG_KO', 'KEGG_Pathway', 'Pfam_Domains', 'GO_Terms',
    'NPFM1', 'NPFM2', 'NPFM3', 'NPFM4', 'Mean_FPKM', 'Protein Description'
]
candidates_df = candidates_df[[c for c in cols_order if c in candidates_df.columns]].sort_values(
    by=['Biological_Process', 'Mean_FPKM'], 
    ascending=[True, False]
)

# Save Master Candidates CSV
candidates_df.to_csv('downstream_results/Targeted_Biological_Processes_Candidate_Genes.csv', index=False)
candidates_df.to_csv('docs/assets/Targeted_Biological_Processes_Candidate_Genes.csv', index=False)

# Multi-sheet Excel workbook builder
def write_multisheet_xlsx(df_dict, xlsx_path):
    sheet_names = list(df_dict.keys())
    
    content_types = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
                     '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">',
                     '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>',
                     '<Default Extension="xml" ContentType="application/xml"/>',
                     '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>',
                     '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>']
    for i in range(1, len(sheet_names) + 1):
        content_types.append(f'<Override PartName="/xl/worksheets/sheet{i}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>')
    content_types.append('</Types>')
    
    root_rels = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
    <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
</Relationships>'''

    wb_sheets = []
    wb_rels = []
    for i, sname in enumerate(sheet_names):
        wb_sheets.append(f'<sheet name="{sname[:31]}" sheetId="{i+1}" r:id="rId{i+1}"/>')
        wb_rels.append(f'<Relationship Id="rId{i+1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i+1}.xml"/>')
    wb_rels.append(f'<Relationship Id="rId{len(sheet_names)+1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>')

    wb = f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
    <sheets>{''.join(wb_sheets)}</sheets>
</workbook>'''

    wb_rels_xml = f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
    {''.join(wb_rels)}
</Relationships>'''

    styles = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
    <fonts count="1"><font><name val="Calibri"/><sz val="11"/></font></fonts>
    <fills count="2"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill></fills>
    <borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>
    <cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>
    <cellXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/></cellXfs>
</styleSheet>'''

    with zipfile.ZipFile(xlsx_path, 'w', compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr('[Content_Types].xml', ''.join(content_types))
        z.writestr('_rels/.rels', root_rels)
        z.writestr('xl/workbook.xml', wb)
        z.writestr('xl/_rels/workbook.xml.rels', wb_rels_xml)
        z.writestr('xl/styles.xml', styles)
        
        for i, (sname, df) in enumerate(df_dict.items()):
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
            z.writestr(f'xl/worksheets/sheet{i+1}.xml', sheet_xml)

# Create multi-sheet dictionary: Overview + 6 individual tabs
sheets_dict = {'All_Candidates': candidates_df}
for proc in process_rules.keys():
    sub_df = candidates_df[candidates_df['Biological_Process'] == proc].copy()
    sheets_dict[proc] = sub_df

write_multisheet_xlsx(sheets_dict, 'downstream_results/Targeted_Biological_Processes_Candidate_Genes.xlsx')
write_multisheet_xlsx(sheets_dict, 'docs/assets/Targeted_Biological_Processes_Candidate_Genes.xlsx')
print("Generated Multi-Sheet Targeted_Biological_Processes_Candidate_Genes.xlsx")

# -------------------------------------------------------------------------
# 2. Process Summary & Representation Figures
# -------------------------------------------------------------------------
proc_sum_rows = []
for proc, rules in process_rules.items():
    sub = candidates_df[candidates_df['Biological_Process'] == proc]
    top_expressed = sub.head(5)['Gene ID'].tolist()
    top_domains = []
    for p in sub['Pfam_Domains'].dropna():
        for d in str(p).split(','):
            if d.strip(): top_domains.append(d.strip())
    top_dom_str = ", ".join(pd.Series(top_domains).value_counts().head(3).index.tolist())
    
    top_keggs = []
    for pw in sub['KEGG_Pathway'].dropna():
        for k in str(pw).split(','):
            if k.strip() and k.strip() != '-': top_keggs.append(k.strip())
    top_kegg_str = ", ".join(pd.Series(top_keggs).value_counts().head(3).index.tolist())
    
    proc_sum_rows.append({
        'Biological_Process': proc,
        'Candidate_Genes_Count': len(sub),
        'Mean_Expression_FPKM': round(sub['Mean_FPKM'].mean(), 2),
        'Top_Pfam_Domains': top_dom_str,
        'Top_KEGG_Pathways': top_kegg_str,
        'Key_Expressed_Candidates': ", ".join(top_expressed)
    })

proc_sum_df = pd.DataFrame(proc_sum_rows)
proc_sum_df.to_csv('downstream_results/Targeted_Processes_Summary_Stats.csv', index=False)
proc_sum_df.to_csv('docs/assets/Targeted_Processes_Summary_Stats.csv', index=False)
write_multisheet_xlsx({'Summary': proc_sum_df}, 'downstream_results/Targeted_Processes_Summary_Stats.xlsx')
write_multisheet_xlsx({'Summary': proc_sum_df}, 'docs/assets/Targeted_Processes_Summary_Stats.xlsx')
print("Saved Targeted_Processes_Summary_Stats.xlsx:\n", proc_sum_df[['Biological_Process', 'Candidate_Genes_Count', 'Mean_Expression_FPKM']])

# -------------------------------------------------------------------------
# 3. High-Quality Multi-Panel / Faceted Visualizations
# -------------------------------------------------------------------------
# Plot 1: Candidate Counts & Mean Expression by Process
proc_colors = {
    'Digestion': '#f59e0b',           # Amber
    'Detoxification': '#10b981',       # Emerald
    'Sex determination': '#ec4899',   # Rose/Pink
    'Immunity': '#ef4444',            # Red
    'Silk production': '#8b5cf6',     # Violet
    'Toxin biogenesis': '#06b6d4'     # Cyan
}

fig_bar = go.Figure()
fig_bar.add_trace(go.Bar(
    x=proc_sum_df['Biological_Process'],
    y=proc_sum_df['Candidate_Genes_Count'],
    marker=dict(
        color=[proc_colors[p] for p in proc_sum_df['Biological_Process']],
        line=dict(color='#0f172a', width=1.5)
    ),
    text=proc_sum_df['Candidate_Genes_Count'],
    textposition='outside'
))

fig_bar.update_layout(
    title=dict(
        text="<b>Candidate Genes Identified across Key Spider Physiological Processes</b><br><sup>Digestion, Detoxification, Sex Determination, Immunity, Silk Production & Toxin Biogenesis</sup>",
        font=dict(size=18, color="#0f172a", family="Plus Jakarta Sans, sans-serif"),
        x=0.5
    ),
    xaxis=dict(title="<b>Biological Process</b>", showline=True, linecolor="#cbd5e1", linewidth=1.5),
    yaxis=dict(title="<b>Identified Candidate Genes</b>", showline=True, linecolor="#cbd5e1", linewidth=1.5, showgrid=True, gridcolor="#f1f5f9"),
    paper_bgcolor="white",
    plot_bgcolor="white",
    width=950,
    height=600,
    margin=dict(l=70, r=40, t=90, b=60)
)

fig_bar.write_html("new plots/html/plotly_target_processes_bar.html")
fig_bar.write_image("visualizations/18_Target_Processes_Candidates_Bar.png", scale=3)
fig_bar.write_image("visualizations/18_Target_Processes_Candidates_Bar.pdf")
fig_bar.write_image("docs/assets/18_Target_Processes_Candidates_Bar.png", scale=3)
fig_bar.write_image("docs/assets/18_Target_Processes_Candidates_Bar.pdf")
print("Saved Target Processes Candidate Bar Chart")

# Plot 2: Expression Profiles of Top Candidate Genes across the 6 Processes
top_candidates = candidates_df.groupby('Biological_Process').head(4).copy()
top_candidates['Label'] = top_candidates['Gene ID'] + " (" + top_candidates['Preferred_Name'].fillna('Unknown') + ")"

fig_heat = go.Figure()
for s in fpkm_cols:
    fig_heat.add_trace(go.Bar(
        name=s,
        x=top_candidates['Label'],
        y=top_candidates[s],
        marker=dict(line=dict(width=0.5))
    ))

fig_heat.update_layout(
    barmode='group',
    title=dict(
        text="<b>Expression (FPKM) of Top Candidate Transcripts per Process</b><br><sup>Comparison across NPFM1 - NPFM4</sup>",
        font=dict(size=18, color="#0f172a", family="Plus Jakarta Sans, sans-serif"),
        x=0.5
    ),
    xaxis=dict(title="<b>Candidate Gene & Preferred Symbol</b>", tickangle=-45, showline=True, linecolor="#cbd5e1"),
    yaxis=dict(title="<b>Expression Level (FPKM)</b>", showline=True, linecolor="#cbd5e1", showgrid=True, gridcolor="#f1f5f9"),
    paper_bgcolor="white",
    plot_bgcolor="white",
    width=1100,
    height=650,
    legend=dict(title="<b>Sample</b>", bgcolor="rgba(255,255,255,0.9)", bordercolor="#cbd5e1", borderwidth=1),
    margin=dict(l=70, r=40, t=90, b=150)
)

fig_heat.write_html("new plots/html/plotly_target_candidates_expression.html")
fig_heat.write_image("visualizations/19_Target_Candidates_Expression_Profile.png", scale=3)
fig_heat.write_image("visualizations/19_Target_Candidates_Expression_Profile.pdf")
fig_heat.write_image("docs/assets/19_Target_Candidates_Expression_Profile.png", scale=3)
fig_heat.write_image("docs/assets/19_Target_Candidates_Expression_Profile.pdf")
print("Saved Target Candidates Expression Profile Bar Chart")

print("All Targeted Process Deliverables and Visualizations built successfully!")
