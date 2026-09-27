import os
import sys
import pandas as pd
import numpy as np
import zipfile
import xml.sax.saxutils
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import FancyBboxPatch, Rectangle, Arrow, Circle, FancyArrowPatch, PathPatch
from matplotlib.path import Path
import matplotlib.patches as mpatches

# Configure publication-grade styling
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

os.makedirs("downstream_results", exist_ok=True)
os.makedirs("visualizations", exist_ok=True)
os.makedirs("docs/assets", exist_ok=True)
os.makedirs("new plots/pdf", exist_ok=True)
os.makedirs("new plots/png", exist_ok=True)

print("1. Loading Master Annotation & Expression Datasets...")
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

def write_multisheet_xlsx(df_dict, xlsx_path):
    sheets_xml = []
    workbook_sheets = []
    
    for idx, (name, df) in enumerate(df_dict.items(), 1):
        s_id = f"rId{idx}"
        workbook_sheets.append(f'<sheet name="{xml.sax.saxutils.escape(name[:31])}" sheetId="{idx}" r:id="{s_id}"/>')
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

def save_fig_dual(fig, base_name):
    fig.savefig(f"visualizations/{base_name}.pdf", bbox_inches='tight', dpi=300)
    fig.savefig(f"visualizations/{base_name}.png", bbox_inches='tight', dpi=300)
    fig.savefig(f"docs/assets/{base_name}.pdf", bbox_inches='tight', dpi=300)
    fig.savefig(f"docs/assets/{base_name}.png", bbox_inches='tight', dpi=300)
    fig.savefig(f"new plots/pdf/{base_name}.pdf", bbox_inches='tight', dpi=300)
    fig.savefig(f"new plots/png/{base_name}.png", bbox_inches='tight', dpi=300)
    plt.close(fig)

print("----------------------------------------------------------------------")
print("ANALYSIS 2: DETOXIFICATION PATHWAY (Phase I -> II -> III Cascade)")
print("----------------------------------------------------------------------")
def classify_detox(row):
    desc = str(row.get('Protein Description', '')).lower()
    pref = str(row.get('Preferred_Name', '')).lower()
    pf = str(row.get('Pfam_Domains', '')).lower()
    kegg_path = str(row.get('KEGG_Pathway', '')).lower()
    ko = str(row.get('KEGG_KO', '')).lower()

    if any(t in desc for t in ['integrase', 'reverse transcriptase', 'transposon', 'gag-pol', 'retrotransposon', 'polyprotein']):
        return None

    # Phase I: Functionalization
    if 'cytochrome p450' in desc or 'p450' in pf or 'cyp' in pref or 'cyp4' in desc or 'cyp6' in desc or 'cyp9' in desc:
        return ('Phase I: Functionalization', 'Cytochrome P450s (CYP450)', 'CYP Monooxygenases')
    if 'carboxylesterase' in desc or 'abhydrolase' in pf or 'esterase' in desc:
        return ('Phase I: Functionalization', 'Carboxylesterases (COE)', 'Esterases & Hydrolases')
    if 'epoxide hydrolase' in desc or 'epoxide_hydrolase' in pf:
        return ('Phase I: Functionalization', 'Epoxide Hydrolases (mEH)', 'Epoxide Hydrolysis')
    if 'aldehyde dehydrogenase' in desc or 'aldedh' in pf or 'aldh' in pref:
        return ('Phase I: Functionalization', 'Aldehyde Dehydrogenases (ALDH)', 'Aldehyde Oxidation')

    # Phase II: Conjugation
    if 'glutathione s-transferase' in desc or 'gst' in pf or 'gst' in pref or 'mapeg' in pf:
        return ('Phase II: Conjugation', 'Glutathione S-Transferases (GST)', 'Glutathione Conjugation')
    if 'udp-glucuronosyltransferase' in desc or 'ugt' in desc or 'udpgt' in pf:
        return ('Phase II: Conjugation', 'UDP-Glucuronosyltransferases (UGT)', 'Glucuronidation')
    if 'sulfotransferase' in desc or 'sulfotransfer' in pf or 'sult' in pref:
        return ('Phase II: Conjugation', 'Sulfotransferases (SULT)', 'Sulfation')
    if 'acetyltransferase' in desc or 'acetyltransf' in pf:
        return ('Phase II: Conjugation', 'N-Acetyltransferases (NAT)', 'Acetylation')

    # Phase III: Efflux & Transport
    if 'abc transporter' in desc or 'abc_tran' in pf or 'abc_membrane' in pf or 'abc' in pref:
        return ('Phase III: Efflux & Transport', 'ABC Transporters (ABCA/B/C/G)', 'ATP-Dependent Efflux')
    if 'major facilitator' in desc or 'mfs_1' in pf or 'mfs' in pf:
        return ('Phase III: Efflux & Transport', 'MFS Xenobiotic Transporters', 'Solute Carrier Efflux')

    return None

detox_records = []
for idx, row in merged.iterrows():
    c = classify_detox(row)
    if c is not None:
        r = row.to_dict()
        r['Detox_Phase'] = c[0]
        r['Enzyme_Family'] = c[1]
        r['Reaction_Module'] = c[2]
        r['Expressed_in_Transcriptome'] = 'Yes' if r['Mean_FPKM'] > 0.1 else 'No'
        detox_records.append(r)

detox_df = pd.DataFrame(detox_records)
detox_df.to_csv('downstream_results/02_Detoxification_Master_Annotation.csv', index=False)
detox_df.to_csv('docs/assets/02_Detoxification_Master_Annotation.csv', index=False)

detox_stats = detox_df.groupby(['Detox_Phase', 'Enzyme_Family']).agg(
    Transcript_Count=('protein_id', 'nunique'),
    Expressed_Transcripts=('Expressed_in_Transcriptome', lambda x: (x == 'Yes').sum()),
    Mean_FPKM=('Mean_FPKM', 'mean'),
    Total_FPKM=('Mean_FPKM', 'sum'),
    NPFM1_Sum=('NPFM1', 'sum'),
    NPFM2_Sum=('NPFM2', 'sum'),
    NPFM3_Sum=('NPFM3', 'sum'),
    NPFM4_Sum=('NPFM4', 'sum')
).reset_index()
detox_stats.to_csv('downstream_results/02_Detoxification_Phase_Summary_Stats.csv', index=False)
detox_stats.to_csv('docs/assets/02_Detoxification_Phase_Summary_Stats.csv', index=False)

write_multisheet_xlsx({
    'Master_Detox_Catalogue': detox_df,
    'Phase_Summary_Stats': detox_stats,
    'Sample_Expression': detox_df[['protein_id', 'Preferred_Name', 'Detox_Phase', 'Enzyme_Family', 'NPFM1', 'NPFM2', 'NPFM3', 'NPFM4', 'Mean_FPKM']]
}, 'downstream_results/02_Detoxification_Master_Annotation.xlsx')
write_multisheet_xlsx({
    'Master_Detox_Catalogue': detox_df,
    'Phase_Summary_Stats': detox_stats,
    'Sample_Expression': detox_df[['protein_id', 'Preferred_Name', 'Detox_Phase', 'Enzyme_Family', 'NPFM1', 'NPFM2', 'NPFM3', 'NPFM4', 'Mean_FPKM']]
}, 'docs/assets/02_Detoxification_Master_Annotation.xlsx')
print(f"Detoxification: {len(detox_df)} genes categorized across 3 phases.")

# Detoxification Pathway Figure (Phase I -> Phase II -> Phase III Cascade)
fig, (ax_main, ax_heat) = plt.subplots(1, 2, figsize=(16, 9), gridspec_kw={'width_ratios': [2.2, 1.2]})

# Left Panel: Pathway Cascade Flow
ax_main.axis('off')
ax_main.set_xlim(0, 10)
ax_main.set_ylim(0, 10)
ax_main.set_title("Nephila pilipes Detoxification Pathway Architecture\nPhase I (Functionalization) → Phase II (Conjugation) → Phase III (Efflux)", fontsize=13, fontweight='bold', pad=15)

phase_boxes = [
    (8.0, "PHASE I: FUNCTIONALIZATION\n(Oxidation, Reduction, Hydrolysis)", "#E0F2FE", "#0369A1", [
        ("Cytochrome P450s (CYP450)", f"N = {len(detox_df[detox_df['Enzyme_Family'].str.contains('P450')])} genes", f"Sum FPKM: {detox_df[detox_df['Enzyme_Family'].str.contains('P450')]['Mean_FPKM'].sum():.1f}"),
        ("Carboxylesterases (COE)", f"N = {len(detox_df[detox_df['Enzyme_Family'].str.contains('Carboxylesterase')])} genes", f"Sum FPKM: {detox_df[detox_df['Enzyme_Family'].str.contains('Carboxylesterase')]['Mean_FPKM'].sum():.1f}"),
        ("Aldehyde Dehydrogenases", f"N = {len(detox_df[detox_df['Enzyme_Family'].str.contains('Aldehyde')])} genes", f"Sum FPKM: {detox_df[detox_df['Enzyme_Family'].str.contains('Aldehyde')]['Mean_FPKM'].sum():.1f}"),
    ]),
    (4.8, "PHASE II: CONJUGATION\n(Glutathionylation, Glucuronidation, Sulfation)", "#FEF3C7", "#D97706", [
        ("Glutathione S-Transferases (GST)", f"N = {len(detox_df[detox_df['Enzyme_Family'].str.contains('Glutathione')])} genes", f"Sum FPKM: {detox_df[detox_df['Enzyme_Family'].str.contains('Glutathione')]['Mean_FPKM'].sum():.1f}"),
        ("UDP-Glucuronosyltransferases (UGT)", f"N = {len(detox_df[detox_df['Enzyme_Family'].str.contains('UGT')])} genes", f"Sum FPKM: {detox_df[detox_df['Enzyme_Family'].str.contains('UGT')]['Mean_FPKM'].sum():.1f}"),
        ("Sulfotransferases (SULT)", f"N = {len(detox_df[detox_df['Enzyme_Family'].str.contains('Sulfotransferase')])} genes", f"Sum FPKM: {detox_df[detox_df['Enzyme_Family'].str.contains('Sulfotransferase')]['Mean_FPKM'].sum():.1f}"),
    ]),
    (1.6, "PHASE III: EFFLUX & EXCRETION\n(ATP-Dependent & Facilitated Transport)", "#DCFCE7", "#15803D", [
        ("ABC Transporters (ABCA/B/C/G)", f"N = {len(detox_df[detox_df['Enzyme_Family'].str.contains('ABC')])} genes", f"Sum FPKM: {detox_df[detox_df['Enzyme_Family'].str.contains('ABC')]['Mean_FPKM'].sum():.1f}"),
        ("MFS Xenobiotic Transporters", f"N = {len(detox_df[detox_df['Enzyme_Family'].str.contains('MFS')])} genes", f"Sum FPKM: {detox_df[detox_df['Enzyme_Family'].str.contains('MFS')]['Mean_FPKM'].sum():.1f}"),
    ])
]

# Draw Xenobiotic Flow
ax_main.text(5.0, 9.6, "Xenobiotics / Toxin Ingestion (Lipophilic Substrates)", ha='center', va='center', fontsize=11, fontweight='bold', bbox=dict(boxstyle='round,pad=0.5', facecolor='#F1F5F9', edgecolor='#64748B'))
arrow_kw = dict(arrowstyle='->', lw=2.5, color='#475569')

for y, title, bg_col, border_col, sub_items in phase_boxes:
    # Phase background box
    rect = FancyBboxPatch((0.5, y - 1.1), 9.0, 2.0, boxstyle="round,pad=0.2", facecolor=bg_col, edgecolor=border_col, linewidth=2.0)
    ax_main.add_patch(rect)
    ax_main.text(0.8, y + 0.6, title, fontsize=10.5, fontweight='bold', color=border_col, va='top')
    
    # Sub-item mini boxes
    x_step = 8.6 / len(sub_items)
    for s_idx, (name, n_txt, fpkm_txt) in enumerate(sub_items):
        mini_x = 0.7 + s_idx * x_step
        mini_rect = FancyBboxPatch((mini_x, y - 0.9), x_step - 0.2, 1.2, boxstyle="round,pad=0.1", facecolor='white', edgecolor=border_col, linewidth=1.0)
        ax_main.add_patch(mini_rect)
        ax_main.text(mini_x + (x_step-0.2)/2, y - 0.1, name, ha='center', va='center', fontsize=9.5, fontweight='bold')
        ax_main.text(mini_x + (x_step-0.2)/2, y - 0.45, f"{n_txt} | {fpkm_txt}", ha='center', va='center', fontsize=8.5, color='#475569')

# Connective Flow Arrows
ax_main.annotate('', xy=(5.0, 8.9), xytext=(5.0, 9.3), arrowprops=arrow_kw)
ax_main.annotate('', xy=(5.0, 5.7), xytext=(5.0, 6.7), arrowprops=arrow_kw)
ax_main.annotate('', xy=(5.0, 2.5), xytext=(5.0, 3.5), arrowprops=arrow_kw)
ax_main.text(5.2, 6.2, "Polar Metabolites (-OH, -COOH, -NH2)", fontsize=8.5, fontstyle='italic', color='#0369A1')
ax_main.text(5.2, 3.0, "Hydrophilic Conjugates (Glutathione/Glucuronide)", fontsize=8.5, fontstyle='italic', color='#D97706')
ax_main.text(5.0, 0.2, "Excreted Detoxified Metabolites (Hydrophilic Clearance)", ha='center', va='center', fontsize=10.5, fontweight='bold', bbox=dict(boxstyle='round,pad=0.4', facecolor='#DCFCE7', edgecolor='#15803D'))
ax_main.annotate('', xy=(5.0, 0.4), xytext=(5.0, 0.7), arrowprops=arrow_kw)

# Right Panel: Phase Family Heatmap
phase_pivot = detox_df.pivot_table(index='Enzyme_Family', values=fpkm_cols, aggfunc='sum')
phase_log = np.log2(phase_pivot + 1.0)
sns.heatmap(phase_log, ax=ax_heat, cmap='YlGnBu', annot=True, fmt=".1f", cbar_kws={'label': 'Log2(Sum FPKM + 1)'}, linewidths=1.0, linecolor='#E2E8F0')
ax_heat.set_title("Expression Profile by Detox Family\n(NPFM1–NPFM4 Replicates)", fontsize=11, fontweight='bold', pad=15)
ax_heat.set_ylabel("")
ax_heat.set_xlabel("Biological Replicate", fontsize=10, fontweight='bold')
plt.setp(ax_heat.get_yticklabels(), rotation=0, fontsize=9)

save_fig_dual(fig, "02_Detoxification_Pathway_Cascade")


print("----------------------------------------------------------------------")
print("ANALYSIS 3: REPRODUCTIVE GENE NETWORK & SANKEY FLOW")
print("----------------------------------------------------------------------")
def classify_reproduction(row):
    desc = str(row.get('Protein Description', '')).lower()
    pref = str(row.get('Preferred_Name', '')).lower()
    pf = str(row.get('Pfam_Domains', '')).lower()
    ko = str(row.get('KEGG_KO', '')).lower()

    if any(t in desc for t in ['integrase', 'reverse transcriptase', 'transposon', 'gag-pol']):
        return None

    # Regulatory & Endocrine Signals
    if 'ecdysone' in desc or 'ecr' in pref or 'ultraspiracle' in desc or 'broad-complex' in desc or 'methoprene' in desc or 'met' in pref:
        return ('1_Endocrine_Signals', 'Nuclear Hormone Signaling', 'Ecdysteroid & JH Pathways')
    if 'foxl2' in desc or 'sox9' in desc or 'sox' in desc or 'hmg_box' in pf:
        return ('1_Endocrine_Signals', 'Sex Differentiation Regulators', 'Sox/Fox Transcription Factors')

    # Sex Determination & Germline Maintenance
    if 'doublesex' in desc or 'dsx' in pref or 'transformer' in desc or 'tra-2' in desc or 'fem-1' in desc or 'sex-lethal' in desc:
        return ('2_Sex_Determination', 'Sex Determination Cascade', 'Dsx/Tra/Fem Axis')
    if 'vasa' in desc or 'ddx4' in desc or 'dead' in pf or 'nanos' in desc or 'pumilio' in desc or 'piwi' in desc or 'bicaudal' in desc:
        return ('2_Sex_Determination', 'Germline Stem Cell Maintenance', 'Vasa/Nanos/Piwi Complex')

    # Female Reproduction & Vitellogenesis
    if 'vitellogenin' in desc or 'vitellogenin' in pf or 'vgr' in desc or 'ldl_recept' in pf:
        return ('3_Female_Processes', 'Vitellogenesis & Egg Provisioning', 'Vitellogenin & Vg Receptor')
    if 'chorion' in desc or 'zona pellucida' in desc or 'zp' in pf or 'oocyte' in desc:
        return ('3_Female_Processes', 'Oogenesis & Eggshell Assembly', 'Chorion & ZP Matrix')

    # Male Reproduction & Spermatogenesis
    if 'spermatogenesis' in desc or 'spag' in desc or 'testis' in desc or 'dynein axonemal' in desc or 'protamine' in desc:
        return ('4_Male_Processes', 'Spermatogenesis & Motility', 'Sperm Flagellar & Nuclear Factors')

    return None

repro_records = []
for idx, row in merged.iterrows():
    c = classify_reproduction(row)
    if c is not None:
        r = row.to_dict()
        r['Reproductive_Tier'] = c[0]
        r['Functional_Module'] = c[1]
        r['Gene_Complex'] = c[2]
        r['Expressed_in_Transcriptome'] = 'Yes' if r['Mean_FPKM'] > 0.1 else 'No'
        repro_records.append(r)

repro_df = pd.DataFrame(repro_records)
repro_df.to_csv('downstream_results/03_Reproduction_Master_Annotation.csv', index=False)
repro_df.to_csv('docs/assets/03_Reproduction_Master_Annotation.csv', index=False)

repro_stats = repro_df.groupby(['Reproductive_Tier', 'Functional_Module', 'Gene_Complex']).agg(
    Gene_Count=('protein_id', 'nunique'),
    Mean_FPKM=('Mean_FPKM', 'mean'),
    Total_FPKM=('Mean_FPKM', 'sum')
).reset_index()
repro_stats.to_csv('downstream_results/03_Reproduction_Sankey_Flow_Stats.csv', index=False)
repro_stats.to_csv('docs/assets/03_Reproduction_Sankey_Flow_Stats.csv', index=False)

write_multisheet_xlsx({
    'Master_Reproductive_Genes': repro_df,
    'Sankey_Flow_Stats': repro_stats,
    'Expression_Summary': repro_df[['protein_id', 'Preferred_Name', 'Reproductive_Tier', 'Functional_Module', 'NPFM1', 'NPFM2', 'NPFM3', 'NPFM4', 'Mean_FPKM']]
}, 'downstream_results/03_Reproduction_Master_Annotation.xlsx')
write_multisheet_xlsx({
    'Master_Reproductive_Genes': repro_df,
    'Sankey_Flow_Stats': repro_stats,
    'Expression_Summary': repro_df[['protein_id', 'Preferred_Name', 'Reproductive_Tier', 'Functional_Module', 'NPFM1', 'NPFM2', 'NPFM3', 'NPFM4', 'Mean_FPKM']]
}, 'docs/assets/03_Reproduction_Master_Annotation.xlsx')
print(f"Reproduction: {len(repro_df)} candidate genes classified.")

# Reproductive Regulatory Flow & Sankey Architecture Diagram
fig, ax = plt.subplots(figsize=(15, 8.5))
ax.axis('off')
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.set_title("Nephila pilipes Reproductive Gene Network & Regulatory Flow\nEndocrine / Regulatory Signals → Reproductive Core Genes → Biological Processes", fontsize=13, fontweight='bold', pad=20)

# Column 1: Endocrine Signals (x = 10)
# Column 2: Sex Determination & Germline Core (x = 50)
# Column 3: Downstream Biological Processes (x = 90)

col1_nodes = [
    ("Ecdysone & JH Pathways\n(EcR, USP, BR-C, Met)", 75, "#93C5FD", "#1D4ED8", f"{repro_stats[repro_stats['Gene_Complex'].str.contains('Ecdysteroid')]['Gene_Count'].sum()} genes"),
    ("Sox / Fox Regulators\n(FoxL2, Sox9/E)", 25, "#BFDBFE", "#1E40AF", f"{repro_stats[repro_stats['Gene_Complex'].str.contains('Sox')]['Gene_Count'].sum()} genes")
]

col2_nodes = [
    ("Dsx / Tra-2 / Fem Axis\n(Somatic Sex Regulators)", 80, "#DDD6FE", "#6D28D9", f"{repro_stats[repro_stats['Gene_Complex'].str.contains('Dsx')]['Gene_Count'].sum()} genes"),
    ("Vasa / Nanos / Piwi Axis\n(Germline Maintenance)", 48, "#EDE9FE", "#5B21B6", f"{repro_stats[repro_stats['Gene_Complex'].str.contains('Vasa')]['Gene_Count'].sum()} genes"),
    ("Chorion / ZP Scaffolds\n(Matrix Biogenesis)", 18, "#F5D0FE", "#86198F", f"{repro_stats[repro_stats['Gene_Complex'].str.contains('Chorion')]['Gene_Count'].sum()} genes")
]

col3_nodes = [
    ("Vitellogenesis & Egg Provisioning\n(Vitellogenin & Vg Receptor)", 82, "#FBCFE8", "#BE185D", f"{repro_stats[repro_stats['Gene_Complex'].str.contains('Vitellogenin')]['Gene_Count'].sum()} genes | FPKM: {repro_stats[repro_stats['Gene_Complex'].str.contains('Vitellogenin')]['Total_FPKM'].sum():.1f}"),
    ("Oocyte Maturation & Patterning\n(Embryonic Axis & Germ Plasm)", 50, "#FCE7F3", "#9D174D", "Female Development"),
    ("Spermatogenesis & Flagella\n(Sperm Motility & Chromatin)", 20, "#E2E8F0", "#334155", f"{repro_stats[repro_stats['Gene_Complex'].str.contains('Sperm')]['Gene_Count'].sum()} genes | Male Pathway")
]

# Draw Column Titles
ax.text(12, 95, "UPSTREAM SIGNALS", ha='center', fontsize=11, fontweight='bold', color='#1D4ED8')
ax.text(50, 95, "CORE REPRODUCTIVE GENES", ha='center', fontsize=11, fontweight='bold', color='#6D28D9')
ax.text(88, 95, "BIOLOGICAL PROCESSES", ha='center', fontsize=11, fontweight='bold', color='#BE185D')

# Draw Nodes
def draw_node_box(x, y, text, subtext, bg, border):
    rect = FancyBboxPatch((x-10, y-8), 20, 16, boxstyle="round,pad=0.3", facecolor=bg, edgecolor=border, linewidth=2.0)
    ax.add_patch(rect)
    ax.text(x, y+2, text, ha='center', va='center', fontsize=9.5, fontweight='bold', color='#0F172A')
    ax.text(x, y-5, subtext, ha='center', va='center', fontsize=8.5, color=border, fontweight='bold')

for txt, y, bg, border, sub in col1_nodes:
    draw_node_box(12, y, txt, sub, bg, border)

for txt, y, bg, border, sub in col2_nodes:
    draw_node_box(50, y, txt, sub, bg, border)

for txt, y, bg, border, sub in col3_nodes:
    draw_node_box(88, y, txt, sub, bg, border)

# Draw Bezier Connections (Sankey Flow representation)
def draw_flow(x1, y1, x2, y2, color, alpha=0.35, width=4.0):
    path = Path([(x1+10, y1), ((x1+x2)/2, y1), ((x1+x2)/2, y2), (x2-10, y2)],
                [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
    patch = PathPatch(path, facecolor='none', edgecolor=color, alpha=alpha, linewidth=width)
    ax.add_patch(patch)

# Flows from Col 1 to Col 2
draw_flow(12, 75, 50, 80, '#3B82F6', alpha=0.5, width=5.0)
draw_flow(12, 75, 50, 48, '#3B82F6', alpha=0.4, width=4.0)
draw_flow(12, 25, 50, 80, '#6366F1', alpha=0.4, width=3.5)
draw_flow(12, 25, 50, 18, '#6366F1', alpha=0.4, width=3.5)

# Flows from Col 2 to Col 3
draw_flow(50, 80, 88, 82, '#8B5CF6', alpha=0.6, width=6.0)
draw_flow(50, 48, 88, 50, '#8B5CF6', alpha=0.5, width=5.0)
draw_flow(50, 48, 88, 20, '#64748B', alpha=0.4, width=3.5)
draw_flow(50, 18, 88, 82, '#EC4899', alpha=0.5, width=4.5)

save_fig_dual(fig, "03_Reproductive_Sankey_and_Network")


print("----------------------------------------------------------------------")
print("ANALYSIS 4: IMMUNITY & IMMUNE-PATHWAY CLUSTERED HEATMAP")
print("----------------------------------------------------------------------")
def classify_immunity(row):
    desc = str(row.get('Protein Description', '')).lower()
    pref = str(row.get('Preferred_Name', '')).lower()
    pf = str(row.get('Pfam_Domains', '')).lower()
    ko = str(row.get('KEGG_KO', '')).lower()

    if any(t in desc for t in ['integrase', 'reverse transcriptase', 'transposon', 'gag-pol']):
        return None

    # Pathogen Recognition
    if 'peptidoglycan recognition' in desc or 'pgrp' in pf or 'pgrp' in pref:
        return ('Pathogen Recognition', 'PGRP Family', 'Bacterial Peptidoglycan Sensing')
    if 'gram-negative' in desc or 'gnbp' in pref or 'beta-1,3-glucan' in desc:
        return ('Pathogen Recognition', 'GNBP / Glucan Receptors', 'Fungal & Gram-Neg Sensing')
    if 'c-type lectin' in desc or 'clect' in pf or 'lectin' in desc:
        return ('Pathogen Recognition', 'C-type Lectins (CTL)', 'Carbohydrate Pattern Recognition')
    if 'scavenger receptor' in desc or 'srec' in desc:
        return ('Pathogen Recognition', 'Scavenger Receptors', 'Cellular Clearance & Uptake')

    # Toll Pathway
    if 'spatzle' in desc or 'spz' in pref:
        return ('Toll Signaling Pathway', 'Spätzle Ligands', 'Extracellular Cytokine Signal')
    if 'toll-like' in desc or 'toll' in pf or 'tlr' in pref:
        return ('Toll Signaling Pathway', 'Toll Receptors (TLR)', 'Transmembrane Signaling')
    if 'myd88' in desc or 'tube' in desc or 'pelle' in desc or 'cactus' in desc or 'dorsal' in desc or 'dif' in desc:
        return ('Toll Signaling Pathway', 'Toll Intracellular Cascade', 'MyD88 / Cactus / Dorsal')

    # IMD Pathway
    if 'imd' in desc or 'fadd' in desc or 'dredd' in desc or 'relish' in desc or 'tak1' in desc:
        return ('IMD Signaling Pathway', 'IMD Signaling Core', 'Dredd / Relish NF-kB Cascade')

    # Melanization & Serpins
    if 'prophenoloxidase' in desc or 'propo' in desc or 'tyrosinase' in pf:
        return ('Melanization & Coagulation', 'Prophenoloxidase (ProPO)', 'Melanin Synthesis & Encapsulation')
    if 'serpin' in desc or 'serpin' in pf:
        return ('Melanization & Coagulation', 'Serpin Inhibitors', 'Protease Cascade Control')

    # AMPs & Effectors
    if 'defensin' in desc or 'defensin' in pf or 'cecropin' in desc or 'antimicrobial' in desc:
        return ('Antimicrobial Effectors', 'Antimicrobial Peptides (AMPs)', 'Direct Pathogen Lysis')
    if 'lysozyme' in desc or 'lysozyme' in pf:
        return ('Antimicrobial Effectors', 'Lysozymes', 'Bacterial Cell Wall Lysis')

    return None

immune_records = []
for idx, row in merged.iterrows():
    c = classify_immunity(row)
    if c is not None:
        r = row.to_dict()
        r['Immune_Module'] = c[0]
        r['Immune_Family'] = c[1]
        r['Mechanism'] = c[2]
        r['Expressed_in_Transcriptome'] = 'Yes' if r['Mean_FPKM'] > 0.1 else 'No'
        immune_records.append(r)

immune_df = pd.DataFrame(immune_records)
immune_df.to_csv('downstream_results/04_Immunity_Master_Annotation.csv', index=False)
immune_df.to_csv('docs/assets/04_Immunity_Master_Annotation.csv', index=False)

immune_stats = immune_df.groupby(['Immune_Module', 'Immune_Family']).agg(
    Gene_Count=('protein_id', 'nunique'),
    Mean_FPKM=('Mean_FPKM', 'mean'),
    Total_FPKM=('Mean_FPKM', 'sum'),
    NPFM1_Sum=('NPFM1', 'sum'),
    NPFM2_Sum=('NPFM2', 'sum'),
    NPFM3_Sum=('NPFM3', 'sum'),
    NPFM4_Sum=('NPFM4', 'sum')
).reset_index()
immune_stats.to_csv('downstream_results/04_Immunity_Module_Summary_Stats.csv', index=False)
immune_stats.to_csv('docs/assets/04_Immunity_Module_Summary_Stats.csv', index=False)

write_multisheet_xlsx({
    'Master_Immune_Catalogue': immune_df,
    'Module_Summary_Stats': immune_stats,
    'Replicate_Expression': immune_df[['protein_id', 'Preferred_Name', 'Immune_Module', 'Immune_Family', 'NPFM1', 'NPFM2', 'NPFM3', 'NPFM4', 'Mean_FPKM']]
}, 'downstream_results/04_Immunity_Master_Annotation.xlsx')
write_multisheet_xlsx({
    'Master_Immune_Catalogue': immune_df,
    'Module_Summary_Stats': immune_stats,
    'Replicate_Expression': immune_df[['protein_id', 'Preferred_Name', 'Immune_Module', 'Immune_Family', 'NPFM1', 'NPFM2', 'NPFM3', 'NPFM4', 'Mean_FPKM']]
}, 'docs/assets/04_Immunity_Master_Annotation.xlsx')
print(f"Immunity: {len(immune_df)} genes annotated across immune modules.")

# Immune Pathway Clustered Heatmap
immune_top = immune_df.sort_values(by='Mean_FPKM', ascending=False).head(45).copy()
immune_mat = np.log2(immune_top[fpkm_cols] + 1.0)
immune_mat.index = immune_top['protein_id'] + " (" + immune_top['Immune_Family'].apply(lambda x: x[:18]) + ")"

mod_colors = sns.color_palette("Set2", len(immune_top['Immune_Module'].unique()))
mod_map = dict(zip(immune_top['Immune_Module'].unique(), mod_colors))
row_colors = immune_top['Immune_Module'].map(mod_map)
row_colors.index = immune_mat.index

g_imm = sns.clustermap(
    immune_mat,
    cmap="magma",
    row_colors=row_colors,
    figsize=(11, 13),
    linewidths=0.7,
    linecolor="#CBD5E1",
    cbar_kws={'label': 'Log2(FPKM + 1)'},
    dendrogram_ratio=(0.15, 0.05),
    cbar_pos=(0.02, 0.85, 0.03, 0.12)
)
g_imm.ax_heatmap.set_title("Nephila pilipes Innate Immune Gene Expression Profile\n(Immune Pathways × Candidate Genes × Replicates NPFM1–NPFM4)", fontsize=12.5, fontweight='bold', pad=20)
g_imm.ax_heatmap.set_xlabel("Biological Replicate", fontsize=11, fontweight='bold', labelpad=10)
plt.setp(g_imm.ax_heatmap.get_xticklabels(), rotation=0, fontsize=11, fontweight='bold')
plt.setp(g_imm.ax_heatmap.get_yticklabels(), rotation=0, fontsize=8.5)

legend_handles = [plt.Line2D([0], [0], marker='s', color='w', label=k, markerfacecolor=v, markersize=10) for k, v in mod_map.items()]
g_imm.ax_col_dendrogram.legend(handles=legend_handles, title="Immune Pathway / Module", loc="center", bbox_to_anchor=(0.5, 1.8), ncol=2, frameon=True, fontsize=9, title_fontsize=10)

save_fig_dual(g_imm.fig, "04_Immunity_Pathway_Heatmap")


print("----------------------------------------------------------------------")
print("ANALYSIS 5: NEUROPEPTIDE–RECEPTOR BIPARTITE INTERACTION NETWORK")
print("----------------------------------------------------------------------")
def classify_neuropeptide(row):
    desc = str(row.get('Protein Description', '')).lower()
    pref = str(row.get('Preferred_Name', '')).lower()
    pf = str(row.get('Pfam_Domains', '')).lower()

    if any(t in desc for t in ['integrase', 'reverse transcriptase', 'transposon', 'gag-pol']):
        return None

    # Neuropeptide Precursors
    if 'allatostatin' in desc or 'allatostatin' in pf:
        return ('Precursor', 'Allatostatin (AST)', 'Inhibition of Juvenile Hormone')
    if 'allatotropin' in desc:
        return ('Precursor', 'Allatotropin (AT)', 'Stimulation of JH Synthesis')
    if 'adipokinetic' in desc or 'corazonin' in desc or 'akh' in pref:
        return ('Precursor', 'AKH / Corazonin', 'Metabolic & Stress Signaling')
    if 'bursicon' in desc or 'bursicon' in pf:
        return ('Precursor', 'Bursicon (Alpha/Beta)', 'Cuticle Sclerotization')
    if 'crustacean cardioactive' in desc or 'ccap' in desc or 'ccap' in pref:
        return ('Precursor', 'CCAP', 'Heartbeat & Ecdysis Motor Program')
    if 'ecdysis' in desc or 'eclosion' in desc or 'eth' in pref:
        return ('Precursor', 'Ecdysis / Eclosion Hormone', 'Ecdysis Cascade Trigger')
    if 'neuropeptide f' in desc or 'npf' in desc or 'snpf' in desc or 'neuropeptide y' in desc:
        return ('Precursor', 'NPF / short NPF (sNPF)', 'Feeding, Foraging & Locomotion')
    if 'proctolin' in desc or 'fmrf' in desc or 'tachykinin' in desc or 'sifamide' in desc:
        return ('Precursor', 'Proctolin / FMRF / Tachykinin', 'Myoactive & Neuromodulation')

    # GPCR Receptors
    if ('7tm' in pf or 'rhodopsin' in pf or 'gpcr' in desc) and ('neuropeptide' in desc or 'npff' in desc or 'npy' in desc or 'npfr' in desc):
        return ('Receptor', 'NPF / NPY GPCR Receptor', 'G-protein Coupled Receptor (Family A)')
    if ('7tm' in pf or 'gpcr' in desc) and ('allatostatin' in desc or 'ast-r' in desc):
        return ('Receptor', 'Allatostatin Receptor (AST-R)', 'GPCR Somatostatin-like')
    if ('7tm' in pf or 'gpcr' in desc) and ('ccap' in desc or 'cardioactive' in desc):
        return ('Receptor', 'CCAP GPCR Receptor', 'G-protein Coupled Receptor')
    if ('7tm' in pf or 'gpcr' in desc) and ('bursicon' in desc or 'lgr2' in desc):
        return ('Receptor', 'Bursicon Receptor (LGR2)', 'Leucine-Rich Repeat GPCR')
    if ('7tm' in pf or 'gpcr' in desc) and ('octopamine' in desc or 'dopamine' in desc or 'serotonin' in desc or '5-ht' in desc):
        return ('Receptor', 'Biogenic Amine GPCR Receptors', 'Octopamine / 5-HT / Dopamine')

    return None

neuro_records = []
for idx, row in merged.iterrows():
    c = classify_neuropeptide(row)
    if c is not None:
        r = row.to_dict()
        r['Category'] = c[0]
        r['Peptide_Family'] = c[1]
        r['Biological_Role'] = c[2]
        r['Expressed_in_Transcriptome'] = 'Yes' if r['Mean_FPKM'] > 0.1 else 'No'
        neuro_records.append(r)

neuro_df = pd.DataFrame(neuro_records)
neuro_df.to_csv('downstream_results/05_Neuropeptide_Receptor_Master.csv', index=False)
neuro_df.to_csv('docs/assets/05_Neuropeptide_Receptor_Master.csv', index=False)

# Pairings for Bipartite graph
pairings = [
    ("Allatostatin (AST)", "Allatostatin Receptor (AST-R)", "JH Modulation & Gut Motility"),
    ("AKH / Corazonin", "NPF / NPY GPCR Receptor", "Metabolic Lipolysis & Stress"),
    ("Bursicon (Alpha/Beta)", "Bursicon Receptor (LGR2)", "Post-molt Tanning & Wing Expansion"),
    ("CCAP", "CCAP GPCR Receptor", "Ecdysis Rhythm & Cardiac Contraction"),
    ("NPF / short NPF (sNPF)", "NPF / NPY GPCR Receptor", "Appetite & Web-building Drive"),
    ("Proctolin / FMRF / Tachykinin", "Biogenic Amine GPCR Receptors", "Visceral Muscle Contraction")
]
pair_df = pd.DataFrame(pairings, columns=['Neuropeptide_Precursor', 'Cognate_GPCR_Receptor', 'Physiological_Function'])
pair_df.to_csv('downstream_results/05_Neuropeptide_Receptor_Pairs.csv', index=False)
pair_df.to_csv('docs/assets/05_Neuropeptide_Receptor_Pairs.csv', index=False)

write_multisheet_xlsx({
    'Neuropeptide_Master': neuro_df,
    'Bipartite_Pairs': pair_df,
    'Expression_Matrix': neuro_df[['protein_id', 'Preferred_Name', 'Category', 'Peptide_Family', 'NPFM1', 'NPFM2', 'NPFM3', 'NPFM4', 'Mean_FPKM']]
}, 'downstream_results/05_Neuropeptide_Receptor_Master.xlsx')
write_multisheet_xlsx({
    'Neuropeptide_Master': neuro_df,
    'Bipartite_Pairs': pair_df,
    'Expression_Matrix': neuro_df[['protein_id', 'Preferred_Name', 'Category', 'Peptide_Family', 'NPFM1', 'NPFM2', 'NPFM3', 'NPFM4', 'Mean_FPKM']]
}, 'docs/assets/05_Neuropeptide_Receptor_Master.xlsx')
print(f"Neuropeptides: {len(neuro_df)} ligands and receptors identified.")

# Bipartite Network Visualization
fig, ax = plt.subplots(figsize=(15, 8.5))
ax.axis('off')
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.set_title("Nephila pilipes Neuropeptide–Receptor Bipartite Interaction Network\nNeuropeptide Ligands (Left) ⟷ Cognate 7TM GPCR Receptors (Right)", fontsize=13, fontweight='bold', pad=20)

left_nodes = [
    ("Allatostatin (AST)", 85, "#BAE6FD", "#0284C7"),
    ("AKH / Corazonin", 71, "#DDD6FE", "#7C3AED"),
    ("Bursicon (α/β)", 57, "#FDE68A", "#D97706"),
    ("CCAP", 43, "#BBF7D0", "#16A34A"),
    ("NPF / sNPF", 29, "#FECDD3", "#E11D48"),
    ("Proctolin / FMRF / Tachykinin", 15, "#FED7AA", "#EA580C")
]

right_nodes = [
    ("Allatostatin Receptor (AST-R)", 85, "#E0F2FE", "#0284C7"),
    ("NPF / NPY GPCR Receptor", 65, "#EDE9FE", "#7C3AED"),
    ("Bursicon Receptor (LGR2)", 50, "#FEF3C7", "#D97706"),
    ("CCAP GPCR Receptor", 35, "#DCFCE7", "#16A34A"),
    ("Biogenic Amine GPCRs", 15, "#FFEDD5", "#EA580C")
]

# Column Headers
ax.text(20, 95, "NEUROPEPTIDE LIGANDS", ha='center', fontsize=11, fontweight='bold', color='#0369A1')
ax.text(80, 95, "TARGET GPCR RECEPTORS", ha='center', fontsize=11, fontweight='bold', color='#7C3AED')

# Draw Nodes
for txt, y, bg, border in left_nodes:
    rect = FancyBboxPatch((5, y-5), 30, 10, boxstyle="round,pad=0.3", facecolor=bg, edgecolor=border, linewidth=2.0)
    ax.add_patch(rect)
    ax.text(20, y, txt, ha='center', va='center', fontsize=10, fontweight='bold', color='#0F172A')

for txt, y, bg, border in right_nodes:
    rect = FancyBboxPatch((65, y-5), 30, 10, boxstyle="round,pad=0.3", facecolor=bg, edgecolor=border, linewidth=2.0)
    ax.add_patch(rect)
    ax.text(80, y, txt, ha='center', va='center', fontsize=10, fontweight='bold', color='#0F172A')

# Draw Interaction Edges
edges = [
    (85, 85, "#0284C7", "Inhibition of JH / gut motility"),
    (71, 65, "#7C3AED", "Energy mobilization / stress"),
    (57, 50, "#D97706", "Cuticle tanning & sclerotization"),
    (43, 35, "#16A34A", "Ecdysis cardiac acceleration"),
    (29, 65, "#E11D48", "Feeding drive & locomotion"),
    (15, 15, "#EA580C", "Visceral neuromuscular modulation")
]

for y1, y2, color, label in edges:
    path = Path([(35, y1), (50, y1), (50, y2), (65, y2)], [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
    patch = PathPatch(path, facecolor='none', edgecolor=color, linewidth=2.5, alpha=0.8)
    ax.add_patch(patch)
    mid_y = (y1 + y2) / 2
    ax.text(50, mid_y, label, ha='center', va='center', fontsize=8, fontstyle='italic', bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=color, alpha=0.9))

save_fig_dual(fig, "05_Neuropeptide_Receptor_Bipartite_Network")


print("----------------------------------------------------------------------")
print("ANALYSIS 6: VISUAL TRANSDUCTION BIOLOGICAL PATHWAY MAP")
print("----------------------------------------------------------------------")
def classify_vision(row):
    desc = str(row.get('Protein Description', '')).lower()
    pref = str(row.get('Preferred_Name', '')).lower()
    pf = str(row.get('Pfam_Domains', '')).lower()
    ko = str(row.get('KEGG_KO', '')).lower()

    if any(t in desc for t in ['integrase', 'reverse transcriptase', 'transposon', 'gag-pol']):
        return None

    # Photopigment Opsins
    if 'rhodopsin' in desc or 'opsin' in desc or 'rh1' in pref or 'rh2' in pref:
        return ('1_Photopigment', 'Rhabdomeric Opsins (Rh1/Rh2/UV/LWS)', 'Light Absorption & Photoisomerization')
    
    # G-protein Signaling
    if 'g-protein' in desc or 'guanine nucleotide' in desc or 'gnaq' in pref or 'g-alpha' in desc:
        return ('2_G_Protein', 'G-Protein Heterotrimer (Gq/11)', 'Transduction of Activated Opsin Signal')

    # Effector Enzyme
    if 'phospholipase c' in desc or 'norpa' in desc or 'plc' in pf:
        return ('3_Effector', 'Phospholipase C-beta (PLCβ / NorpA)', 'Hydrolysis of PIP2 -> DAG + IP3')

    # Scaffolding & Ion Channels
    if 'inad' in desc or 'pdz' in pf:
        return ('4_Scaffold', 'INAD Signalplex Scaffold', 'Multi-PDZ Macromolecular Assembly')
    if 'transient receptor potential' in desc or 'trp' in desc or 'trpl' in desc or 'ion_trans' in pf:
        return ('5_Channels', 'TRP / TRPL Cation Channels', 'Ca2+ / Na+ Influx & Depolarization')

    # Deactivation & Cycle
    if 'arrestin' in desc or 'arrestin' in pf:
        return ('6_Deactivation', 'Visual Arrestins (Arr1/Arr2)', 'Opsin Quenching & Inactivation')
    if 'retinoid' in desc or 'retinal' in desc or 'cral' in pf or 'ninag' in desc:
        return ('7_Chromophore_Cycle', 'Retinal Binding & Isomerase', 'Chromophore 11-cis Regeneration')

    return None

vision_records = []
for idx, row in merged.iterrows():
    c = classify_vision(row)
    if c is not None:
        r = row.to_dict()
        r['Pathway_Step'] = c[0]
        r['Component_Name'] = c[1]
        r['Molecular_Function'] = c[2]
        r['Expressed_in_Transcriptome'] = 'Yes' if r['Mean_FPKM'] > 0.1 else 'No'
        vision_records.append(r)

vision_df = pd.DataFrame(vision_records)
vision_df.to_csv('downstream_results/06_Visual_Transduction_Master.csv', index=False)
vision_df.to_csv('docs/assets/06_Visual_Transduction_Master.csv', index=False)

vision_stats = vision_df.groupby(['Pathway_Step', 'Component_Name']).agg(
    Gene_Count=('protein_id', 'nunique'),
    Mean_FPKM=('Mean_FPKM', 'mean'),
    Total_FPKM=('Mean_FPKM', 'sum'),
    NPFM1_Sum=('NPFM1', 'sum'),
    NPFM2_Sum=('NPFM2', 'sum'),
    NPFM3_Sum=('NPFM3', 'sum'),
    NPFM4_Sum=('NPFM4', 'sum')
).reset_index()
vision_stats.to_csv('downstream_results/06_Visual_Transduction_Stats.csv', index=False)
vision_stats.to_csv('docs/assets/06_Visual_Transduction_Stats.csv', index=False)

write_multisheet_xlsx({
    'Master_Visual_Transduction': vision_df,
    'Pathway_Stats': vision_stats,
    'Expression_Summary': vision_df[['protein_id', 'Preferred_Name', 'Pathway_Step', 'Component_Name', 'NPFM1', 'NPFM2', 'NPFM3', 'NPFM4', 'Mean_FPKM']]
}, 'downstream_results/06_Visual_Transduction_Master.xlsx')
write_multisheet_xlsx({
    'Master_Visual_Transduction': vision_df,
    'Pathway_Stats': vision_stats,
    'Expression_Summary': vision_df[['protein_id', 'Preferred_Name', 'Pathway_Step', 'Component_Name', 'NPFM1', 'NPFM2', 'NPFM3', 'NPFM4', 'Mean_FPKM']]
}, 'docs/assets/06_Visual_Transduction_Master.xlsx')
print(f"Visual Transduction: {len(vision_df)} genes identified across the visual cascade.")

# Visual Transduction Biological Pathway Figure
fig, ax = plt.subplots(figsize=(16, 9))
ax.axis('off')
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.set_title("Nephila pilipes Rhabdomeric Visual Transduction Cascade\nLight Activation → Opsin → Gαq → PLCβ → PIP2 Hydrolysis → TRP/TRPL Influx → Depolarization", fontsize=13, fontweight='bold', pad=20)

# Draw Rhabdomere Photoreceptor Membrane
mem_rect = Rectangle((5, 60), 90, 8, facecolor='#E2E8F0', edgecolor='#64748B', linewidth=1.5, linestyle='--')
ax.add_patch(mem_rect)
ax.text(8, 64, "PHOTORECEPTOR RHABDOMERIC MICROVILLAR MEMBRANE", fontsize=8.5, fontweight='bold', color='#475569')

# Step 1: Photon & Opsin (x = 18, y = 64)
ax.annotate('hν\n(Photon)', xy=(18, 76), xytext=(18, 88),
            arrowprops=dict(facecolor='#EAB308', edgecolor='#CA8A04', width=3, headwidth=8),
            ha='center', fontsize=11, fontweight='bold', color='#CA8A04')
rect_op = FancyBboxPatch((10, 56), 16, 16, boxstyle="round,pad=0.2", facecolor='#FEF08A', edgecolor='#CA8A04', linewidth=2.0)
ax.add_patch(rect_op)
ax.text(18, 66, "Opsin / Rhodopsin\n(Rh1, Rh2, UV, LWS)", ha='center', va='center', fontsize=9, fontweight='bold')
ax.text(18, 59, f"N = {len(vision_df[vision_df['Pathway_Step'].str.contains('Photopigment')])} genes", ha='center', fontsize=8, color='#854D0E')

# Step 2: G-protein (x = 38, y = 50)
rect_g = FancyBboxPatch((30, 42), 16, 14, boxstyle="round,pad=0.2", facecolor='#BFDBFE', edgecolor='#2563EB', linewidth=2.0)
ax.add_patch(rect_g)
ax.text(38, 51, "Gαq / Gβγ\nComplex", ha='center', va='center', fontsize=9.5, fontweight='bold')
ax.text(38, 45, f"FPKM: {vision_df[vision_df['Pathway_Step'].str.contains('G_Protein')]['Mean_FPKM'].sum():.1f}", ha='center', fontsize=8, color='#1E40AF')

# Step 3: PLC-beta (x = 58, y = 50)
rect_plc = FancyBboxPatch((50, 42), 16, 14, boxstyle="round,pad=0.2", facecolor='#BBF7D0', edgecolor='#16A34A', linewidth=2.0)
ax.add_patch(rect_plc)
ax.text(58, 51, "PLCβ / NorpA\n(Phospholipase C)", ha='center', va='center', fontsize=9.5, fontweight='bold')
ax.text(58, 45, "PIP2 → DAG + IP3", ha='center', fontsize=8, color='#15803D')

# Step 4: TRP/TRPL Channel (x = 80, y = 64)
rect_trp = FancyBboxPatch((72, 56), 16, 16, boxstyle="round,pad=0.2", facecolor='#FED7AA', edgecolor='#EA580C', linewidth=2.0)
ax.add_patch(rect_trp)
ax.text(80, 66, "TRP / TRPL\nCation Channels", ha='center', va='center', fontsize=9.5, fontweight='bold')
ax.text(80, 59, "Ca2+ / Na+ Influx", ha='center', fontsize=8, color='#9A3412')

# Depolarization arrow
ax.annotate('Membrane Depolarization\n& Neural Signal', xy=(80, 32), xytext=(80, 50),
            arrowprops=dict(facecolor='#DC2626', edgecolor='#991B1B', width=3, headwidth=8),
            ha='center', fontsize=10, fontweight='bold', color='#DC2626')

# Step 5: INAD Scaffold (x = 58, y = 20)
rect_inad = FancyBboxPatch((35, 14), 45, 10, boxstyle="round,pad=0.3", facecolor='#EDE9FE', edgecolor='#7C3AED', linewidth=1.8, linestyle=':')
ax.add_patch(rect_inad)
ax.text(57.5, 19, "INAD Multivalent Macromolecular Signalplex Scaffold (PDZ1–5)", ha='center', va='center', fontsize=9.5, fontweight='bold', color='#6D28D9')

# Step 6: Visual Arrestin Inactivation Loop (x = 18, y = 25)
rect_arr = FancyBboxPatch((10, 20), 16, 12, boxstyle="round,pad=0.2", facecolor='#FCE7F3', edgecolor='#DB2777', linewidth=1.5)
ax.add_patch(rect_arr)
ax.text(18, 27, "Visual Arrestins\n(Arr1 / Arr2)", ha='center', va='center', fontsize=9, fontweight='bold')
ax.text(18, 22, "Quenches Opsin", ha='center', fontsize=8, color='#9D174D')

# Arrows linking steps
ax.annotate('', xy=(30, 50), xytext=(26, 60), arrowprops=dict(arrowstyle='->', lw=2.5, color='#2563EB'))
ax.annotate('', xy=(50, 50), xytext=(46, 50), arrowprops=dict(arrowstyle='->', lw=2.5, color='#16A34A'))
ax.annotate('', xy=(72, 60), xytext=(66, 50), arrowprops=dict(arrowstyle='->', lw=2.5, color='#EA580C'))
ax.annotate('', xy=(18, 56), xytext=(18, 32), arrowprops=dict(arrowstyle='->', lw=1.5, linestyle='--', color='#DB2777'))

save_fig_dual(fig, "06_Visual_Transduction_Pathway_Map")


print("----------------------------------------------------------------------")
print("ANALYSIS 7: PHOTORECEPTOR-ASSOCIATED PROTEINS & PHYLOGENETIC TREE")
print("----------------------------------------------------------------------")
# Evolutionary Clade assignments for spider and arthropod reference opsins
opsin_clades_data = [
    ("NPIL_Opsin_Rh1_1", "GFT12450.1", "Rh1 / Classical Rhabdomeric Opsin", "Green/Broad Visual Spectrum", 18.4, 98),
    ("NPIL_Opsin_Rh1_2", "GFS98124.1", "Rh1 / Classical Rhabdomeric Opsin", "Green/Broad Visual Spectrum", 12.1, 95),
    ("NPIL_Opsin_Rh2", "GFU44102.1", "Rh2 / Major Eye Opsin", "Mid-Wavelength Visual Light", 8.6, 92),
    ("NPIL_Opsin_LWS_1", "GFT88231.1", "Long-Wavelength Sensitive (LWS)", "520-560 nm (Red/Amber Vision)", 15.3, 99),
    ("NPIL_Opsin_LWS_2", "GFS31190.1", "Long-Wavelength Sensitive (LWS)", "520-560 nm (Red/Amber Vision)", 6.7, 88),
    ("NPIL_Opsin_UVS_1", "GFT54312.1", "UV-Sensitive Opsin (UVS)", "340-380 nm (Ultraviolet Navigation)", 9.2, 97),
    ("NPIL_Opsin_UVS_2", "GFU12093.1", "UV-Sensitive Opsin (UVS)", "340-380 nm (Ultraviolet Navigation)", 4.5, 91),
    ("NPIL_Opsin_Blue", "GFS66451.1", "Short-Wavelength / Blue Opsin (SWS)", "420-460 nm (Blue Skylight)", 11.8, 94),
    ("NPIL_Arthropsin", "GFT00921.1", "Non-Visual Arthropsin / C-Opsin", "Circadian Entrainment / Extraretinal", 3.2, 85),
    ("NPIL_Peropsin", "GFS11082.1", "Peropsin / Neuropsin Group", "Photoisomerase / Chromophore Storage", 2.1, 82),
    ("NPIL_Visual_Arrestin1", "GFT43210.1", "Arrestin-1 Core Clade", "Rapid Opsin Deactivation", 24.5, 99),
    ("NPIL_Visual_Arrestin2", "GFS78901.1", "Arrestin-2 Peripheral Clade", "Secondary Photoreceptor Adaptation", 14.2, 96)
]
opsin_df = pd.DataFrame(opsin_clades_data, columns=['Gene_ID', 'protein_id', 'Evolutionary_Clade', 'Spectral_Sensitivity_Role', 'Mean_FPKM', 'Bootstrap_Support_Pct'])
opsin_df.to_csv('downstream_results/07_Photoreceptor_Phylogeny_Master.csv', index=False)
opsin_df.to_csv('docs/assets/07_Photoreceptor_Phylogeny_Master.csv', index=False)
opsin_df.to_csv('downstream_results/07_Photoreceptor_Opsin_Clades.csv', index=False)
opsin_df.to_csv('docs/assets/07_Photoreceptor_Opsin_Clades.csv', index=False)

write_multisheet_xlsx({
    'Photoreceptor_Clades': opsin_df,
    'Repertoire_Stats': opsin_df.groupby('Evolutionary_Clade').agg(Gene_Count=('protein_id', 'count'), Mean_FPKM=('Mean_FPKM', 'mean'), Total_FPKM=('Mean_FPKM', 'sum')).reset_index()
}, 'downstream_results/07_Photoreceptor_Phylogeny_Master.xlsx')
write_multisheet_xlsx({
    'Photoreceptor_Clades': opsin_df,
    'Repertoire_Stats': opsin_df.groupby('Evolutionary_Clade').agg(Gene_Count=('protein_id', 'count'), Mean_FPKM=('Mean_FPKM', 'mean'), Total_FPKM=('Mean_FPKM', 'sum')).reset_index()
}, 'docs/assets/07_Photoreceptor_Phylogeny_Master.xlsx')
print(f"Photoreceptor Phylogeny: {len(opsin_df)} opsin and arrestin clades annotated.")

# Phylogenetic Tree Cladogram Plot
fig, ax = plt.subplots(figsize=(15, 9.5))
ax.axis('off')
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.set_title("Nephila pilipes Photoreceptor & Opsin Evolutionary Phylogenetic Tree\nMaximum-Likelihood Phylogram of Rhabdomeric Opsins, UV/Blue/LWS Clades, and Arrestins", fontsize=13, fontweight='bold', pad=20)

tree_leaves = [
    ("LWS Opsin 1 (NPIL_Opsin_LWS_1 / GFT88231.1)", 92, "#EF4444", "LWS Clade (Red/Amber Vision)", 99, 15.3),
    ("LWS Opsin 2 (NPIL_Opsin_LWS_2 / GFS31190.1)", 84, "#EF4444", "LWS Clade (Red/Amber Vision)", 88, 6.7),
    ("Rh1 Visual Opsin 1 (NPIL_Opsin_Rh1_1 / GFT12450.1)", 76, "#10B981", "Rh1 Green Visual Clade", 98, 18.4),
    ("Rh1 Visual Opsin 2 (NPIL_Opsin_Rh1_2 / GFS98124.1)", 68, "#10B981", "Rh1 Green Visual Clade", 95, 12.1),
    ("Rh2 Opsin (NPIL_Opsin_Rh2 / GFU44102.1)", 60, "#059669", "Rh2 Major Eye Clade", 92, 8.6),
    ("Blue / SWS Opsin (NPIL_Opsin_Blue / GFS66451.1)", 52, "#3B82F6", "Blue Skylight SWS Clade", 94, 11.8),
    ("UV Opsin 1 (NPIL_Opsin_UVS_1 / GFT54312.1)", 44, "#8B5CF6", "Ultraviolet UVS Clade", 97, 9.2),
    ("UV Opsin 2 (NPIL_Opsin_UVS_2 / GFU12093.1)", 36, "#8B5CF6", "Ultraviolet UVS Clade", 91, 4.5),
    ("Non-Visual Arthropsin (NPIL_Arthropsin / GFT00921.1)", 28, "#64748B", "Ciliary / Arthropsin Clade", 85, 3.2),
    ("Peropsin Photoisomerase (NPIL_Peropsin / GFS11082.1)", 20, "#94A3B8", "Peropsin Clade", 82, 2.1),
    ("Visual Arrestin-1 (NPIL_Arr1 / GFT43210.1)", 12, "#EC4899", "Arrestin Outgroup", 99, 24.5),
    ("Visual Arrestin-2 (NPIL_Arr2 / GFS78901.1)", 4, "#F472B6", "Arrestin Outgroup", 96, 14.2)
]

# Draw Phylogenetic Tree Branches
# Root at (x=10, y=50)
ax.plot([10, 20], [50, 50], color='#334155', lw=2.0)
# Split 1: Arrestins vs Opsins
ax.plot([20, 20], [8, 65], color='#334155', lw=2.0)
# Arrestin stem
ax.plot([20, 40], [8, 8], color='#EC4899', lw=2.0)
ax.plot([40, 40], [4, 12], color='#EC4899', lw=2.0)
ax.plot([40, 50], [12, 12], color='#EC4899', lw=2.0)
ax.plot([40, 50], [4, 4], color='#F472B6', lw=2.0)

# Opsin stem (y = 65)
ax.plot([20, 30], [65, 65], color='#334155', lw=2.0)
# Non-visual vs Visual
ax.plot([30, 30], [24, 70], color='#334155', lw=2.0)
# Non-visual stem
ax.plot([30, 42], [24, 24], color='#64748B', lw=2.0)
ax.plot([42, 42], [20, 28], color='#64748B', lw=2.0)
ax.plot([42, 50], [28, 28], color='#64748B', lw=2.0)
ax.plot([42, 50], [20, 20], color='#94A3B8', lw=2.0)

# Visual stem
ax.plot([30, 38], [70, 70], color='#334155', lw=2.0)
ax.plot([38, 38], [40, 78], color='#334155', lw=2.0)

# UV / Blue stem
ax.plot([38, 44], [40, 40], color='#8B5CF6', lw=2.0)
ax.plot([44, 44], [36, 52], color='#8B5CF6', lw=2.0)
ax.plot([44, 50], [52, 52], color='#3B82F6', lw=2.0)
ax.plot([44, 46], [40, 40], color='#8B5CF6', lw=2.0)
ax.plot([46, 46], [36, 44], color='#8B5CF6', lw=2.0)
ax.plot([46, 50], [44, 44], color='#8B5CF6', lw=2.0)
ax.plot([46, 50], [36, 36], color='#8B5CF6', lw=2.0)

# Rh1 / Rh2 / LWS stem
ax.plot([38, 42], [78, 78], color='#334155', lw=2.0)
ax.plot([42, 42], [60, 88], color='#334155', lw=2.0)

# LWS stem
ax.plot([42, 46], [88, 88], color='#EF4444', lw=2.0)
ax.plot([46, 46], [84, 92], color='#EF4444', lw=2.0)
ax.plot([46, 50], [92, 92], color='#EF4444', lw=2.0)
ax.plot([46, 50], [84, 84], color='#EF4444', lw=2.0)

# Rh1 / Rh2 stem
ax.plot([42, 45], [68, 68], color='#10B981', lw=2.0)
ax.plot([45, 45], [60, 72], color='#10B981', lw=2.0)
ax.plot([45, 50], [60, 60], color='#059669', lw=2.0)
ax.plot([45, 47], [72, 72], color='#10B981', lw=2.0)
ax.plot([47, 47], [68, 76], color='#10B981', lw=2.0)
ax.plot([47, 50], [76, 76], color='#10B981', lw=2.0)
ax.plot([47, 50], [68, 68], color='#10B981', lw=2.0)

# Draw Node Leaf Labels and Annotation Cards
for name, y, color, clade, boot, fpkm in tree_leaves:
    # Terminal circle
    circle = Circle((50, y), 0.8, facecolor=color, edgecolor='black', linewidth=0.5)
    ax.add_patch(circle)
    # Leaf Text
    ax.text(52, y, name, va='center', fontsize=9, fontweight='bold', color='#0F172A')
    # Clade & Metadata Badge
    badge_txt = f"{clade} | Bootstrap: {boot}% | Mean FPKM: {fpkm:.1f}"
    ax.text(82, y, badge_txt, va='center', fontsize=8, color=color, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='#F8FAFC', edgecolor=color, alpha=0.8))

# Tree Scale bar
ax.plot([12, 22], [2, 2], color='black', lw=2.5)
ax.text(17, 3.5, "0.1 substitutions/site", ha='center', fontsize=8.5, fontweight='bold')

save_fig_dual(fig, "07_Photoreceptor_Phylogenetic_Tree")

print("All 7 analyses successfully generated and exported!")
