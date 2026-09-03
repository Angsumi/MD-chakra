import pandas as pd
import numpy as np
import plotly.graph_objects as go
import os

print("Generating pure classic Boxplot and standalone Violin plot...")

fpkm_df = pd.read_csv('downstream_results/6_FPKM_normalized_counts_individual.csv')
fpkm_df.rename(columns={'Unnamed: 0': 'Gene ID'}, inplace=True)
samples = ['NPFM1', 'NPFM2', 'NPFM3', 'NPFM4']

log2_fpkm = fpkm_df.copy()
for s in samples:
    log2_fpkm[f'log2_{s}'] = np.log2(log2_fpkm[s] + 1)

melted_log2 = log2_fpkm.melt(
    id_vars=['Gene ID'], 
    value_vars=[f'log2_{s}' for s in samples],
    var_name='Sample_col', 
    value_name='Log2_FPKM'
)
melted_log2['Sample'] = melted_log2['Sample_col'].str.replace('log2_', '')
melted_expressed = melted_log2[melted_log2['Log2_FPKM'] > 0.1].copy()

palette = {
    'NPFM1': '#6366f1', # Indigo
    'NPFM2': '#06b6d4', # Cyan
    'NPFM3': '#10b981', # Emerald
    'NPFM4': '#f59e0b'  # Amber
}

# -------------------------------------------------------------------------
# 1. Pure Classic Publication Boxplot (Tukey style)
# -------------------------------------------------------------------------
fig_pure_box = go.Figure()

for s in samples:
    sub = melted_expressed[melted_expressed['Sample'] == s]['Log2_FPKM']
    fig_pure_box.add_trace(go.Box(
        y=sub,
        name=s,
        marker_color=palette[s],
        fillcolor=f"rgba{tuple(list(int(palette[s].lstrip('#')[i:i+2], 16) for i in (0, 2, 4)) + [0.35])}",
        line=dict(color=palette[s], width=2),
        boxmean=True, # Shows both mean (dashed) and median (solid)
        boxpoints=False # Outliers suppressed for clean aesthetic
    ))

fig_pure_box.update_layout(
    title=dict(
        text="<b>Log2-Expression Boxplot (NPFM1 - NPFM4)</b><br><sup>Classic boxplot showing median (solid line), mean (dashed line), IQR (25th-75th), and whiskers</sup>",
        font=dict(size=20, color="#0f172a", family="Plus Jakarta Sans, sans-serif"),
        x=0.5
    ),
    xaxis=dict(title="<b>Sample</b>", showline=True, linecolor="#cbd5e1", linewidth=1.5),
    yaxis=dict(title="<b>log2(FPKM + 1)</b>", showline=True, linecolor="#cbd5e1", linewidth=1.5, showgrid=True, gridcolor="#f1f5f9"),
    paper_bgcolor="white",
    plot_bgcolor="white",
    width=900,
    height=600,
    showlegend=False,
    margin=dict(l=70, r=40, t=90, b=60)
)

fig_pure_box.write_html("new plots/html/plotly_expression_boxplot_pure.html")
fig_pure_box.write_image("visualizations/12_Expression_Distribution_Boxplot.png", scale=3)
fig_pure_box.write_image("visualizations/12_Expression_Distribution_Boxplot.pdf")
fig_pure_box.write_image("docs/assets/12_Expression_Distribution_Boxplot.png", scale=3)
fig_pure_box.write_image("docs/assets/12_Expression_Distribution_Boxplot.pdf")
print("Saved Pure Classic Boxplot to 12_Expression_Distribution_Boxplot")

# -------------------------------------------------------------------------
# 2. Standalone Violin Plot
# -------------------------------------------------------------------------
fig_violin = go.Figure()
for s in samples:
    sub = melted_expressed[melted_expressed['Sample'] == s]['Log2_FPKM']
    fig_violin.add_trace(go.Violin(
        x=[s]*len(sub),
        y=sub,
        name=s,
        box_visible=False,
        meanline_visible=True,
        fillcolor=f"rgba{tuple(list(int(palette[s].lstrip('#')[i:i+2], 16) for i in (0, 2, 4)) + [0.45])}",
        line_color=palette[s],
        points=False
    ))

fig_violin.update_layout(
    title=dict(
        text="<b>Log2-Expression Violin Distributions</b><br><sup>Violin kernel density profiles showing continuous expression distribution per sample</sup>",
        font=dict(size=20, color="#0f172a", family="Plus Jakarta Sans, sans-serif"),
        x=0.5
    ),
    xaxis=dict(title="<b>Sample</b>", showline=True, linecolor="#cbd5e1", linewidth=1.5),
    yaxis=dict(title="<b>log2(FPKM + 1)</b>", showline=True, linecolor="#cbd5e1", linewidth=1.5, showgrid=True, gridcolor="#f1f5f9"),
    paper_bgcolor="white",
    plot_bgcolor="white",
    width=900,
    height=600,
    showlegend=False,
    margin=dict(l=70, r=40, t=90, b=60)
)

fig_violin.write_html("new plots/html/plotly_expression_violin.html")
fig_violin.write_image("visualizations/12b_Expression_Distribution_Violin.png", scale=3)
fig_violin.write_image("visualizations/12b_Expression_Distribution_Violin.pdf")
fig_violin.write_image("docs/assets/12b_Expression_Distribution_Violin.png", scale=3)
fig_violin.write_image("docs/assets/12b_Expression_Distribution_Violin.pdf")
print("Saved Standalone Violin Plot to 12b_Expression_Distribution_Violin")
