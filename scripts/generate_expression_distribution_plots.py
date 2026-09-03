import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from scipy import stats
import os

print("Generating Expression Distribution Plots...")

os.makedirs("visualizations", exist_ok=True)
os.makedirs("new plots/html", exist_ok=True)
os.makedirs("new plots/pdf", exist_ok=True)
os.makedirs("new plots/png", exist_ok=True)
os.makedirs("docs/assets", exist_ok=True)

# 1. Load FPKM
fpkm_df = pd.read_csv('downstream_results/6_FPKM_normalized_counts_individual.csv')
fpkm_df.rename(columns={'Unnamed: 0': 'Gene ID'}, inplace=True)
samples = ['NPFM1', 'NPFM2', 'NPFM3', 'NPFM4']

# Filter expressed genes (at least one sample > 0.1 FPKM) for distribution plots, or all
# In RNA-seq QC, distribution of expressed genes (log2(FPKM + 1)) is standard:
log2_fpkm = fpkm_df.copy()
for s in samples:
    log2_fpkm[f'log2_{s}'] = np.log2(log2_fpkm[s] + 1)

# Melt for plotting
melted_log2 = log2_fpkm.melt(
    id_vars=['Gene ID'], 
    value_vars=[f'log2_{s}' for s in samples],
    var_name='Sample_col', 
    value_name='Log2_FPKM'
)
melted_log2['Sample'] = melted_log2['Sample_col'].str.replace('log2_', '')

# Filter out unexpressed (log2_FPKM > 0.1) for boxplot to show active distribution clearly
melted_expressed = melted_log2[melted_log2['Log2_FPKM'] > 0.1].copy()

palette = {
    'NPFM1': '#6366f1', # Indigo
    'NPFM2': '#06b6d4', # Cyan
    'NPFM3': '#10b981', # Emerald
    'NPFM4': '#f59e0b'  # Amber
}

# -------------------------------------------------------------------------
# Plot 1: Boxplots & Violin Plots (Log2 Expression Distributions)
# -------------------------------------------------------------------------
fig_box = go.Figure()
for s in samples:
    sub = melted_expressed[melted_expressed['Sample'] == s]['Log2_FPKM']
    fig_box.add_trace(go.Violin(
        x=[s]*len(sub),
        y=sub,
        name=s,
        box_visible=True,
        meanline_visible=True,
        fillcolor=palette[s],
        opacity=0.6,
        line_color=palette[s],
        points=False
    ))

fig_box.update_layout(
    title=dict(
        text="<b>Log2-Expression Distribution across Samples</b><br><sup>Distribution of log2(FPKM + 1) for expressed genes</sup>",
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

fig_box.write_html("new plots/html/plotly_expression_boxplot.html")
fig_box.write_image("visualizations/12_Expression_Distribution_Boxplot.png", scale=3)
fig_box.write_image("visualizations/12_Expression_Distribution_Boxplot.pdf")
fig_box.write_image("docs/assets/12_Expression_Distribution_Boxplot.png", scale=3)
fig_box.write_image("docs/assets/12_Expression_Distribution_Boxplot.pdf")
print("Saved Expression Distribution Boxplot (HTML, PNG, PDF)")

# -------------------------------------------------------------------------
# Plot 2: Overlaid Density Plot (Kernel Density Estimation - KDE)
# -------------------------------------------------------------------------
fig_density = go.Figure()

x_eval = np.linspace(0, 14, 400)
for s in samples:
    vals = melted_expressed[melted_expressed['Sample'] == s]['Log2_FPKM'].values
    kde = stats.gaussian_kde(vals)
    density = kde(x_eval)
    
    fig_density.add_trace(go.Scatter(
        x=x_eval,
        y=density,
        mode='lines',
        name=s,
        line=dict(color=palette[s], width=3),
        fill='tozeroy',
        fillcolor=f"rgba{tuple(list(int(palette[s].lstrip('#')[i:i+2], 16) for i in (0, 2, 4)) + [0.15])}"
    ))

fig_density.update_layout(
    title=dict(
        text="<b>Expression Density Distributions (KDE)</b><br><sup>Density profile of log2(FPKM + 1) showing library comparability</sup>",
        font=dict(size=20, color="#0f172a", family="Plus Jakarta Sans, sans-serif"),
        x=0.5
    ),
    xaxis=dict(title="<b>log2(FPKM + 1)</b>", showline=True, linecolor="#cbd5e1", linewidth=1.5, showgrid=True, gridcolor="#f1f5f9"),
    yaxis=dict(title="<b>Density</b>", showline=True, linecolor="#cbd5e1", linewidth=1.5, showgrid=True, gridcolor="#f1f5f9"),
    paper_bgcolor="white",
    plot_bgcolor="white",
    width=900,
    height=600,
    legend=dict(title="<b>Sample</b>", bgcolor="rgba(255,255,255,0.8)", bordercolor="#cbd5e1", borderwidth=1),
    margin=dict(l=70, r=40, t=90, b=60)
)

fig_density.write_html("new plots/html/plotly_expression_density.html")
fig_density.write_image("visualizations/13_Expression_Density_Plot.png", scale=3)
fig_density.write_image("visualizations/13_Expression_Density_Plot.pdf")
fig_density.write_image("docs/assets/13_Expression_Density_Plot.png", scale=3)
fig_density.write_image("docs/assets/13_Expression_Density_Plot.pdf")
print("Saved Expression Density Plot (HTML, PNG, PDF)")

# -------------------------------------------------------------------------
# Plot 3: Counts-Per-Million (CPM) & Cumulative Distribution (ECDF)
# -------------------------------------------------------------------------
counts_file = "counts/gene_counts.txt"
counts_df = pd.read_csv(counts_file, sep='\t', comment='#')
sample_cols = [c for c in counts_df.columns if 'NPFM' in c or 'sorted.bam' in c]

# Compute CPM
cpm_dict = {}
for i, s in enumerate(samples):
    col = sample_cols[i]
    total_counts = counts_df[col].sum()
    cpm_dict[s] = (counts_df[col] / total_counts) * 1e6

cpm_df = pd.DataFrame(cpm_dict)
cpm_df['Gene ID'] = counts_df['Geneid']

# Save CPM matrix
cpm_df[['Gene ID'] + samples].to_csv('downstream_results/CPM_normalized_counts_individual.csv', index=False)
cpm_df[['Gene ID'] + samples].to_csv('docs/assets/CPM_normalized_counts_individual.csv', index=False)
print("Saved CPM_normalized_counts_individual.csv")

fig_cpm = go.Figure()
for s in samples:
    # Sorted CPM in descending order to plot cumulative percentage
    sorted_cpm = np.sort(cpm_df[s].values)[::-1]
    cum_pct = (np.cumsum(sorted_cpm) / sorted_cpm.sum()) * 100
    x_rank = np.arange(1, len(sorted_cpm) + 1)
    
    # Sample every 50 points to keep plot lightweight
    step = 50
    fig_cpm.add_trace(go.Scatter(
        x=x_rank[::step],
        y=cum_pct[::step],
        mode='lines',
        name=s,
        line=dict(color=palette[s], width=2.5)
    ))

fig_cpm.update_layout(
    title=dict(
        text="<b>Cumulative CPM Read Fraction (Library Diversity)</b><br><sup>Percentage of total reads consumed as genes are cumulatively ranked</sup>",
        font=dict(size=20, color="#0f172a", family="Plus Jakarta Sans, sans-serif"),
        x=0.5
    ),
    xaxis=dict(title="<b>Cumulative Gene Rank (Ranked by Expression)</b>", type="log", showline=True, linecolor="#cbd5e1", linewidth=1.5, showgrid=True, gridcolor="#f1f5f9"),
    yaxis=dict(title="<b>Cumulative Reads (%)</b>", range=[0, 105], showline=True, linecolor="#cbd5e1", linewidth=1.5, showgrid=True, gridcolor="#f1f5f9"),
    paper_bgcolor="white",
    plot_bgcolor="white",
    width=900,
    height=600,
    legend=dict(title="<b>Sample</b>", bgcolor="rgba(255,255,255,0.8)", bordercolor="#cbd5e1", borderwidth=1),
    margin=dict(l=70, r=40, t=90, b=60)
)

fig_cpm.write_html("new plots/html/plotly_cpm_cumulative.html")
fig_cpm.write_image("visualizations/14_Cumulative_CPM_Distribution.png", scale=3)
fig_cpm.write_image("visualizations/14_Cumulative_CPM_Distribution.pdf")
fig_cpm.write_image("docs/assets/14_Cumulative_CPM_Distribution.png", scale=3)
fig_cpm.write_image("docs/assets/14_Cumulative_CPM_Distribution.pdf")
print("Saved Cumulative CPM Distribution Plot (HTML, PNG, PDF)")

# -------------------------------------------------------------------------
# Summary Table: Gene Detection Breakdown across Expression Thresholds
# -------------------------------------------------------------------------
cutoffs = [0.1, 1.0, 5.0, 10.0, 50.0]
stats_rows = []
total_genes = len(fpkm_df)

for cut in cutoffs:
    row = {'FPKM Cutoff': f'> {cut} FPKM'}
    for s in samples:
        cnt = (fpkm_df[s] > cut).sum()
        pct = (cnt / total_genes) * 100
        row[s] = f"{cnt:,} ({pct:.1f}%)"
    stats_rows.append(row)

dist_stats_df = pd.DataFrame(stats_rows)
dist_stats_df.to_csv('downstream_results/Expression_Distribution_Threshold_Stats.csv', index=False)
dist_stats_df.to_csv('docs/assets/Expression_Distribution_Threshold_Stats.csv', index=False)
print("Saved Expression_Distribution_Threshold_Stats.csv:\n", dist_stats_df)

print("All distribution analyses and plots generated successfully!")
