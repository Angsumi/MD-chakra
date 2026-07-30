import pandas as pd
import plotly.express as px
import plotly.io as pio

pio.templates.default = "plotly_dark"
table1_csv = 'processed_data/Table 1.csv'
df_table1 = pd.read_csv(table1_csv)

if 'Reference Genome (All Genes)' in df_table1.columns:
    df_table1 = df_table1.drop(columns=['Reference Genome (All Genes)'])

for col in ['NPFM1', 'NPFM2', 'NPFM3', 'NPFM4']:
    if df_table1[col].dtype == object:
        df_table1[col] = df_table1[col].astype(str).str.replace('%', '').astype(float)

df_melt = df_table1.melt(id_vars=['Metric'], value_vars=['NPFM1', 'NPFM2', 'NPFM3', 'NPFM4'], 
                         var_name='Sample', value_name='Value')

# Vibrant but soothing color palette
soothing_vibrant_colors = ['#FF8FA3', '#FFD166', '#06D6A0', '#118AB2']

fig_table1 = px.bar(df_melt, x='Metric', y='Value', color='Sample', barmode='stack',
                    title='RNA-seq Data and Expressed Genes Statistics (Stacked)',
                    color_discrete_sequence=soothing_vibrant_colors,
                    log_y=True,
                    height=700)

fig_table1.update_layout(title_font_size=24, title_x=0.5, xaxis_tickangle=-45,
                         xaxis_title="Metric", yaxis_title="Value (Log Scale)")

fig_table1.write_html("new plots/html/plotly_table1_stats.html")
fig_table1.write_image("new plots/png/plotly_table1_stats.png", scale=3)
fig_table1.write_image("new plots/pdf/plotly_table1_stats.pdf")
print("Regenerated Table 1 stats as a single unified stacked bar plot.")
