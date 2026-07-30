import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio

pio.templates.default = "plotly_white" # White background like the screenshot
go_csv = 'processed_data/Figure2_GO_Classification_Data.csv'
df_go = pd.read_csv(go_csv)
df_go_melt = df_go.melt(id_vars=['Ontology', 'Description'], value_vars=['All Unigene', 'DEG Unigene'], var_name='Type', value_name='Count')

# We want them grouped by Ontology on the x-axis, so we sort by Ontology, then by Count descending
# We map Ontology to ensure the order: BP, CC, MF
ontology_order = {'Biological Process': 1, 'Cellular Component': 2, 'Molecular Function': 3}
df_go_melt['Ontology_Order'] = df_go_melt['Ontology'].map(ontology_order)
df_go_melt = df_go_melt.sort_values(by=['Ontology_Order', 'Count'], ascending=[True, False])

# Create a combined category for coloring
df_go_melt['Color_Category'] = df_go_melt['Ontology'] + ' - ' + df_go_melt['Type']

# Vibrant and soothing colors matching the grouped style
# Reds for BP, Greens for CC, Blues for MF
color_map = {
    'Biological Process - All Unigene': '#FCA5A5', # Light Red
    'Biological Process - DEG Unigene': '#EF4444', # Dark Red
    'Cellular Component - All Unigene': '#86EFAC', # Light Green
    'Cellular Component - DEG Unigene': '#22C55E', # Dark Green
    'Molecular Function - All Unigene': '#93C5FD', # Light Blue
    'Molecular Function - DEG Unigene': '#3B82F6'  # Dark Blue
}

# Add a tiny pseudo-count to handle 0 on log scale gracefully
df_go_melt['Count_Plot'] = df_go_melt['Count'].replace(0, 0.1)

fig_go = px.bar(df_go_melt, x='Description', y='Count_Plot', color='Color_Category', barmode='group',
                color_discrete_map=color_map,
                log_y=True,
                height=800, width=1400)

fig_go.update_layout(
    title='GO Classification (Grouped by Ontology)',
    title_font_size=24, 
    title_x=0.5,
    xaxis_tickangle=-60,
    xaxis_title="",
    yaxis_title="Number of Genes (Log Scale)",
    plot_bgcolor='white',
    paper_bgcolor='white',
    legend_title_text='',
    bargap=0.02,       # Reduces gap between different GO categories (makes bars thicker)
    bargroupgap=0.0    # Reduces gap between All Unigene and DEG Unigene bars
)

# Add lines or text to separate the three ontologies on the x-axis
# This gives the "grouped" look at the bottom
fig_go.write_html("new plots/html/plotly_go_classification.html")
fig_go.write_image("new plots/png/plotly_go_classification.png", scale=3)
fig_go.write_image("new plots/pdf/plotly_go_classification.pdf")
print("Regenerated GO classification plot to match the screenshot style.")
