import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import os

pio.templates.default = "plotly_dark"
pio.templates["plotly_dark"].layout.paper_bgcolor = "white"
pio.templates["plotly_dark"].layout.plot_bgcolor = "white"
pio.templates["plotly_dark"].layout.font.color = "black"

go_csv = 'processed_data/Figure2_GO_Classification_Data.csv'
df_go = pd.read_csv(go_csv)

total_all = df_go['All Unigene'].sum()
total_deg = df_go['DEG Unigene'].sum()
df_go['All Unigene %'] = (df_go['All Unigene'] / total_all) * 100
df_go['DEG Unigene %'] = (df_go['DEG Unigene'] / total_deg) * 100

df_go_melt = df_go.melt(id_vars=['Ontology', 'Description'], value_vars=['All Unigene %', 'DEG Unigene %'], 
                        var_name='Type', value_name='Percentage')

df_go_counts = df_go.melt(id_vars=['Ontology', 'Description'], value_vars=['All Unigene', 'DEG Unigene'],
                          var_name='Type_Count', value_name='Count')
df_go_melt['Count'] = df_go_counts['Count']

ontology_order = {'Biological Process': 1, 'Cellular Component': 2, 'Molecular Function': 3}
df_go_melt['Ontology_Order'] = df_go_melt['Ontology'].map(ontology_order)
df_go_melt = df_go_melt.sort_values(by=['Ontology_Order', 'Percentage'], ascending=[True, False])

df_go_melt['Category'] = df_go_melt['Ontology'] + ' - ' + df_go_melt['Type'].str.replace(' %', '')

color_map = {
    'Biological Process - All Unigene': '#FCA5A5',
    'Biological Process - DEG Unigene': '#EF4444',
    'Cellular Component - All Unigene': '#86EFAC',
    'Cellular Component - DEG Unigene': '#22C55E',
    'Molecular Function - All Unigene': '#93C5FD',
    'Molecular Function - DEG Unigene': '#3B82F6'
}

df_go_melt['Percentage_Plot'] = df_go_melt['Percentage'].replace(0, 0.05)

fig_go = px.bar(df_go_melt, x='Description', y='Percentage_Plot', color='Category', barmode='group',
                color_discrete_map=color_map,
                log_y=True,
                hover_data=['Percentage', 'Count'],
                height=900, width=2600)

fig_go.update_traces(marker_line_color='black', marker_line_width=0.8)

fig_go.update_layout(
    title='',
    xaxis_tickangle=-60,
    xaxis_title="",
    yaxis_title="<b>Percentage of genes</b>",
    plot_bgcolor='white',
    paper_bgcolor='white',
    font_color='black',
    showlegend=False,
    bargap=0.0,
    bargroupgap=0.0,
    margin=dict(b=200, r=150, t=50, l=100)
)

fig_go.update_xaxes(showline=True, linewidth=1, linecolor='black', gridcolor='white', tickfont=dict(size=12, color='black'))
fig_go.update_yaxes(showline=True, linewidth=1, linecolor='black', gridcolor='lightgray',
                    tickvals=[0.1, 1, 10, 100], ticktext=["0.1", "1", "10", "100"], range=[-1.5, 2.05],
                    title_font=dict(size=16, color='black'), tickfont=dict(size=14, color='black'))

# Add right y-axis
fig_go.add_trace(go.Scatter(x=[df_go_melt['Description'].iloc[0]], y=[100], yaxis='y2', mode='markers', marker=dict(color='white', size=1), showlegend=False))

fig_go.update_layout(
    yaxis2=dict(
        title="<b>Number of genes</b>",
        title_font=dict(size=16, color='black'),
        overlaying='y',
        side='right',
        type='log',
        range=[-1.5, 2.05],
        tickvals=[0.1, 1, 10, 100],
        ticktext=[
            f"{int(total_deg * 0.001)}<br>{int(total_all * 0.001)}",
            f"{int(total_deg * 0.01)}<br>{int(total_all * 0.01)}",
            f"{int(total_deg * 0.1)}<br>{int(total_all * 0.1)}",
            f"{int(total_deg)}<br>{int(total_all)}"
        ],
        showgrid=False,
        showline=True, linewidth=1, linecolor='black',
        tickfont=dict(size=14, color='black')
    )
)

# Custom legend (using annotations and shapes)
legend_x = 0.85
legend_y_deg = 0.95
legend_y_all = 0.91

# DEG
fig_go.add_shape(type="rect", xref="paper", yref="paper", x0=legend_x, y0=legend_y_deg-0.015, x1=legend_x+0.01, y1=legend_y_deg+0.015, fillcolor="#EF4444", line_width=0)
fig_go.add_shape(type="rect", xref="paper", yref="paper", x0=legend_x+0.012, y0=legend_y_deg-0.015, x1=legend_x+0.022, y1=legend_y_deg+0.015, fillcolor="#22C55E", line_width=0)
fig_go.add_shape(type="rect", xref="paper", yref="paper", x0=legend_x+0.024, y0=legend_y_deg-0.015, x1=legend_x+0.034, y1=legend_y_deg+0.015, fillcolor="#3B82F6", line_width=0)
fig_go.add_annotation(x=legend_x+0.04, y=legend_y_deg, xref="paper", yref="paper", text="<b>DEG Unigene</b>", showarrow=False, font=dict(size=14, color="black"), xanchor="left")

# All
fig_go.add_shape(type="rect", xref="paper", yref="paper", x0=legend_x, y0=legend_y_all-0.015, x1=legend_x+0.01, y1=legend_y_all+0.015, fillcolor="#FCA5A5", line_width=0)
fig_go.add_shape(type="rect", xref="paper", yref="paper", x0=legend_x+0.012, y0=legend_y_all-0.015, x1=legend_x+0.022, y1=legend_y_all+0.015, fillcolor="#86EFAC", line_width=0)
fig_go.add_shape(type="rect", xref="paper", yref="paper", x0=legend_x+0.024, y0=legend_y_all-0.015, x1=legend_x+0.034, y1=legend_y_all+0.015, fillcolor="#93C5FD", line_width=0)
fig_go.add_annotation(x=legend_x+0.04, y=legend_y_all, xref="paper", yref="paper", text="<b>All Unigene</b>", showarrow=False, font=dict(size=14, color="black"), xanchor="left")

# X-axis brackets
def add_bracket(fig, df, ontology, text, color, y_offset=-1.0, y_depth=0.03):
    cats = df[df['Ontology'] == ontology]['Description'].unique()
    if len(cats) == 0: return
    x0 = cats[0]
    x1 = cats[-1]
    
    # Horizontal line
    fig.add_shape(type="line", xref="x", yref="paper", x0=x0, y0=y_offset, x1=x1, y1=y_offset, line=dict(color="black", width=1.5))
    # Vertical ticks at ends
    fig.add_shape(type="line", xref="x", yref="paper", x0=x0, y0=y_offset, x1=x0, y1=y_offset+y_depth, line=dict(color="black", width=1.5))
    fig.add_shape(type="line", xref="x", yref="paper", x0=x1, y0=y_offset, x1=x1, y1=y_offset+y_depth, line=dict(color="black", width=1.5))
    
    # Label
    # To place it in the center, we need the paper coordinate of the center, but xref="x" doesn't easily let us find the midpoint of categorical axis for annotations unless we pass the exact categorical value. We can just pick the middle element.
    mid_idx = len(cats) // 2
    x_mid = cats[mid_idx]
    fig.add_annotation(x=x_mid, y=y_offset - 0.05, xref="x", yref="paper", text=f"<b>{text}</b>", showarrow=False, font=dict(size=16, color="black"), xanchor="center")

fig_go.write_image("GO_Classification_without_ontology_labels.png", scale=2)

# Now add extra margin for the brackets
fig_go.update_layout(margin=dict(b=500, r=150, t=50, l=100))

add_bracket(fig_go, df_go, 'Biological Process', 'biological process', 'black')
add_bracket(fig_go, df_go, 'Cellular Component', 'cellular component', 'black')
add_bracket(fig_go, df_go, 'Molecular Function', 'molecular function', 'black')

fig_go.write_image("GO_Classification_with_ontology_labels.png", scale=2)
fig_go.write_image("test_go.png", scale=2)
print("Done")
