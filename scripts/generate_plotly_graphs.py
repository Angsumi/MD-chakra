import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from sklearn.decomposition import PCA
from scipy.spatial.distance import pdist, squareform
import os

# Set a modern, sleek default template
pio.templates.default = "plotly_dark"
pio.templates["plotly_dark"].layout.paper_bgcolor = "white"
pio.templates["plotly_dark"].layout.plot_bgcolor = "white"
pio.templates["plotly_dark"].layout.font.color = "black"

print("Generating Plotly graphs...")

# ---------------------------------------------------------
# 1. Table 1 Statistics Bar Charts
# ---------------------------------------------------------
table1_csv = 'Table 1.csv'
if os.path.exists(table1_csv):
    df_table1 = pd.read_csv(table1_csv)
    # Drop reference genome column
    if 'Reference Genome (All Genes)' in df_table1.columns:
        df_table1 = df_table1.drop(columns=['Reference Genome (All Genes)'])
    
    for col in ['NPFM1', 'NPFM2', 'NPFM3', 'NPFM4']:
        if df_table1[col].dtype == object:
            df_table1[col] = df_table1[col].astype(str).str.replace('%', '').astype(float)
    
    # Save processed CSV for this specific graph
    plotly_table1_csv = 'plotly_table1_data.csv'
    df_table1.to_csv(plotly_table1_csv, index=False)
    
    # Melt for plotly
    df_melt = df_table1.melt(id_vars=['Metric'], value_vars=['NPFM1', 'NPFM2', 'NPFM3', 'NPFM4'], 
                             var_name='Sample', value_name='Value')
    
    # --- Original Faceted Table 1 Plot ---
    fig_table1 = px.bar(df_melt, x='Sample', y='Value', color='Sample', facet_col='Metric', facet_col_wrap=3,
                        title='RNA-seq Data and Expressed Genes Statistics',
                        color_discrete_sequence=['#00e1d9', '#6366f1', '#ec4899', '#f59e0b'],
                        height=1000, width=1600,
                        facet_col_spacing=0.08, facet_row_spacing=0.12)
    
    fig_table1.update_xaxes(showgrid=False, zeroline=False)
    fig_table1.update_yaxes(matches=None, showticklabels=True, showgrid=False, zeroline=False)
    fig_table1.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    fig_table1.update_layout(title_font_size=28, title_x=0.5, margin=dict(t=120))
    fig_table1.write_html("plotly_table1_stats.html")
    fig_table1.write_image("plotly_table1_stats.png", scale=3)
    fig_table1.write_image("plotly_table1_stats.pdf")
    
    # --- New Table 1 Plot (GO Classification Style) ---
    fig_table1_style2 = px.bar(df_melt, x='Metric', y='Value', color='Sample', barmode='group',
                               title='RNA-seq Data and Expressed Genes Statistics (Grouped Style)',
                               color_discrete_sequence=['#E63946', '#F4A261', '#2A9D8F', '#264653'],
                               log_y=True,
                               height=800, width=1600)
                               
    fig_table1_style2.update_layout(
        title_font_size=28, title_x=0.5,
        plot_bgcolor='white', paper_bgcolor='white', font_color='black',
        xaxis_title="", yaxis_title="<b>Value (log scale)</b>",
        xaxis_tickangle=-45, bargap=0.1, bargroupgap=0.0,
        legend=dict(title="<b>Sample</b>", x=1.01, y=1, bgcolor='white', bordercolor='black', borderwidth=1),
        margin=dict(b=200, r=150, t=100, l=100)
    )
    
    fig_table1_style2.update_xaxes(showline=True, linewidth=1, linecolor='black', showgrid=False, zeroline=False, tickfont=dict(size=16, color='black'))
    fig_table1_style2.update_yaxes(showline=True, linewidth=1, linecolor='black', showgrid=False, zeroline=False, tickfont=dict(size=14, color='black'))
    fig_table1_style2.update_traces(marker_line_color='black', marker_line_width=1)
    
    fig_table1_style2.write_html("plotly_table1_stats_style2.html")
    fig_table1_style2.write_image("plotly_table1_stats_style2.png", scale=3)
    fig_table1_style2.write_image("plotly_table1_stats_style2.pdf")
    print("Saved plotly_table1_stats_style2.html, .png, .pdf")
    print("Saved plotly_table1_stats.html, .png, .pdf")


# ---------------------------------------------------------
# 2. GO Classification Data
# ---------------------------------------------------------
go_csv = 'Figure2_GO_Classification_Data.csv'
if os.path.exists(go_csv):
    df_go = pd.read_csv(go_csv)
    
    # Calculate percentages
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
        margin=dict(b=600, r=150, t=50, l=100)
    )
    
    fig_go.update_xaxes(showline=True, linewidth=1, linecolor='black', showgrid=False, zeroline=False, tickfont=dict(size=12, color='black'))
    fig_go.update_yaxes(showline=True, linewidth=1, linecolor='black', showgrid=False, zeroline=False,
                        tickvals=[0.1, 1, 10, 100], ticktext=["0.1", "1", "10", "100"], range=[-1, 2.05],
                        title_font=dict(size=16, color='black'), tickfont=dict(size=14, color='black'))

    # Add dummy trace for right y-axis
    fig_go.add_trace(go.Scatter(x=[df_go_melt['Description'].iloc[0]], y=[100], yaxis='y2', mode='markers', marker=dict(color='white', size=1), showlegend=False))

    fig_go.update_layout(
        yaxis2=dict(
            title="<b>Number of genes</b>",
            title_font=dict(size=16, color='black'),
            overlaying='y',
            side='right',
            type='log',
            range=[-1, 2.05],
            tickvals=[0.1, 1, 10, 100],
            ticktext=[
                f"{int(total_deg * 0.001)}<br>{int(total_all * 0.001)}",
                f"{int(total_deg * 0.01)}<br>{int(total_all * 0.01)}",
                f"{int(total_deg * 0.1)}<br>{int(total_all * 0.1)}",
                f"{int(total_deg)}<br>{int(total_all)}"
            ],
            showgrid=False, zeroline=False,
            showline=True, linewidth=1, linecolor='black',
            tickfont=dict(size=14, color='black')
        )
    )

    # Custom legend (using annotations and shapes)
    legend_x = 0.85
    legend_y_deg = 0.98
    legend_y_all = 0.93
    w = 0.016
    h = 0.025

    # DEG Unigene
    fig_go.add_shape(type="rect", xref="paper", yref="paper", x0=legend_x, y0=legend_y_deg-h, x1=legend_x+w, y1=legend_y_deg+h, fillcolor="#EF4444", line_width=0)
    fig_go.add_shape(type="rect", xref="paper", yref="paper", x0=legend_x+w, y0=legend_y_deg-h, x1=legend_x+2*w, y1=legend_y_deg+h, fillcolor="#22C55E", line_width=0)
    fig_go.add_shape(type="rect", xref="paper", yref="paper", x0=legend_x+2*w, y0=legend_y_deg-h, x1=legend_x+3*w, y1=legend_y_deg+h, fillcolor="#3B82F6", line_width=0)
    fig_go.add_annotation(x=legend_x+3.5*w, y=legend_y_deg, xref="paper", yref="paper", text="<b>DEG Unigene</b>", showarrow=False, font=dict(size=18, color="black"), xanchor="left")

    # All Unigene
    fig_go.add_shape(type="rect", xref="paper", yref="paper", x0=legend_x, y0=legend_y_all-h, x1=legend_x+w, y1=legend_y_all+h, fillcolor="#FCA5A5", line_width=0)
    fig_go.add_shape(type="rect", xref="paper", yref="paper", x0=legend_x+w, y0=legend_y_all-h, x1=legend_x+2*w, y1=legend_y_all+h, fillcolor="#86EFAC", line_width=0)
    fig_go.add_shape(type="rect", xref="paper", yref="paper", x0=legend_x+2*w, y0=legend_y_all-h, x1=legend_x+3*w, y1=legend_y_all+h, fillcolor="#93C5FD", line_width=0)
    fig_go.add_annotation(x=legend_x+3.5*w, y=legend_y_all, xref="paper", yref="paper", text="<b>All Unigene</b>", showarrow=False, font=dict(size=18, color="black"), xanchor="left")

    # X-axis brackets
    def add_bracket(fig, df, ontology, text, color, y_offset=-0.65, y_depth=0.02):
        cats = df[df['Ontology'] == ontology]['Description'].unique()
        if len(cats) == 0: return
        x0 = cats[0]
        x1 = cats[-1]
        
        # Horizontal line
        fig.add_shape(type="line", xref="x", yref="paper", x0=x0, y0=y_offset, x1=x1, y1=y_offset, line=dict(color="black", width=2.0))
        # Vertical ticks at ends
        fig.add_shape(type="line", xref="x", yref="paper", x0=x0, y0=y_offset, x1=x0, y1=y_offset+y_depth, line=dict(color="black", width=2.0))
        fig.add_shape(type="line", xref="x", yref="paper", x0=x1, y0=y_offset, x1=x1, y1=y_offset+y_depth, line=dict(color="black", width=2.0))
        
        mid_idx = len(cats) // 2
        x_mid = cats[mid_idx]
        fig.add_annotation(x=x_mid, y=y_offset - 0.01, xref="x", yref="paper", text=f"<b>{text}</b>", showarrow=False, font=dict(size=24, color="black"), xanchor="center", yanchor="top")

    add_bracket(fig_go, df_go, 'Biological Process', 'biological process', 'black')
    add_bracket(fig_go, df_go, 'Cellular Component', 'cellular component', 'black')
    add_bracket(fig_go, df_go, 'Molecular Function', 'molecular function', 'black')

    fig_go.write_html("plotly_go_classification.html")
    fig_go.write_image("plotly_go_classification.png", scale=3)
    fig_go.write_image("plotly_go_classification.pdf")
    print("Saved plotly_go_classification.html, .png, .pdf")

# ---------------------------------------------------------
# 3. PCA, Heatmap, and Distance Matrix from VST Counts
# ---------------------------------------------------------
vst_csv = "vst_normalized_counts.csv"
if os.path.exists(vst_csv):
    df_vst = pd.read_csv(vst_csv, index_col=0)
    
    # --- PCA ---
    pca_csv = "plotly_pca_data.csv"
    if os.path.exists(pca_csv):
        print(f"Using existing {pca_csv} for PCA Plot")
        pca_df = pd.read_csv(pca_csv)
    else:
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(df_vst.T)
        pca_df = pd.DataFrame({
            'PC1': pca_result[:, 0],
            'PC2': pca_result[:, 1],
            'Sample': df_vst.columns
        })
        pca_df.to_csv(pca_csv, index=False)
        print(f"Generated and saved {pca_csv}")
    
    fig_pca = px.scatter(pca_df, x='PC1', y='PC2', color='Sample', 
                         title='Principal Component Analysis (PCA)', 
                         size_max=25, size=[15]*len(pca_df),
                         color_discrete_sequence=['#E63946', '#F4A261', '#2A9D8F', '#264653'],
                         height=600)
    fig_pca.update_traces(marker=dict(size=20, line=dict(width=2, color='white')))
    fig_pca.update_layout(title_font_size=24, title_x=0.5)
    fig_pca.update_xaxes(showgrid=False, zeroline=False, showline=True, linewidth=1, linecolor='black', mirror=True)
    fig_pca.update_yaxes(showgrid=False, zeroline=False, showline=True, linewidth=1, linecolor='black', mirror=True)
    fig_pca.write_html("plotly_pca_plot.html")
    fig_pca.write_image("plotly_pca_plot.png", scale=3)
    fig_pca.write_image("plotly_pca_plot.pdf")
    print("Saved plotly_pca_plot.html, .png, .pdf")

    # --- Heatmap Top 50 Variable Genes ---
    heatmap_csv = "plotly_heatmap_top50_data.csv"
    if os.path.exists(heatmap_csv):
        print(f"Using existing {heatmap_csv} for Heatmap")
        df_top50_z = pd.read_csv(heatmap_csv, index_col=0)
    else:
        gene_vars = df_vst.var(axis=1)
        top50_genes = gene_vars.nlargest(50).index
        df_top50 = df_vst.loc[top50_genes]
        # Z-score normalization
        df_top50_z = df_top50.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
        df_top50_z.to_csv(heatmap_csv)
        print(f"Generated and saved {heatmap_csv}")

    fig_heat = px.imshow(df_top50_z, text_auto=False, aspect="auto", 
                         title="Expression Heatmap (Top 50 Most Variable Genes)",
                         color_continuous_scale="viridis", height=800)
    fig_heat.update_layout(title_font_size=24, title_x=0.5)
    fig_heat.update_xaxes(showgrid=False, zeroline=False)
    fig_heat.update_yaxes(showgrid=False, zeroline=False)
    fig_heat.write_html("plotly_heatmap.html")
    fig_heat.write_image("plotly_heatmap.png", scale=3)
    fig_heat.write_image("plotly_heatmap.pdf")
    print("Saved plotly_heatmap.html, .png, .pdf")

    # --- Sample Distance Matrix ---
    dist_csv = "plotly_distance_matrix_data.csv"
    if os.path.exists(dist_csv):
        print(f"Using existing {dist_csv} for Distance Matrix")
        dist_df = pd.read_csv(dist_csv, index_col=0)
    else:
        dist_matrix = pdist(df_vst.T, metric='euclidean')
        sq_dist = squareform(dist_matrix)
        dist_df = pd.DataFrame(sq_dist, index=df_vst.columns, columns=df_vst.columns)
        dist_df.to_csv(dist_csv)
        print(f"Generated and saved {dist_csv}")
    
    fig_dist = px.imshow(dist_df, text_auto=".1f", aspect="auto",
                         title="Sample Distance Matrix (Euclidean)",
                         color_continuous_scale="plasma", height=600)
    fig_dist.update_layout(title_font_size=24, title_x=0.5)
    fig_dist.update_xaxes(showgrid=False, zeroline=False)
    fig_dist.update_yaxes(showgrid=False, zeroline=False)
    fig_dist.write_html("plotly_distance_matrix.html")
    fig_dist.write_image("plotly_distance_matrix.png", scale=3)
    fig_dist.write_image("plotly_distance_matrix.pdf")
    print("Saved plotly_distance_matrix.html, .png, .pdf")

print("All tasks completed.")
