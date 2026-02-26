import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import math

def create_topology_plot(mapping, sequences_data, available_transcripts, title, x_label):
    print("Creating topology plot...")
    
    # Define colors and labels
    color_map = {
        'S': '#FF0000', 
        'O': '#E69F00', 
        'I': '#95B8C8', 
        '-': "#808080",
        'M': '#006D6F', 
    }

    letter_to_label = {
        'S': 'Signal peptide', 
        'O': 'Extracellular', 
        'I': 'Intracellular',
        '-': 'Alignment gap', 
        'M': 'Transmembrane', 
    }

    fig = go.Figure()
    isoforms = list(dict.fromkeys(mapping.values()))
    y_labels = [transcript_id for isoform in isoforms for transcript_id in isoform.split("<br>")]

    y_label_available = [y_label in available_transcripts for y_label in y_labels]

    added_to_legend = set()

    # Loop through sequences (isoforms)
    for i, seq_data in enumerate(sequences_data):
        if not y_label_available[i]:
            continue
        y_val = y_labels[i]
        
        for feature, ranges in seq_data.items():
            for start, width in ranges:
                # Add a bar for each segment
                fig.add_trace(go.Bar(
                    name=letter_to_label.get(feature, feature),
                    x=[width],
                    base=[start],
                    y=[y_val],
                    orientation='h',
                    marker_color=color_map.get(feature, '#000000'),
                    showlegend=True if feature not in added_to_legend else False,
                    legendgroup=feature, # Groups identical features in the legend
                    hovertemplate=f"<b>{letter_to_label.get(feature)}</b><br>" +
                                  f"Range: {start} - {start + width}<br>" +
                                  f"Length: {width}<extra></extra>"
                ))
                added_to_legend.add(feature)

    # Calculate height to match your original logic
    num_isoforms = len(y_label_available)
    calculated_height = num_isoforms * 50 + 200 

    fig.update_layout(
        title={
            'text': title,
            'x': 0.5,
            'y': 1,
            'xanchor': 'center',
            'yanchor': 'top',
            'font': {'size': 20}
        },
        xaxis_title=x_label,
        barmode='stack', # This ensures bars with the same 'base' don't overlap vertically
        bargap=0.02,
        height=calculated_height,
        plot_bgcolor='white',
        legend=dict(
            title="Features",
            orientation="v",
            yanchor="top",
            y=1,
            xanchor="left",
            x=1.02
        ),
        yaxis=dict(
            autorange="reversed", # Invert y-axis like your Matplotlib code
            showgrid=False,
            zeroline=False
        ),
        margin=dict(l=10, r=10, t=50, b=50)
    )

    return fig

def plot_expression_data(expression_df):
    print("Creating expression plot...")

    cancer_types = expression_df.loc[:,"cancer_type"].unique().tolist()

    nb_cancer_types = len(cancer_types)
    rows_count = math.ceil(nb_cancer_types / 2)
    pixels_per_row = 600

    fig = make_subplots(rows=math.ceil(nb_cancer_types/2), 
                        cols=2 if nb_cancer_types > 1 else 1,
                        shared_yaxes="all",
                        subplot_titles=cancer_types,
                        vertical_spacing=150 / (rows_count * pixels_per_row), # 150 pixels vertical spacing between plots
                        horizontal_spacing= 0.03, # 3% horizontal spacing between plots
                        )

    colors = fig.layout.template.layout.colorway

    for i, cancer_type in enumerate(cancer_types):
        data_for_each_cancer_type = expression_df.loc[(expression_df.loc[:,"cancer_type"] == cancer_type), :]
        current_color = colors[i % len(colors)]

        fig.add_trace(
            go.Box(
                x=data_for_each_cancer_type["protein"].iloc[0],
                q1=data_for_each_cancer_type["q1"].iloc[0],
                q3=data_for_each_cancer_type["q3"].iloc[0],
                median=data_for_each_cancer_type["median"].iloc[0],
                lowerfence=data_for_each_cancer_type["lowerfence"].iloc[0],
                upperfence=data_for_each_cancer_type["upperfence"].iloc[0],
                showlegend=False,
                marker_color=current_color,
            ), 
            row=i//2+1, 
            col=i%2+1,
        )

        proteins = data_for_each_cancer_type["protein"].iloc[0]
        outliers_lists = data_for_each_cancer_type["y"].iloc[0]

        # We combine all X and Y coordinates into single flat lists for this trace
        all_x = []
        all_y = []

        for protein, values in zip(proteins, outliers_lists):
            all_x.extend([protein] * len(values))
            all_y.extend(values)

        # Same color as the box plot
        fig.add_trace(
            go.Scatter(
                x=all_x, 
                y=all_y,
                mode='markers',
                marker=dict(
                    size=5,
                    symbol='circle-open',
                    color=current_color,
                ),
                showlegend=False,
            ),
            row=i//2+1,
            col=i%2+1
        )

    fig.update_layout(
        height=rows_count * pixels_per_row, 
        yaxis=dict(
            title=dict(
                text="log2(TPM+1)",
            )
        ),
    )

    return fig