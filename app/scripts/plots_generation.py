import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import math

def create_topology_plot(mapping, sequences_data, title, x_label):
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
    y_labels = list(mapping.values())
    added_to_legend = set()

    # Loop through sequences (isoforms)
    for i, seq_data in enumerate(sequences_data):
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
    num_isoforms = len(sequences_data)
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

    indexes = expression_df.index.droplevel(level="transcript").droplevel(level="protein").unique()
    expression_data_to_plot = []

    for index in indexes:
        expression_data = expression_df.loc[index].values
        expression_data_to_plot.append(expression_data.flatten())

    cancer_types = expression_df.index.get_level_values('cancer_type').unique().tolist()
    nb_cancer_types = len(cancer_types)
    rows_count = math.ceil(nb_cancer_types / 2)
    pixels_per_row = 300

    isoforms_list = expression_df.index.get_level_values('protein').tolist()

    fig = make_subplots(rows=math.ceil(nb_cancer_types/2), 
                        cols=2 if nb_cancer_types > 1 else 1,
                        shared_yaxes="all",
                        subplot_titles=cancer_types,
                        )

    isoforms_seen = []
    for i, data_for_each_cancer_type in enumerate(expression_data_to_plot):
        for j, data_for_each_transcript in enumerate(data_for_each_cancer_type):
            data_items = data_for_each_transcript.copy()

            data_items.pop("y")
            data_items.pop("boxpoints",)

            if isoforms_list[j] not in isoforms_seen:
                isoforms_seen.append(isoforms_list[j])
                fig.add_trace(
                    go.Box(
                        data_items,
                        name=isoforms_list[j],
                        showlegend=True,
                        legendgroup=isoforms_list[j],
                    ), 
                    row=i//2+1, 
                    col=i%2+1,
                )

            else:
                fig.add_trace(
                    go.Box(
                        data_items,
                        name=isoforms_list[j],
                        showlegend=False,
                        legendgroup=isoforms_list[j],
                    ), 
                    row=i//2+1, 
                    col=i%2+1,
                )

            if len(data_for_each_transcript["y"]) > 0:
                    fig.add_trace(
                        go.Scatter(
                            # Repeat the X name for every outlier point
                            x=[data_for_each_transcript["x"][0]] * len(data_for_each_transcript["y"]), 
                            y=data_for_each_transcript["y"],
                            mode='markers',
                            marker=dict(
                                color=data_for_each_transcript["marker_color"],
                                size=5, # Adjust outlier size as needed
                                symbol='circle-open' # Optional: make them hollow circles
                            ),
                            showlegend=False,
                            name=isoforms_list[j],
                            legendgroup=isoforms_list[j],
                        ),
                        row=i//2+1,
                        col=i%2+1
                    )



    fig.update_layout(
        height=rows_count * pixels_per_row, 
        showlegend=True,
        legend=dict(
            orientation="v",
            yanchor="top",
            y=1,
            xanchor="left",
            x=1.02
        ),
    )

    return fig