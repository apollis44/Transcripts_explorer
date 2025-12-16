import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import pandas as pd
import math

def create_topology_plot(mapping, sequences_data, title, x_label):
    print("Creating topology plot...")
    
    # Calculate height based on number of isoforms
    num_isoforms = len(sequences_data)
    calculated_height = num_isoforms * 0.5 + 2 

    # Define colors
    color_map = {
        'S': '#FF0000', 'O': '#E69F00', 'I': '#95B8C8', '-': "#808080",
        'M': '#006D6F', 'E': "#970005", 'N': '#95B8C8', 'X': '#E69F00',
    }

    # Create the plot with calculated height
    fig, ax = plt.subplots(figsize=(10, calculated_height))

    y_labels = mapping.values()
    unique_letters = set()

    # Loop through sequences
    for i, seq_data in enumerate(sequences_data):
        for feature, ranges in seq_data.items():
            unique_letters.add(feature)
            ax.broken_barh(
                xranges=ranges,
                yrange=(i + 0.1, 0.8), # (y_start, height) -> bar covers i to i+1
                facecolors=color_map[feature],
            )

    # Center ticks in the middle of each bar (i + 0.5)
    ax.set_yticks([i + 0.5 for i in range(len(y_labels))])
    ax.set_yticklabels(y_labels, fontsize=10)
    
    # Remove the tick marks themselves (the little lines), but keep the labels
    ax.tick_params(axis='y', length=0) 

    # Style the plot
    ax.set_ylim(0, len(sequences_data))
    
    # Ensure x-limit covers the full protein length
    # (Assuming membrane_topology has the sequence length info correct)
    max_length = max([seq_data[-1][0] + seq_data[-1][1] for seq_data in sequences_data[0].values()])
    ax.set_xlim(0, max_length)
    
    ax.set_xlabel(x_label, fontsize=14)
    ax.spines[['top', 'right', 'left']].set_visible(False)

    # Invert y-axis so the first isoform is at the top
    ax.invert_yaxis()

    # Legend
    letter_to_label = {
        'S': 'Signal peptide', 'O': 'Extracellular', 'I': 'Intracellular',
        '-': 'Alignment gap', 'M': 'Transmembrane', 'E': 'Epitope',
        'N': 'Intron', 'X': 'Exon',
    }
    legend_patches = [mpatches.Patch(color=color, label=letter_to_label[label])
                      for label, color in color_map.items() if label in unique_letters]
    
    # Move legend outside
    ax.legend(handles=legend_patches, bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0., fontsize=7)

    plt.title(title, fontsize=16)
    return fig

def plot_expression_data(expression_df):
    print("Creating expression plot...")

    cancer_types = expression_df.index.get_level_values('cancer_type').unique()
    nb_plots = len(cancer_types)

    fig, axes = plt.subplots(math.ceil(nb_plots/2), 2 if nb_plots > 1 else 1, figsize=(10, nb_plots * 5))

    for i, cancer_type in enumerate(cancer_types):
        current_ax = None
        if nb_plots > 2:
            current_ax = axes[i//2, i%2]
        elif nb_plots == 2:
            current_ax = axes[i%2]
        else:
            current_ax = axes
        stats = expression_df.xs(cancer_type, level="cancer_type")

        boxes = current_ax.bxp(stats.loc[:,"expression"], patch_artist=True)
        
        for box, color in zip(boxes['boxes'], [stat['color'] for stat in stats.loc[:,"expression"]]):
            box.set_facecolor(color)

        # Add title
        current_ax.set_title(cancer_type)

        # Rotate x-axis labels
        current_ax.tick_params(axis='x', rotation=90)

    # Add legend
    handles = []
    labels = []
    for unique_isoform in stats.index.get_level_values('protein').unique():
        color = stats.xs(unique_isoform, level="protein")["expression"].iloc[0]['color']
        handles.append(mpatches.Patch(facecolor=color))
        labels.append(unique_isoform)

    fig.legend(handles, labels, loc='upper right')

    if nb_plots%2 != 0 and nb_plots > 1:
        axes[-1, -1].axis('off')  # Hide the last subplot if odd number of cancer types

    return fig