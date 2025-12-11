from shiny import App, ui, render, module
import os
import matplotlib.pyplot as plt
import pandas as pd
import json
from plots_generation import create_topology_plot

base_dir = os.path.dirname(os.path.abspath(__file__))

# ==========================================
# 1. Define the Module
# ==========================================

# The Module UI: Defines what ONE protein tab looks like
@module.ui
def protein_tab_ui(protein_name):
    return ui.nav_panel(
        protein_name, 
        ui.navset_pill_list(
            ui.nav_panel("Description", 
                ui.markdown(f"## Description of {protein_name}"),
            ),
            ui.nav_panel("Topology", 
                ui.card(
                    ui.card_header("Isoform Topology"),
                    # Note: We just use simple ID "topology" here. 
                    # The module handles the unique ID prefixes automatically.
                    ui.output_plot("topology_plot", fill=True),
                    full_screen=True,
                )
            ),
            ui.nav_panel("Expression",
                ui.markdown(f"## Expression of {protein_name}"),
                ui.card(
                    ui.output_plot("expression_plot", fill=True),
                    full_screen=True,
                )
            ),
            widths=(3, 9)
        ),
        value=protein_name
    )

# The Module Server: Defines the logic for ONE protein
@module.server
def protein_tab_server(input, output, session, protein_name):
    
    @render.plot
    def topology_plot():
        title = "Membrane topology"
        x_label = "Amino acid position in MSA"
        membrane_topology = pd.read_csv(base_dir + "/files_for_plots/" + protein_name + "/membrane_topology.csv", index_col=0)
        sequences_data = json.load(open(base_dir + "/files_for_plots/" + protein_name + "/sequences_data.json"))
        fig = create_topology_plot(membrane_topology, sequences_data, title, x_label)
        return fig

    @render.plot
    def expression_plot():        
        protein = input.main_navbar()
        mapping = json.load(open(base_dir + "/files_for_plots/" + protein + "/transcript_to_isoform_mapping.json"))
        expression_df = pd.read_csv(base_dir + "/files_for_plots/" + protein + "/TCGA_GTEx_expression_data.csv")
        fig = plot_expression(expression_df, mapping, type_data="GTEX") # Using GTEX as default or example
        return fig

# ==========================================
# 2. Main Application
# ==========================================

# Get protein list (mocked for this example)
proteins = os.listdir(os.path.join(base_dir, "files_for_plots"))

# Generate the UI List
# We assign a unique 'id' to each module call (e.g., "tab_Protein_A")
my_tabs = [
    protein_tab_ui(id=f"tab_{p}", protein_name=p) for p in proteins
]

app_ui = ui.page_navbar(
    *my_tabs,  # Unpack the list of tabs
    id="main_navbar",
    title="Page title",
    fillable=True
)

def server(input, output, session):
    # Call the server logic for each protein
    for p in proteins:
        protein_tab_server(id=f"tab_{p}", protein_name=p)

app = App(app_ui, server)
