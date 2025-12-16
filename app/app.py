from shiny import App, ui, render, module, reactive
import os
import matplotlib.pyplot as plt
import pandas as pd
import json
from plots_generation import create_topology_plot, plot_expression_data


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
                    ui.output_plot("topology_plot", fill=True),
                    full_screen=True,
                )
            ),
            ui.nav_panel("Expression",
                ui.popover(
                    ui.input_action_button("btn", "A button", class_="mt-3"),
                    ui.input_selectize(
                        "isoforms", 
                        "Isoforms:", 
                        choices=[], 
                        selected=[],
                        multiple=True,
                        width="100%"
                    ),
                    ui.input_checkbox_group(
                        "study",
                        "Study:",
                        choices=[],
                        selected=[],
                        width="100%"
                    ),
                    ui.input_selectize(
                        "cancer_type", 
                        "Cancer type:", 
                        choices=[], 
                        selected=[],
                        multiple=True,
                        width="100%"
                    ),
                    ui.input_action_button("btn_update", "Update Plot"),
                    ui.input_action_button("btn_cancel", "Cancel"),
                    title="Parameters for " + protein_name,
                    id="popover_id"
                ),
                ui.card(
                    ui.card_header("Expression of the isoforms"),
                    ui.output_plot("expression_plot", height = "10000px"),
                )
            ),
            widths=(3, 9)
        ),
        value=protein_name
    )

# The Module Server: Defines the logic for ONE protein
@module.server
def protein_tab_server(input, output, session, expression_rvs, protein_name):

    @reactive.calc
    def get_isoform_list():
        expression_df = import_expression_data()
        return expression_df.index.get_level_values('protein').unique().tolist()

    @reactive.calc
    def get_study_list():
        expression_df = import_expression_data()
        return expression_df.index.get_level_values('study').unique().tolist()

    @reactive.calc
    def get_cancer_type_list():
        expression_df = import_expression_data()
        return expression_df.index.get_level_values('cancer_type').unique().tolist()
    
    @reactive.effect
    def initiliase_expression_UIs():
        isoforms = get_isoform_list()
        studies = get_study_list()
        cancer_types = get_cancer_type_list()

        # We want the default values to be all cancers and all isoforms for TCGA
        expression_rvs["isoforms"].set([])
        expression_rvs["study"].set(["TCGA"])
        expression_rvs["cancer_type"].set([])
        
        # Update the UI element with the loaded list
        ui.update_selectize(
            "isoforms", 
            choices=isoforms, 
            selected=[]
        )

        ui.update_checkbox_group(
            "study", 
            choices=studies, 
            selected=["TCGA"]
        )

        ui.update_selectize(
            "cancer_type", 
            choices=cancer_types, 
            selected=[]
        )

    @reactive.effect
    @reactive.event(input.btn_update)
    def update_expression_rvs():
        ui.update_popover("popover_id", show=False)
        expression_rvs["isoforms"].set(input.isoforms())
        expression_rvs["study"].set(input.study())
        expression_rvs["cancer_type"].set(input.cancer_type())
        
    @reactive.effect
    @reactive.event(input.btn_cancel)
    def cancel_expression_rvs():
        ui.update_popover("popover_id", show=False)
        isoforms = get_isoform_list()
        cancer_types = get_cancer_type_list()

        # We reset the reactive variables
        expression_rvs["isoforms"].set(isoforms)
        expression_rvs["study"].set(["TCGA"])
        expression_rvs["cancer_type"].set(cancer_types)

        # We reset the popover UI
        ui.update_selectize(
            "isoforms", 
            choices=isoforms, 
            selected=[]
        )

        ui.update_checkbox_group(
            "study", 
            choices=studies, 
            selected=["TCGA"]
        )

        ui.update_selectize(
            "cancer_type", 
            choices=cancer_types, 
            selected=[]
        )

    @render.plot
    def topology_plot():
        title = "Membrane topology"
        x_label = "Amino acid position in MSA"
        # Ensure paths exist in your actual implementation
        mapping = json.load(open(base_dir + "/files_for_plots/" + protein_name + "/transcript_to_isoform_mapping.json"))
        sequences_data = json.load(open(base_dir + "/files_for_plots/" + protein_name + "/sequences_data.json"))
        fig = create_topology_plot(mapping, sequences_data, title, x_label)
        return fig

    @reactive.calc
    def import_expression_data():
        print("Importing expression data...")
        path = base_dir + "/files_for_plots/" + protein_name + "/TCGA_GTEx_plotting_data.csv"
        expression_df = pd.read_csv(path, index_col=[0,1,2,3])

        expression_df['expression'] = expression_df['expression'].apply(json.loads)

        return expression_df

    @reactive.calc
    def compute_expression_fig():
        """
        Calculates the figure. This only runs when the data changes. It does NOT run on window resize.
        """
        expression_df = import_expression_data()
        
        ## Filter logic

        idx = pd.IndexSlice

        # Filter isoforms
        selected_isoforms = expression_rvs["isoforms"]()

        if len(selected_isoforms) > 0:
            expression_df = expression_df.loc[idx[:,:,list(selected_isoforms),:]]

        # Filter study
        selected_study = expression_rvs["study"]()

        if len(selected_study) == 0:
            empty_fig = plt.figure()
            plt.text(0.5, 1, "No study selected", fontsize=16, ha="center", va="center")
            plt.axis("off")
            return empty_fig

        expression_df = expression_df.loc[idx[list(selected_study),:,:,:]]

        # Filter cancer type
        selected_cancer_type = expression_rvs["cancer_type"]()
        
        if len(selected_cancer_type) > 0:
            expression_df = expression_df.loc[idx[:,list(selected_cancer_type),:,:]]

        # Generate the heavy plot
        fig = plot_expression_data(expression_df)
        return fig

    @render.plot
    def expression_plot():
        """
        This just serves the pre-calculated figure. 
        If the modal opens and resizes the page, this runs, but it instantly 
        retrieves the result from compute_expression_fig() without re-calculating.
        """
        return compute_expression_fig()

# ==========================================
# 2. Main Application
# ==========================================

# Get protein list (mocked for this example)
# Ensure this directory exists or use a dummy list for testing
files_path = os.path.join(base_dir, "files_for_plots")
if os.path.exists(files_path):
    proteins = os.listdir(files_path)
else:
    proteins = ["Protein_A", "Protein_B"] # Fallback for testing

my_tabs = [
    protein_tab_ui(id=f"tab_{p}", protein_name=p) for p in proteins
]

custom_css = ui.tags.style("""
    .popover {
        max-width: 800px !important; /* Adjust this value as needed */
        width: auto;                 /* Allows it to fit content up to max-width */
    }
""")

app_ui = ui.page_navbar(
    *my_tabs,
    id="main_navbar",
    title="Page title",
    fillable=True,
    header=custom_css
)

def server(input, output, session):
    for p in proteins:
        # Create a separate reactive value for each protein loop
        # We wrap it in a function/closure context if needed, but passing it as arg works
        isoforms_rv = reactive.Value([])
        study_rv = reactive.Value(["TCGA"])
        cancer_type_rv = reactive.Value([])
        expression_rvs = {"isoforms": isoforms_rv, 
                          "study": study_rv,
                          "cancer_type": cancer_type_rv}
        protein_tab_server(id=f"tab_{p}", protein_name=p, expression_rvs=expression_rvs)

app = App(app_ui, server)