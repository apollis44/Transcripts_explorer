import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, dcc, html, State, ctx, ALL, no_update
from scripts.plots_generation import *
import pandas as pd
import os
from functools import lru_cache
import shelve
import numpy as np
import urllib.parse

# Initial values for the dropdowns
tissue_types_inital_value = []

app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server
base_dir = os.path.dirname(os.path.abspath(__file__))

# Get the protein options from the shelve file
db_names = shelve.open(f"{base_dir}/files_for_plots/genes")
protein_options = list(db_names.keys())

def getting_gene_names(gene_name):
    if gene_name not in db_names:
        return None
    gene_names = db_names[gene_name]
    return gene_names

header = html.Div(
    [
        html.H3("Transcript Explorer", className="me-4"),
        html.Div(
            [
                dbc.Input(id="protein-input", 
                          placeholder="Type a protein name...",
                          autoComplete="off"),
                dbc.ListGroup(
                    id="protein-options",
                    style={
                        "position": "absolute", 
                        "top": "100%", 
                        "left": 0, 
                        "width": "100%", 
                        "zIndex": 1000,
                        "maxHeight": "200px",
                        "overflowY": "auto",
                        "display": "none" # Hidden by default
                    }
                ),
            ], 
            style={"position": "relative"},
        ),
        dbc.Button("Submit", id="protein-submit", n_clicks=0),
    ],
    className="d-flex align-items-center border-bottom mb-4 pt-2 ps-3",
    style={"height": "60px"}
)

sidebar = html.Div(
    [
        dbc.Nav(
            [
                dbc.NavLink("Description", id="link-description", href="/", active="exact"),
                dbc.NavLink("Localization", id="link-localization", href="/Localization", active="exact"),
                dbc.NavLink("Topology", id="link-topology", href="/Topology", active="exact"),
                dbc.NavLink("Expression", id="link-expression", href="/Expression", active="exact"),
            ],
            vertical=True,
            pills=True,
        ),
    ],
    style={
            "backgroundColor": "#f8f9fa",
            "padding": "10px",
            "border": "1px solid #dee2e6",
            "borderRadius": "4px",
            "height": "100%"
        },
)

content = html.Div(id="page-content")

app.layout = dbc.Container(
    [   
        dcc.Location(id="url", refresh=False),
        header,
        dbc.Row(
            [
                # Left Column (Sidebar) - Width 3/12
                dbc.Col(
                    sidebar,
                    width=3
                ),
                
                # Right Column (Content) - Width 9/12
                dbc.Col(
                    [
                        content,
                    ],
                    width=9
                ),
            ]
        ),
    ],
    fluid=True, # Uses full width of the screen
    className="p-0"
)

@app.callback(
    Output("protein-options", "children"),
    Output("protein-options", "style"),
    Input("protein-input", "value"),
    State({"type": "result-item", "index": ALL}, "n_clicks"),
)
def update_protein_options(value, n_clicks):
    if value is None:
        return [], {"display": "none"}
    
    if n_clicks != [] and (np.array(n_clicks) != None).any():
        return no_update, no_update
    
    options = []

    valid_protein_options = [protein for protein in sorted(protein_options) if protein.upper().startswith(value.upper())]
    for protein in valid_protein_options[:10]:
            options.append(
                dbc.ListGroupItem(
                    protein,
                    id={"type": "result-item", "index": protein},
                    action=True
                )
            )

    if len(valid_protein_options) > 10:
        options.append(
            dbc.ListGroupItem(
                f"{len(valid_protein_options) - 10} more options",
                disabled=True,
                color="light"
            )
        )

    style = {
        "position": "absolute", "top": "100%", "left": 0, 
        "width": "100%", "zIndex": 1000, "display": "block"
    }
    return options, style

@app.callback(
    Output("protein-input", "value"),
    Output("protein-options", "style", allow_duplicate=True),
    Output("protein-submit", "n_clicks", allow_duplicate=True),
    Input({"type": "result-item", "index": ALL}, "n_clicks"),
    State("protein-submit", "n_clicks"),
    prevent_initial_call=True
)
def select_item(n_clicks, protein_submit_n_clicks):
    if not ctx.triggered:
        return no_update, no_update, no_update

    if (np.array(n_clicks) == None).all():
        return no_update, no_update, no_update

    triggered_id = ctx.triggered_id
    selected_value = triggered_id['index']

    return selected_value, {"display": "none"}, protein_submit_n_clicks + 1

@lru_cache(maxsize=10)
def get_localization_data(protein):
    db = shelve.open(f"{base_dir}/files_for_plots/deeploc2_output")
    df = db[protein]
    db.close()
    df = df.iloc[:,:-4] # We exclude the last 4 columns that contains data we don't use
    return df

@lru_cache(maxsize=10)
def get_expression_data(protein):
    db = shelve.open(f"{base_dir}/files_for_plots/TCGA_GTEx_plotting_data")
    df = db[protein]
    db.close()
    if df is None:
        return None
    df = df.groupby(['study', 'tissue_type'], sort=False).agg(list).reset_index()
    return df

@lru_cache(maxsize=10)
def get_topology_data(protein):
    db = shelve.open(f"{base_dir}/files_for_plots/transcripts_to_isoforms_mapping")
    mapping = db[protein]
    db.close()

    db = shelve.open(f"{base_dir}/files_for_plots/membrane_topology_objects")
    sequences_data = db[protein]
    db.close()

    db = shelve.open(f"{base_dir}/files_for_plots/TCGA_GTEx_plotting_data")
    df = db[protein]
    available_transcripts = []
    for available_transcript in df.loc[:,"protein"].unique():
        available_transcripts.extend(available_transcript.split("<br>"))

    return mapping, sequences_data, available_transcripts

@lru_cache(maxsize=10)
def get_query_data(search):
    search_str = search.lstrip("?") if "?" in search else ""
    parsed_search = urllib.parse.parse_qs(search_str)
    if "protein" not in parsed_search:
        return None
    return parsed_search["protein"][0]

@app.callback(
    Output("url", "search"),
    Input("protein-input", "n_submit"),
    Input("protein-submit", "n_clicks"),
    State("url", "search"),
    State("protein-input", "value"),
    prevent_initial_call=True,
)
def update_query(n_submit, n_clicks, search, protein):
    if not ctx.triggered:
        return no_update
    
    search_str = search.lstrip("?") if "?" in search else ""
    parsed_search = urllib.parse.parse_qs(search_str)
    parsed_search["protein"] = [protein]
    search_str = urllib.parse.urlencode(parsed_search, doseq=True)
    return "?" + search_str

@app.callback(
    [
        Output("link-description", "href"),
        Output("link-localization", "href"),
        Output("link-topology", "href"),
        Output("link-expression", "href"),
    ],
    Input("url", "search")
)
def update_nav_links(search):
    # If search is None or empty, just return the base paths
    query = search if search else ""
    
    # We return the base path + the current query string for each link
    return f"/{query}", f"/Localization{query}",f"/Topology{query}", f"/Expression{query}"

@app.callback(
    Output("page-content", "children"),
    Output("protein-options", "style", allow_duplicate=True),
    Input("url", "search"),
    Input("url", "pathname"),
    prevent_initial_call=True,
)
def render_page_content(query, pathname):

    protein = get_query_data(query)

    if protein is None:
        return html.P("Please select a protein"), {"display": "none"}

    protein = getting_gene_names(protein.upper())

    if protein is None:
        return html.P("This protein is not in the database"), {"display": "none"}

    if pathname == "/":
        return html.P("This is the content of the home page of " + protein + "!"), {"display": "none"}

    elif pathname == "/Localization":
        return html.Div(
            [
                dbc.Spinner(
                    children=dcc.Graph(id="localization-plot"),
                    size="lg",
                    color="primary",
                    type="border",
                    fullscreen=False,
                    id="localization-spinner",
                    spinner_style={
                        "position": "absolute", 
                        "top": "20px", 
                        "left": "50%", 
                    },
                ),
            ]
        ), {"display": "none"}

    elif pathname == "/Topology":
        return html.Div(
            [
                dbc.Checklist(
                    options=[
                        {"label": "All transcripts", "value": "all"},
                    ],
                    value=["all"],
                    switch=True,
                    id="topology-checklist",
                ),
                dbc.Spinner(
                    children=dcc.Graph(id="topology-plot"),
                    size="lg",
                    color="primary",
                    type="border",
                    fullscreen=False,
                    id="topology-spinner",
                    spinner_style={
                        "position": "absolute", 
                        "top": "20px", 
                        "left": "50%", 
                    },
                ),
            ]
        ), {"display": "none"}

    elif pathname == "/Expression":
        return html.Div(id="expression-container"), {"display": "none"}

@app.callback(
    Output("expression-container", "children"),
    Input("expression-container", "id"),
    State("url", "search"),
    optional=True,
)
def manage_expression_page(expression_container_id, search):
    protein = get_query_data(search)
    protein = getting_gene_names(protein.upper())
    expression_df = get_expression_data(protein)

    if expression_df is None:
        return html.P("Expression data is not available for this protein's transcripts")

    return html.Div([
        html.Div([
                html.P("Select the cancer types to plot:"),
                dcc.Dropdown(
                    options=expression_df.loc[:,"tissue_type"].unique(),
                    value=tissue_types_inital_value,
                    multi=True,
                    id="expression-cancer-type-dropdown",
                    placeholder="Select or leave empty to plot all cancer types",
                ),
                html.Br(),
                # We want the buttons next to each other
                dbc.Button("Load plot", id="expression-load-button", className="me-1")
            ],
            id="expression-parameters-container"
        ),
        html.Div([
            dbc.Button("Reset parameters", id="expression-reset-button", className="me-1"),
            dbc.Spinner(
                children=dcc.Graph(id="expression-plot"),
                size="lg",
                color="primary",
                type="border",
                fullscreen=False,
                id="expression-spinner",
                spinner_style={
                    "position": "absolute", 
                    "top": "20px", 
                    "left": "50%", 
                },
            ),
            ],
            id="expression-container",
            style={'display': 'none'}
        ),
        dcc.Store(id="expression-parameters")
    ])

@app.callback(
    Output("expression-container", "style"),
    Output("expression-parameters-container", "style"),
    Output("expression-parameters", "data"),
    Output("expression-cancer-type-dropdown", "value"),
    Input("expression-load-button", "n_clicks"),
    Input("expression-reset-button", "n_clicks"),
    State("expression-cancer-type-dropdown", "value"),
    prevent_initial_call=True,
    optional=True,
)
def expression_container_style(_1, _2, tissue_types):
    ctx = dash.callback_context
    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
    if trigger_id == "expression-load-button":
        return {'display': 'block'}, {'display': 'none'}, tissue_types, tissue_types
    elif trigger_id == "expression-reset-button":
        return {'display': 'none'}, {'display': 'block'}, tissue_types_inital_value, tissue_types_inital_value

@app.callback(
    Output("localization-plot", "figure"),
    Input("url", "search"),
    optional=True,
)
def localization_plot(search):
    protein = get_query_data(search)
    protein = getting_gene_names(protein.upper())
    localization_data = get_localization_data(protein)
    fig = create_localization_plot(localization_data)
    return fig

@app.callback(
    Output("topology-plot", "figure"),
    Input("url", "search"),
    Input("topology-checklist", "value"),
    optional=True,
)
def topology_plot(search, all_transcripts):
    protein = get_query_data(search)
    protein = getting_gene_names(protein.upper())
    title = "Membrane topology"
    x_label = "Amino acid position in MSA"
    mapping, sequences_data, available_transcripts = get_topology_data(protein)
    fig = create_topology_plot(mapping, sequences_data, available_transcripts, title, x_label, all_transcripts)
    return fig

@app.callback(
    Output('expression-plot', 'figure'),
    Input("expression-parameters", "data"),
    State("url", "search"),
    prevent_initial_call=True,
    optional=True,
)
def update_expression_plot(parameters, search):
    protein = get_query_data(search)
    protein = getting_gene_names(protein.upper())
    expression_df = get_expression_data(protein)
    tissue_types = parameters

    # Filter cancer type
    if len(tissue_types) > 0:
        expression_df = expression_df.loc[expression_df.loc[:,"tissue_type"].isin(tissue_types),:]
    
    # Generate the plot
    fig = plot_expression_data(expression_df)

    # Reset the buttons and close the popover
    return fig

if __name__ == "__main__":
    app.run(debug=True)