import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, dcc, html, State, ctx
from scripts.plots_generation import create_topology_plot, plot_expression_data
import pandas as pd
import os
from functools import lru_cache
import shelve

# Initial values for the dropdowns
isoforms_inital_value = []
cancer_types_inital_value = []

app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server
base_dir = os.path.dirname(os.path.abspath(__file__))

# Get the protein options from the shelve file
db = shelve.open(f"{base_dir}/files_for_plots/membrane_topology_objects")
protein_options = list(db.keys())
db.close()

header = html.Div(
    [
        html.H3("Transcript Explorer", className="me-4"),
        html.Div(
            [
                dbc.Input(id="protein-input", 
                          list="protein-options", 
                          placeholder="Type a protein name...",
                          autoComplete="off"),
                html.Datalist(id="protein-options")
            ],
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
                dbc.NavLink("Description", href="/", active="exact"),
                dbc.NavLink("Topology", href="/Topology", active="exact"),
                dbc.NavLink("Expression", href="/Expression", active="exact"),
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
        dcc.Location(id="url"),
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
    Input("protein-input", "value")
)
def update_protein_options(value):
    if value is None:
        return []
    options = [html.Option(value=protein) for protein in protein_options if protein.upper().startswith(value.upper())]
    if len(options) > 10:
        options = options[:10] + [html.Option(value="...", disabled=True)]
    return options

@lru_cache(maxsize=10)
def get_expression_data(protein):
    db = shelve.open(f"{base_dir}/files_for_plots/TCGA_GTEx_plotting_data")
    df = db[protein]
    db.close()
    return df

@lru_cache(maxsize=10)
def get_topology_data(protein):
    db = shelve.open(f"{base_dir}/files_for_plots/transcripts_to_isoforms_mapping")
    mapping = db[protein]
    db.close()

    db = shelve.open(f"{base_dir}/files_for_plots/membrane_topology_objects")
    sequences_data = db[protein]
    db.close()

    return mapping, sequences_data

@app.callback(
    Output("page-content", "children"),
    Input("url", "pathname"),
    Input("protein-input", "n_submit"),
    Input("protein-submit", "n_clicks"),
    State("protein-input", "value")
)
def render_page_content(pathname, n_submit, n_clicks, protein):
    if protein is None:
        return html.P("Please select a protein")

    protein = protein.upper()

    if protein not in protein_options:
        return html.P("This protein is not in the database")

    expression_df = get_expression_data(protein)

    if pathname == "/":
        return html.P("This is the content of the home page of " + protein + "!")

    elif pathname == "/Topology":
        return html.Div(
            [
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
        )

    elif pathname == "/Expression":
        return html.Div(id="expression-container")

@app.callback(
    Output("expression-container", "children"),
    Input("expression-container", "id"),
    State("protein-input", "value"),
    optional=True,
)
def manage_expression_page(expression_container_id, protein):
    protein = protein.upper()
    expression_df = get_expression_data(protein)

    return html.Div([
        html.Div([
                html.P("Select the isoforms to plot:"),
                dcc.Dropdown(
                    options=expression_df.loc[:,"protein"].unique(),
                    value=isoforms_inital_value,
                    multi=True,
                    id="expression-isoforms-dropdown",
                    placeholder="Select or leave empty to plot all isoforms",
                ),
                html.Br(),
                html.P("Select the cancer types to plot:"),
                dcc.Dropdown(
                    options=expression_df.loc[:,"cancer_type"].unique(),
                    value=cancer_types_inital_value,
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
            dbc.Button("Reset parameters", id="expression-reset-button", className="me-1"),
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
    Output("expression-isoforms-dropdown", "value"),
    Output("expression-cancer-type-dropdown", "value"),
    Input("expression-load-button", "n_clicks"),
    Input("expression-reset-button", "n_clicks"),
    State("expression-isoforms-dropdown", "value"),
    State("expression-cancer-type-dropdown", "value"),
    prevent_initial_call=True,
    optional=True,
)
def expression_container_style(_1, _2, isoforms, cancer_types):
    ctx = dash.callback_context
    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
    if trigger_id == "expression-load-button":
        return {'display': 'block'}, {'display': 'none'}, [isoforms, cancer_types], isoforms, cancer_types
    elif trigger_id == "expression-reset-button":
        return {'display': 'none'}, {'display': 'block'}, dash.no_update, isoforms_inital_value, cancer_types_inital_value

@app.callback(
    Output("topology-plot", "figure"),
    Input("protein-input", "value"),
    optional=True,
)
def topology_plot(protein):
    protein = protein.upper()
    title = "Membrane topology"
    x_label = "Amino acid position in MSA"
    mapping, sequences_data = get_topology_data(protein)
    fig = create_topology_plot(mapping, sequences_data, title, x_label)
    return fig

@app.callback(
    Output('expression-plot', 'figure'),
    Input("expression-parameters", "data"),
    State("protein-input", "value"),
    prevent_initial_call=True,
    optional=True,
)
def update_expression_plot(parameters, protein):
    protein = protein.upper()
    expression_df = get_expression_data(protein)
    isoforms, cancer_types = parameters

    ## Filter logic
    idx = pd.IndexSlice

    # Filter isoforms
    if len(isoforms) > 0:
        expression_df = expression_df.loc[expression_df.loc[:,"protein"].isin(isoforms),:]

    # Filter cancer type
    if len(cancer_types) > 0:
        expression_df = expression_df.loc[expression_df.loc[:,"cancer_type"].isin(cancer_types),:]

    # Generate the plot
    fig = plot_expression_data(expression_df)

    # Reset the buttons and close the popover
    return fig

if __name__ == "__main__":
    app.run(debug=True)