import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, dcc, html, State, ctx
from scripts.plots_generation import create_topology_plot, plot_expression_data
import pandas as pd
import os
import json
from functools import lru_cache
import shelve

# Initial values for the dropdowns
isoforms_inital_value = []
studies_inital_value = ["TCGA"]
cancer_types_inital_value = []

app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])
base_dir = os.path.dirname(os.path.abspath(__file__))

header = html.Div(
    [
        html.H3("Transcript Explorer", className="me-4"),
        dbc.Tabs(
            [
                dbc.Tab(label="HER2", tab_id="HER2", label_style={"cursor": "pointer"}),
                dbc.Tab(label="CD20", tab_id="CD20", label_style={"cursor": "pointer"}),
            ],
            id="top-tabs",
            active_tab="HER2",
            className="border-0",
        ),
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
            "backgroundColor": "#f8f9fa", # Light gray background
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

@lru_cache(maxsize=10)
def get_expression_data(active_tab):
    db = shelve.open(f"{base_dir}/files_for_plots/TCGA_GTEx_plotting_data")
    df = db[active_tab]
    db.close()
    return df

@lru_cache(maxsize=10)
def get_topology_data(active_tab):
    db = shelve.open(f"{base_dir}/files_for_plots/transcripts_to_isoforms_mapping")
    mapping = db[active_tab]
    db.close()

    db = shelve.open(f"{base_dir}/files_for_plots/membrane_topology_objects")
    sequences_data = db[active_tab]
    db.close()
    # mapping = json.load(open(base_dir + "/files_for_plots/" + active_tab + "/transcript_to_isoform_mapping.json"))
    # sequences_data = json.load(open(base_dir + "/files_for_plots/" + active_tab + "/sequences_data.json"))
    return mapping, sequences_data

@app.callback(
    Output("page-content", "children"),
    Input("url", "pathname"),
    Input("top-tabs", "active_tab")
)
def render_page_content(pathname, active_tab):
    expression_df = get_expression_data(active_tab)

    if pathname == "/":
        return html.P("This is the content of the home page of " + active_tab + "!")
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
        return html.Div(
            [
                dbc.Button(
                    "Choose parameters",
                    id="expression-popover-button",
                    className="me-1",
                ),
                dbc.Popover(
                    dbc.PopoverBody(
                        [
                            html.P("Select the isoforms to plot:"),
                            dcc.Dropdown(
                                options=expression_df.index.get_level_values("protein").unique(),
                                value=isoforms_inital_value,
                                multi=True,
                                id="expression-isoforms-dropdown",
                            ),
                            html.Br(),
                            html.P("Select the studies to plot:"),
                            dcc.Checklist(
                                options=expression_df.index.get_level_values("study").unique(),
                                value=studies_inital_value,
                                id="expression-study-checklist",
                            ),
                            html.Br(),
                            html.P("Select the cancer types to plot:"),
                            dcc.Dropdown(
                                options=expression_df.index.get_level_values("cancer_type").unique(),
                                value=cancer_types_inital_value,
                                multi=True,
                                id="expression-cancer-type-dropdown",
                            ),
                            html.Br(),
                            # We want the buttons next to each other
                            dbc.Button("Update Plot", id="expression-update-button", className="me-1"),
                            dbc.Button("Cancel", id="expression-cancel-button", className="me-1"),
                        ]
                    ),
                    target="expression-popover-button",
                    trigger="click",
                    id="expression-popover",
                    style={"width": "500px"}
                ),
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
        )
         
    # If the user tries to reach a different page, return a 404 message
    return html.Div(
        [
            html.H1("404: Not found", className="text-danger"),
            html.Hr(),
            html.P(f"The pathname {pathname} was not recognised..."),
        ],
        className="p-3 bg-light rounded-3",
    )

@app.callback(
    Output("topology-plot", "figure"),
    Input("top-tabs", "active_tab"),
    optional=True,
)
def topology_plot(active_tab):
    title = "Membrane topology"
    x_label = "Amino acid position in MSA"
    mapping, sequences_data = get_topology_data(active_tab)

    fig = create_topology_plot(mapping, sequences_data, title, x_label)
    return fig

@app.callback(
    Output('expression-plot', 'figure'),
    Output("expression-isoforms-dropdown", "value"),
    Output("expression-study-checklist", "value"),
    Output("expression-cancer-type-dropdown", "value"),
    Input('expression-update-button', 'n_clicks'),
    Input('expression-cancel-button', 'n_clicks'),
    State("expression-isoforms-dropdown", "value"),
    State("expression-study-checklist", "value"),
    State("expression-cancer-type-dropdown", "value"),
    State("top-tabs", "active_tab"),
    optional=True,
    running=[
        (Output("expression-popover", "is_open"), False, False),
    ]
    )
def update_expression_plot(_1, _2, isoforms, studies, cancer_types, active_tab):
    button_id = ctx.triggered_id

    expression_df = get_expression_data(active_tab)

    # If cancel button was the one pushed we reset the values
    if button_id == "expression-cancel-button":
        isoforms = isoforms_inital_value
        studies = studies_inital_value
        cancer_types = cancer_types_inital_value

    # For the first call, the values are not initialized yet
    if isoforms is None or studies is None or cancer_types is None:
        isoforms = isoforms_inital_value
        studies = studies_inital_value
        cancer_types = cancer_types_inital_value
    
    ## Filter logic
    idx = pd.IndexSlice

    # Filter isoforms
    if len(isoforms) > 0:
        expression_df = expression_df.loc[idx[:,:,list(isoforms),:]]

    # Filter study
    # If no study is selected, return an empty figure
    if len(studies) == 0:
        empty_fig = go.Figure()
        empty_fig.update_layout(
            xaxis =  { "visible": False },
            yaxis = { "visible": False },
            annotations = [
                {   
                    "text": "No study selected",
                    "xref": "paper",
                    "yref": "paper",
                    "showarrow": False,
                    "font": {
                        "size": 28
                    }
                }
            ]
        )
        return empty_fig

    expression_df = expression_df.loc[idx[list(studies),:,:,:]]

    # Filter cancer type
    if len(cancer_types) > 0:
        expression_df = expression_df.loc[idx[:,list(cancer_types),:,:]]

    # Generate the plot
    fig = plot_expression_data(expression_df)

    # Reset the buttons and close the popover
    return fig, isoforms, studies, cancer_types

@app.callback(
    Output("expression-cancer-type-dropdown", "options"),
    Input("expression-study-checklist", "value"),
    State("top-tabs", "active_tab"),
    optional=True,
)
def update_cancer_type_dropdown(studies, active_tab):
    expression_df = get_expression_data(active_tab)
    idx = pd.IndexSlice
    if studies:
        return expression_df.loc[idx[list(studies),:,:,:]].index.get_level_values("cancer_type").unique().tolist()
    return dash.no_update

if __name__ == "__main__":
    app.run(debug=True)