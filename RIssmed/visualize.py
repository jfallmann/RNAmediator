import gzip
from typing import Dict, List
import csv
import sqlite3

import dash  # (version 1.12.0) pip install dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.express as px  # (version 4.7.0)
import plotly.graph_objects as go
import plotly.io as pio
from dash.dependencies import Input, Output

PLOTLY_COLORS = px.colors.qualitative.Light24
pio.templates["plotly_white"].update({"layout": {
    # e.g. you want to change the background to transparent
    'paper_bgcolor': 'rgba(0,0,0,0)',
    'plot_bgcolor': ' rgba(0,0,0,0)',
    "font": dict(color="white")
}})


def read_data(file):
    return pd.read_csv(file, delimiter='\t',
                       names=['Chr', 'Start', 'End', 'Constraint', 'Accessibility_difference', 'Strand',
                              'Distance_to_constraint', 'Accessibility_no_constraint', 'Accessibility_constraint',
                              'Energy_Difference', 'Kd_change',
                              'Zscore'])  # ,'ChrBS','StartBS','EndBS','NameBS','ScoreBS','StrandBS'])


def main():
    data = []
    with gzip.open("/home/rabsch/Documents/egg_mountpoint/Collection_unpaired.bed.gz", "rb") as handle:
        for x, line in enumerate(handle):
            if x > 20000:
                line = line.decode()
                data.append(line)
            if x == 300000:
                break
    with open("testfile.bed", "w") as handle:
        handle.write("".join(data))


def dropdown_menues(gene_id_options: List[Dict], constraint_options: List[Dict]):
    dropdown_style = {'width': "100%", "margin": "auto", }
    dropdowns = [
        dcc.Dropdown(id="slct_chrom",
                     options=gene_id_options,
                     multi=False,
                     value=gene_id_options[0]["value"],
                     style=dropdown_style),
        dcc.Dropdown(id="slct_constraint",
                     options=constraint_options,
                     multi=False,
                     value=constraint_options[0]["value"],
                     style=dropdown_style
                     ),
    ]
    menu = html.Table(
        html.Tr(
            [html.Td(html.Div(col, style={"width": "100%"}), className="selection-table_interactive") for col in
             dropdowns]
        ),
        id="selection-table", style={"width": "80%", "margin": "auto"}

    )

    return menu


def get_app_layout(app: dash.Dash, df: pd.DataFrame):
    sel = selectors(df)
    app.layout = html.Div([

        html.Div(html.H1("RIssmed Dasboard", style={'text-align': 'center'}), className="page-header"),
        html.Div(dropdown_menues(*sel), className="databox"),

        html.Div([html.H3(id='header', children=[], style={"text-align": "center"}),
                  dcc.Graph(id='plotly_graph', figure={"layout": {"height": 700}})], className="plotly-graph"),

    ], id="wrapper"
    )


def selectors(df: pd.DataFrame):
    gene_id_opts = [{"label": entry, "value": entry} for entry in sorted(df["Chr"].unique())]
    constraint_opts = [{"label": entry, "value": entry} for entry in sorted(df["Constraint"].unique())]
    return gene_id_opts, constraint_opts


app = dash.Dash(__name__, external_stylesheets=[dbc.themes.DARKLY])
df = read_data("testfile.bed")
get_app_layout(app, df)


@app.callback(
    [Output(component_id='header', component_property='children'),
     Output(component_id='plotly_graph', component_property='figure')],
    [Input(component_id='slct_chrom', component_property='value'),
     Input(component_id='slct_constraint', component_property='value')]

)
def update_graph(slct_chrom, slct_constraint):
    header = f"Changes in probabilities of being unpaired (Accessibility) for Gene {slct_chrom} at constraint {slct_constraint}"
    dff = df.copy()
    dff = dff[dff["Chr"] == slct_chrom]
    fig = go.Figure()
    if type(slct_constraint) == str:
        slct_constraint = [slct_constraint]
    for col_idx, element in enumerate(slct_constraint):
        pdf = dff[dff["Constraint"] == element]
        x = pdf["Distance_to_constraint"]
        y = pdf["Accessibility_difference"]
        fig.add_trace(go.Scatter(x=x, y=y, line={"width": 4, "color": PLOTLY_COLORS[col_idx]}))
    fig.layout.template = "plotly_white"
    fig.update_yaxes(range=[-1, 1])
    return [header, fig]


if __name__ == '__main__':
    app.run_server(debug=True, port=8080, host="0.0.0.0")
