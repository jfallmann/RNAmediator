import csv
import gzip
import os
import sqlite3
from typing import Dict, List, Union

import dash  # (version 1.12.0) pip install dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import numpy as np
import pandas as pd
import plotly.express as px  # (version 4.7.0)
import plotly.graph_objects as go
import plotly.io as pio
from dash import callback_context
from dash.dependencies import Input, Output, State, ALL

app = dash.Dash("FOO", external_stylesheets=[dbc.themes.DARKLY])

MODE = "db"
NUMBER_OF_INTERESTING = 10
TABLE_COLUMNS = 1
PLOTLY_COLORS = px.colors.qualitative.Light24
COLUMN_NAMES = {'Chr': "VARCHAR(20)",
                'Start': "INT",
                'End': "INT",
                'Fold_Constraint': "VARCHAR(40)",
                'Accessibility_difference': "REAL",
                'Strand': "VARCHAR(2)",
                'Distance_to_constraint': "INT",
                'Accessibility_no_constraint': "REAL",
                'Accessibility_constraint': "REAL",
                'Energy_Difference': "REAL",
                'Kd_change': "REAL",
                'Zscore': "REAL"}
pio.templates["plotly_white"].update({"layout": {
    # e.g. you want to change the background to transparent
    'paper_bgcolor': 'rgba(0,0,0,0)',
    'plot_bgcolor': ' rgba(0,0,0,0)',
    "font": dict(color="white")
}})


def read_data(file: str):
    return pd.read_csv(file, delimiter='\t',
                       names=list(COLUMN_NAMES))  # ,'ChrBS','StartBS','EndBS','NameBS','ScoreBS','StrandBS'])


def csv_to_sqlite(file: str, db_path: str):
    assert not os.path.exists(db_path), "database already exists."
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    names = [f"{key} {COLUMN_NAMES[key]}" for key in COLUMN_NAMES]
    columns = ", ".join(names)
    call = f'CREATE TABLE IF NOT EXISTS test ( {columns} ) '
    cur.execute(call)
    cur.execute("CREATE INDEX constraint_index ON test (Fold_Constraint)")
    if ".gz" in file:
        file_handle = gzip.open(file, "rt")
    else:
        file_handle = open(file)
    csv_reader = csv.reader(file_handle, delimiter="\t", )
    for row in csv_reader:
        cur.execute("INSERT INTO test VALUES (?,?,?,?,?,?,?,?,?,?,?,?);", row)
    file_handle.close()
    con.commit()
    con.close()
    insert_interesting_table(db_path)


def insert_interesting_table(db_path: str):
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    cur.execute(f'CREATE TABLE IF NOT EXISTS importance '
                f'( Chr VARCHAR(5), '
                f'Fold_Constraint VARCHAR(40) PRIMARY KEY, '
                f'Mean_Value REAL, '
                f'Max_Value REAL, '
                f' FOREIGN KEY (Fold_Constraint) REFERENCES test(Fold_Constraint)) ')
    cur.execute("SELECT DISTINCT Fold_Constraint, Chr FROM test")
    constraints = cur.fetchall()
    for entry in constraints:
        cur.execute("SELECT Distance_to_constraint, Accessibility_difference, Chr FROM test "
                    "WHERE Fold_Constraint=? AND Chr=?",
                    entry)
        values = cur.fetchall()
        constraint_max = np.max([abs(diff[1]) for diff in values])
        constraint_mean = np.mean([abs(diff[1]) for diff in values])
        cur.execute("INSERT INTO importance VALUES (?, ?, ?, ?)", [entry[1], entry[0], constraint_mean, constraint_max])
    con.commit()
    con.close()


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


def interesting_table(interesting: List[sqlite3.Row], prev_clicks: int = 0, next_clicks: int = 0):
    page = next_clicks - prev_clicks
    header = ["index"] + interesting[0].keys()
    interesting = [(x + page*NUMBER_OF_INTERESTING, *element) for x, element in enumerate(interesting)]
    menu = html.Div([html.Div(html.Button("", id="prev-button", n_clicks=prev_clicks,),
                              style={"display": "table-cell", "margin": "auto", 'align-items': 'center',
                                     "justify-content": "center", "vertical-align": "middle"}),
            html.Div(html.Table(
                [html.Tr(
                    [html.Td(element,
                             className="interesting-table-header",
                             style={"font-weight": "bold"}) for element in header]
                )] +
                [html.Tr(
                    list(tablerow_generator(entry, y))
                ) for y, entry in enumerate(interesting)],
                className="interesting-table-tablerows", style={"width": "100%", "margin": "auto", 'text-align': 'center'}

            ), style={"display": "table-cell", "width": "80%"}, className="interesting-table"),
            html.Div(html.Button("", id="next-button", n_clicks=next_clicks,
                                 style={"margin": "auto"}),
                     style={"display": "table-cell", "margin": "auto", 'align-items': 'center',
                            "justify-content": "center", "vertical-align": "middle"})
            ], style={"display": "table", "width": "100%", 'text-align': 'center', "vertical-align": "middle"},
                    className="interesting-table-all")
    return menu


def tablerow_generator(row, row_idx: int):
    for x, col in enumerate(row):
        if x == 0:
            style = {"width": "7%"}
            classname = "interesting-table-index-column"
        else:
            style = {}
            classname = "interesting-table-column"
        column = html.Td(html.Button(col, id={"index": f"{x}-{row_idx}", "type": "interesting-table-button",
                                              "chrom": f"{row[1]}", "constraint": row[2]},
                                     n_clicks=0, className=f"{row[2]}",
                                     style={"white-space": "nowrap", "margin": "auto", "width": "100%"}),

                         className=classname, style=style)
        yield column


def get_app_layout(app: dash.Dash, df: Union[pd.DataFrame, str]):
    if MODE == "db":
        sel = selectors_sql(df)
        interesting = get_interesting(df)
    else:
        sel = selectors_pandas(df)
        interesting = []

    app.layout = html.Div([

        html.Div(html.H1("RIssmed Dasboard", style={'text-align': 'center'}), className="page-header"),
        html.Div(dropdown_menues(*sel), className="databox"),

        html.Div([html.H3(id='header', children=[], style={"text-align": "center"}),
                  dcc.Graph(id='plotly_graph', figure={"layout": {"height": 700}})], className="plotly-graph"),
        html.Div([html.H3(id='header-2', children=[], style={"text-align": "center"})]),
        html.Div(interesting_table(interesting), className="databox", id="interesting-table-div")],
        id="wrapper"
    )


def get_interesting(db_path: str, page: int = 0):
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()
    cur.execute(f"SELECT * FROM importance "
                f"ORDER BY Max_Value DESC "
                f"LIMIT {NUMBER_OF_INTERESTING} OFFSET {NUMBER_OF_INTERESTING * page}")
    return_list = cur.fetchall()
    conn.close()
    global TABLE_COLUMNS
    TABLE_COLUMNS = len(return_list[0])
    return return_list


def selectors_pandas(df: pd.DataFrame):
    gene_id_opts = [{"label": entry, "value": entry} for entry in sorted(df["Chr"].unique())]
    constraint_opts = [{"label": entry, "value": entry} for entry in sorted(df["Fold_Constraint"].unique())]
    return gene_id_opts, constraint_opts


def selectors_sql(db_path: str):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute("SELECT DISTINCT Fold_Constraint FROM test")
    constraint_opts = cur.fetchall()
    cur.execute("SELECT DISTINCT Chr FROM test")
    gene_id_opts = cur.fetchall()
    conn.close()
    gene_id_opts = [{"label": entry[0], "value": entry[0]} for entry in gene_id_opts]
    constraint_opts = [{"label": entry[0], "value": entry[0]} for entry in constraint_opts]
    return gene_id_opts, constraint_opts


def update_graph_via_sql(slct_chrom, slct_constraint):
    conn = sqlite3.connect("whole_bed.db")
    cur = conn.cursor()
    cur.execute("SELECT Distance_to_constraint, Accessibility_difference FROM test WHERE Fold_Constraint=? AND Chr=?",
                (slct_constraint, slct_chrom))
    rows = cur.fetchall()
    x, y = zip(*rows)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x, y=y, line={"width": 4, "color": PLOTLY_COLORS[0]}))
    fig.layout.template = "plotly_white"
    fig.update_yaxes(range=[-1, 1])
    conn.close()
    return fig


def update_graph_via_pandas(slct_chrom, slct_constraint):
    dff = df.copy()
    dff = dff[dff["Chr"] == slct_chrom]
    fig = go.Figure()
    if type(slct_constraint) == str:
        slct_constraint = [slct_constraint]
    for col_idx, element in enumerate(slct_constraint):
        pdf = dff[dff["Fold_Constraint"] == element]
        x = pdf["Distance_to_constraint"]
        y = pdf["Accessibility_difference"]
        fig.add_trace(go.Scatter(x=x, y=y, line={"width": 4, "color": PLOTLY_COLORS[col_idx]}))
    fig.layout.template = "plotly_white"
    fig.update_yaxes(range=[-1, 1])
    return fig


@app.callback(
    [Output(component_id='header', component_property='children'),
     Output(component_id='plotly_graph', component_property='figure'),
     ],
    [Input(component_id='slct_chrom', component_property='value'),
     Input(component_id='slct_constraint', component_property='value'),
     ]

)
def update_graph(slct_chrom, slct_constraint):
    header = f"Changes in probabilities of being unpaired (Accessibility) for Gene {slct_chrom} at constraint {slct_constraint}"
    if MODE == "db":
        fig = update_graph_via_sql(slct_chrom, slct_constraint)
    else:
        fig = update_graph_via_pandas(slct_chrom, slct_constraint)
    return [header, fig]


@app.callback(
    [

        Output(component_id="slct_constraint", component_property="value"),
        Output(component_id="slct_chrom", component_property="value"),

    ],
    [
        Input({"type": "interesting-table-button", "index": ALL, "constraint": ALL, "chrom": ALL}, "n_clicks")
    ],
    [
        State(component_id="slct_constraint", component_property="value"),
        State(component_id="slct_chrom", component_property="value"),
     ]

)
def table_click_callback(*args):
    callback_dict = callback_context.triggered[0]["prop_id"].split(".")[0]
    if callback_dict != "":
        callback_dict = eval(callback_dict)
        chrom = callback_dict["chrom"]
        constraint = callback_dict["constraint"]
    else:
        constraint, chrom = args[-2:]
    return [constraint, chrom]



@app.callback(
    [
        Output(component_id="interesting-table-div", component_property="children")
    ],
    [
        Input(component_id="prev-button", component_property="n_clicks"),
        Input(component_id="next-button", component_property="n_clicks"),
    ]

)
def table_switch_callback(prev_clicks, next_clicks):
    page = next_clicks - prev_clicks
    if page < 0:
        next_clicks, prev_clicks, page = 0, 0, 0
    interesting = get_interesting(df, page)
    html_table = interesting_table(interesting, prev_clicks, next_clicks)
    return [html_table]


if __name__ == '__main__':
    from tempfile import TemporaryDirectory
    with TemporaryDirectory() as handle:
        csv_to_sqlite("testfile.bed", os.path.join(handle, "test.db"))
        global df
        #df = read_data("testfile.bed")
        df = os.path.join(handle, "test.db")
        get_app_layout(app, df)
        get_interesting(df)
        app.run_server(debug=True, port=8080, host="0.0.0.0")

