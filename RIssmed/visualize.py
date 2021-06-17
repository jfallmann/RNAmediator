from __future__ import annotations
from tempfile import TemporaryDirectory
import base64
import csv
import gzip
import os
import sqlite3
from typing import Dict, List, Union
from RNAtweaks.RIssmedArgparsers import visualiziation_parser
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

FILEDIR = os.path.dirname(os.path.abspath(__file__))
ASSETS_DIR = os.path.join(FILEDIR, "assets")
app = dash.Dash("FOO", external_stylesheets=[dbc.themes.DARKLY], assets_folder=ASSETS_DIR)



MODE = "db"
NUMBER_OF_INTERESTING = 10
TABLE_COLUMNS = 1
PLOTLY_COLORS = px.colors.qualitative.Light24
COLUMN_NAMES = {'Chr': "VARCHAR(20)",
                'Start': "INT",
                'End': "INT",
                'Accessibility_difference': "REAL",
                'Strand': "VARCHAR(2)",
                'Distance_to_constraint': "INT",
                'Accessibility_no_constraint': "REAL",
                'Accessibility_constraint': "REAL",
                'Energy_Difference': "REAL",
                'Kd_change': "REAL",
                'Zscore': "REAL",
                'Genomic_Start': 'INT',
                'Genomic_END': 'INT',
                'Gene_of_interest': 'VARCHAR(20)',
                }
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
    cur.execute("CREATE INDEX goi_index ON test (Gene_of_interest)")
    cur.execute("CREATE INDEX chr_index ON test (Chr)")
    cur.execute("CREATE INDEX gstart_index ON test (Genomic_Start)")
    cur.execute("CREATE INDEX gend_index ON test (Genomic_End)")
    if ".gz" in file:
        file_handle = gzip.open(file, "rt")
    else:
        file_handle = open(file)
    csv_reader = csv.reader(file_handle, delimiter="\t", )
    for row in csv_reader:
        goi, pos, genomic_pos = row[3].split("|")
        gstart, gend = genomic_pos.split("-")
        row = row + [int(gstart), int(gend), goi]
        row.pop(3)
        cur.execute("INSERT INTO test VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?);", row)
    file_handle.close()
    con.commit()
    con.close()
    insert_interesting_table(db_path)


def insert_interesting_table(db_path: str):
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    cur.execute(f'CREATE TABLE IF NOT EXISTS importance '
                f'( Chr VARCHAR(5), '
                f'Gene_of_interest VARCHAR(20), '
                f'Genomic_Start INT,'
                f'Genomic_End INT, '
                f'Mean_Value REAL, '
                f'Max_Value REAL, '
                f' FOREIGN KEY (Gene_of_interest) REFERENCES test(Gene_of_interest)) ')
    cur.execute("SELECT DISTINCT Chr, Gene_of_interest, Genomic_Start, Genomic_End FROM test")
    constraints = cur.fetchall()
    for entry in constraints:
        cur.execute("SELECT Distance_to_constraint, Accessibility_difference, Chr FROM test "
                    "WHERE Chr=? AND Gene_of_interest=? AND Genomic_Start=? AND Genomic_End=?",
                    entry)
        values = cur.fetchall()
        constraint_max = np.max([abs(diff[1]) for diff in values])
        constraint_mean = np.mean([abs(diff[1]) for diff in values])
        cur.execute("INSERT INTO importance VALUES (?, ?, ?, ?, ?, ?)", [entry[0], entry[1], entry[2], entry[3],
                                                                         constraint_mean, constraint_max])
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


def search_inputs():
    menu = html.Div(
        [
            html.Nobr([
                html.Label(
                    "Chr",
                    htmlFor="search-chr-input",
                    className="search-label"

                ),
                dcc.Input(id="search-chr-input",
                          className="search-input",
                          placeholder="Search Chromosome",
                          )
            ]),
            html.Nobr([
                html.Label(
                    "Fold Constraint",
                    htmlFor="search-chr-input",
                    className="search-label"
                ),
                dcc.Input(id="search-input",
                          placeholder="Search Constraint",
                          className="search-input",
                          )
            ]),
            html.Nobr([
                html.Label(
                    "Gene of Interest",
                    htmlFor="search-goi-input",
                    className="search-label"
                ),
                dcc.Input(id="search-goi-input",
                          placeholder="Search Gene",
                          className="search-input",
                          )
            ]),
            html.Nobr([
                html.Label(
                    "Span Start",
                    htmlFor="search-span-start-input",
                    className="search-label"
                ),
                dcc.Input(id="search-span-start-input",
                          placeholder="Search Gene",
                          className="search-input",
                          )
            ]),
            html.Nobr([
                html.Label(
                    "Span End",
                    htmlFor="search-span-end-input",
                    className="search-label"
                ),
                dcc.Input(id="search-span-end-input",
                          placeholder="Search Gene",
                          className="search-input",
                          )
            ]),
         ],
        className="search-wrapper",
        style={

            "display": "table",
            'text-align': 'center'

        }
    )
    return menu


def interesting_table(interesting: List[sqlite3.Row], prev_clicks: int = 0, next_clicks: int = 0,
                      sorting: str = "Max_Value", sorting_clicks: int = 0):
    page = next_clicks - prev_clicks
    if len(interesting) == 0:
        header = ["WARNING", "No matching entries found"]
        interesting = [(0, "-", "-")]
        clickable = False
    else:
        clickable = True
        header = interesting[0].keys()
        interesting = [(x + page*NUMBER_OF_INTERESTING, *element) for x, element in enumerate(interesting)]
    menu = [html.Div([html.Div(html.Button("", id="prev-button", n_clicks=prev_clicks,),
                              style={"display": "table-cell", "margin": "auto", 'align-items': 'center',
                                     "justify-content": "center", "vertical-align": "middle"}),
            html.Div(html.Table(
                [html.Tr([html.Td(html.Button("#", style={"width": "100%", "pointer-events": "None"}), 
                                  className="interesting-table-header")] +
                         list(table_header_generator(header, sort=sorting, sorting_clicks=sorting_clicks,
                                                     clickable=clickable))
                         )] +
                [html.Tr(
                    list(tablerow_generator(entry))
                ) for y, entry in enumerate(interesting)],
                className="interesting-table-tablerows", style={"width": "100%", "margin": "auto", 'text-align': 'center'}

            ), style={"display": "table-cell", "width": "80%"}, className="interesting-table"),
            html.Div(html.Button("", id="next-button", n_clicks=next_clicks,
                                 style={"margin": "auto"}),
                     style={"display": "table-cell", "margin": "auto", 'align-items': 'center',
                            "justify-content": "center", "vertical-align": "middle"})
            ], style={"display": "table", "width": "100%", 'text-align': 'center', "vertical-align": "middle"},
                    className="interesting-table-all", id="interesting-table-all")]
    return menu


def table_header_generator(row, sort, sorting_clicks, clickable: bool = True):
    for element in row:
        if clickable:
            style = {"white-space": "nowrap", "margin": "auto", "width": "100%",
                                            "font-weight": "bold"}
        else:
            style = {"white-space": "nowrap", "margin": "auto", "width": "100%",
                     "font-weight": "bold", "pointer-events": "None", "color": "#FF7676"}
        if element == sort:
            if sorting_clicks % 2:
                order = "\U000025B2"
            else:
                order = "\U000025BC"

        else:
            order = ""
        column = html.Td(html.Button(f"{element.replace('_', ' ')}{order}", n_clicks=sorting_clicks,
                                     id={"index": "sorting",
                                         "type": "interesting-table-header-button",
                                         "name": f"{element}"},
                                     style=style),
                         className="interesting-table-header",
                         style={"font-weight": "bold"})
        yield column


def tablerow_generator(row):
    for x, col in enumerate(row):
        if x == 0:
            style = {"width": "7%"}
            classname = "interesting-table-index-column"
        elif x == 1:
            style = {"width": "10%"}
            classname = "interesting-table-column"
        else:
            style = {"width": f"{(100-17)/(len(row)-2)}%"}
            classname = "interesting-table-column"
        if type(col) == float:
            col_to_show = round(col, 3)
        else:
            col_to_show = col
            col = ""
        column = html.Td(html.Button(col_to_show, id={"index": f"{x}-{row[1]}-{row[2]}", "type": "interesting-table-button",
                                              "chrom": f"{row[1]}", "goi": row[2], "start": row[3], "end": row[4]},
                                     n_clicks=0, className="interesting-table-button", title=f"{col}",
                                     style={"white-space": "nowrap", "margin": "auto", "width": "100%"}),

                         className=classname, style=style)
        yield column


def get_app_layout(app: dash.Dash, df: Union[pd.DataFrame, str]):
    if MODE == "db":
        interesting = get_interesting(df)
    else:
        interesting = []

    app.layout = html.Div([
        dcc.Location(id="url", refresh=False),
        html.Div(html.H3("RIssmed Dasboard"), className="page-header"),

        html.Div([html.H4(id='header', children=[], style={"text-align": "center"}),
                  dcc.Graph(id='plotly_graph', style={"height": "375px"})], className="databox", id="graph-box"),
        html.Div([search_inputs()] + interesting_table(interesting),
                 className="databox", id="interesting-table-div"),
        html.Div([],
                className="", id="ingo")],
        id="wrapper"
    )


def get_ingo():
    img_path = os.path.join(ASSETS_DIR, "animation")
    encoded_img = base64.b64encode(open(img_path, "rb").read())
    svg = 'data:image/svg+xml;base64,{}'.format(encoded_img.decode())
    ingo = html.Div(
        [
            html.H3("Say Hello to Ingo", style={"text-align": "center"}),
            html.Img(src=svg)
            ], style={"position": "relative", "height": "300"}, className="databox")
    return ingo


def get_interesting(db_path: str, page: int = 0, ordering: str = "Max_Value", sorting_clicks: int = 0,
                    substrings: SearchSettings = None):
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    if sorting_clicks % 2:
        deasc = "ASC"
    else:
        deasc = "DESC"
    if substrings is None:
        substrings = SearchSettings()
    print(substrings)
    cur = conn.cursor()
    cur.execute(f"SELECT * FROM importance "
                f"WHERE Chr like '{substrings.chr}%' "
                f"AND Gene_of_interest like '{substrings.goi}%' "
                f"AND Genomic_Start >= ? "
                f"AND Genomic_End <= ? "
                f"ORDER BY {ordering} {deasc} "
                f"LIMIT {NUMBER_OF_INTERESTING} OFFSET {NUMBER_OF_INTERESTING * page}",
                (substrings.span_start, substrings.span_end))
    return_list = cur.fetchall()
    conn.close()
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


def update_graph_via_sql(slct_chrom, slct_goi, slct_start, slct_end):
    conn = sqlite3.connect(database)
    cur = conn.cursor()
    cur.execute("SELECT Distance_to_constraint, Accessibility_difference, Accessibility_no_constraint, "
                "Accessibility_constraint FROM test WHERE Gene_of_interest=? AND Genomic_Start=? AND "
                "Genomic_End=? AND Chr=?",
                (slct_goi, slct_start, slct_end, slct_chrom))
    rows = cur.fetchall()
    if len(rows) > 0:
        distance, acc_diff, acc_no_const, acc_cons = zip(*rows)
    else:
        distance = acc_diff = acc_no_const = acc_cons = []
    test = set(distance)
    distance = list(distance)
    acc_diff = list(acc_diff)
    acc_no_const = list(acc_no_const)
    acc_cons = list(acc_cons)
    for x in range(len(test)):
        if x not in test:
            distance.append(x)
            acc_diff.append("")
            acc_no_const.append("")
            acc_cons.append("")
        if -x not in test:
            distance.append(-x)
            acc_diff.append("")
            acc_no_const.append("")
            acc_cons.append("")
        if x in test and -x in test:
            break
    if len(test) > 0:
        sorted_data = sorted(zip(distance, acc_diff, acc_no_const, acc_cons))
        distance, acc_diff, acc_no_const, acc_cons = zip(*sorted_data)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=distance, y=acc_diff, line={"width": 4, "color": PLOTLY_COLORS[7]},
                             name="Accessibility Difference", visible="legendonly", connectgaps=False))
    fig.add_trace(go.Scatter(x=distance, y=acc_no_const, line={"width": 4, "color": PLOTLY_COLORS[1]},
                             name="Accessibility no constraint", connectgaps=False))
    fig.add_trace(go.Scatter(x=distance, y=acc_cons, line={"width": 4, "color": PLOTLY_COLORS[0]},
                             name="Accessibility with constraint", connectgaps=False))
    fig.layout.template = "plotly_white"
    fig.update_yaxes(range=[-1, 1], title="Probability of being unpaired")
    x_range = np.max(np.abs(distance)) if len(distance) > 0 else 120
    fig.update_xaxes(range=[-x_range, x_range], title="Distance to constraint")
    conn.close()
    return fig


def update_graph_via_pandas(slct_chrom, slct_constraint):
    dff = database.copy()
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
    [
     Output(component_id='plotly_graph', component_property='figure'),
     ],
    [Input(component_id='header', component_property='children'),
     ]

)
def update_graph(header):
    slct_chrom, slct_goi, slct_start, slct_end = header.split(" ")
    if MODE == "db":
        fig = update_graph_via_sql(slct_chrom, slct_goi, slct_start, slct_end)
    else:
        raise NotImplementedError
    return [fig]


@app.callback(
    [

        Output(component_id="header", component_property="children"),


    ],
    [
        Input({"type": "interesting-table-button", "index": ALL, "goi": ALL, "chrom": ALL,
               "start": ALL, "end": ALL}, "n_clicks")
    ],


)
def table_click_callback(*args):
    callback_dict = callback_context.triggered[0]["prop_id"].split(".")[0]
    x = callback_context
    if callback_dict != "":
        callback_dict = eval(callback_dict)
        chrom = callback_dict["chrom"]
        goi = callback_dict["goi"]
        start = callback_dict["start"]
        end = callback_dict["end"]
    else:
        goi = start = end = chrom = ""
    return [f"{chrom} {goi} {start} {end}"]



@app.callback(
    [
        Output(component_id="interesting-table-all", component_property="children"),
        Output(component_id="url", component_property="pathname"),
        Output(component_id="ingo", component_property="children")
    ],
    [
        Input(component_id="prev-button", component_property="n_clicks"),
        Input(component_id="next-button", component_property="n_clicks"),
        Input({"type": "interesting-table-header-button", "index": ALL, "name": ALL}, "n_clicks"),
        Input(component_id="search-input", component_property="value"),
        Input(component_id="search-chr-input", component_property="value"),
        Input(component_id="search-goi-input", component_property="value"),
        Input(component_id="search-span-start-input", component_property="value"),
        Input(component_id="search-span-end-input", component_property="value"),
    ],
    [
        State(component_id="url", component_property="pathname"),
        State(component_id={"index": "sorting", "name": ALL, "type": ALL}, component_property="n_clicks"),
    ]

)
def table_switch_callback(prev_clicks, next_clicks, inputs, search_input, search_chr_input,
                          search_goi, search_span_start, search_span_end, url, last_sort): # inputs is necessary for callback_context
    trigger = callback_context.triggered[0]["prop_id"].split(".")[0]
    if trigger in ["next-button", "prev-button", "", "url"]:
        page = next_clicks - prev_clicks
        url = url[1:] if url.startswith("/") else url
        if page < 0:
            next_clicks = prev_clicks = page = 0
        if url == "" or url == "/":
            sorting = "Max_Value"
        else:
            sorting = url
        sorting_clicks = last_sort[0]
    elif trigger in ["search-input", "search-chr-input", "search-goi-input", "search-span-start-input", "search-span-end-input"]:
        if url == "" or url == "/":
            sorting = "Max_Value"
        else:
            sorting = url
        next_clicks = prev_clicks = page = 0
        sorting_clicks = 0

    else:
        sorting_clicks = callback_context.triggered[0]["value"]
        trigger_dict = eval(trigger)
        next_clicks = prev_clicks = page = 0
        sorting = trigger_dict["name"]
    search_settings = SearchSettings(search_input, search_chr_input, search_goi, search_span_start, search_span_end)
    interesting = get_interesting(database, page, sorting, sorting_clicks, search_settings)
    if search_chr_input in ["ingo", "Ingo"] and search_goi in ["Flamingo", "flamingo"]:
        ingo = get_ingo()
    else:
        ingo = []
    html_table = interesting_table(interesting, prev_clicks, next_clicks, sorting, sorting_clicks=sorting_clicks)
    return [html_table, sorting, ingo]


class SearchSettings:
    def __init__(self, fold_constraint="", chromosome="", goi="", span_start=0, span_end="Infinity"):
        self.fold_constraint = fold_constraint if fold_constraint is not None else ""
        self.chr = chromosome if chromosome is not None else ""
        self.goi = goi if goi is not None else ""
        try:
            self.span_start = int(span_start)
        except (ValueError, TypeError):
            self.span_start = 0
        try:
            self.span_end = int(span_end)
        except (ValueError, TypeError):
            self.span_end = "Infinity"

    def __str__(self):
        return str(self.__dict__)


if __name__ == '__main__':
    args = visualiziation_parser()
    bed_file = args.file
    global database
    if args.tmp is False:
        database = args.database
        tmpdir = None
    else:
        tmpdir = TemporaryDirectory()
        database = os.path.join(tmpdir.name, "rissmed.db")

    if not os.path.exists(database):
        csv_to_sqlite(bed_file, database)
    get_app_layout(app, database)
    app.run_server(debug=True, port=8080, host="0.0.0.0")
    if tmpdir:
        tmpdir.cleanup()

