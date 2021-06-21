from __future__ import annotations

import base64
import csv
import gzip
import os
import sqlite3
import zipfile
from tempfile import TemporaryDirectory
from typing import List, Union

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

from RIssmed.RNAtweaks.RIssmedArgparsers import visualiziation_parser

FILEDIR = os.path.dirname(os.path.abspath(__file__))
ASSETS_DIR = os.path.join(FILEDIR, "assets")
app = dash.Dash("FOO", external_stylesheets=[dbc.themes.DARKLY], assets_folder=ASSETS_DIR,
                index_string=open(os.path.join(ASSETS_DIR, "index.html")).read())

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
        cur.execute('INSERT INTO test VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?);', row)
    file_handle.close()
    con.commit()
    con.close()
    insert_interesting_table(db_path)


def insert_interesting_table(db_path: str):
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    cur.execute(f"CREATE TABLE IF NOT EXISTS importance "
                f"( Chr VARCHAR(5), "
                f"Gene_of_interest VARCHAR(20), "
                f"Genomic_Start INT,"
                f"Genomic_End INT, "
                f"Mean_Value REAL, "
                f"Max_Value REAL, "
                f" FOREIGN KEY (Gene_of_interest) REFERENCES test(Gene_of_interest)) ")
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
    with open("../testfile.bed", "w") as handle:
        handle.write("".join(data))


def search_inputs():
    columns = "p-1 col-md-3 col-10 d-flex justify-content-center align-items-center text-center"
    menu = html.Div(
        [
            html.Div(
                html.Button("", className="prev-button", id="prev-button", n_clicks=0),
                className="col-2 col-md-1 text-center order-2 order-md-1 m-1 m-md-0"),
            html.Div(html.Div([
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
                ], className=columns),
                html.Nobr([
                    html.Label(
                        "GOI",
                        htmlFor="search-goi-input",
                        className="search-label"
                    ),
                    dcc.Input(id="search-goi-input",
                              placeholder="Search Gene",
                              className="search-input",
                              )
                ], className=columns),
                html.Nobr([
                    html.Label(
                        "Start",
                        htmlFor="search-span-start-input",
                        className="search-label"
                    ),
                    dcc.Input(id="search-span-start-input",
                              placeholder="Start",
                              className="search-input",
                              )
                ], className=columns),
                html.Nobr([
                    html.Label(
                        "End",
                        htmlFor="search-span-end-input",
                        className="search-label"
                    ),
                    dcc.Input(id="search-span-end-input",
                              placeholder="End",
                              className="search-input",
                              )
                ], className=columns),
            ], className="row justify-content-center justify-content-md-between"), className="col-12 col-md-10 order-md-2 order-1"),
            html.Div(html.Button("", className="next-button", id="next-button", n_clicks=0), className="col-2 col-md-1 text-center order-3 order-md-3 m-1 m-md-0")
        ],
        className="row justify-content-around align-items-center p-2",

    )
    return menu


def interesting_table(interesting: List[sqlite3.Row], prev_clicks: int = 0, next_clicks: int = 0,
                      sorting: str = "Max_Value", sorting_clicks: int = 0):
    page = next_clicks - prev_clicks
    if len(interesting) == 0:
        header = ["WARNING", "No", "matching", "entries"]
        interesting = [(0, "-", "-", "-", "-")]
        clickable = False
    else:
        clickable = True
        header = interesting[0].keys()
        interesting = [(x + 1 + page * NUMBER_OF_INTERESTING, *element) for x, element in enumerate(interesting)]
    menu = [html.Div([
                      html.Div(html.Table(
                          [html.Tr([html.Td(html.Button("#", style={"width": "100%", "pointer-events": "None"}),
                                            className="interesting-table-header")] +
                                   list(table_header_generator(header, sort=sorting, sorting_clicks=sorting_clicks,
                                                               clickable=clickable))
                                   )] +
                          [html.Tr(
                              list(tablerow_generator(entry))
                          ) for y, entry in enumerate(interesting)],
                          className="interesting-table-tablerows",
                          style={"width": "100%", "margin": "auto", 'text-align': 'center'}

                      ), style={"width": "100%", "overflow": "auto"}, className="interesting-table col-12 m-0 p-0"),
                      ],
                     className="col-11 justify-content-center")]
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
            style = {"width": f"{(100 - 17) / (len(row) - 2)}%"}
            classname = "interesting-table-column"
        if type(col) == float:
            col_to_show = round(col, 3)
        else:
            col_to_show = col
            col = ""
        column = html.Td(html.Button(col_to_show, id={"index": f"{x}-{row[1]}-{row[2]}",
                                                      "type": "interesting-table-button",
                                                      "chrom": f"{row[1]}", "goi": row[2],
                                                      "start": row[3],
                                                      "end": row[4]},
                                     n_clicks=0,
                                     className="interesting-table-button",
                                     title=f"{col}",
                                     style={"white-space": "nowrap", "margin": "auto", "width": "100%"}),

                         className=classname, style=style)
        yield column


def modal_image_download():
    modal = dbc.Modal([
                              dbc.ModalHeader("Select Format"),
                              dbc.ModalBody([
                                  html.Div([
                                      html.Button("SVG", id="svg-download", className="btn btn-primary col-3"),
                                      html.Button("PNG", id="png-download", className="btn btn-primary col-3")
                                  ], className="row justify-content-around")
                              ]),
                              dbc.ModalFooter(
                                  dbc.Button(
                                      "Close", id="close", className="ml-auto", n_clicks=0
                                  )
                              ),
                          ], id="modal")
    return modal


def get_app_layout(dash_app: dash.Dash, df: Union[pd.DataFrame, str]):
    if MODE == "db":
        interesting = get_interesting(df)
    else:
        interesting = []

    dash_app.layout = html.Div([
        dcc.Location(id="url", refresh=False),

        html.Div([
            html.Div(html.Div(html.H3("RIssmed Dasboard"), className="databox", style={"text-align": "center"}),
                     className="col-12 p-1 justify-content-center"),
            html.Div([
                html.Div([
                    html.Div([
                        html.H4(id='header', children=[], style={"text-align": "center"}),
                        dcc.Graph(id='plotly_graph', style={"height": "375px"}, className="col-12"),
                        dcc.Download(id="download-image"),
                        html.Button("Download", type="button", id="open", className="btn btn-primary m-2 col-5"),
                        modal_image_download(),
                    ], className="row justify-content-center"),

                ], className="databox",
                         id="graph-box"),
            ], className="col-12 p-1"),

            html.Div([
                html.Div([search_inputs(), html.Div(interesting_table(interesting), id="interesting-table-all", className="row justify-content-center m-1")],
                         className="databox", id="interesting-table-div",  style={"height": "100%"}),

            ], className="col-12 p-1 justify-content-center",  style={"max-height": "50%"}),
        ], className="row justify-content-between", id="100-box"),
        html.Div([],
                 className="row", id="ingo"),
    ],
        id="wrapper", className="container-fluid"
    )


def get_ingo(first_name, name):
    if first_name in ["ingo", "Ingo"] and name in ["Flamingo", "flamingo"]:
        img_path = os.path.join(ASSETS_DIR, "animation.zip")
        with zipfile.ZipFile(img_path) as file:
            img = file.read("animation_2.svg", pwd=b"ingo")
        encoded_img = base64.b64encode(img)
        svg = 'data:image/svg+xml;base64,{}'.format(encoded_img.decode())
        ingo = html.Div(html.Div(
            [
                html.H3("Say Hello to Ingo", style={"text-align": "center"}),
                html.Div(html.Img(src=svg, style={"margin": "auto"}), style={"margin": "auto", "margin-top": "20px",
                                                                             "text-align": "center"})
            ], className="databox"), className="col-12 p-1 justify-content-center", id="ingo-box")
    else:
        ingo = []
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

    fig.update_layout(legend=dict(
        orientation="h",
        yanchor="bottom",
        y=1.02,
        xanchor="right",
        x=1
    ))
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
def table_click_callback(*click_args):
    assert click_args is not None
    callback_dict = callback_context.triggered[0]["prop_id"].split(".")[0]
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
        Output(component_id="ingo", component_property="children"),
        Output(component_id="prev-button", component_property="n_clicks"),
        Output(component_id="next-button", component_property="n_clicks"),
    ],
    [
        Input(component_id="prev-button", component_property="n_clicks"),
        Input(component_id="next-button", component_property="n_clicks"),
        Input({"type": "interesting-table-header-button", "index": ALL, "name": ALL}, "n_clicks"),
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
def table_switch_callback(prev_clicks, next_clicks, inputs, search_chr_input,
                          search_goi, search_span_start, search_span_end, url,
                          last_sort):  # inputs is necessary for callback_context
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
    elif trigger in ["search-input", "search-chr-input", "search-goi-input", "search-span-start-input",
                     "search-span-end-input"]:
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
    search_settings = SearchSettings(search_chr_input, search_goi, search_span_start, search_span_end)
    interesting = get_interesting(database, page, sorting, sorting_clicks, search_settings)
    ingo = get_ingo(search_chr_input, search_goi)
    html_table = interesting_table(interesting, prev_clicks, next_clicks, sorting, sorting_clicks=sorting_clicks)
    return [html_table, sorting, ingo, prev_clicks, next_clicks]


class SearchSettings:
    def __init__(self, chromosome="", goi="", span_start=0, span_end="Infinity"):
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

@app.callback(
    Output("download-image", "data"),
    [Input("svg-download", "n_clicks"),
     Input("png-download", "n_clicks")],
    State(component_id='plotly_graph', component_property='figure'),
    prevent_initial_call=True
)
def download_image(svg, pdf, fig):
    trigger = callback_context.triggered[0]["prop_id"].split(".")[0]
    format = trigger.split("-")[0]
    fig = go.Figure(fig)
    fig.layout.template = "plotly"
    if format == "png":
        with TemporaryDirectory() as handle:
            file_name = os.path.join(handle, "foo.png")
            fig.write_image(file_name, format="png")
            return dcc.send_file(file_name)


    else:
        image = fig.to_image(format=format)
        image = image.decode()
        return dict(content=image, filename=f"foo.{format}")


@app.callback(
    Output("modal", "is_open"),
    [Input("open", "n_clicks"), Input("close", "n_clicks"), Input("png-download", "n_clicks")],
    [State("modal", "is_open")],
)
def toggle_modal(n1, n2, n3, is_open):
    if n1 or n2 or n3:
        return not is_open
    return is_open


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
